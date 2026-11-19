# Figures 1 to 3 and supplemental figures S1 and S2.

blue <- "#1F5C8B"
orange <- "#C0562E"
green <- "#2E7D64"

theme_paper <- theme_classic(base_size = 9) +
  theme(axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3),
        plot.title = element_text(face = "bold", size = 11, hjust = 0),
        legend.title = element_blank(),
        legend.key.height = unit(10, "pt"),
        plot.margin = margin(4, 6, 4, 4))

save_figure <- function(plot, name, width, height) {
  ggsave(file.path("output", paste0(name, ".pdf")), plot,
         width = width, height = height)
  ggsave(file.path("output", paste0(name, ".png")), plot,
         width = width, height = height, dpi = 300)
}

# Figures 2 and S1 show change from each patient's first inlet concentration.
# The models in 03 use the measured log concentrations before this rescaling.
change <- samples |>
  group_by(patient) |>
  mutate(change = 100 * (c_pre / first(c_pre) - 1)) |>
  ungroup()

# Figure 1. Mean concentrations and extraction ratio with 95% t confidence
# intervals at each sampling point, plotted at the mean recorded time.
mean_by_sample <- function(column) {
  samples |>
    group_by(sample) |>
    summarise(time = mean(time), n = n(),
              mean = mean(.data[[column]]), sd = sd(.data[[column]])) |>
    mutate(lower = mean - qt(0.975, n - 1) * sd / sqrt(n),
           upper = mean + qt(0.975, n - 1) * sd / sqrt(n))
}

# Offset the two series slightly on the time axis so both intervals are visible.
concentrations <- bind_rows(
  mutate(mean_by_sample("c_pre"), series = "Pre-dialyzer",
         time = time - 0.035),
  mutate(mean_by_sample("c_post"), series = "Post-dialyzer",
         time = time + 0.035)
) |>
  mutate(series = factor(series, levels = c("Pre-dialyzer", "Post-dialyzer")))

panel_a <- ggplot(concentrations, aes(time, mean, colour = series)) +
  geom_line(linewidth = 0.65) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.06,
                linewidth = 0.45) +
  geom_point(size = 1.7) +
  scale_colour_manual(values = c("Pre-dialyzer" = blue,
                                 "Post-dialyzer" = orange)) +
  scale_x_continuous(breaks = 0:4, expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 50, 10), expand = c(0, 0)) +
  # The upper end of the first pre-dialyzer interval runs past the panel.
  coord_cartesian(xlim = c(-0.15, 4.15), ylim = c(-2, 52)) +
  labs(title = "A", x = "Time from dialysis start (h)",
       y = "Plasma TXA concentration (mg/L)") +
  theme_paper +
  theme(legend.position = "inside", legend.position.inside = c(0.78, 0.93))

panel_b <- ggplot(mean_by_sample("er"), aes(time, mean)) +
  geom_hline(yintercept = 0, linetype = "dotted", linewidth = 0.3) +
  geom_line(colour = green, linewidth = 0.65) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.06,
                linewidth = 0.45, colour = green) +
  geom_point(colour = green, size = 1.7) +
  scale_x_continuous(breaks = 0:4, expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2), expand = c(0, 0)) +
  coord_cartesian(xlim = c(-0.15, 4.15), ylim = c(0, 1)) +
  labs(title = "B", x = "Time from dialysis start (h)",
       y = "Dialyzer extraction ratio") +
  theme_paper

save_figure(panel_a + panel_b, "figure_1", width = 7.4, height = 3.3)

# Figure 2. Observed change with the selected power-time curve, its 95%
# confidence band and the 90% prediction band for patient curves.
markers <- data.frame(
  time = c(halving_time, 3.5),
  change = -100 * c(0.5, rr_3.5h[["power"]]),
  label = c(sprintf("50%% reduction at %.1f h", halving_time),
            sprintf("%.0f%% reduction at 3.5 h", 100 * rr_3.5h[["power"]])),
  label_x = c(halving_time + 0.08, 3.5 - 0.08),
  hjust = c(0, 1)
)
curve_label <- bquote(Delta * C(t) == 100 * group("{",
  exp(-.(round(k, 3)) * t^.(round(p, 3))) - 1, "}"))

figure_2 <- ggplot(change, aes(time, change, group = patient)) +
  geom_ribbon(data = curve_band, aes(time, ymin = pi_low, ymax = pi_high,
                                     fill = "prediction"),
              inherit.aes = FALSE, alpha = 0.70) +
  geom_ribbon(data = curve_band, aes(time, ymin = ci_low, ymax = ci_high,
                                     fill = "confidence"),
              inherit.aes = FALSE, alpha = 0.82) +
  geom_line(colour = "grey48", linewidth = 0.55, alpha = 0.85) +
  geom_point(shape = 21, colour = "grey38", fill = "white", stroke = 0.5,
             size = 1.8, alpha = 0.95) +
  geom_line(data = curve_band, aes(time, change, colour = "model"),
            inherit.aes = FALSE, linewidth = 1.05) +
  geom_segment(data = markers,
               aes(x = 0, xend = time, y = change, yend = change),
               inherit.aes = FALSE, linetype = "dotted", colour = "grey42",
               linewidth = 0.4) +
  geom_segment(data = markers,
               aes(x = time, xend = time, y = change, yend = -100),
               inherit.aes = FALSE, linetype = "dotted", colour = "grey42",
               linewidth = 0.4) +
  geom_point(data = markers, aes(time, change), inherit.aes = FALSE,
             shape = 21, size = 2.2, stroke = 0.55, colour = blue,
             fill = "white") +
  geom_text(data = markers, aes(label_x, -90, label = label, hjust = hjust),
            inherit.aes = FALSE, size = 3.0, fontface = "bold",
            colour = "grey22") +
  scale_colour_manual(values = c(model = blue), labels = list(curve_label)) +
  scale_fill_manual(
    values = c(confidence = "#8FB1C8", prediction = "#D9E6EF"),
    breaks = c("confidence", "prediction"),
    labels = c("95% model-curve confidence interval",
               "90% patient prediction interval")
  ) +
  scale_x_continuous(breaks = 0:4, expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(-100, 0, 20), expand = c(0, 0)) +
  coord_cartesian(xlim = c(0, 4.12), ylim = c(-100, 0)) +
  labs(x = "Time from dialysis start (h)",
       y = "Pre-dialyzer TXA change from baseline (%)") +
  theme_paper +
  theme(axis.text = element_text(size = 10.5),
        legend.position = "bottom",
        legend.text = element_text(size = 8.5),
        legend.key.width = unit(22, "pt"),
        legend.box = "vertical") +
  guides(colour = guide_legend(order = 1), fill = guide_legend(order = 2))

save_figure(figure_2, "figure_2", width = 7, height = 4.8)

# Figure 3. Each point represents one resample of 10 patients. Compare their
# mean observed reduction with the prediction at their mean sampling duration.
calibration_lines <- bind_rows(
  data.frame(curve = "Bootstrap calibration", x = calibration_curve$predicted,
             y = calibration_curve$observed),
  data.frame(curve = "Optimism-corrected", x = calibration_curve$predicted,
             y = calibration_curve$corrected)
)
points_range <- 100 * range(selected$predicted_rr, selected$observed_rr)
limits <- c(floor((points_range[1] - 1) / 2) * 2,
            ceiling((points_range[2] + 1) / 2) * 2)

figure_3 <- ggplot() +
  geom_point(data = selected, aes(100 * predicted_rr, 100 * observed_rr),
             colour = "grey55", size = 0.9, alpha = 0.22) +
  geom_abline(data = data.frame(curve = "Ideal"),
              aes(intercept = 0, slope = 1, colour = curve, linetype = curve),
              linewidth = 0.8) +
  # The dashed bootstrap curve is drawn over the solid corrected curve so
  # that both stay visible.
  geom_line(data = filter(calibration_lines, curve == "Optimism-corrected"),
            aes(100 * x, 100 * y, colour = curve, linetype = curve),
            linewidth = 0.95) +
  geom_line(data = filter(calibration_lines, curve == "Bootstrap calibration"),
            aes(100 * x, 100 * y, colour = curve, linetype = curve),
            linewidth = 1.05) +
  scale_colour_manual(
    breaks = c("Ideal", "Bootstrap calibration", "Optimism-corrected"),
    values = c(Ideal = "grey35", "Bootstrap calibration" = orange,
               "Optimism-corrected" = blue)
  ) +
  scale_linetype_manual(
    breaks = c("Ideal", "Bootstrap calibration", "Optimism-corrected"),
    values = c(Ideal = "dotted", "Bootstrap calibration" = "longdash",
               "Optimism-corrected" = "solid")
  ) +
  scale_x_continuous(breaks = pretty(limits, n = 6), expand = c(0, 0)) +
  scale_y_continuous(breaks = pretty(limits, n = 6), expand = c(0, 0)) +
  coord_equal(xlim = limits, ylim = limits) +
  labs(x = "Model-predicted mean first-to-final RR (%)",
       y = "Observed mean first-to-final RR (%)") +
  theme_paper +
  theme(legend.position = "bottom", legend.text = element_text(size = 8),
        legend.key.width = unit(24, "pt"))

save_figure(figure_3, "figure_3", width = 5.8, height = 5.4)

# Figure S1. The 3 candidate concentration-time curves.
candidate_names <- c(constant = "Constant-rate exponential",
                     power = "Power-time (selected)",
                     asymptotic = "Asymptotic exponential")
candidates <- expand.grid(time = seq(0, 4, length.out = 201),
                          model = c("power", "asymptotic", "constant"),
                          stringsAsFactors = FALSE)
candidates$change <- mapply(function(t, m) {
  100 * (exp(log_ratio(m, fixef(fits[[m]]), t)) - 1)
}, candidates$time, candidates$model)
candidates$model <- factor(candidates$model,
                           levels = c("power", "asymptotic", "constant"))
candidate_labels <- sprintf("%s: RR at 3.5 h %.1f%%",
                            candidate_names[levels(candidates$model)],
                            100 * rr_3.5h[levels(candidates$model)])

figure_s1 <- ggplot(change, aes(time, change, group = patient)) +
  geom_line(colour = "grey78", linewidth = 0.45, alpha = 0.8) +
  geom_point(shape = 21, colour = "grey58", fill = "white", stroke = 0.4,
             size = 1.5, alpha = 0.85) +
  geom_vline(xintercept = 3.5, linetype = "dotted", colour = "grey50",
             linewidth = 0.4) +
  geom_line(data = candidates,
            aes(time, change, colour = model, linetype = model),
            inherit.aes = FALSE, linewidth = 1.0) +
  annotate("text", x = 3.5, y = -2, label = "3.5 h", hjust = -0.12, vjust = 0,
           size = 2.7, colour = "grey35") +
  scale_colour_manual(values = c(power = blue, asymptotic = green,
                                 constant = orange),
                      labels = candidate_labels) +
  scale_linetype_manual(values = c(power = "solid", asymptotic = "longdash",
                                   constant = "dotdash"),
                        labels = candidate_labels) +
  scale_x_continuous(breaks = 0:4, expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(-100, 0, 20), expand = c(0, 0)) +
  coord_cartesian(xlim = c(0, 4.12), ylim = c(-100, 0)) +
  labs(x = "Time from dialysis start (h)",
       y = "Pre-dialyzer TXA change from baseline (%)") +
  theme_paper +
  theme(legend.position = "bottom", legend.text = element_text(size = 8),
        legend.key.width = unit(24, "pt")) +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE),
         linetype = guide_legend(nrow = 2, byrow = TRUE))

save_figure(figure_s1, "figure_s1", width = 7.4, height = 4.8)

# Figure S2. Published models from the intravenous dose to the first sample,
# for the representative patient.
pk_plot <- mutate(pk_curves,
                  renal_scenario = factor(renal_scenario,
                                          levels = c("low", "20", "40", "90")))
# Labels: the observed mean on the left of the diamond, and the implied
# renal input or the zero-elimination boundary above the low curve.
observed <- data.frame(
  model = names(pk_models),
  label = sprintf("atop('Observed mean at %.1f h', '%.1f mg/L')",
                  representative$hours, representative$c_first)
)
implied_labels <- pk_plot |>
  filter(renal_scenario == "low") |>
  group_by(model) |>
  slice_min(abs(hours - 9), n = 1, with_ties = FALSE) |>
  ungroup() |>
  left_join(select(renal_fit, model, boundary_at_zero), by = "model") |>
  mutate(label = if_else(
    boundary_at_zero, "atop('Zero systemic', 'TXA clearance')",
    sprintf("atop('Model-implied %s', '%.1f')",
            if_else(model == "model_1", "CrCl", "eGFR"), input_used)
  ))

figure_s2 <- ggplot(pk_plot,
                    aes(hours, concentration, colour = renal_scenario)) +
  geom_line(linewidth = 1.0) +
  annotate("segment", x = 0, xend = representative$hours,
           y = representative$c_first, yend = representative$c_first,
           linetype = "dotted", colour = "grey45", linewidth = 0.45) +
  annotate("point", x = representative$hours, y = representative$c_first,
           shape = 23, size = 3.2, stroke = 0.65, colour = "black",
           fill = "black") +
  geom_text(data = observed,
            aes(x = representative$hours * 0.98, y = 70, label = label),
            inherit.aes = FALSE, hjust = 1, size = 2.7, lineheight = 0.95,
            colour = "grey20", parse = TRUE) +
  geom_text(data = implied_labels,
            aes(hours, concentration + 20, label = label),
            inherit.aes = FALSE, size = 2.7, lineheight = 0.95,
            fontface = "bold", colour = "#6A3D9A", parse = TRUE) +
  scale_colour_manual(
    name = "Renal input",
    values = c("#6A3D9A", "#D55E00", "#009E73", "#0072B2"),
    labels = c("Low elimination", "20", "40", "90")
  ) +
  facet_wrap(vars(model), nrow = 1,
             labeller = as_labeller(c(
               model_1 = "PK model 1 · CrCl (mL/min)",
               model_2 = "PK model 2 · eGFR (mL/min/1.73 m²)"))) +
  scale_x_continuous(breaks = seq(0, 24, 4),
                     expand = expansion(mult = c(0, 0.02))) +
  scale_y_continuous(breaks = seq(0, 250, 50),
                     expand = expansion(mult = c(0, 0.03))) +
  coord_cartesian(ylim = c(0, 250)) +
  labs(x = "Time since intravenous dose (h)",
       y = "Plasma TXA concentration (mg/L)") +
  theme_paper +
  theme(legend.position = "bottom",
        legend.title = element_text(face = "bold", size = 8.5),
        legend.text = element_text(size = 7.8),
        legend.key.width = unit(34, "pt"),
        panel.spacing = unit(14, "pt"),
        strip.background = element_rect(fill = "grey93", colour = NA),
        strip.text = element_text(face = "bold", size = 9.5))

save_figure(figure_s2, "figure_s2", width = 7.4, height = 4.8)
