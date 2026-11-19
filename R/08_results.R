# Format the analysis results for comparison with the manuscript and tables.
# Rounding happens here, after the estimates and intervals are calculated.

fmt <- function(x, digits = 0) {
  x[abs(x) < 0.5 * 10^-digits] <- 0   # print 0 rather than -0
  sprintf(paste0("%.", digits, "f"), x)
}
ci_text <- function(ci, digits = 0) {
  sprintf("95%% CI, %s to %s", fmt(ci[1], digits), fmt(ci[2], digits))
}
fmt_ci <- function(estimate, ci, digits = 0) {
  sprintf("%s (%s)", fmt(estimate, digits), ci_text(ci, digits))
}
# For the output of t_summary().
fmt_mean <- function(s, digits = 0, scale = 1) {
  fmt_ci(scale * s[["mean"]], scale * s[c("lower", "upper")], digits)
}
fmt_sd <- function(x, digits = 0) {
  paste(fmt(mean(x), digits), "\u00b1", fmt(sd(x), digits))
}
fmt_range <- function(x, digits = 0) {
  paste(fmt(min(x), digits), "to", fmt(max(x), digits))
}
fmt_pair <- function(x, digits = 0) paste(fmt(x, digits), collapse = " and ")

md_table <- function(df) {
  c(paste("|", paste(names(df), collapse = " | "), "|"),
    paste0("|", strrep(" --- |", ncol(df))),
    apply(df, 1, function(row) paste("|", paste(row, collapse = " | "), "|")),
    "")
}
section <- function(title, ...) {
  values <- c(...)
  c(paste("##", title), "",
    md_table(data.frame(Result = names(values), Value = unname(values))))
}

rr_ci <- function(model, digits) {
  fmt_ci(100 * rr_3.5h[[model]], 100 * rr_3.5h_ci[, model], digits)
}
boot_range <- function(x, digits, scale = 1) {
  fmt_range(scale * percentile_ci(x), digits)
}
both_77 <- filter(table_s2, scenario == "Both exchange rates reduced by 77%")

participants <- with(by_patient, section(
  "Participants and dialysis (Results, figure 1 legend)",
  "Patients" = nrow(by_patient),
  "Paired samples" = nrow(samples),
  "Patients with 7, 6 and 5 pairs" =
    paste(table(factor(pairs, levels = 7:5)), collapse = ", "),
  "Patients contributing at each sampling point" =
    paste(table(samples$sample), collapse = ", "),
  "Pairs with a lower post-dialyzer than pre-dialyzer concentration" =
    sum(samples$c_post < samples$c_pre),
  "Weight, kg" = fmt_sd(weight),
  "Intravenous tranexamic acid, mg" = fmt_sd(dose_iv),
  "Intravenous tranexamic acid, range, mg" = fmt_range(dose_iv),
  "Topical tranexamic acid, mg (patient)" =
    sprintf("%s (%d)", fmt(sum(dose_topical)), patient[dose_topical > 0]),
  "Oral tranexamic acid, mg (patient)" =
    sprintf("%s (%d)", fmt(sum(dose_oral)), patient[dose_oral > 0]),
  "Intravenous dose to first sample, h" =
    fmt_sd(dose_to_first_sample_min / 60, 1),
  "Baseline pre-dialyzer concentration, mg/L" = fmt_sd(c_first, 1),
  "Baseline pre-dialyzer concentration, range, mg/L" = fmt_range(c_first, 1),
  "Final pre-dialyzer concentration, mg/L" = fmt_sd(c_last, 1),
  "Final pre-dialyzer concentration, range, mg/L" = fmt_range(c_last, 1),
  "Pre-dialysis hematocrit, %" = fmt_sd(hct, 1),
  "Blood flow at sampling, mL/min" = fmt_sd(qb),
  "Patients with dialysate flow 500 and 700 mL/min" =
    paste(sum(qd == 500), "and", sum(qd == 700)),
  "Ultrafiltration rate at sampling, mL/h" = fmt_sd(quf),
  "First to final pair, h" = fmt_sd(duration, 1),
  "First to final pair, range, h" = fmt_range(duration, 1),
  "Dialyzer models" = n_distinct(dialyzer),
  "Patients with FX-class and Theranova dialyzers" =
    paste(sum(grepl("^FX", dialyzer)), "and",
          sum(grepl("^Theranova", dialyzer)))
))

outcomes <- section(
  "Clearance, extraction and reduction (Abstract, Results)",
  "Dialytic clearance, mL/min" = fmt_mean(clearance),
  "One-sided p value for mean clearance above 0" =
    format.pval(primary_test$p.value, eps = 0.001),
  "Extraction ratio" = fmt_mean(extraction, 2),
  "Extraction ratio, range of patient means" = fmt_range(by_patient$er, 2),
  "Observed first-to-final reduction ratio, %" = fmt_mean(reduction, 0, 100),
  "Observed reduction ratio, range, %" = fmt_range(100 * by_patient$rr)
)

concentration_models <- section(
  "Concentration-time models (Results, figure 2, figure S1)",
  "AIC, power time" = fmt(aic[["power"]], 1),
  "AIC, asymptotic" = fmt(aic[["asymptotic"]], 1),
  "AIC, constant rate" = fmt(aic[["constant"]], 1),
  "Bootstrap samples selecting power time, asymptotic, constant rate" =
    paste(times_selected[c("power", "asymptotic", "constant")],
          collapse = ", "),
  "Bootstrap fits that did not converge, of 3000" = failed_fits,
  "Power time: k and p" = fmt_pair(c(k, p), 3),
  "Constant rate: k" = fmt(exp(fixef(fits$constant)[["log_k"]]), 3),
  "Asymptotic: d and r" =
    fmt_pair(exp(fixef(fits$asymptotic)[c("log_d", "log_r")]), 3),
  "Halving time, h" = fmt_ci(halving_time, halving_time_ci, 1),
  "RR at 3.5 h, power time, %" = rr_ci("power", 1),
  "RR at 3.5 h, constant rate, %" = rr_ci("constant", 1),
  "RR at 3.5 h, asymptotic, %" = rr_ci("asymptotic", 1)
)

validation <- section(
  "Internal validation (Results, figure 3 legend)",
  "Optimism-corrected mean prediction error" =
    fmt(corrected[["mean_error"]], 2),
  "Optimism-corrected root mean square error" = fmt(corrected[["rmse"]], 2),
  "Optimism-corrected calibration intercept" = fmt(corrected[["intercept"]], 2),
  "Optimism-corrected calibration slope" =
    fmt_ci(corrected[["slope"]], corrected_slope_ci, 2),
  "Mean first-to-final duration, h" = fmt(original[["duration"]], 2),
  "Mean observed first-to-final RR, %" = fmt(100 * original[["rr"]], 1),
  "Predicted RR at that duration, %" = fmt(100 * original_predicted_rr, 1),
  "Mean optimism in calibration error, percentage points" =
    fmt(100 * calibration_optimism, 2),
  "Bootstrap 95% interval, mean duration, h" = boot_range(selected$duration, 2),
  "Bootstrap 95% interval, observed mean RR, %" =
    boot_range(selected$observed_rr, 1, 100),
  "Bootstrap 95% interval, predicted mean RR, %" =
    boot_range(selected$predicted_rr, 1, 100)
)

sensitivity <- section(
  "Sensitivity analyses (Results)",
  "Clearance without the first pair, mL/min" = fmt_mean(clearance_later),
  "Extraction ratio without the first pair" = fmt_mean(extraction_later, 2),
  "Rebound 2 h after a 3.5 h session, %" = fmt(mean(rebound_pct), 1),
  "RR at 3.5 h before rebound, %" = rr_ci("power", 0),
  "RR 2 h after rebound, %" =
    fmt_ci(100 * mean(rr_after_rebound), 100 * rr_after_rebound_ci),
  "Clearance with rising hematocrit, mL/min" = fmt_mean(clearance_hct),
  "Patients without oral or topical dosing" = nrow(iv_only),
  "Their mean first concentration, mg/L" = fmt(representative$c_first, 1),
  "Model 1 implied creatinine clearance, mL/min" =
    fmt(renal_fit$implied_input[renal_fit$model == "model_1"], 2),
  "Model 2 implied indexed 2021 CKD-EPI creatinine eGFR, mL/min/1.73 m²" =
    if (oracle_fit$boundary_at_zero)
      "No nonnegative eGFR reproduced the observed mean"
    else fmt(oracle_fit$implied_input, 2),
  "Cumulative loss from plasma, mg" = fmt_mean(loss_summary[, "primary"]),
  "Observed interval, h" = fmt(mean(plasma_loss$duration), 1),
  "Observed interval, range, h" = fmt_range(plasma_loss$duration, 1),
  "Loss with fitted curves, mg" = fmt_mean(loss_summary[, "fitted_curve"]),
  "Loss with rising hematocrit, mg" = fmt_mean(loss_summary[, "hct"])
)

discussion <- section(
  "Discussion",
  "Model 2, both exchange rates reduced by 77%: reduction at 2 and 24 h, %" =
    fmt_pair(both_77$reduction[both_77$hours %in% c(2, 24)], 1)
)

# Table 1. Mark additional topical doses with a dagger and oral doses
# with an asterisk.
mark <- with(by_patient, ifelse(dose_topical > 0, "\u2020",
                                ifelse(dose_oral > 0, "*", "")))
table_1 <- with(by_patient, data.frame(
  "Patient" = as.character(patient),
  "Weight, kg" = fmt(weight, 1),
  "Dialyzer" = dialyzer,
  "Hct, %" = fmt(hct),
  "Mean QB, mL/min" = fmt(qb),
  "QD, mL/min" = fmt(qd),
  "Mean QUF, mL/h" = fmt(quf),
  "Duration to final sample, h" = fmt(duration, 1),
  "IV TXA, mg" = paste0(fmt(dose_iv), mark),
  "IV dose to baseline, h" = fmt(dose_to_first_sample_min / 60, 1),
  "Pairs" = as.character(pairs),
  "Baseline TXA, mg/L" = fmt(c_first, 1),
  "Final TXA, mg/L" = fmt(c_last, 1),
  "Clearance, mL/min" = fmt(cl),
  "ER" = fmt(er, 2),
  "RR, %" = fmt(100 * rr),
  check.names = FALSE
))
table_1[nrow(table_1) + 1, ] <- with(by_patient, c(
  "Overall", fmt_sd(weight, 1), "", fmt_sd(hct, 1), fmt_sd(qb), "",
  fmt_sd(quf), fmt_sd(duration, 1), fmt_sd(dose_iv),
  fmt_sd(dose_to_first_sample_min / 60, 1), paste(sum(pairs), "total"),
  fmt_sd(c_first, 1), fmt_sd(c_last, 1),
  sprintf("%s (%s)", fmt_sd(cl), ci_text(clearance[c("lower", "upper")])),
  sprintf("%s (%s)", fmt_sd(er, 2),
          ci_text(extraction[c("lower", "upper")], 2)),
  sprintf("%s (%s)", fmt_sd(100 * rr),
          ci_text(100 * reduction[c("lower", "upper")]))
))

supplement <- section(
  "Supplemental Digital Content 1",
  "Plasma loss, range, mg" = fmt_range(plasma_loss$primary),
  "Plasma loss, mean observed interval, h" = fmt(mean(plasma_loss$duration), 2),
  "Plasma loss, observed interval, range, h" =
    fmt_range(plasma_loss$duration, 2),
  "Change in mean loss with rising hematocrit, mg" =
    fmt_mean(loss_summary[, "hct_change"]),
  "Change in mean loss with fitted curves, mg" =
    fmt_mean(loss_summary[, "curve_change"]),
  "Change in mean loss using the first or the last tied sample, mg" =
    fmt_pair(loss_summary["mean", c("tied_first", "tied_last")] -
               loss_summary["mean", "primary"], 2),
  "Model 1 rebound: concentration at the end of dialysis and 2 h later, mg/L" =
    fmt_pair(c(mean(c_end_dialysis), mean(c_after_rebound)), 1),
  "Representative patient: dose, mg" = fmt(representative$dose, 1),
  "Representative patient: weight, kg" = fmt(representative$weight, 2),
  "Representative patient: dose to first sample, h" =
    fmt(representative$hours, 2),
  "Representative patient: first concentration, mg/L" =
    fmt(representative$c_first, 2),
  "Model 2: pre-dialysis systemic TXA clearance used for rebound, L/h" =
    fmt(m2$cl, 2),
  "Model 2: predicted minus observed first concentration, mg/L" =
    fmt(oracle_fit$predicted_first - representative$c_first, 2),
  "Model 2: central, fast and slow peripheral concentrations, mg/L" =
    paste(fmt(c_start, 2), collapse = ", "),
  "Figure S2 legend: dose (mg), weight (kg) and dose to first sample (h)" =
    paste(fmt(representative$dose), fmt(representative$weight, 1),
          fmt(representative$hours, 1), sep = ", ")
)

model_names <- c(constant = "Constant rate", power = "Power time",
                 asymptotic = "Asymptotic")
variance_table <- data.frame(
  "Model" = model_names[rownames(variance_components)],
  "SD, log baseline" = fmt(variance_components[, "sd_baseline"], 3),
  "SD, log decline parameter" = fmt(variance_components[, "sd_decline"], 3),
  "Correlation" = fmt(variance_components[, "correlation"], 3),
  "Residual SD" = fmt(variance_components[, "sd_residual"], 3),
  check.names = FALSE
)

table_s1 <- with(plasma_loss, data.frame(
  "Patient" = c(patient, "Mean"),
  "Observed interval, h" = fmt(c(duration, mean(duration)), 2),
  "Primary, mg" = c(fmt(primary), fmt_mean(loss_summary[, "primary"])),
  "Hematocrit sensitivity, mg" = c(fmt(hct), fmt_mean(loss_summary[, "hct"])),
  "Latent curve, mg" =
    c(fmt(fitted_curve), fmt_mean(loss_summary[, "fitted_curve"])),
  check.names = FALSE
))

cells <- paste(fmt(table_s2$rebound, 1), "/", fmt(table_s2$reduction, 1))
table_s2_wide <- data.frame(
  "Exchange during dialysis and recovery" = unique(table_s2$scenario),
  "2 h: rebound / reduction, %" = cells[table_s2$hours == 2],
  "6 h: rebound / reduction, %" = cells[table_s2$hours == 6],
  "24 h: rebound / reduction, %" = cells[table_s2$hours == 24],
  check.names = FALSE
)

table_s3_text <- with(table_s3, data.frame(
  "Exchange retained during dialysis, %" = fmt(100 * retained),
  "Exchange after dialysis" =
    ifelse(persists, "Reduction persists", "Immediately restored"),
  "Rebound at 2 h, %" = fmt(rebound, 1),
  "Reduction from baseline at 2 h, %" = fmt(reduction, 1),
  check.names = FALSE
))

writeLines(c(
  "# Results",
  "",
  sprintf("Generated by run_all.R with R %s (%s) and nlme %s.",
          getRversion(), R.version$nickname,
          packageDescription("nlme")$Version),
  "",
  participants, outcomes, concentration_models, validation, sensitivity,
  discussion,
  "## Table 1", "", md_table(table_1),
  supplement,
  "### Variance components", "", md_table(variance_table),
  "### Table S1", "", md_table(table_s1),
  "### Table S2", "", md_table(table_s2_wide),
  "### Table S3", "", md_table(table_s3_text)
), "output/results.md", useBytes = TRUE)
