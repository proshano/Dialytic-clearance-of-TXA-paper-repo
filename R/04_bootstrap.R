# Internal validation in 1000 patient-cluster bootstrap samples.
#
# Each bootstrap sample draws 10 patients with replacement and keeps all of
# their measurements together. Refit all three models, starting from their
# full-data estimates, and select the lowest AIC again in each resample.
# This repeats the choice of model when assessing predictive performance.

n_boot <- 1000
percentile_ci <- function(x) quantile(x, c(0.025, 0.975))

set.seed(20260805)
boot_ids <- replicate(n_boot, sample(by_patient$patient, replace = TRUE),
                      simplify = FALSE)

# Give each draw a new patient ID. Otherwise nlme would treat two draws of
# the same patient as a single patient with twice as many measurements.
resample <- function(ids) {
  pieces <- lapply(seq_along(ids), function(i) {
    model_data |> filter(patient == ids[i]) |> mutate(patient = i)
  })
  bind_rows(pieces) |> mutate(patient = factor(patient))
}

# For validation, compare each later log concentration with that patient's
# first. The baseline log ratio is zero by definition, so omit it here.
# A well-calibrated curve has intercept 0 and slope 1.
curve_performance <- function(data, model, fixed) {
  comparisons <- data |>
    group_by(patient) |>
    mutate(observed = log_c - first(log_c)) |>
    ungroup() |>
    filter(sample > 1) |>
    mutate(predicted = log_ratio(model, fixed, time))
  error <- comparisons$observed - comparisons$predicted
  calibration <- coef(lm(observed ~ predicted, data = comparisons))
  c(mean_error = mean(error),
    rmse = sqrt(mean(error^2)),
    intercept = calibration[[1]],
    slope = calibration[[2]])
}

# Mean first-to-final duration and reduction ratio over a set of patients.
first_to_final <- function(data) {
  ends <- data |>
    group_by(patient) |>
    summarise(duration = last(time),
              rr = (first(c_pre) - last(c_pre)) / first(c_pre))
  c(duration = mean(ends$duration), rr = mean(ends$rr))
}
original <- first_to_final(model_data)

fit_replicate <- function(b) {
  data <- resample(boot_ids[[b]])
  resampled <- first_to_final(data)
  rows <- lapply(models, function(model) {
    fit <- tryCatch(
      suppressWarnings(fit_model(model, data, fixef(fits[[model]]))),
      error = function(e) NULL
    )
    # A failed fit has no AIC. Keep the other fits from that resample.
    if (is.null(fit)) return(NULL)
    fixed <- fixef(fit)
    is_power <- model == "power"
    sd_log_k <- NA
    if (is_power) sd_log_k <- as.numeric(VarCorr(fit)["log_k", "StdDev"])
    # Performance in the bootstrap sample (apparent) and in the original data
    # (test). The difference is the optimism of the bootstrap fit.
    apparent <- curve_performance(data, model, fixed)
    test <- curve_performance(model_data, model, fixed)
    data.frame(
      replicate = b,
      model = model,
      aic = AIC(fit),
      rr_3.5h = 1 - exp(log_ratio(model, fixed, 3.5)),
      k = if (is_power) exp(fixed[["log_k"]]) else NA,
      p = if (is_power) exp(fixed[["log_p"]]) else NA,
      sd_log_k = sd_log_k,
      optimism_error = apparent[["mean_error"]] - test[["mean_error"]],
      optimism_rmse = apparent[["rmse"]] - test[["rmse"]],
      optimism_intercept = apparent[["intercept"]] - test[["intercept"]],
      optimism_slope = apparent[["slope"]] - test[["slope"]],
      apparent_slope = apparent[["slope"]],
      # Mean reduction ratio of the bootstrap sample, and the model's
      # predictions at the mean durations of the bootstrap and original data.
      duration = resampled[["duration"]],
      observed_rr = resampled[["rr"]],
      predicted_rr = 1 - exp(log_ratio(model, fixed, resampled[["duration"]])),
      predicted_rr_original =
        1 - exp(log_ratio(model, fixed, original[["duration"]]))
    )
  })
  bind_rows(rows)
}

# The patient draws above are fixed before work is divided between cores.
cores <- if (.Platform$OS.type == "windows") 1L else getOption("mc.cores", 2L)
message("Refitting three models in each of ", n_boot, " patient resamples.")
boot <- bind_rows(parallel::mclapply(seq_len(n_boot), fit_replicate,
                                     mc.cores = cores))
failed_fits <- n_boot * length(models) - nrow(boot)

selected <- boot |>
  group_by(replicate) |>
  slice_min(aic, n = 1, with_ties = FALSE) |>
  ungroup()
times_selected <- table(factor(selected$model, levels = models))

# For each metric, subtract the mean optimism of the selected bootstrap
# fits from the full-data performance of the selected power-time model.
optimism <- colMeans(select(selected,
  mean_error = optimism_error, rmse = optimism_rmse,
  intercept = optimism_intercept, slope = optimism_slope
))
full_data <- curve_performance(model_data, "power", fixef(fits$power))
corrected <- full_data - optimism[names(full_data)]
corrected_slope_ci <- percentile_ci(selected$apparent_slope -
                                      optimism[["slope"]])

# Calibration of the mean observed first-to-final reduction ratio (figure 3).
# Each bootstrap sample contributes its mean observed ratio and the selected
# model's prediction at that sample's mean duration. The calibration curve is
# a lowess smooth without robustness iterations.
original_predicted_rr <- 1 - exp(log_ratio("power", fixef(fits$power),
                                           original[["duration"]]))
calibration_optimism <- mean(
  (selected$observed_rr - selected$predicted_rr) -
    (original[["rr"]] - selected$predicted_rr_original)
)
smooth <- lowess(selected$predicted_rr, selected$observed_rr, f = 2 / 3,
                 iter = 0)
calibration_curve <- data.frame(
  predicted = seq(min(selected$predicted_rr), max(selected$predicted_rr),
                  length.out = 101)
)
calibration_curve$observed <- approx(smooth$x, smooth$y,
                                     xout = calibration_curve$predicted,
                                     rule = 2, ties = mean)$y
calibration_curve$corrected <- calibration_curve$observed - calibration_optimism

# Percentile intervals. Those for the reduction ratio at 3.5 h come from the
# bootstrap fits of each model, so they are conditional on its form.
power_boot <- filter(boot, model == "power")
halving_time_ci <- percentile_ci((log(2) / power_boot$k)^(1 / power_boot$p))
rr_3.5h_ci <- sapply(models, function(m) {
  percentile_ci(boot$rr_3.5h[boot$model == m])
})

# Pointwise bands around the typical-patient curve for figure 2, as percentage
# change from baseline. The 95% confidence band spans the bootstrap curves.
# For the 90% prediction band, each bootstrap curve is paired with one patient
# drawn from its estimated between-patient distribution of log k. The band
# describes underlying patient curves. The random baseline cancels in a
# within-patient ratio. We do not add residual error from individual samples.
set.seed(20260811)
patient_k <- power_boot$k * exp(power_boot$sd_log_k * rnorm(nrow(power_boot)))
curve_band <- bind_rows(lapply(seq(0, 4, length.out = 201), function(t) {
  typical <- 100 * (exp(-power_boot$k * t^power_boot$p) - 1)
  patient_change <- 100 * (exp(-patient_k * t^power_boot$p) - 1)
  data.frame(
    time = t,
    change = 100 * (exp(-k * t^p) - 1),
    ci_low = unname(quantile(typical, 0.025)),
    ci_high = unname(quantile(typical, 0.975)),
    pi_low = unname(quantile(patient_change, 0.05)),
    pi_high = unname(quantile(patient_change, 0.95))
  )
}))
