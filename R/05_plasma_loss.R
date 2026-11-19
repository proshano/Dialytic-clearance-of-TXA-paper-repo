# Loss of tranexamic acid from plasma (post hoc analysis, table S1)
#
# At each paired draw, clearance (mL/min) times inlet concentration (mg/L)
# divided by 1000 gives the loss rate in mg/min. This is the amount entering
# the dialyzer in plasma each minute minus the amount leaving in plasma.
# We estimate loss only from the first to the last draw in all 10 patients.

# For each interval, multiply its duration by the average of its two rates.
trapezoid <- function(time, rate) {
  mean_rate <- (head(rate, -1) + tail(rate, -1)) / 2
  sum(diff(time) * mean_rate)
}

samples <- samples |>
  mutate(loss_rate = cl * c_pre / 1000,
         loss_rate_hct = cl_hct * c_pre / 1000)

# If draws share a recorded time, use their mean loss rate, then repeat
# the calculation with the first and last draws at that time separately.
by_time <- samples |>
  group_by(patient, time) |>
  summarise(
    rate = mean(loss_rate),
    rate_first = first(loss_rate),
    rate_last = last(loss_rate),
    rate_hct = mean(loss_rate_hct),
    cl = mean(cl),
    .groups = "drop"
  )

# In the sensitivity analysis, coef() gives the fitted curve for each
# patient, including their random effects. Exponentiate that curve without
# a lognormal mean correction. Interpolate clearance between paired draws.
curve_parameters <- coef(fits$power)

curve_loss <- function(id) {
  draws <- filter(by_time, patient == id)
  fit <- curve_parameters[as.character(id), ]
  concentration <- function(t) exp(fit$b0 - exp(fit$log_k) * t^exp(fit$log_p))
  clearance <- approxfun(draws$time, draws$cl)
  rate <- function(t) clearance(t) * concentration(t) / 1000
  intervals <- sapply(seq_len(nrow(draws) - 1), function(j) {
    integrate(rate, draws$time[j], draws$time[j + 1],
              rel.tol = 1e-9, abs.tol = 1e-9)$value
  })
  60 * sum(intervals)
}

# The time column is in hours. Multiplying by 60 gives the amount in mg.
plasma_loss <- by_time |>
  group_by(patient) |>
  summarise(
    duration = max(time),
    primary = 60 * trapezoid(time, rate),
    hct = 60 * trapezoid(time, rate_hct),
    tied_first = 60 * trapezoid(time, rate_first),
    tied_last = 60 * trapezoid(time, rate_last)
  ) |>
  mutate(fitted_curve = sapply(patient, curve_loss))

loss_summary <- sapply(
  list(primary = plasma_loss$primary,
       hct = plasma_loss$hct,
       fitted_curve = plasma_loss$fitted_curve,
       hct_change = plasma_loss$hct - plasma_loss$primary,
       curve_change = plasma_loss$fitted_curve - plasma_loss$primary,
       tied_first = plasma_loss$tied_first,
       tied_last = plasma_loss$tied_last),
  t_summary
)
