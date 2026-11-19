# Rebound after dialysis and renal function before dialysis
#
# Model 1 (Lanoiselée et al.) uses two compartments and creatinine clearance
# (mL/min). Model 2 (Gurunathan et al.) uses three compartments. In model 2 both
# intercompartmental clearances scale with weight to the power of 0.75, and
# all three volumes scale linearly with weight. The renal covariate is
# 2021 CKD-EPI creatinine eGFR indexed to body surface area (mL/min/1.73 m2)
# (Dr. Usha Gurunathan, personal communication, November 18, 2026).
# Both models use volumes (v) in liters and systemic and intercompartmental
# clearances (cl and q) in L/h. We use their published parameter estimates.
model_1 <- function(weight, crcl = 0) {
  list(cl = 5.32 * (crcl / 85)^0.548,
       v = c(6.09 * (weight / 70), 9.55),
       q = 23.3)
}

model_2 <- function(weight, egfr_ml_min_1_73m2 = 0) {
  weight_ratio <- weight / 70
  list(cl = 5.0 * weight_ratio^0.75 * (egfr_ml_min_1_73m2 / 89.5)^1.65,
       v = c(6.8, 9.5, 56.9) * weight_ratio,
       q = c(17.6, 1.0) * weight_ratio^0.75)
}

# Each column describes where drug goes from one compartment. Column 1
# represents plasma, which exchanges drug with each peripheral compartment.
# Only plasma has systemic elimination.
rate_matrix <- function(cl, v, q) {
  n <- length(v)
  rates <- matrix(0, n, n)
  rates[1, 1] <- -(cl + sum(q)) / v[1]
  for (j in 2:n) {
    rates[j, 1] <- q[j - 1] / v[1]
    rates[1, j] <- q[j - 1] / v[j]
    rates[j, j] <- -q[j - 1] / v[j]
  }
  rates
}

# Solve the constant-rate compartment equations over t hours. Each term in
# the solution follows an exponential with one eigenvalue of the rate matrix.
propagate <- function(rates, initial_amounts, t) {
  modes <- eigen(rates)
  initial_weights <- solve(modes$vectors, initial_amounts)
  weights_at_t <- exp(modes$values * t) * initial_weights
  as.vector(Re(modes$vectors %*% weights_at_t))
}

# Central concentration (mg/L) t hours after an intravenous bolus.
bolus_concentration <- function(model, dose, t) {
  a0 <- c(dose, rep(0, length(model$v) - 1))
  propagate(rate_matrix(model$cl, model$v, model$q), a0, t)[1] / model$v[1]
}

# No samples were taken after dialysis. For these projections, impose the
# same fitted relative plasma decline throughout the 3.5 h session. Tissue
# exchanges drug with plasma at q_during. After dialysis, allow redistribution
# at q_after with zero systemic clearance and no further dose.
after_dialysis <- function(c0, v, q_during, q_after, k, p, hours_after,
                           duration = 3.5) {
  central <- function(t) c0[1] * exp(-k * t^p)
  c_end <- central(duration)
  # End-of-session tissue concentration has two parts. Some drug was already
  # in tissue at the start, and some entered from plasma during dialysis.
  for (j in seq_along(q_during)) {
    rate <- q_during[j] / v[j + 1]
    remaining <- c0[j + 1] * exp(-rate * duration)
    from_plasma <- function(t) rate * exp(-rate * (duration - t)) * central(t)
    inflow <- integrate(from_plasma, 0, duration, rel.tol = 1e-10)$value
    c_end[j + 1] <- remaining + inflow
  }
  rates <- rate_matrix(0, v, q_after)
  later <- sapply(hours_after, function(u) propagate(rates, c_end * v, u)[1])
  c(end = c_end[1], later / v[1])
}

# Two-compartment rebound at 2 h (post hoc, all 10 patients)
# Start plasma and tissue at the same concentration, scaled to 1. Use each
# patient's weight for the central volume and the fixed power-time curve for
# the fall in plasma. Repeat with each bootstrap curve for the interval.
rebound_by_patient <- function(k, p, retained = 1, persists = FALSE) {
  t(sapply(by_patient$weight, function(weight) {
    m <- model_1(weight)
    q_after <- if (persists) m$q * retained else m$q
    after_dialysis(c(1, 1), m$v, m$q * retained, q_after, k, p, 2)
  }))
}

rebound <- rebound_by_patient(k, p)
rebound_pct <- 100 * (rebound[, 2] / rebound[, 1] - 1)
rr_after_rebound <- 1 - rebound[, 2]
c_end_dialysis <- by_patient$c_first * rebound[, 1]
c_after_rebound <- by_patient$c_first * rebound[, 2]

rr_after_rebound_ci <- percentile_ci(mapply(function(k, p) {
  mean(1 - rebound_by_patient(k, p)[, 2])
}, power_boot$k, power_boot$p))

# Predialysis concentrations after intravenous dosing (figure S2)
# Other routes have unknown systemic input. Use the eight patients who
# received only intravenous tranexamic acid to define one representative
# patient. The clearance and plasma-loss analyses above use all 10.
iv_only <- filter(by_patient, dose_oral == 0, dose_topical == 0)
representative <- list(
  dose = mean(iv_only$dose_iv),
  weight = mean(iv_only$weight),
  hours = mean(iv_only$dose_to_first_sample_min) / 60,
  c_first = mean(iv_only$c_first)
)

pk_models <- list(model_1 = model_1, model_2 = model_2)

# Serum creatinine was unavailable. Find the model input that reproduces the
# group's mean first concentration. Zero elimination gives the highest
# prediction. If even that falls short, use zero for the curve and leave the
# implied kidney function missing.
renal_fit <- bind_rows(lapply(pk_models, function(model) {
  predict_first <- function(input) {
    bolus_concentration(model(representative$weight, input),
                        representative$dose, representative$hours)
  }
  target <- representative$c_first
  at_zero <- predict_first(0)
  boundary <- at_zero < target
  stopifnot(predict_first(200) <= target)
  input <- 0
  if (!boundary) {
    difference <- function(x) predict_first(x) - target
    input <- uniroot(difference, c(0, 200), tol = 1e-10)$root
  }
  data.frame(input_used = input, boundary_at_zero = boundary,
             implied_input = if (boundary) NA_real_ else input,
             predicted_first = predict_first(input))
}), .id = "model")

# The same comparison numbers refer to CrCl for model 1 and indexed eGFR
# for model 2. The low-elimination curve uses each model's selected input.
pk_curves <- expand.grid(
  hours = seq(0, representative$hours, length.out = 301),
  renal_scenario = c("low", "20", "40", "90"),
  model = names(pk_models), stringsAsFactors = FALSE
) |>
  left_join(select(renal_fit, model, input_used), by = "model") |>
  mutate(renal_input = coalesce(as.numeric(na_if(renal_scenario, "low")),
                                input_used))
pk_curves$concentration <- mapply(function(hours, renal_input, model) {
  m <- pk_models[[model]](representative$weight, renal_input)
  bolus_concentration(m, representative$dose, hours)
}, pk_curves$hours, pk_curves$renal_input, pk_curves$model)

# Slower exchange during dialysis (tables S2 and S3)
# Model 2 starts with the amounts predicted after the intravenous dose. At
# zero systemic clearance its first plasma concentration falls slightly below
# the observed mean. Use that predicted starting point in every scenario.
oracle_fit <- filter(renal_fit, model == "model_2")
m2 <- model_2(representative$weight, oracle_fit$input_used)
c_start <- propagate(rate_matrix(m2$cl, m2$v, m2$q),
                     c(representative$dose, 0, 0), representative$hours) / m2$v

stopifnot(abs(c_start[1] - oracle_fit$predicted_first) < 1e-8)

# Fractions of usual fast and slow exchange during dialysis. The 77%
# reduction comes from the study of N-acetylprocainamide by Stec et al.
# Restore usual exchange at dialysis end except in the final row.
scenarios_2 <- tribble(
  ~scenario,                                  ~fast, ~slow, ~persists,
  "Exchange unchanged",                           1,     1, FALSE,
  "Slow exchange halved",                         1,   0.5, FALSE,
  "Both exchange rates halved",                 0.5,   0.5, FALSE,
  "Slow exchange reduced by 77%",                 1,  0.23, FALSE,
  "Both exchange rates reduced by 77%",        0.23,  0.23, FALSE,
  "Slow exchange halved, reduction persists",     1,   0.5, TRUE
)
# Rebound is relative to plasma concentration at the end of dialysis.
# Remaining reduction is relative to the concentration just before dialysis.
table_s2 <- bind_rows(lapply(seq_len(nrow(scenarios_2)), function(i) {
  s <- scenarios_2[i, ]
  q_during <- m2$q * c(s$fast, s$slow)
  q_after <- if (s$persists) q_during else m2$q
  conc <- after_dialysis(c_start, m2$v, q_during, q_after, k, p, c(2, 6, 24))
  data.frame(scenario = s$scenario, hours = c(2, 6, 24),
             rebound = 100 * (conc[-1] / conc[1] - 1),
             reduction = 100 * (1 - conc[-1] / c_start[1]))
}))

# Table S3 applies the same exchange fractions to the two-compartment
# projection with initial equilibrium in all 10 patients.
table_s3 <- expand.grid(persists = c(FALSE, TRUE), retained = c(1, 0.5, 0.23))
table_s3[c("rebound", "reduction")] <- t(mapply(function(retained, persists) {
  conc <- rebound_by_patient(k, p, retained, persists)
  c(mean(100 * (conc[, 2] / conc[, 1] - 1)), mean(100 * (1 - conc[, 2])))
}, table_s3$retained, table_s3$persists))
