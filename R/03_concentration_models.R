# Concentration during dialysis (figure 2 and supplemental figure S1)
#
# We fit all 66 inlet concentrations on the natural log scale. Time is in
# hours from each patient's first sample. Each curve gives the expected log
# concentration for a patient with that baseline and decline parameter.
#
#   constant rate   log C = b0 - k t
#   power time      log C = b0 - k t^p
#   asymptotic      log C = b0 - d (1 - exp(-r t))
#
# In the asymptotic model d is the eventual drop in log concentration.
# Between patients, b0 (log baseline) and log_k (or log_d) each have a fixed
# mean and a patient-specific deviation. The two deviations follow a bivariate
# normal distribution with mean zero. Each fit estimates both SDs and their
# correlation. The power p and asymptotic rate r are shared across patients.
#
# Within a patient, each observed log concentration equals the curve plus an
# independent normal error with mean zero and a common SD (sigma). These
# errors are independent of the patient effects. Sigma includes variation
# from sampling, the assay and other departures from the fitted curve.
# In the equations in Supplemental Digital Content 1, u0i and u1i are the
# patient deviations and epsilon_ij is the error for draw j in patient i.

models <- c("constant", "power", "asymptotic")

model_data <- samples |>
  mutate(patient = factor(patient), log_c = log(c_pre))

fit_model <- function(model, data, start) {
  terms <- switch(model,
    constant = list(
      curve = log_c ~ b0 - exp(log_k) * time,
      fixed = b0 + log_k ~ 1,
      random = b0 + log_k ~ 1 | patient
    ),
    power = list(
      curve = log_c ~ b0 - exp(log_k) * time^exp(log_p),
      fixed = b0 + log_k + log_p ~ 1,
      random = b0 + log_k ~ 1 | patient
    ),
    asymptotic = list(
      curve = log_c ~ b0 - exp(log_d) * (1 - exp(-exp(log_r) * time)),
      fixed = b0 + log_d + log_r ~ 1,
      random = b0 + log_d ~ 1 | patient
    )
  )

  nlme(
    model = terms$curve, data = data, fixed = terms$fixed,
    random = terms$random,   # Unstructured covariance for the patient effects
    correlation = NULL,      # Independent residuals given the patient effects
    weights = NULL,          # One residual SD on the log scale
    start = start, method = "ML",
    control = nlmeControl(maxIter = 200, pnlsMaxIter = 50, msMaxIter = 200)
  )
}

# For a typical patient, set the random effects to zero. Taking
# log(C(t) / C(0)) cancels the baseline and leaves just the decline.
log_ratio <- function(model, fixed, t) {
  switch(model,
    constant = -exp(fixed[["log_k"]]) * t,
    power = -exp(fixed[["log_k"]]) * t^exp(fixed[["log_p"]]),
    asymptotic = -exp(fixed[["log_d"]]) * (1 - exp(-exp(fixed[["log_r"]]) * t))
  )
}

# Start with the mean log concentration at the first sample. Use the baseline
# from the power-time fit to initialize the two other fits.
start_b0 <- mean(model_data$log_c[model_data$sample == 1])
fits <- list()
fits$power <- fit_model("power", model_data,
                        c(b0 = start_b0, log_k = log(0.54), log_p = log(0.71)))
b0_hat <- fixef(fits$power)[["b0"]]
fits$constant <- fit_model("constant", model_data,
                           c(b0 = b0_hat, log_k = log(0.37)))
fits$asymptotic <- fit_model("asymptotic", model_data,
                             c(b0 = b0_hat, log_d = log(2), log_r = log(0.5)))
fits <- fits[models]

# Select by Akaike information criterion (AIC). The figures and rebound
# calculations below use the power-time curve selected in the study.
aic <- sapply(fits, AIC)
if (names(which.min(aic)) != "power") {
  stop("A different model has the lowest AIC. The rebound calculations ",
       "and figure 2 require the power-time model selected in the study.")
}

# Solve exp(-k t^p) = 0.5 for the first time concentration halves.
k <- exp(fixef(fits$power)[["log_k"]])
p <- exp(fixef(fits$power)[["log_p"]])
halving_time <- (log(2) / k)^(1 / p)

# Reduction ratio after a typical 3.5 h session, under each candidate model.
rr_3.5h <- sapply(models, function(m) {
  1 - exp(log_ratio(m, fixef(fits[[m]]), 3.5))
})

# VarCorr gives the between-patient SDs (log baseline, then log decline)
# and their correlation. sigma gives the within-patient residual SD.
# The Variance components table in output/results.md reports all four.
variance_components <- t(sapply(fits, function(fit) {
  vc <- VarCorr(fit)
  c(sd_baseline = as.numeric(vc["b0", "StdDev"]),
    sd_decline = as.numeric(vc[2, "StdDev"]),
    correlation = as.numeric(vc[2, "Corr"]),
    sd_residual = sigma(fit))
}))
