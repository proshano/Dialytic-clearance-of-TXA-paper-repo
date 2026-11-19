# Dialytic clearance (primary), extraction ratio and reduction ratio.
# Some patients gave more samples than others. Using one mean clearance per
# patient gives each of the 10 patients equal weight in the cohort estimate.

# Mean with a two-sided 95% t confidence interval, SD and range.
t_summary <- function(x) {
  ci <- t.test(x, alternative = "two.sided", conf.level = 0.95)$conf.int
  c(mean = mean(x), lower = ci[1], upper = ci[2], sd = sd(x),
    min = min(x), max = max(x))
}

clearance <- t_summary(by_patient$cl)
extraction <- t_summary(by_patient$er)
reduction <- t_summary(by_patient$rr)

# The interval above uses variation across the 10 patients (9 degrees of
# freedom). The protocol also specified a one-sided test against zero.
primary_test <- t.test(by_patient$cl, mu = 0, alternative = "greater")

# Prespecified sensitivity analysis. Blood flow may still have been rising
# at the first sample, so recalculate each patient's mean without that pair.
later_pairs <- samples |>
  filter(sample > 1) |>
  group_by(patient) |>
  summarise(cl = mean(cl), er = mean(er))
clearance_later <- t_summary(later_pairs$cl)
extraction_later <- t_summary(later_pairs$er)

# Post hoc sensitivity analysis. Hematocrit was measured before dialysis.
# Assume it rises linearly to 110% of that value at the final sample.
samples <- samples |>
  group_by(patient) |>
  mutate(
    hct_rising = (hct / 100) * (1 + 0.10 * (time / last(time))),
    cl_hct = qb * (1 - hct_rising) * er + (quf / 60) * (c_post / c_pre)
  ) |>
  ungroup()

clearance_hct <- samples |>
  group_by(patient) |>
  summarise(cl = mean(cl_hct)) |>
  pull(cl) |>
  t_summary()
