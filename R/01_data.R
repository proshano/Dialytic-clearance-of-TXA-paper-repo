# Read the study data and calculate clearance at each paired sample.

input_files <- c("data/patients.csv", "data/samples.csv")
missing_files <- input_files[!file.exists(input_files)]
if (length(missing_files)) {
  stop("Missing ", paste(missing_files, collapse = " and "), ". ",
       "Study data are available from the corresponding author on ",
       "reasonable request. Run from the folder containing run_all.R.",
       call. = FALSE)
}

patients <- read.csv("data/patients.csv")
samples <- read.csv("data/samples.csv")

# Each row pairs inlet (c_pre) and outlet (c_post) concentrations at one
# sampling time. Less plasma leaves the dialyzer than enters it because of
# ultrafiltration. Convert its recorded rate from mL/h to mL/min for clearance.
samples <- samples |>
  left_join(select(patients, patient, hct), by = "patient") |>
  mutate(
    time = time_min / 60,                          # h since the first sample
    qp = qb * (1 - hct / 100),                     # plasma flow, mL/min
    er = (c_pre - c_post) / c_pre,                 # extraction ratio
    cl = qp * er + (quf / 60) * (c_post / c_pre)   # dialytic clearance, mL/min
  ) |>
  arrange(patient, sample)

# For each patient, average clearance and extraction across pairs. The
# reduction ratio compares the first and last inlet concentrations. Sample
# number gives the order of draws with the same recorded time.
by_patient <- samples |>
  group_by(patient) |>
  summarise(
    pairs = n(),
    duration = last(time),
    c_first = first(c_pre),
    c_last = last(c_pre),
    rr = (c_first - c_last) / c_first,
    er = mean(er),
    cl = mean(cl),
    qb = mean(qb),
    quf = mean(quf)
  ) |>
  left_join(patients, by = "patient")
