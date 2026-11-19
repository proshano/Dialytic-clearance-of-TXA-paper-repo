# Clearance of tranexamic acid by hemodialysis
#
# Place patients.csv and samples.csv in data/ and run Rscript run_all.R
# from this folder. The scripts share one R session. Open output/results.md
# for the numerical results, including the variance components and tables.
# The five figures are in the same folder.

library(dplyr)
library(ggplot2)
library(nlme)
library(patchwork)

source("R/01_data.R")
dir.create("output", showWarnings = FALSE)

source("R/02_clearance.R")
source("R/03_concentration_models.R")
source("R/04_bootstrap.R")
source("R/05_plasma_loss.R")
source("R/06_pharmacokinetic_models.R")
source("R/07_figures.R")
source("R/08_results.R")

writeLines(capture.output(sessionInfo()), "output/session_info.txt")
message("Finished. See output/results.md and the figures in output/.")
