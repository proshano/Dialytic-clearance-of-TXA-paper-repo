# Clearance of tranexamic acid by hemodialysis

Analysis code for the prospective observational study by Roshanov et al. We measured tranexamic acid at the dialyzer inlet and outlet during one postoperative hemodialysis session in each of 10 patients. The scripts start with the 66 paired samples and produce the clearance estimates, concentration models, tables and figures.

## Reading the analysis

For the mixed models, start with `R/03_concentration_models.R`. The `random` terms let patients differ in baseline concentration and in the rate or extent of decline. The residual standard deviation (SD) describes the scatter of measurements around each patient's curve on the log scale. The last block extracts both between-patient SDs, their correlation and the residual SD for each model. Supplemental Digital Content 1 reports the estimates.

## Running the analysis

Place `patients.csv` and `samples.csv` in the `data` folder beside `run_all.R`. The corresponding author can provide the study data on reasonable request. The column definitions are in `data/README.md`. The demographic summaries in the paper come from the case report forms.

We used R 4.6.1 with dplyr, ggplot2, nlme and patchwork. With those packages installed, open R in the folder containing `run_all.R` and run:

```r
source("run_all.R")
```

From a terminal in the same folder, use `Rscript run_all.R`. Open `output/results.md` for the numerical results and tables. The same folder contains all five figures in PDF and PNG and the session information.

To use the package versions from the study, start R in the analysis folder and run the following once. The lockfile also records ragg, which ggplot2 uses for PNG output when available.

```r
install.packages("renv")
renv::load(".")
renv::restore()
source("run_all.R")
```

Most of the run time goes to 3,000 bootstrap fits. The script uses two cores by default on macOS and Linux and one on Windows. For a single-core run on any platform, set `options(mc.cores = 1L)` before sourcing the script.

## Contents

| File | What to look for |
| --- | --- |
| `run_all.R` | Run the scripts in order in one R session |
| `R/01_data.R` | From paired inlet and outlet samples to a clearance for each patient |
| `R/02_clearance.R` | Mean clearance and its two-sided 95% CI across patients, extraction and reduction ratios, and sensitivity analyses |
| `R/03_concentration_models.R` | The three mixed models, their variance components and the time to a 50% fall in concentration |
| `R/04_bootstrap.R` | Resample patients, repeat model selection, assess calibration and calculate the bands in figure 2 |
| `R/05_plasma_loss.R` | Table S1. Integrate measured loss rates and repeat with each patient's fitted concentration curve |
| `R/06_pharmacokinetic_models.R` | Figure S2 and tables S2 and S3. Compare renal inputs and project rebound under different rates of tissue exchange |
| `R/07_figures.R` | Draw figures 1 to 3 and supplemental figures S1 and S2 |
| `R/08_results.R` | Write estimates, intervals and tables to `output/results.md` |

## Notes

Small differences in numerical convergence can occur across computers. The bootstrap keeps successful fits from each resample. The concentration-model results in `output/results.md` include the count of failed fits.
