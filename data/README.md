# Data

The corresponding author can provide the two data files on reasonable request. The times are elapsed minutes. All fields are complete for the recorded draws. The files omit age, sex, height and type of surgery. Their summaries in the paper come from the case report forms.

## patients.csv

One row per patient.

| Column | Description | Unit |
| --- | --- | --- |
| `patient` | Study number, as in Table 1 | |
| `weight` | Pre-dialysis weight | kg |
| `dialyzer` | FX 600, 800 or 1000 (Fresenius Medical Care) or Theranova 400 or 500 (Baxter) | |
| `hct` | Most recent pre-dialysis hematocrit, measured on the day of the session | % |
| `qd` | Prescribed dialysate flow for the session | mL/min |
| `dose_iv` | Intravenous tranexamic acid given during surgery | mg |
| `dose_oral` | Total oral tranexamic acid dose across administrations | mg |
| `dose_topical` | Topical tranexamic acid irrigation during surgery. The amount absorbed was not measured | mg |
| `dose_to_first_sample_min` | Time from intravenous dosing to the first pair. The models treat the whole intravenous dose as given at the recorded time | min |

## samples.csv

One row per paired sample, 66 in all.

| Column | Description | Unit |
| --- | --- | --- |
| `patient` | Study number | |
| `sample` | Scheduled draw. Samples 1 to 6 correspond to 0, 0.5, 1, 1.5, 2 and 3 h. Sample 7 marks the end of dialysis | |
| `time_min` | Recorded time since the first paired sample | min |
| `qb` | Blood flow at the time of the sample | mL/min |
| `quf` | Ultrafiltration rate at the time of the sample | mL/h |
| `c_pre` | Plasma tranexamic acid at the dialyzer inlet (pre-dialyzer) | mg/L |
| `c_post` | Plasma tranexamic acid at the dialyzer outlet (post-dialyzer), drawn at the same time | mg/L |

Sessions may end before all scheduled draws. When draws have the same recorded time, the sample number gives their order.
