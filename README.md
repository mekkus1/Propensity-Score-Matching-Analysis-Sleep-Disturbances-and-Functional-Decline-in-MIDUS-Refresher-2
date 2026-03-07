# Propensity Score–Matched Analysis of Sleep Disturbances and Functional Limitation in MIDUS Refresher 2

This repository contains R code for a propensity score–matched analysis examining the association between **sleep disturbances** and **functional limitation** in midlife using data from the **Midlife in the United States (MIDUS) Refresher 2** study.

The analysis reproduces the results for the manuscript:

> **Dose‑Response Association Between Sleep Disturbances and Functional Limitation in Midlife: A Propensity Score–Matched Study**  
> Philips N. Okeagu, Godisgreat O. Okeke, Yabets T. Kebede

---

## Overview

The main script implements the following steps:

1. **Load and clean MIDUS Refresher 2 data** (SPSS `.sav` file).
2. **Create outcome variables**:
   - Primary outcome: functional limitation based on validated ADL/IADL items.
   - Multimorbidity and secondary outcomes for sensitivity analyses.
3. **Define a 3‑level sleep exposure** from *feeling unrested during the day* (`rb1sa57d`):
   - **Good**: never or rarely unrested,
   - **Intermediate**: sometimes unrested,
   - **Poor**: often or almost always unrested.
4. **Estimate multinomial propensity scores** for the 3‑level sleep exposure using 15 covariates:
   - Demographic, socioeconomic, behavioral, mental health, and physical health variables.
5. **Perform pairwise 1:3 nearest‑neighbor matching** with a caliper of 0.2 (logit scale) for:
   - Poor vs Good,
   - Poor vs Intermediate,
   - Intermediate vs Good.
6. **Assess covariate balance** using standardized mean differences and Love plots.
7. **Estimate matched odds ratios** (ORs) for functional limitation using logistic regression with cluster‑robust standard errors.
8. **Conduct dose–response analyses**:
   - Ordinal sleep severity score (0 = Good, 1 = Intermediate, 2 = Poor),
   - Cumulative sleep problem count (0–3: trouble falling asleep, night waking, feeling unrested).
9. **Run subgroup analyses** (age <50 vs ≥50; sex) and **sensitivity analyses** (1:1 matching, stricter caliper).
10. **Export tables and figures** for manuscript use.

---

## Data

This project uses restricted MIDUS Refresher 2 survey data (SPSS format), which are **not included** in this repository.

Update the path in the script to point to your local `.sav` file, for example:

```r
midus_raw <- read_sav("path/to/MR2_P1_SURVEY_N2154_20251003.sav")
psm_sleep_function_midus_refresher2.R
Key sections:

Load required libraries
Load and clean MIDUS data (haven, janitor)
Create functional limitation outcomes (ADL/IADL‑based)
Create 3‑level sleep exposure and original dichotomous sleep variable
Recode smoking and other covariates
Assemble the analysis dataset
Check and impute missing values (median/mode)
Estimate multinomial propensity scores (nnet::multinom)
Perform pairwise 1:3 matching (MatchIt::matchit)
Assess balance (cobalt::bal.tab, Love plot)
Estimate matched ORs with cluster‑robust SEs (sandwich, lmtest)
Dose–response models (sleep_level, sleep_score)
Subgroup analyses (age group, sex)
Sensitivity analyses (1:1 matching, caliper 0.1)
Visualizations (forest plot, dose–response plot)
Export results (CSVs, RDS, summary text)
Requirements
R version
R ≥ 4.2 (tested with R 4.4.2)
R packages
install.packages(c(
  "tidyverse",
  "nnet",
  "MatchIt",
  "cobalt",
  "sandwich",
  "lmtest",
  "broom",
  "ggplot2",
  "haven",
  "janitor"
))
library(tidyverse)
library(nnet)
library(MatchIt)
library(cobalt)
library(sandwich)
library(lmtest)
library(broom)
library(ggplot2)
library(haven)
library(janitor)
How to Run
Clone or download this repository.
Open psm_sleep_function_midus_refresher2.R in R or RStudio.
Edit the read_sav() path to point to your local MIDUS .sav file.
(Optional) Set a different random seed via set.seed().
Source or run the script from top to bottom.
On completion you should see:

“=== ANALYSIS COMPLETE ===”
A summary of created files and output directory.
Outputs
The script creates a directory:
results_matched_3level/
and saves:

CSV tables
primary_results.csv – ORs and 95% CIs for:

Poor vs Good,
Poor vs Intermediate,
Intermediate vs Good.
dose_response.csv – ORs and 95% CIs for intermediate and poor sleep vs good sleep.

linear_trend.csv – OR per one‑level increase in sleep severity.

subgroup_results.csv – ORs by age group (<50, ≥50) and sex (if available).

sensitivity_results.csv – ORs for:

Primary 1:3 matching (caliper 0.2),
1:1 matching,
1:3 matching with caliper 0.1.
Figures (PNG)
love_plot.png – Covariate balance (Poor vs Good).
forest_plot.png – Forest plot of primary matched comparisons.
dose_response.png – Dose–response plot by sleep category.
R objects and summary
analysis_results.rds – List containing datasets and model results:

data, matched_pg, matched_pi, matched_ig,
combined, primary, dose, trend, subgroup, sensitivity.
manuscript_summary.txt – Human‑readable summary of:

Sample characteristics,
Primary matched results,
Dose–response models,
Subgroup and sensitivity analyses.
Reproducing Manuscript Results
With the official MIDUS Refresher 2 dataset and the script as provided, you should reproduce the key results reported in the manuscript, including:

Poor vs Good sleep (primary PSM analysis):

OR ≈ 2.26 (95% CI ≈ 1.61–3.18) for functional limitation.

Dose–response models:

OR ≈ 1.45 per one‑level increase in sleep severity (Good → Intermediate → Poor).
OR ≈ 1.43 per additional sleep problem (0–3).
Subgroup and sensitivity estimates should match manuscript tables within rounding error.

Notes and Limitations
Cross‑sectional design; associations should not be interpreted as causal.
Sleep and functional outcomes are self‑reported MIDUS survey measures.
Missing data are handled via simple median/mode imputation.
Matching and PS estimation adjust for measured covariates only; unmeasured confounding may remain.


## Citation
If you use or adapt this code, please cite:

Okeagu PN, Okeke GO, Kebede YT. Dose‑Response Association Between Sleep Disturbances and Functional Limitation in Midlife: A Propensity Score–Matched Study. [Manuscript].

and:

Williams DR, Lachman ME, Krueger RF, et al. Midlife in the United States (MIDUS Refresher 2), 2022–2024. Inter‑university Consortium for Political and Social Research (ICPSR).

Contact
Philips N. Okeagu

Department of Oral Health Policy and Epidemiology

Harvard School of Dental Medicine

188 Longwood Avenue, Boston, MA 02115, USA

Email: philipsokeagu@hsph.harvard.edu

188 Longwood Avenue, Boston, MA 02115, USA

Email: philipsokeagu@hsph.harvard.edu
---

