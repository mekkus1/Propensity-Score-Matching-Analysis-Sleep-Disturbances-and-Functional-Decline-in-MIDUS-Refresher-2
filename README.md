# Propensity Score–Matched Analysis of Sleep Disturbances and Functional Limitation in MIDUS Refresher 2

This repository contains R code for a propensity score–matched analysis examining the association between **sleep disturbances** and **functional limitation** in midlife using data from the **Midlife in the United States (MIDUS) Refresher 2** study.

The analysis reproduces the results for the manuscript:

> **Dose‑Response Association Between Sleep Disturbances and Functional Limitation in Midlife: A Propensity Score–Matched Study**  
> Philips N. Okeagu, Godisgreat O. Okeke, Yabets T. Kebede

---

## Overview

The main script (`psm_sleep_function_midus_refresher2.R`) implements the following steps:

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
7. **Estimate matched odds ratios (ORs)** for functional limitation using logistic regression with cluster‑robust standard errors.
8. **Conduct dose–response analyses**:
   - Ordinal sleep severity score (0 = Good, 1 = Intermediate, 2 = Poor),
   - Cumulative sleep problem count (0–3: trouble falling asleep, night waking, feeling unrested).
9. **Run subgroup analyses** (age <50 vs ≥50; sex) and **sensitivity analyses** (1:1 matching, stricter caliper).
10. **Export tables and figures** for manuscript use.

---

## Data

This project uses **restricted** MIDUS Refresher 2 survey data (SPSS format), which are **not included** in this repository.

Update the path in the script to point to your local `.sav` file, for example:

```r
midus_raw <- read_sav("path/to/MR2_P1_SURVEY_N2154_20251003.sav")

## Citation

If you use or adapt this code, please cite:

Okeagu PN, Okeke GO, Kebede YT. *Dose‑Response Association Between Sleep Disturbances and Functional Limitation in Midlife: A Propensity Score–Matched Study.* [Manuscript].

and:

Williams DR, Lachman ME, Krueger RF, et al. *Midlife in the United States (MIDUS Refresher 2), 2022–2024.* Inter‑university Consortium for Political and Social Research (ICPSR).

## Contact

Philips N. Okeagu  
Department of Oral Health Policy and Epidemiology  
Harvard School of Dental Medicine  
188 Longwood Avenue, Boston, MA 02115, USA  
Email: philipsokeagu@hsph.harvard.edu
