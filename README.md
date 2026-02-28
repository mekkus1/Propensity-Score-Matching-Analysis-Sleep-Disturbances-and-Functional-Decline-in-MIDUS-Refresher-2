# Propensity Score Matching Analysis: Sleep Disturbances and Functional Decline in MIDUS Refresher 2

## About This Repository

This repository contains R code for a propensity score matched analysis examining the association between poor sleep and functional decline using data from the Midlife in the United States (MIDUS) Refresher 2 study.

## 📄 Related Publication

This code accompanies the manuscript:

> **Dose-Response Association Between Poor Sleep and Functional Limitation in Midlife: A Propensity Score–Matched Study**  
> Okeagu PN, Okeke GO, Odu EC, Apreku A, Kebede YT  
> *Under review*

## 🔍 Re-Analysis Context

**Important note:** This repository contains a **re-analysis** of the data. The exact original analysis pipeline that produced the published results is no longer available. This code represents our best effort to reconstruct the analytic approach based on the methods described in the manuscript.

| Analysis | Odds Ratio (95% CI) | Sample Size |
|----------|---------------------|--------------|
| **Published manuscript** | 1.81 (1.32–2.48) | N = 866 (matched) |
| **This re-analysis** | 0.99 (0.68–1.44) | N = 870 (matched) |


## 📊 What This Code Does

1. Loads a frozen dataset (`data/ml_data_psm_final.rds`)
2. Defines exposure (poor sleep: feeling unrested often/always vs. never/rarely)
3. Selects 15 confounders (demographics, socioeconomic factors, health behaviors, mental health)
4. Performs 1:3 nearest neighbor propensity score matching with caliper 0.2
5. Estimates the association between poor sleep and functional decline in the matched sample

## 🖥️ Requirements

```r
# Install required packages
install.packages(c("tidyverse", "MatchIt", "cobalt", "sandwich", 
                   "lmtest", "tableone"))
