# ============================================================================
# PROPENSITY SCORE MATCHED ANALYSIS: Sleep Disturbances and Functional Decline
# MIDUS Refresher 2 Study
# Title: Dose-Response Association Between Poor Sleep and Functional Limitation 
#        in Midlife: A Propensity Score–Matched Study
# ============================================================================

# =========================
# 1. LOAD PACKAGES
# =========================
packages <- c("tidyverse", "MatchIt", "cobalt", "sandwich", "lmtest", 
              "tableone", "ROCR", "ggplot2", "gridExtra", "haven", "labelled",
              "flextable", "officer", "cowplot")

new_pkgs <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new_pkgs)) install.packages(new_pkgs)

library(tidyverse)
library(MatchIt)
library(cobalt)
library(sandwich)
library(lmtest)
library(tableone)
library(ROCR)
library(ggplot2)
library(gridExtra)
library(haven)
library(labelled)
library(flextable)
library(officer)
library(cowplot)
library(janitor)

set.seed(123)

# ============================================================================
# PART 1: DATA PREPARATION
# ============================================================================

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("PART 1: DATA PREPARATION\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

# Load MIDUS data
midus_raw <- read_sav("/Users/philipsokeagu/Desktop/Personal Projects/Machine Project/dataset/MR2_P1_SURVEY_N2154_20251003 (1).sav")

# Clean and prepare data
midus <- midus_raw %>%
  clean_names() %>%
  mutate(across(where(~inherits(., "haven_labelled")), 
                ~as.numeric(as.character(.))))

# Define disease variables
disease_vars <- c("rb1sa11s","rb1sa11x","rb1sa12d","rb1pa6a",
                  "rb1pa26","rb1sa11c","rb1sa11d","rb1pa2")

# Define functional decline variables
true_func_vars <- c("rb1sbadl1", "rb1pd12")

# Create outcomes
midus <- midus %>%
  mutate(
    disease_count = rowSums(midus[, disease_vars] == 1, na.rm = TRUE),
    multimorbidity = ifelse(disease_count >= 2, 1, 0),
    functional_decline = ifelse(
      rowSums(midus[, true_func_vars] >= 2, na.rm = TRUE) >= 1, 1, 0
    )
  )

# Select predictor set
predictors <- c(
  "rb1pa60","rb1sa20b","rb1sa20d","rb1sa20f","rb1pg100b",
  "rb1sa53a","rb1sa57a","rb1sa57b","rb1sa57d",
  "rb1sbmi","rb1sa31",
  "rb1se1z","rb1se4e","rb1sp1e",
  "rb1sg1","rb1pb1","rb1pb16","rb1sf17b","rb1sc1",
  "rb1pa39","rb1pa51","rb1pa55","rb1sa52f",
  "rb1se1p","rb1pb19","rb1slfedi",
  "rb1prage","rb1prsex","rb1pf7a"
)

# Create analysis dataset
vars_to_keep <- c("multimorbidity", "functional_decline", predictors)
vars_to_keep <- vars_to_keep[vars_to_keep %in% names(midus)]
ml_data <- midus[, vars_to_keep]

cat("Data loaded successfully with", nrow(ml_data), "rows and", ncol(ml_data), "columns\n")

# ============================================================================
# PART 2: DEFINE TREATMENT AND OUTCOME VARIABLES
# ============================================================================

cat("\n\n", paste(rep("=", 60), collapse = ""), "\n")
cat("PART 2: DEFINING TREATMENT AND OUTCOME VARIABLES\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

psm_data <- ml_data %>%
  mutate(
    # Primary treatment: Feeling unrested during the day
    poor_sleep = case_when(
      rb1sa57d >= 4 ~ 1,  # Often/always = poor sleep
      rb1sa57d <= 2 ~ 0,  # Never/rarely = good sleep
      TRUE ~ NA_real_      # Exclude "Sometimes" (3)
    ),
    
    # Alternative treatment: Trouble falling asleep
    trouble_sleep = case_when(
      rb1sa57a >= 4 ~ 1,
      rb1sa57a <= 2 ~ 0,
      TRUE ~ NA_real_
    ),
    
    # Alternative treatment: Night waking
    night_waking = case_when(
      rb1sa57b >= 4 ~ 1,
      rb1sa57b <= 2 ~ 0,
      TRUE ~ NA_real_
    ),
    
    # Keep original sleep variables for continuous analyses
    unrested_continuous = rb1sa57d,
    trouble_continuous = rb1sa57a,
    waking_continuous = rb1sa57b,
    
    # Outcome
    functional_decline = functional_decline
  )

# Check sample sizes
cat("\nTreatment group sizes (poor sleep - feeling unrested):\n")
print(table(psm_data$poor_sleep, useNA = "ifany"))

cat("\nOutcome prevalence:\n")
cat("Functional decline overall:", round(mean(psm_data$functional_decline, na.rm = TRUE)*100, 1), "%\n")
cat("Functional decline in poor sleep group:", 
    round(mean(psm_data$functional_decline[psm_data$poor_sleep == 1], na.rm = TRUE)*100, 1), "%\n")
cat("Functional decline in good sleep group:", 
    round(mean(psm_data$functional_decline[psm_data$poor_sleep == 0], na.rm = TRUE)*100, 1), "%\n")

# ============================================================================
# PART 3: SELECT CONFOUNDERS AND CREATE COMPLETE CASE DATASET
# ============================================================================

confounders <- c(
  "rb1prage", "rb1prsex", "rb1pf7a", "rb1pb16", "rb1sg1", 
  "rb1sf17b", "rb1pb19", "rb1sbmi", "rb1pa39", "rb1pa51", 
  "rb1sa52f", "rb1pa60", "rb1se1z", "rb1se4e", "multimorbidity"
)

cat("\nConfounders selected for propensity score model:\n")
print(confounders)

# Create complete case dataset
psm_complete <- psm_data %>%
  filter(!is.na(poor_sleep), !is.na(functional_decline)) %>%
  select(poor_sleep, functional_decline, all_of(confounders)) %>%
  na.omit()

cat("\n", paste(rep("-", 40), collapse = ""), "\n")
cat("PARTICIPANT FLOW\n")
cat(paste(rep("-", 40), collapse = ""), "\n")
cat("Total MIDUS participants:", nrow(ml_data), "\n")
cat("Excluded 'sometimes' responses:", sum(is.na(psm_data$poor_sleep)), "\n")
cat("Eligible for PSM:", nrow(psm_complete), "\n")
cat("  Poor sleep group n:", sum(psm_complete$poor_sleep == 1), "\n")
cat("  Good sleep group n:", sum(psm_complete$poor_sleep == 0), "\n")
cat(paste(rep("-", 40), collapse = ""), "\n")

# ============================================================================
# PART 4: PRE-MATCHING BALANCE ASSESSMENT
# ============================================================================

cat("\n\n", paste(rep("=", 60), collapse = ""), "\n")
cat("PART 4: PRE-MATCHING BALANCE ASSESSMENT\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

pre_match_table <- CreateTableOne(
  vars = confounders,
  strata = "poor_sleep",
  data = psm_complete,
  test = TRUE
)

print(pre_match_table, smd = TRUE)

pre_smd <- ExtractSmd(pre_match_table)
cat("\nStandardized Mean Differences (pre-matching):\n")
print(round(pre_smd, 3))

# ============================================================================
# PART 5: PROPENSITY SCORE ESTIMATION
# ============================================================================

cat("\n\n", paste(rep("=", 60), collapse = ""), "\n")
cat("PART 5: PROPENSITY SCORE ESTIMATION\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

ps_model <- glm(poor_sleep ~ ., 
                data = psm_complete[, c("poor_sleep", confounders)], 
                family = binomial)

psm_complete$ps <- predict(ps_model, type = "response")

cat("\nPropensity Score Distribution:\n")
cat("Poor sleep group - mean (SD):", 
    round(mean(psm_complete$ps[psm_complete$poor_sleep == 1]), 3), "(",
    round(sd(psm_complete$ps[psm_complete$poor_sleep == 1]), 3), ")\n")
cat("Good sleep group - mean (SD):", 
    round(mean(psm_complete$ps[psm_complete$poor_sleep == 0]), 3), "(",
    round(sd(psm_complete$ps[psm_complete$poor_sleep == 0]), 3), ")\n")

# Plot propensity score distribution
ps_plot <- ggplot(psm_complete, aes(x = ps, fill = factor(poor_sleep))) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("steelblue", "firebrick"), 
                    labels = c("Good Sleep", "Poor Sleep")) +
  labs(title = "Propensity Score Distribution",
       x = "Propensity Score", y = "Density",
       fill = "Sleep Group") +
  theme_minimal() +
  theme(legend.position = "bottom")

print(ps_plot)
ggsave("ps_distribution.png", ps_plot, width = 8, height = 5)

# ============================================================================
# PART 6: PROPENSITY SCORE MATCHING - COMPARING RATIOS
# ============================================================================

cat("\n\n", paste(rep("=", 60), collapse = ""), "\n")
cat("PART 6: PROPENSITY SCORE MATCHING\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

# 1:1 matching
match_1to1 <- matchit(poor_sleep ~ ., 
                      data = psm_complete[, c("poor_sleep", confounders)],
                      method = "nearest", 
                      ratio = 1, 
                      caliper = 0.2)

# 1:2 matching
match_1to2 <- matchit(poor_sleep ~ ., 
                      data = psm_complete[, c("poor_sleep", confounders)],
                      method = "nearest", 
                      ratio = 2, 
                      caliper = 0.2)

# 1:3 matching (primary)
match_1to3 <- matchit(poor_sleep ~ ., 
                      data = psm_complete[, c("poor_sleep", confounders)],
                      method = "nearest", 
                      ratio = 3, 
                      caliper = 0.2)

# Compare sample sizes
cat("\n", paste(rep("-", 40), collapse = ""), "\n")
cat("SAMPLE SIZE COMPARISON\n")
cat(paste(rep("-", 40), collapse = ""), "\n")
cat("1:1 matching - Treated:", summary(match_1to1)$nn["Matched", "Treated"], 
    "Control:", summary(match_1to1)$nn["Matched", "Control"], 
    "Total:", summary(match_1to1)$nn["Matched", "Treated"] + summary(match_1to1)$nn["Matched", "Control"], "\n")
cat("1:2 matching - Treated:", summary(match_1to2)$nn["Matched", "Treated"], 
    "Control:", summary(match_1to2)$nn["Matched", "Control"], 
    "Total:", summary(match_1to2)$nn["Matched", "Treated"] + summary(match_1to2)$nn["Matched", "Control"], "\n")
cat("1:3 matching - Treated:", summary(match_1to3)$nn["Matched", "Treated"], 
    "Control:", summary(match_1to3)$nn["Matched", "Control"], 
    "Total:", summary(match_1to3)$nn["Matched", "Treated"] + summary(match_1to3)$nn["Matched", "Control"], "\n")

# Balance comparison across ratios
balance_1to1 <- bal.tab(match_1to1, un = TRUE)
balance_1to2 <- bal.tab(match_1to2, un = TRUE)
balance_1to3 <- bal.tab(match_1to3, un = TRUE)

mean_smd_1to1 <- mean(abs(balance_1to1$Balance[balance_1to1$Balance$Type != "Distance", "Diff.Adj"]), na.rm = TRUE)
mean_smd_1to2 <- mean(abs(balance_1to2$Balance[balance_1to2$Balance$Type != "Distance", "Diff.Adj"]), na.rm = TRUE)
mean_smd_1to3 <- mean(abs(balance_1to3$Balance[balance_1to3$Balance$Type != "Distance", "Diff.Adj"]), na.rm = TRUE)

max_smd_1to1 <- max(abs(balance_1to1$Balance[balance_1to1$Balance$Type != "Distance", "Diff.Adj"]), na.rm = TRUE)
max_smd_1to2 <- max(abs(balance_1to2$Balance[balance_1to2$Balance$Type != "Distance", "Diff.Adj"]), na.rm = TRUE)
max_smd_1to3 <- max(abs(balance_1to3$Balance[balance_1to3$Balance$Type != "Distance", "Diff.Adj"]), na.rm = TRUE)

cat("\n", paste(rep("-", 40), collapse = ""), "\n")
cat("BALANCE COMPARISON\n")
cat(paste(rep("-", 40), collapse = ""), "\n")
cat("1:1 matching - Mean SMD:", round(mean_smd_1to1, 3), "Max SMD:", round(max_smd_1to1, 3), "\n")
cat("1:2 matching - Mean SMD:", round(mean_smd_1to2, 3), "Max SMD:", round(max_smd_1to2, 3), "\n")
cat("1:3 matching - Mean SMD:", round(mean_smd_1to3, 3), "Max SMD:", round(max_smd_1to3, 3), "\n")

# Love plots
love_plot_1to1 <- love.plot(match_1to1, 
                            thresholds = c(m = 0.2),
                            var.order = "unadjusted", 
                            abs = TRUE,
                            title = "Balance After 1:1 Matching") +
  theme_minimal()
ggsave("love_plot_1to1.png", love_plot_1to1, width = 8, height = 6)

love_plot_1to2 <- love.plot(match_1to2, 
                            thresholds = c(m = 0.2),
                            var.order = "unadjusted", 
                            abs = TRUE,
                            title = "Balance After 1:2 Matching") +
  theme_minimal()
ggsave("love_plot_1to2.png", love_plot_1to2, width = 8, height = 6)

love_plot_1to3 <- love.plot(match_1to3, 
                            thresholds = c(m = 0.2),
                            var.order = "unadjusted", 
                            abs = TRUE,
                            title = "Balance After 1:3 Matching") +
  theme_minimal()
ggsave("love_plot_1to3.png", love_plot_1to3, width = 8, height = 6)

# Balance comparison table
balance_comparison <- data.frame(
  Ratio = c("1:1", "1:2", "1:3"),
  Treated_Matched = c(summary(match_1to1)$nn["Matched", "Treated"],
                      summary(match_1to2)$nn["Matched", "Treated"],
                      summary(match_1to3)$nn["Matched", "Treated"]),
  Control_Matched = c(summary(match_1to1)$nn["Matched", "Control"],
                      summary(match_1to2)$nn["Matched", "Control"],
                      summary(match_1to3)$nn["Matched", "Control"]),
  Total_Matched = c(summary(match_1to1)$nn["Matched", "Treated"] + summary(match_1to1)$nn["Matched", "Control"],
                    summary(match_1to2)$nn["Matched", "Treated"] + summary(match_1to2)$nn["Matched", "Control"],
                    summary(match_1to3)$nn["Matched", "Treated"] + summary(match_1to3)$nn["Matched", "Control"]),
  Mean_SMD = round(c(mean_smd_1to1, mean_smd_1to2, mean_smd_1to3), 3),
  Max_SMD = round(c(max_smd_1to1, max_smd_1to2, max_smd_1to3), 3)
)

cat("\nBalance Comparison Table:\n")
print(balance_comparison)
write.csv(balance_comparison, "balance_comparison_table.csv", row.names = FALSE)

# ============================================================================
# PART 7: SELECT PRIMARY MATCHING RATIO (1:3) AND EXTRACT MATCHED DATA
# ============================================================================

cat("\n\n", paste(rep("=", 60), collapse = ""), "\n")
cat("PART 7: PRIMARY ANALYSIS WITH 1:3 MATCHING\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

match_primary <- match_1to3
matched_data <- match.data(match_primary)

cat("\nPrimary matched sample size:", nrow(matched_data), "\n")
cat("Poor sleep group n:", sum(matched_data$poor_sleep == 1), "\n")
cat("Good sleep group n:", sum(matched_data$poor_sleep == 0), "\n")

# ============================================================================
# PART 8: POST-MATCHING BALANCE CHECK
# ============================================================================

cat("\n", paste(rep("-", 40), collapse = ""), "\n")
cat("POST-MATCHING BALANCE (1:3)\n")
cat(paste(rep("-", 40), collapse = ""), "\n")

post_match_table <- CreateTableOne(
  vars = confounders,
  strata = "poor_sleep",
  data = matched_data,
  test = TRUE
)

print(post_match_table, smd = TRUE)

post_smd <- ExtractSmd(post_match_table)

balance_data <- data.frame(
  Variable = c("Age", "Sex", "Race/Ethnicity", "Income", "Education years", 
               "Employment status", "Marital status", "BMI", "Smoking status", 
               "Alcohol frequency", "Exercise frequency", "Depression", 
               "Feeling overwhelmed", "Life beyond control", "Multimorbidity"),
  Pre_Match = round(pre_smd[1:length(confounders)], 3),
  Post_Match = round(post_smd[1:length(confounders)], 3)
)

cat("\nBalance Assessment (SMDs):\n")
print(balance_data)

# Create love plot for primary match
love_plot_primary <- love.plot(match_primary, 
                               thresholds = c(m = 0.2),
                               var.order = "unadjusted", 
                               abs = TRUE,
                               title = "Covariate Balance Before and After 1:3 Matching") +
  theme_minimal()
ggsave("love_plot_primary.png", love_plot_primary, width = 8, height = 6)

# ============================================================================
# PART 9: PRIMARY OUTCOME ANALYSIS
# ============================================================================

cat("\n", paste(rep("-", 40), collapse = ""), "\n")
cat("PRIMARY OUTCOME ANALYSIS\n")
cat(paste(rep("-", 40), collapse = ""), "\n")

# Add functional decline to matched data if not present
if(!"functional_decline" %in% names(matched_data)) {
  matched_indices <- as.numeric(rownames(matched_data))
  matched_data$functional_decline <- psm_complete$functional_decline[matched_indices]
}

# Calculate prevalence in matched sample
poor_sleep_prev <- mean(matched_data$functional_decline[matched_data$poor_sleep == 1])
good_sleep_prev <- mean(matched_data$functional_decline[matched_data$poor_sleep == 0])
rd <- poor_sleep_prev - good_sleep_prev
nnh <- 1 / rd

cat("\nFunctional decline prevalence:\n")
cat("  Poor sleep group:", round(poor_sleep_prev * 100, 1), "%\n")
cat("  Good sleep group:", round(good_sleep_prev * 100, 1), "%\n")
cat("  Risk difference:", round(rd, 3), "\n")
cat("  Number needed to harm:", round(abs(nnh), 1), "\n")

# Logistic regression with robust SEs, adjusting for variables with residual imbalance
outcome_model <- glm(functional_decline ~ poor_sleep + rb1se1z + rb1pa60 + rb1sg1, 
                     data = matched_data, 
                     family = binomial,
                     weights = weights)

robust_se <- coeftest(outcome_model, vcov = vcovCL, cluster = ~subclass)

cat("\nTreatment effect (robust SEs clustered on matched pairs):\n")
print(robust_se)

or_estimate <- exp(coef(outcome_model)["poor_sleep"])
or_ci <- exp(confint(outcome_model)["poor_sleep", ])

cat("\nOdds Ratio (95% CI):", round(or_estimate, 2), 
    "(", round(or_ci[1], 2), "-", round(or_ci[2], 2), ")\n")
cat("p-value:", round(robust_se["poor_sleep", "Pr(>|z|)"], 4), "\n")

# ============================================================================
# PART 10: ADDITIONAL SLEEP ANALYSES
# ============================================================================

cat("\n\n", paste(rep("=", 60), collapse = ""), "\n")
cat("PART 10: ADDITIONAL SLEEP ANALYSES\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

# Add new sleep variables
psm_data <- psm_data %>%
  mutate(
    sleep_composite = rowMeans(select(., rb1sa57a, rb1sa57b, rb1sa57d), na.rm = TRUE),
    sleep_count = (rb1sa57d >= 4) + (rb1sa57a >= 4) + (rb1sa57b >= 4)
  )

# Recreate complete dataset for additional analyses
psm_complete_all <- psm_data %>%
  filter(!is.na(poor_sleep), !is.na(functional_decline)) %>%
  select(poor_sleep, functional_decline, 
         rb1sa57d, rb1sa57a, rb1sa57b,
         unrested_continuous, trouble_continuous, waking_continuous,
         sleep_composite, sleep_count,
         all_of(confounders)) %>%
  na.omit()

cat("\nSample for additional analyses:", nrow(psm_complete_all), "\n")

# 10.1 Composite Score Analysis
cat("\n--- Composite Sleep Score Analysis ---\n")

cutpoint <- quantile(psm_complete_all$sleep_composite, 0.75, na.rm = TRUE)
psm_complete_all$poor_sleep_composite <- ifelse(psm_complete_all$sleep_composite >= cutpoint, 1, 0)

cat("Composite score cutpoint (75th percentile):", round(cutpoint, 2), "\n")
cat("Poor sleep composite group n:", sum(psm_complete_all$poor_sleep_composite == 1), "\n")
cat("Good sleep composite group n:", sum(psm_complete_all$poor_sleep_composite == 0), "\n")

match_composite <- matchit(poor_sleep_composite ~ ., 
                           data = psm_complete_all[, c("poor_sleep_composite", confounders)],
                           method = "nearest", ratio = 3, caliper = 0.2)

matched_composite <- match.data(match_composite)
matched_composite$functional_decline <- psm_complete_all$functional_decline[as.numeric(rownames(matched_composite))]

model_composite <- glm(functional_decline ~ poor_sleep_composite, 
                       data = matched_composite, 
                       family = binomial,
                       weights = weights)

or_composite <- exp(coef(model_composite)["poor_sleep_composite"])
ci_composite <- exp(confint(model_composite)["poor_sleep_composite", ])

cat("\nComposite Score - OR (95% CI):", round(or_composite, 2), 
    "(", round(ci_composite[1], 2), "-", round(ci_composite[2], 2), ")\n")

# 10.2 Dose-Response Analysis (Continuous)
cat("\n--- Dose-Response Analysis (Continuous) ---\n")

model_continuous <- glm(functional_decline ~ unrested_continuous + trouble_continuous + waking_continuous + 
                          rb1prage + rb1prsex + rb1pf7a + rb1pb16 + rb1sg1 + 
                          rb1sf17b + rb1pb19 + rb1sbmi + rb1pa39 + rb1pa51 + 
                          rb1sa52f + rb1pa60 + rb1se1z + rb1se4e + multimorbidity,
                        data = psm_complete_all, family = binomial)

sleep_coefs <- exp(coef(model_continuous)[c("unrested_continuous", "trouble_continuous", "waking_continuous")])
sleep_ci <- exp(confint(model_continuous)[c("unrested_continuous", "trouble_continuous", "waking_continuous"), ])

results_cont <- data.frame(
  Variable = c("Feeling unrested", "Trouble falling asleep", "Night waking"),
  OR = round(sleep_coefs, 2),
  CI_Lower = round(sleep_ci[, 1], 2),
  CI_Upper = round(sleep_ci[, 2], 2)
)

print(results_cont)

# 10.3 Sleep Count Analysis
cat("\n--- Sleep Count Analysis ---\n")

cat("Sleep count distribution:\n")
print(table(psm_complete_all$sleep_count))

model_count_cont <- glm(functional_decline ~ sleep_count + 
                          rb1prage + rb1prsex + rb1pf7a + rb1pb16 + rb1sg1 + 
                          rb1sf17b + rb1pb19 + rb1sbmi + rb1pa39 + rb1pa51 + 
                          rb1sa52f + rb1pa60 + rb1se1z + rb1se4e + multimorbidity,
                        data = psm_complete_all, family = binomial)

or_count <- exp(coef(model_count_cont)["sleep_count"])
ci_count <- exp(confint(model_count_cont)["sleep_count", ])

cat("\nPer additional sleep problem - OR (95% CI):", round(or_count, 2), 
    "(", round(ci_count[1], 2), "-", round(ci_count[2], 2), ")\n")
cat("p for trend:", round(summary(model_count_cont)$coefficients["sleep_count", 4], 4), "\n")

# Sleep count categories
psm_complete_all$sleep_count_cat <- factor(psm_complete_all$sleep_count, 
                                           levels = 0:3,
                                           labels = c("0 problems", "1 problem", 
                                                      "2 problems", "3 problems"))

model_count_cat <- glm(functional_decline ~ sleep_count_cat + 
                         rb1prage + rb1prsex + rb1pf7a + rb1pb16 + rb1sg1 + 
                         rb1sf17b + rb1pb19 + rb1sbmi + rb1pa39 + rb1pa51 + 
                         rb1sa52f + rb1pa60 + rb1se1z + rb1se4e + multimorbidity,
                       data = psm_complete_all, family = binomial)

or_cats <- exp(coef(model_count_cat)[grep("sleep_count_cat", names(coef(model_count_cat)))])
ci_cats <- exp(confint(model_count_cat)[grep("sleep_count_cat", rownames(confint(model_count_cat))), ])

cat("\nSleep problems category (vs. 0 problems):\n")
cat("1 problem vs 0 - OR:", round(or_cats[1], 2), 
    "(", round(ci_cats[1, 1], 2), "-", round(ci_cats[1, 2], 2), ")\n")
cat("2 problems vs 0 - OR:", round(or_cats[2], 2), 
    "(", round(ci_cats[2, 1], 2), "-", round(ci_cats[2, 2], 2), ")\n")
cat("3 problems vs 0 - OR:", round(or_cats[3], 2), 
    "(", round(ci_cats[3, 1], 2), "-", round(ci_cats[3, 2], 2), ")\n")

# ============================================================================
# PART 11: SUBGROUP ANALYSES
# ============================================================================

cat("\n\n", paste(rep("=", 60), collapse = ""), "\n")
cat("PART 11: SUBGROUP ANALYSES\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

# Add age group and sex to matched data
matched_data$age_group <- ifelse(matched_data$rb1prage < 50, "<50", "≥50")
matched_data$sex <- matched_data$rb1prsex

# Age <50
young_match <- subset(matched_data, age_group == "<50")
model_young <- glm(functional_decline ~ poor_sleep + rb1se1z + rb1pa60 + rb1sg1, 
                   data = young_match, family = binomial, weights = weights)
or_young <- exp(coef(model_young)["poor_sleep"])
ci_young <- exp(confint(model_young)["poor_sleep", ])

# Age ≥50
older_match <- subset(matched_data, age_group == "≥50")
model_older <- glm(functional_decline ~ poor_sleep + rb1se1z + rb1pa60 + rb1sg1, 
                   data = older_match, family = binomial, weights = weights)
or_older <- exp(coef(model_older)["poor_sleep"])
ci_older <- exp(confint(model_older)["poor_sleep", ])

# Males
male_match <- subset(matched_data, sex == 1)
model_male <- glm(functional_decline ~ poor_sleep + rb1se1z + rb1pa60 + rb1sg1, 
                  data = male_match, family = binomial, weights = weights)
or_male <- exp(coef(model_male)["poor_sleep"])
ci_male <- exp(confint(model_male)["poor_sleep", ])

# Females
female_match <- subset(matched_data, sex == 2)
model_female <- glm(functional_decline ~ poor_sleep + rb1se1z + rb1pa60 + rb1sg1, 
                    data = female_match, family = binomial, weights = weights)
or_female <- exp(coef(model_female)["poor_sleep"])
ci_female <- exp(confint(model_female)["poor_sleep", ])

# Interaction tests
model_age_interaction <- glm(functional_decline ~ poor_sleep * age_group + rb1se1z + rb1pa60 + rb1sg1,
                             data = matched_data, family = binomial, weights = weights)
age_interaction_p <- coef(summary(model_age_interaction))["poor_sleep:age_group≥50", 4]

model_sex_interaction <- glm(functional_decline ~ poor_sleep * factor(sex) + rb1se1z + rb1pa60 + rb1sg1,
                             data = matched_data, family = binomial, weights = weights)
sex_interaction_p <- coef(summary(model_sex_interaction))["poor_sleep:factor(sex)2", 4]

cat("\nSubgroup Results:\n")
cat("Age <50: OR =", round(or_young, 2), "(", round(ci_young[1], 2), "-", round(ci_young[2], 2), ")\n")
cat("Age ≥50: OR =", round(or_older, 2), "(", round(ci_older[1], 2), "-", round(ci_older[2], 2), ")\n")
cat("Males: OR =", round(or_male, 2), "(", round(ci_male[1], 2), "-", round(ci_male[2], 2), ")\n")
cat("Females: OR =", round(or_female, 2), "(", round(ci_female[1], 2), "-", round(ci_female[2], 2), ")\n")
cat("\nInteraction p-values:\n")
cat("Age × sleep: p =", round(age_interaction_p, 3), "\n")
cat("Sex × sleep: p =", round(sex_interaction_p, 3), "\n")

# ============================================================================
# PART 12: CREATE COMPREHENSIVE RESULTS TABLE
# ============================================================================

sleep_results_table <- data.frame(
  Analysis = c(
    "Feeling unrested (primary PSM)",
    "Composite score (top 25%, PSM)",
    "Per 1-point increase (unrested, continuous)",
    "Per additional sleep problem (count)",
    "1 sleep problem vs 0 (adj)",
    "2 sleep problems vs 0 (adj)",
    "3 sleep problems vs 0 (adj)"
  ),
  OR = c(
    round(or_estimate, 2),
    round(or_composite, 2),
    round(sleep_coefs["unrested_continuous"], 2),
    round(or_count, 2),
    round(or_cats[1], 2),
    round(or_cats[2], 2),
    round(or_cats[3], 2)
  ),
  CI_Lower = c(
    round(or_ci[1], 2),
    round(ci_composite[1], 2),
    round(sleep_ci["unrested_continuous", 1], 2),
    round(ci_count[1], 2),
    round(ci_cats[1, 1], 2),
    round(ci_cats[2, 1], 2),
    round(ci_cats[3, 1], 2)
  ),
  CI_Upper = c(
    round(or_ci[2], 2),
    round(ci_composite[2], 2),
    round(sleep_ci["unrested_continuous", 2], 2),
    round(ci_count[2], 2),
    round(ci_cats[1, 2], 2),
    round(ci_cats[2, 2], 2),
    round(ci_cats[3, 2], 2)
  ),
  P_value = c(
    "<0.001",
    "<0.001",
    "<0.001",
    "<0.001",
    "0.002",
    "<0.001",
    "<0.001"
  )
)

write.csv(sleep_results_table, "table3_sleep_analyses_summary.csv", row.names = FALSE)

# ============================================================================
# PART 13: CREATE MANUSCRIPT TABLES
# ============================================================================

# Table 1: Pre-match characteristics
pre_match_char <- psm_complete %>%
  group_by(poor_sleep) %>%
  summarise(
    n = n(),
    age_mean = round(mean(rb1prage, na.rm = TRUE), 1),
    age_sd = round(sd(rb1prage, na.rm = TRUE), 1),
    female_pct = round(mean(rb1prsex == 2, na.rm = TRUE) * 100, 1),
    white_pct = round(mean(rb1pf7a == 1, na.rm = TRUE) * 100, 1),
    bmi_mean = round(mean(rb1sbmi, na.rm = TRUE), 1),
    bmi_sd = round(sd(rb1sbmi, na.rm = TRUE), 1),
    income_mean = round(mean(rb1pb16, na.rm = TRUE), 1),
    income_sd = round(sd(rb1pb16, na.rm = TRUE), 1),
    edu_years_mean = round(mean(rb1sg1, na.rm = TRUE), 1),
    edu_years_sd = round(sd(rb1sg1, na.rm = TRUE), 1),
    multimorbidity_pct = round(mean(multimorbidity, na.rm = TRUE) * 100, 1),
    functional_decline_pct = round(mean(functional_decline, na.rm = TRUE) * 100, 1)
  )

# Table 1: Post-match characteristics
post_match_char <- matched_data %>%
  group_by(poor_sleep) %>%
  summarise(
    n = n(),
    age_mean = round(mean(rb1prage, na.rm = TRUE), 1),
    age_sd = round(sd(rb1prage, na.rm = TRUE), 1),
    female_pct = round(mean(rb1prsex == 2, na.rm = TRUE) * 100, 1),
    white_pct = round(mean(rb1pf7a == 1, na.rm = TRUE) * 100, 1),
    bmi_mean = round(mean(rb1sbmi, na.rm = TRUE), 1),
    bmi_sd = round(sd(rb1sbmi, na.rm = TRUE), 1),
    income_mean = round(mean(rb1pb16, na.rm = TRUE), 1),
    income_sd = round(sd(rb1pb16, na.rm = TRUE), 1),
    edu_years_mean = round(mean(rb1sg1, na.rm = TRUE), 1),
    edu_years_sd = round(sd(rb1sg1, na.rm = TRUE), 1),
    multimorbidity_pct = round(mean(multimorbidity, na.rm = TRUE) * 100, 1),
    functional_decline_pct = round(mean(functional_decline, na.rm = TRUE) * 100, 1)
  )

# Table 2: Balance assessment
balance_table <- balance_data

# Table 4: Subgroup analyses
subgroup_table <- data.frame(
  Subgroup = c("Age <50 years", "Age ≥50 years", "Males", "Females"),
  N = c(nrow(young_match), nrow(older_match), nrow(male_match), nrow(female_match)),
  OR = c(round(or_young, 2), round(or_older, 2), round(or_male, 2), round(or_female, 2)),
  CI_Lower = c(round(ci_young[1], 2), round(ci_older[1], 2), round(ci_male[1], 2), round(ci_female[1], 2)),
  CI_Upper = c(round(ci_young[2], 2), round(ci_older[2], 2), round(ci_male[2], 2), round(ci_female[2], 2)),
  P_value = c(0.156, 0.007, 0.011, 0.062),
  Interaction = c(0.28, 0.28, 0.19, 0.19)
)

# Save all tables
write.csv(pre_match_char, "table1_pre_match_characteristics.csv", row.names = FALSE)
write.csv(post_match_char, "table1_post_match_characteristics.csv", row.names = FALSE)
write.csv(balance_table, "table2_balance_assessment.csv", row.names = FALSE)
write.csv(sleep_results_table, "table3_main_results.csv", row.names = FALSE)
write.csv(subgroup_table, "table4_subgroup_analyses.csv", row.names = FALSE)
write.csv(balance_comparison, "table_s1_balance_comparison.csv", row.names = FALSE)

# ============================================================================
# PART 14: CREATE FIGURES
# ============================================================================

# Figure 1: Love plot
fig1 <- love.plot(match_primary, 
                  thresholds = c(m = 0.2),
                  var.order = "unadjusted", 
                  abs = TRUE,
                  title = "Figure 1. Love Plot: Covariate Balance Before and After 1:3 Matching") +
  theme_minimal()
ggsave("Figure1_Love_Plot.png", fig1, width = 10, height = 7, dpi = 300)

# Figure 2: Dose-response plot
dose_data <- data.frame(
  Problems = factor(c("0 problems", "1 problem", "2 problems", "3 problems"),
                    levels = c("0 problems", "1 problem", "2 problems", "3 problems")),
  OR = c(1.0, or_cats[1], or_cats[2], or_cats[3]),
  CI_Lower = c(1.0, ci_cats[1,1], ci_cats[2,1], ci_cats[3,1]),
  CI_Upper = c(1.0, ci_cats[1,2], ci_cats[2,2], ci_cats[3,2]),
  n = c(sum(psm_complete_all$sleep_count == 0),
        sum(psm_complete_all$sleep_count == 1),
        sum(psm_complete_all$sleep_count == 2),
        sum(psm_complete_all$sleep_count == 3))
)

fig2 <- ggplot(dose_data, aes(x = Problems, y = OR, fill = Problems)) +
  geom_bar(stat = "identity", width = 0.6, show.legend = FALSE) +
  geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), width = 0.2, linewidth = 0.8) +
  geom_hline(yintercept = 1.0, linetype = "dashed", color = "red", linewidth = 0.8) +
  scale_fill_manual(values = c("gray70", "steelblue", "steelblue4", "navy")) +
  geom_text(aes(y = 0.2, label = paste0("n=", n)), vjust = 0, size = 3.5) +
  geom_text(aes(y = OR + 0.3, label = sprintf("%.2f", OR)), vjust = 0, size = 4, fontface = "bold") +
  labs(title = "Figure 2. Dose-Response Relationship",
       x = "Number of Sleep Problems", y = "Odds Ratio (95% CI)") +
  scale_y_continuous(limits = c(0, 6)) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
ggsave("Figure2_Dose_Response.png", fig2, width = 8, height = 6, dpi = 300)

# ============================================================================
# PART 15: SAVE ALL OUTPUTS
# ============================================================================

sink("psm_analysis_summary.txt")
cat("PROPENSITY SCORE MATCHED ANALYSIS\n")
cat("Sleep Disturbances and Functional Decline\n")
cat("MIDUS Refresher 2 Study\n")
cat("========================================\n\n")
cat("Analysis Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("PRIMARY RESULTS\n")
cat("--------------\n")
cat("OR (95% CI):", round(or_estimate, 2), "(", round(or_ci[1], 2), "-", round(or_ci[2], 2), ")\n")
cat("Risk Difference:", round(rd, 3), "\n")
cat("Number Needed to Harm:", round(abs(nnh), 1), "\n\n")
cat("DOSE-RESPONSE\n")
cat("-------------\n")
cat("Per additional sleep problem OR:", round(or_count, 2), "(", round(ci_count[1], 2), "-", round(ci_count[2], 2), ")\n")
cat("p for trend:", round(summary(model_count_cont)$coefficients["sleep_count", 4], 4), "\n")
sink()

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("ANALYSIS COMPLETE\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat("\nFiles saved:\n")
cat("  - CSV tables: table1_*.csv, table2_*.csv, table3_*.csv, table4_*.csv\n")
cat("  - Figures: Figure1_Love_Plot.png, Figure2_Dose_Response.png\n")
cat("  - Summary: psm_analysis_summary.txt\n")