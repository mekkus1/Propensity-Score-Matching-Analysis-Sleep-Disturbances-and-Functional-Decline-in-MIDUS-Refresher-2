################################################################################
# PROPENSITY SCORE MATCHED ANALYSIS WITH 3-LEVEL SLEEP EXPOSURE
# Using your actual MIDUS Refresher 2 data 
################################################################################

# ==============================================================================
# 1. LOAD REQUIRED LIBRARIES
# ==============================================================================

library(tidyverse)
library(nnet)
library(MatchIt)
library(cobalt)
library(sandwich)
library(lmtest)
library(broom)
library(ggplot2)
library(haven)      # For reading SPSS files
library(janitor)    # For cleaning names

# Set seed for reproducibility
set.seed(2024)

# ==============================================================================
# 2. LOADING DATA
# ==============================================================================

cat("\n=== Loading MIDUS Refresher 2 Data ===\n")

midus_raw <- read_sav("")

# Convert to dataframe and clean names
midus <- midus_raw %>%
  clean_names() %>%
  # Convert haven_labelled to numeric
  mutate(across(where(~inherits(., "haven_labelled")), 
                ~as.numeric(as.character(.))))

cat("Data loaded successfully\n")
cat("Sample size:", nrow(midus), "participants\n")
cat("Number of variables:", ncol(midus), "\n")

# ==============================================================================
# 3. RECREATE THE OUTCOME VARIABLE (func_decline_true) FROM YOUR ORIGINAL CODE
# ==============================================================================

cat("\n=== Creating Functional Limitation Outcome ===\n")

# Define the disease variables for multimorbidity (from your original code)
disease_vars <- c("rb1sa11s","rb1sa11x","rb1sa12d","rb1pa6a",
                  "rb1pa26","rb1sa11c","rb1sa11d","rb1pa2")

# TRUE functional impairment variables (validated ADL/IADL)
true_func_vars <- c("rb1sbadl1", "rb1pd12")

# Recreate the outcome variables exactly as in your original code
midus <- midus %>%
  mutate(
    # Multimorbidity (>=2 chronic conditions)
    disease_count = rowSums(midus[, disease_vars] == 1, na.rm = TRUE),
    multimorbidity = ifelse(disease_count >= 2, 1, 0),
    
    # TRUE FUNCTIONAL DECLINE - Using validated ADL/IADL items
    # At least ONE limitation in ADL or IADL (score ≥2)
    func_decline_true = ifelse(
      rowSums(midus[, true_func_vars] >= 2, na.rm = TRUE) >= 1,
      1, 0
    ),
    
    # Secondary outcomes for sensitivity analysis
    any_adl = ifelse(rb1sbadl1 > 1, 1, 0),
    any_iadl = ifelse(rb1pd12 > 1, 1, 0),
    
    # Keep original broad definition for comparison
    func_limit_count = rowSums(midus[, c("rb1sa24b","rb1sa24f","rb1sa24c",
                                         "rb1sa24af","rb1sa24ag")] > 1, na.rm = TRUE),
    poor_health = ifelse(rb1pa1 %in% c(4,5), 1, 0)
  )

# Create the functional limitation factor variable
midus <- midus %>%
  mutate(
    functional_limitation = factor(
      ifelse(func_decline_true == 1, "Yes", "No"),
      levels = c("No", "Yes")
    )
  )

# Check outcome prevalence
cat("\n=== Outcome Prevalence ===\n")
cat("Multimorbidity:", round(mean(midus$multimorbidity, na.rm=TRUE)*100, 1), "%\n")
cat("Functional Decline (True ADL/IADL):", 
    round(mean(midus$func_decline_true, na.rm=TRUE)*100, 1), "%\n")
cat("  - Any ADL limitation:", round(mean(midus$any_adl, na.rm=TRUE)*100, 1), "%\n")
cat("  - Any IADL limitation:", round(mean(midus$any_iadl, na.rm=TRUE)*100, 1), "%\n")

# ==============================================================================
# 4. CREATE THE 3-LEVEL SLEEP EXPOSURE VARIABLE
# ==============================================================================

cat("\n=== Creating 3-Level Sleep Exposure ===\n")

midus <- midus %>%
  mutate(
    # Create 3-level sleep variable from rb1sa57d (Unrested during day)
    # Original coding: 1=Never, 2=Rarely, 3=Sometimes, 4=Often, 5=Almost Always
    sleep_3level = factor(case_when(
      rb1sa57d %in% c(1, 2) ~ "Good",        # Never/Rarely
      rb1sa57d == 3 ~ "Intermediate",         # Sometimes
      rb1sa57d %in% c(4, 5) ~ "Poor"          # Often/Almost Always
    ), levels = c("Good", "Intermediate", "Poor")),
    
    # Original dichotomous version (for comparison)
    poor_sleep_original = factor(case_when(
      rb1sa57d %in% c(1, 2) ~ "No",
      rb1sa57d %in% c(4, 5) ~ "Yes",
      TRUE ~ NA_character_
    ), levels = c("No", "Yes"))
  )

# Check distribution
cat("\nSleep distribution (3-level):\n")
print(table(midus$sleep_3level))
print(prop.table(table(midus$sleep_3level)))

cat("\nOriginal dichotomous (excluding 'sometimes'):\n")
print(table(midus$poor_sleep_original, useNA = "ifany"))
cat(sprintf("Excluded 'sometimes': %d participants\n", 
            sum(is.na(midus$poor_sleep_original))))

# ==============================================================================
# 5. HANDLE SMOKING VARIABLE (as in your original code)
# ==============================================================================

cat("\n=== Handling Smoking Variable ===\n")

midus <- midus %>%
  mutate(
    rb1pa39 = case_when(
      is.na(rb1pa39) ~ 3,  # Unknown category
      rb1pa39 == 1 ~ 1,     # Yes
      rb1pa39 == 2 ~ 2,     # No
      TRUE ~ 3               # Default to Unknown
    )
  )

cat("Smoking status (rb1pa39) recoded: 1=Yes, 2=No, 3=Unknown\n")

# ==============================================================================
# 6. CREATE ANALYSIS DATASET
# ==============================================================================

cat("\n=== Creating Analysis Dataset ===\n")

analysis_data <- midus %>%
  dplyr::select(
    # ID (use row number as ID if no specific ID column)
    id = 1,  # Using first column as ID
    
    # Outcomes and exposure
    sleep_3level,
    functional_limitation,
    func_decline_true,
    multimorbidity,
    
    # Demographics
    age = rb1prage,
    sex = rb1prsex,
    race = rb1pf7a,
    
    # Socioeconomic
    education = rb1pb1,
    income = rb1pb16,
    marital = rb1pb19,
    employment = rb1sf17b,
    
    # Health behaviors
    bmi = rb1sbmi,
    smoking = rb1pa39,
    alcohol_freq = rb1pa51,
    exercise_freq = rb1sa52f,
    
    # Mental health
    depression = rb1pa60,
    overwhelmed = rb1se1z,
    beyond_control = rb1se4e,
    
    # Additional sleep variables for dose-response
    trouble_falling = rb1sa57a,
    night_waking = rb1sa57b
  ) %>%
  # Remove rows with missing sleep_3level
  filter(!is.na(sleep_3level)) %>%
  # Convert categorical variables to factors
  mutate(
    sex = factor(sex, labels = c("Male", "Female")),
    race = factor(race),
    marital = factor(marital),
    employment = factor(employment),
    smoking = factor(smoking, labels = c("Yes", "No", "Unknown")),
    alcohol_freq = factor(alcohol_freq),
    exercise_freq = factor(exercise_freq),
    depression = factor(ifelse(depression == 1, "Yes", "No"), 
                        levels = c("No", "Yes")),
    overwhelmed = factor(overwhelmed),
    beyond_control = factor(beyond_control),
    multimorbidity = factor(ifelse(multimorbidity == 1, "Yes", "No"),
                            levels = c("No", "Yes"))
  )

cat(sprintf("Analysis dataset: %d participants\n", nrow(analysis_data)))
cat("Sleep distribution in analysis dataset:\n")
print(table(analysis_data$sleep_3level))

# ==============================================================================
# 7. CHECK FOR MISSING VALUES AND IMPUTE IF NEEDED
# ==============================================================================

cat("\n=== Checking Missing Values ===\n")

missing_summary <- analysis_data %>%
  summarise(across(everything(), ~sum(is.na(.)))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "n_missing") %>%
  filter(n_missing > 0)

if(nrow(missing_summary) > 0) {
  print(missing_summary)
  
  # Simple median/mode imputation
  analysis_data <- analysis_data %>%
    mutate(across(where(is.numeric), ~ifelse(is.na(.), median(., na.rm = TRUE), .))) %>%
    mutate(across(where(is.factor), ~{
      if(any(is.na(.))) {
        mode_val <- names(sort(table(.), decreasing = TRUE))[1]
        fct_explicit_na(., na_level = mode_val)
      } else {
        .
      }
    }))
  
  cat("\nMissing values imputed\n")
} else {
  cat("No missing values found!\n")
}

# ==============================================================================
# 8. MULTINOMIAL PROPENSITY SCORE MODEL
# ==============================================================================

cat("\n=== Estimating Propensity Scores ===\n")

ps_model <- multinom(
  sleep_3level ~ age + sex + race + income + education + employment + marital +
    bmi + smoking + alcohol_freq + exercise_freq + depression + 
    overwhelmed + beyond_control + multimorbidity,
  data = analysis_data,
  trace = FALSE
)

# Get predicted probabilities
probs <- predict(ps_model, type = "probs")
analysis_data$ps_good <- probs[, "Good"]
analysis_data$ps_intermediate <- probs[, "Intermediate"]
analysis_data$ps_poor <- probs[, "Poor"]

cat("Propensity score ranges:\n")
cat(sprintf("  Good: %.3f - %.3f\n", min(analysis_data$ps_good), max(analysis_data$ps_good)))
cat(sprintf("  Intermediate: %.3f - %.3f\n", min(analysis_data$ps_intermediate), max(analysis_data$ps_intermediate)))
cat(sprintf("  Poor: %.3f - %.3f\n", min(analysis_data$ps_poor), max(analysis_data$ps_poor)))

# ==============================================================================
# 9. PAIRWISE PROPENSITY SCORE MATCHING
# ==============================================================================

cat("\n=== Performing Pairwise Propensity Score Matching ===\n")

# 9.1 Poor vs Good matching
data_poor_good <- analysis_data %>%
  filter(sleep_3level %in% c("Good", "Poor")) %>%
  mutate(treatment = ifelse(sleep_3level == "Poor", 1, 0))

cat("\n--- Poor vs Good ---\n")
cat(sprintf("Sample size: %d (Good: %d, Poor: %d)\n", 
            nrow(data_poor_good),
            sum(data_poor_good$sleep_3level == "Good"),
            sum(data_poor_good$sleep_3level == "Poor")))

match_pg <- matchit(
  treatment ~ age + sex + race + income + education + employment + marital +
    bmi + smoking + alcohol_freq + exercise_freq + depression + 
    overwhelmed + beyond_control + multimorbidity,
  data = data_poor_good,
  method = "nearest",
  ratio = 3,
  caliper = 0.2,
  std.caliper = TRUE
)

matched_pg <- match.data(match_pg)
cat(sprintf("Matched sample: %d (Good: %d, Poor: %d)\n", 
            nrow(matched_pg),
            sum(matched_pg$sleep_3level == "Good"),
            sum(matched_pg$sleep_3level == "Poor")))

# 9.2 Poor vs Intermediate matching
data_poor_int <- analysis_data %>%
  filter(sleep_3level %in% c("Intermediate", "Poor")) %>%
  mutate(treatment = ifelse(sleep_3level == "Poor", 1, 0))

cat("\n--- Poor vs Intermediate ---\n")
cat(sprintf("Sample size: %d (Intermediate: %d, Poor: %d)\n", 
            nrow(data_poor_int),
            sum(data_poor_int$sleep_3level == "Intermediate"),
            sum(data_poor_int$sleep_3level == "Poor")))

match_pi <- matchit(
  treatment ~ age + sex + race + income + education + employment + marital +
    bmi + smoking + alcohol_freq + exercise_freq + depression + 
    overwhelmed + beyond_control + multimorbidity,
  data = data_poor_int,
  method = "nearest",
  ratio = 3,
  caliper = 0.2,
  std.caliper = TRUE
)

matched_pi <- match.data(match_pi)
cat(sprintf("Matched sample: %d (Intermediate: %d, Poor: %d)\n", 
            nrow(matched_pi),
            sum(matched_pi$sleep_3level == "Intermediate"),
            sum(matched_pi$sleep_3level == "Poor")))

# 9.3 Intermediate vs Good matching
data_int_good <- analysis_data %>%
  filter(sleep_3level %in% c("Good", "Intermediate")) %>%
  mutate(treatment = ifelse(sleep_3level == "Intermediate", 1, 0))

cat("\n--- Intermediate vs Good ---\n")
cat(sprintf("Sample size: %d (Good: %d, Intermediate: %d)\n", 
            nrow(data_int_good),
            sum(data_int_good$sleep_3level == "Good"),
            sum(data_int_good$sleep_3level == "Intermediate")))

match_ig <- matchit(
  treatment ~ age + sex + race + income + education + employment + marital +
    bmi + smoking + alcohol_freq + exercise_freq + depression + 
    overwhelmed + beyond_control + multimorbidity,
  data = data_int_good,
  method = "nearest",
  ratio = 3,
  caliper = 0.2,
  std.caliper = TRUE
)

matched_ig <- match.data(match_ig)
cat(sprintf("Matched sample: %d (Good: %d, Intermediate: %d)\n", 
            nrow(matched_ig),
            sum(matched_ig$sleep_3level == "Good"),
            sum(matched_ig$sleep_3level == "Intermediate")))

# ==============================================================================
# 10. BALANCE ASSESSMENT
# ==============================================================================

cat("\n=== Balance Assessment ===\n")

# Check balance for each comparison
balance_pg <- bal.tab(match_pg, data = data_poor_good)
balance_pi <- bal.tab(match_pi, data = data_poor_int)
balance_ig <- bal.tab(match_ig, data = data_int_good)

cat("\nPoor vs Good - Maximum SMD after matching:\n")
print(max(abs(balance_pg$Balance$Diff.Adj), na.rm = TRUE))

cat("\nPoor vs Intermediate - Maximum SMD after matching:\n")
print(max(abs(balance_pi$Balance$Diff.Adj), na.rm = TRUE))

cat("\nIntermediate vs Good - Maximum SMD after matching:\n")
print(max(abs(balance_ig$Balance$Diff.Adj), na.rm = TRUE))

# Create Love plot for Poor vs Good
balance_df <- as.data.frame(balance_pg$Balance)
balance_df$Variable <- rownames(balance_df)

love_data <- balance_df %>%
  dplyr::select(Variable, Diff.Un, Diff.Adj) %>%
  pivot_longer(cols = c(Diff.Un, Diff.Adj), 
               names_to = "Type", 
               values_to = "SMD") %>%
  mutate(Type = recode(Type, "Diff.Un" = "Unmatched", "Diff.Adj" = "Matched"))

love_plot <- ggplot(love_data, aes(x = abs(SMD), y = reorder(Variable, abs(SMD)), 
                                   color = Type, shape = Type)) +
  geom_point(size = 2.5, alpha = 0.7) +
  geom_vline(xintercept = 0.1, linetype = "dashed", color = "gray50", alpha = 0.5) +
  geom_vline(xintercept = 0.2, linetype = "dashed", color = "gray30", alpha = 0.5) +
  scale_x_continuous(limits = c(0, max(abs(love_data$SMD), na.rm = TRUE) + 0.1)) +
  scale_color_manual(values = c("Unmatched" = "#E41A1C", "Matched" = "#377EB8")) +
  theme_minimal() +
  labs(title = "Covariate Balance: Poor vs Good Sleep",
       x = "Absolute Standardized Mean Difference", y = "") +
  theme(legend.position = "bottom")

print(love_plot)
ggsave("love_plot.png", love_plot, width = 10, height = 8, dpi = 300)

# ==============================================================================
# 11. PRIMARY OUTCOME ANALYSES
# ==============================================================================

cat("\n=== Primary Matched Analyses ===\n")

# Function to run matched analysis
run_matched_analysis <- function(data) {
  model <- glm(
    functional_limitation ~ treatment + overwhelmed + beyond_control + education,
    data = data,
    family = binomial()
  )
  
  vcov_mat <- vcovCL(model, cluster = data$subclass, type = "HC0")
  coef_test <- coeftest(model, vcov = vcov_mat)
  
  or <- exp(coef(model)[2])
  ci <- exp(confint(model, parm = 2, vcov. = vcov_mat))
  p_val <- coef_test[2, 4]
  
  return(list(or = or, ci_lower = ci[1], ci_upper = ci[2], p = p_val))
}

# Run analyses
pg_results <- run_matched_analysis(matched_pg)
pi_results <- run_matched_analysis(matched_pi)
ig_results <- run_matched_analysis(matched_ig)

# Compile results
primary_results <- data.frame(
  Comparison = c("Poor vs Good", "Poor vs Intermediate", "Intermediate vs Good"),
  OR = c(pg_results$or, pi_results$or, ig_results$or),
  CI_lower = c(pg_results$ci_lower, pi_results$ci_lower, ig_results$ci_lower),
  CI_upper = c(pg_results$ci_upper, pi_results$ci_upper, ig_results$ci_upper),
  p_value = c(pg_results$p, pi_results$p, ig_results$p)
)

print(primary_results)

# ==============================================================================
# 12. DOSE-RESPONSE ANALYSIS
# ==============================================================================

cat("\n=== Dose-Response Analysis ===\n")

# Create sleep scores
matched_pg <- matched_pg %>%
  mutate(sleep_score = ifelse(sleep_3level == "Poor", 2, 0))

matched_pi <- matched_pi %>%
  mutate(sleep_score = case_when(
    sleep_3level == "Poor" ~ 2,
    sleep_3level == "Intermediate" ~ 1
  ))

matched_ig <- matched_ig %>%
  mutate(sleep_score = case_when(
    sleep_3level == "Intermediate" ~ 1,
    sleep_3level == "Good" ~ 0
  ))

# Combine matched samples
combined_matched <- bind_rows(matched_pg, matched_pi, matched_ig) %>%
  distinct(id, .keep_all = TRUE)  # Remove duplicates

cat(sprintf("Combined matched sample size: %d\n", nrow(combined_matched)))
cat("Sleep score distribution:\n")
print(table(combined_matched$sleep_score))

# Create sleep level factor
combined_matched <- combined_matched %>%
  mutate(sleep_level = factor(case_when(
    sleep_score == 0 ~ "Good",
    sleep_score == 1 ~ "Intermediate",
    sleep_score == 2 ~ "Poor"
  ), levels = c("Good", "Intermediate", "Poor")))

# Dose-response model
dose_model <- glm(
  functional_limitation ~ sleep_level + overwhelmed + beyond_control + education,
  data = combined_matched,
  family = binomial()
)

dose_results <- tidy(dose_model, conf.int = TRUE, exponentiate = TRUE) %>%
  filter(grepl("sleep_level", term))

print(dose_results)

# Test for linear trend
trend_model <- glm(
  functional_limitation ~ sleep_score + overwhelmed + beyond_control + education,
  data = combined_matched,
  family = binomial()
)

trend_result <- tidy(trend_model, conf.int = TRUE, exponentiate = TRUE) %>%
  filter(term == "sleep_score")

cat("\nLinear trend per level increase:\n")
print(trend_result)

# ==============================================================================
# 13. SUBGROUP ANALYSES (Poor vs Good)
# ==============================================================================

cat("\n=== Subgroup Analyses ===\n")

matched_pg <- matched_pg %>%
  mutate(age_group = ifelse(age < 50, "<50", "≥50"))

# Function for subgroup analysis
run_subgroup <- function(data, subset_condition) {
  subset_data <- tryCatch({
    data %>% filter(!!rlang::parse_expr(subset_condition))
  }, error = function(e) NULL)
  
  if(is.null(subset_data) || nrow(subset_data) < 50) return(NULL)
  
  model <- glm(
    functional_limitation ~ treatment + overwhelmed + beyond_control + education,
    data = subset_data,
    family = binomial()
  )
  
  if(!model$converged) return(NULL)
  
  or <- exp(coef(model)["treatment"])
  ci <- exp(confint(model, parm = "treatment"))
  
  return(data.frame(OR = or, CI_lower = ci[1], CI_upper = ci[2]))
}

# Run subgroup analyses
subgroup_results <- bind_rows(
  data.frame(Subgroup = "Age <50", run_subgroup(matched_pg, "age < 50")),
  data.frame(Subgroup = "Age ≥50", run_subgroup(matched_pg, "age >= 50")),
  data.frame(Subgroup = "Male", run_subgroup(matched_pg, "sex == 'Male'")),
  data.frame(Subgroup = "Female", run_subgroup(matched_pg, "sex == 'Female'"))
)

print(subgroup_results)

# ==============================================================================
# 14. SENSITIVITY ANALYSES - CORRECTED VERSION
# ==============================================================================

cat("\n=== Sensitivity Analyses ===\n")

# 14.1 1:1 matching
match_pg_1to1 <- matchit(
  treatment ~ age + sex + race + income + education + employment + marital +
    bmi + smoking + alcohol_freq + exercise_freq + depression + 
    overwhelmed + beyond_control + multimorbidity,
  data = data_poor_good,
  method = "nearest",
  ratio = 1,
  caliper = 0.2,
  std.caliper = TRUE
)

matched_pg_1to1 <- match.data(match_pg_1to1)
model_1to1 <- glm(
  functional_limitation ~ treatment + overwhelmed + beyond_control + education,
  data = matched_pg_1to1,
  family = binomial()
)

# 14.2 Stricter caliper (0.1)
match_pg_strict <- matchit(
  treatment ~ age + sex + race + income + education + employment + marital +
    bmi + smoking + alcohol_freq + exercise_freq + depression + 
    overwhelmed + beyond_control + multimorbidity,
  data = data_poor_good,
  method = "nearest",
  ratio = 3,
  caliper = 0.1,
  std.caliper = TRUE
)

matched_pg_strict <- match.data(match_pg_strict)
model_strict <- glm(
  functional_limitation ~ treatment + overwhelmed + beyond_control + education,
  data = matched_pg_strict,
  family = binomial()
)

# Compile sensitivity results - USING PRIMARY_RESULTS FOR THE PRIMARY ANALYSIS
sensitivity_results <- data.frame(
  Analysis = c("Primary analysis (1:3 matching, caliper=0.2)", 
               "1:1 matching", 
               "Stricter caliper (0.1)"),
  OR = c(
    primary_results$OR[primary_results$Comparison == "Poor vs Good"],  # Using the correct 2.26
    exp(coef(model_1to1)["treatment"]),
    exp(coef(model_strict)["treatment"])
  ),
  CI_lower = c(
    primary_results$CI_lower[primary_results$Comparison == "Poor vs Good"],
    exp(confint(model_1to1, parm = "treatment"))[1],
    exp(confint(model_strict, parm = "treatment"))[1]
  ),
  CI_upper = c(
    primary_results$CI_upper[primary_results$Comparison == "Poor vs Good"],
    exp(confint(model_1to1, parm = "treatment"))[2],
    exp(confint(model_strict, parm = "treatment"))[2]
  ),
  N = c(nrow(matched_pg), nrow(matched_pg_1to1), nrow(matched_pg_strict))
)

print(sensitivity_results)

# ==============================================================================
# 15. VISUALIZATIONS
# ==============================================================================

cat("\n=== Generating Visualizations ===\n")

# Forest plot
forest_data <- primary_results %>%
  mutate(Comparison = factor(Comparison, 
                             levels = c("Poor vs Good", "Poor vs Intermediate", "Intermediate vs Good")))

forest_plot <- ggplot(forest_data, aes(x = OR, y = Comparison, 
                                       xmin = CI_lower, xmax = CI_upper)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
  geom_point(size = 3, color = "#377EB8") +
  geom_errorbarh(height = 0.2, color = "#377EB8") +
  scale_x_log10(breaks = c(0.5, 0.75, 1, 1.5, 2, 3)) +
  theme_minimal() +
  labs(title = "Propensity Score Matched Analysis",
       subtitle = "Association between sleep disturbances and functional limitation",
       x = "Odds Ratio (95% CI)", y = "")

print(forest_plot)
ggsave("forest_plot.png", forest_plot, width = 8, height = 4, dpi = 300)

# Dose-response plot
dose_plot_data <- dose_results %>%
  mutate(sleep_level = case_when(
    grepl("Intermediate", term) ~ "Intermediate",
    grepl("Poor", term) ~ "Poor"
  )) %>%
  bind_rows(data.frame(
    term = "Good (ref)",
    estimate = 1,
    conf.low = 1,
    conf.high = 1,
    sleep_level = "Good"
  )) %>%
  mutate(sleep_level = factor(sleep_level, levels = c("Good", "Intermediate", "Poor")))

dose_plot <- ggplot(dose_plot_data, aes(x = sleep_level, y = estimate, 
                                        ymin = conf.low, ymax = conf.high)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
  geom_point(size = 4, color = "#E41A1C") +
  geom_errorbar(width = 0.2, color = "#E41A1C") +
  scale_y_log10(breaks = c(0.8, 1, 1.2, 1.4, 1.6, 1.8, 2.0)) +
  theme_minimal() +
  labs(title = "Dose-Response Relationship",
       subtitle = "Odds of functional limitation by sleep group",
       x = "Sleep Group", y = "Odds Ratio (95% CI)")

print(dose_plot)
ggsave("dose_response.png", dose_plot, width = 6, height = 5, dpi = 300)

# ==============================================================================
# 16. EXPORT RESULTS
# ==============================================================================

cat("\n=== Exporting Results ===\n")

# Create results directory
results_dir <- "results_matched_3level"
if(!dir.exists(results_dir)) dir.create(results_dir)

# Save results as CSV
write.csv(primary_results, file.path(results_dir, "primary_results.csv"), row.names = FALSE)
write.csv(dose_results, file.path(results_dir, "dose_response.csv"), row.names = FALSE)
write.csv(trend_result, file.path(results_dir, "linear_trend.csv"), row.names = FALSE)
if(exists("subgroup_results") && nrow(subgroup_results) > 0) {
  write.csv(subgroup_results, file.path(results_dir, "subgroup_results.csv"), row.names = FALSE)
}
write.csv(sensitivity_results, file.path(results_dir, "sensitivity_results.csv"), row.names = FALSE)

# Save R objects
saveRDS(list(
  data = analysis_data,
  matched_pg = matched_pg,
  matched_pi = matched_pi,
  matched_ig = matched_ig,
  combined = combined_matched,
  primary = primary_results,
  dose = dose_results,
  trend = trend_result,
  subgroup = subgroup_results,
  sensitivity = sensitivity_results
), file = file.path(results_dir, "analysis_results.rds"))

# Create manuscript summary
sink(file.path(results_dir, "manuscript_summary.txt"))

cat("==========================================================\n")
cat("PROPENSITY SCORE MATCHED ANALYSIS WITH 3-LEVEL EXPOSURE\n")
cat("MIDUS Refresher 2 Data\n")
cat("==========================================================\n\n")

cat("SAMPLE CHARACTERISTICS:\n")
cat(sprintf("Total N: %d\n", nrow(analysis_data)))
cat(sprintf("Good Sleep: %d (%.1f%%)\n", 
            sum(analysis_data$sleep_3level == "Good"),
            100 * mean(analysis_data$sleep_3level == "Good")))
cat(sprintf("Intermediate Sleep: %d (%.1f%%)\n", 
            sum(analysis_data$sleep_3level == "Intermediate"),
            100 * mean(analysis_data$sleep_3level == "Intermediate")))
cat(sprintf("Poor Sleep: %d (%.1f%%)\n\n", 
            sum(analysis_data$sleep_3level == "Poor"),
            100 * mean(analysis_data$sleep_3level == "Poor")))

cat("PRIMARY MATCHED RESULTS:\n")
print(primary_results)

cat("\nDOSE-RESPONSE ANALYSIS:\n")
print(dose_results)
cat(sprintf("\nLinear trend: OR = %.2f per level (95%% CI: %.2f-%.2f), p = %.4f\n",
            trend_result$estimate, trend_result$conf.low, 
            trend_result$conf.high, trend_result$p.value))

if(exists("subgroup_results") && nrow(subgroup_results) > 0) {
  cat("\nSUBGROUP ANALYSES (Poor vs Good):\n")
  print(subgroup_results)
}

cat("\nSENSITIVITY ANALYSES (Poor vs Good):\n")
print(sensitivity_results)

sink()

cat("\n=== ANALYSIS COMPLETE ===\n")
cat(sprintf("All results saved to '%s' directory\n", results_dir))
cat("Files created:\n")
cat("  - primary_results.csv\n")
cat("  - dose_response.csv\n")
cat("  - linear_trend.csv\n")
cat("  - subgroup_results.csv\n")
cat("  - sensitivity_results.csv\n")
cat("  - love_plot.png\n")
cat("  - forest_plot.png\n")
cat("  - dose_response.png\n")
cat("  - analysis_results.rds\n")
cat("  - manuscript_summary.txt\n")
