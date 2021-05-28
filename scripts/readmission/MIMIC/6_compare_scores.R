# Packages
require(dplyr)
require(readr)
require(forcats)
require(magrittr)
require(lubridate)
require(ROCR)
require(ggplot2)
require(tidyr)
require(ResourceSelection)
require(s2dverification)
require(scales)
require(caret)
require(glmnet)
require(mice)
require(doParallel)
require(tools)

# Load data
patients <- read_csv("data/predictors.csv")

# Load functions
source("functions/NA_count.R")
source("functions/inverse_logit.R")
source("functions/nomogram_convert.R")
source("functions/calibration.R")
source("functions/brier_extraction.R")

crosstab <- function(varname, .df = patients)
{
  options(dplyr.summarise.inform=F)
  
  # Find
  .df  %>% 
    group_by_at(c("readmission", varname)) %>% 
    summarise(n = n()) %>% 
    mutate(`%` = (n/sum(n) * 100)) %>% 
    mutate(`n (%)` = paste(n, " (", round(`%`,1),
                           "%)", sep = "")) %>% 
    select(any_of(c("readmission", varname, "n (%)"))) %>% 
    pivot_wider(names_from = "readmission",
                values_from = "n (%)") %>% 
    as.data.frame()
}

# Filter missing data----------

# Assess NA
sapply(patients, function(x) sum(is.na(x)))

# Count and exclude missing chart data
sum(patients$chart_missing)
patients %<>% filter(chart_missing == FALSE)

# Count and exclude missing events data
sum(is.na(patients$fluid_balance_5L))
patients %<>% filter(!is.na(patients$fluid_balance_5L))

# Filter patients missing respiratory rate assessments
sum(is.na(patients$respiratory_rate_initial))
patients %<>% filter(!is.na(respiratory_rate_initial))

# Filter patients missing functional mobility assessments
sum(is.na(patients$ambulation))
patients %<>% filter(!is.na(patients$ambulation))

# Filter patients with no initial or final blood labs
sum(is.na(patients$hyperglycemia_initial) | is.na(patients$anaemia_initial) |
      is.na(patients$hyperglycemia_final) | is.na(patients$anaemia_final) |
      is.na(patients$blood_urea_nitrogen_initial) | is.na(patients$blood_urea_nitrogen_final) |
      is.na(patients$serum_chloride_initial) | is.na(patients$serum_chloride_final))
patients %<>% filter(!is.na(patients$hyperglycemia_initial) & !is.na(patients$anaemia_initial) &
                       !is.na(patients$hyperglycemia_final) & !is.na(patients$anaemia_final) &
                       !is.na(patients$blood_urea_nitrogen_initial) & !is.na(patients$blood_urea_nitrogen_final) &
                       !is.na(patients$serum_chloride_initial) & !is.na(patients$serum_chloride_final))

# Confirm no remaining NA
sapply(patients, function(x) sum(is.na(x)))

# Data cleanup-------

# Clear up outcome
patients %<>% 
  mutate(readmission = fct_recode(as.character(readmission),
                                  `No Readmission` = "FALSE",
                                  `Readmitted to ICU` = "TRUE"))
# Readmission
patients %>% 
  group_by(readmission) %>% 
  summarise(n = n()) %>% 
  mutate(`%` = (n/sum(n) * 100))

# Fix age
patients %<>% 
  mutate(age = ifelse(age == ">90", 90, age)) %>% 
  mutate(age = floor(as.numeric(age)))

# Fix admission source
patients %<>% 
  mutate(admission_source = fct_relevel(admission_source, "OT"))

# Hammer: Aligned----------

# Cross-tabulate
crosstab("general_surgery")
crosstab("cardiac_surgery")
crosstab("hyperglycemia_initial")
crosstab("anaemia_initial")
crosstab("high_apache")
crosstab("fluid_balance_5L")
crosstab("ambulation")
crosstab("los_5")

# Score data - coefficients
coefficients_hammer <- log(0.005 / (1 - 0.005)) +
  ifelse(patients$sex == "M", 0.43, 0) +
  ifelse(patients$general_surgery, 0.64, 0) +
  ifelse(patients$cardiac_surgery, -1.01, 0) +
  ifelse(patients$hyperglycemia_initial, 0.36, 0) +
  ifelse(patients$anaemia_initial, 1.11, 0) +
  ifelse(patients$high_apache, 0.44, 0) +
  ifelse(patients$fluid_balance_5L, 0.52, 0) +
  ifelse(patients$ambulation == FALSE, 0.7, 0) +
  ifelse(patients$los_5, 0.88, 0)

# Convert coefficients to probs
probs_coefficients_hammer <- inverse_logit(coefficients_hammer)
probs_hammer_A <- probs_coefficients_hammer

# Martin: Aligned------------

patients %>% 
  group_by(readmission) %>% 
  summarise(mean = mean(respiratory_rate_initial),
            sd = sd(respiratory_rate_initial))

patients %>% 
  group_by(readmission) %>% 
  summarise(mean = mean(blood_urea_nitrogen_initial),
            sd = sd(blood_urea_nitrogen_initial))

patients %>% 
  group_by(readmission) %>% 
  summarise(mean = mean(serum_glucose_initial),
            sd = sd(serum_glucose_initial)) %>% 
  as.data.frame()

patients %>% 
  group_by(readmission) %>% 
  summarise(mean = mean(serum_chloride_initial),
            sd = sd(serum_chloride_initial))%>% 
  as.data.frame()

patients %>% 
  group_by(readmission) %>% 
  summarise(mean = mean(age),
            sd = sd(age))

crosstab("atrial_fibrillation")
crosstab("renal_insufficiency")


# Score data - coefficients
coefficients_martin <- -9.284491 +
  (0.04883 * patients$respiratory_rate_initial) +
  (0.011588 * patients$age) +
  (0.036104 * patients$serum_chloride_initial) +
  (0.004967 * patients$blood_urea_nitrogen_initial) +
  (0.580153 * patients$atrial_fibrillation) +
  (0.458202 * patients$renal_insufficiency) +
  (0.003519 * patients$serum_glucose_initial)
  
# Convert to probabilities
probs_martin_A <- inverse_logit(coefficients_martin)

# Frost: Aligned---------------

patients %>% 
  group_by(readmission) %>% 
  summarise(mean = mean(age),
            sd = sd(age))

crosstab("sex")
crosstab("admission_source")

patients %>% 
  group_by(readmission) %>% 
  summarise(mean = mean(apache_II),
            sd = sd(apache_II))

crosstab("los_7")
crosstab("after_hours_discharge")
crosstab("acute_renal_failure_initial")

# Create scoring system from nomogram for age
age_score_system_input <- data.frame(
  input = seq(15, 90, 5),
  pix = seq(2, 287, length.out = 16)
)
age_score_system_output <- data.frame(
  output = seq(0, 8, 1),
  pix = seq(2, 290, length.out = 9)
)

# Create scoring system from nomogram for APACHE-II
apache_score_system_input <- data.frame(
  input = seq(0, 45, 5),
  pix = seq(12, 733, length.out = 10)
)
apache_score_system_output <- data.frame(
  output = seq(0, 20, 1),
  pix = seq(12, 733, length.out = 21)
)

# Score data - points
scores_frost <- nomogram_convert(patients$age, age_score_system_input,
                 age_score_system_output) +
  ifelse(patients$sex == "M", 2, 0) +
  ifelse(patients$elective_admission == FALSE, 12.25, 0) +
  case_when(
    patients$admission_source == "OT" ~ 0,
    patients$admission_source == "Emergency" ~ 9.5,
    patients$admission_source == "OtherHosp" ~ 10.25,
    patients$admission_source == "Ward" ~ 14.75
  ) +
  nomogram_convert(patients$apache_II, apache_score_system_input,
                   apache_score_system_output) +
  ifelse(patients$los_7, 17, 0) +
  ifelse(patients$after_hours_discharge, 4, 0) +
  ifelse(patients$acute_renal_failure_initial, 10, 0)

# Create nomogram system for total points
points_system_input <- data.frame(
  input = seq(0, 100, 5),
  pix = seq(7, 727, length.out = 21)
)
points_system_output <- data.frame(
  output = seq(0.05, 0.5, 0.05),
  pix = c(212, 330, 403, 460, 506,
          547, 583, 618, 652, 683)
)

# Convert using logarithmic transform
probs_frost_A <- nomogram_convert(scores_frost, points_system_input,
                 points_system_output, log = TRUE)

# Hammer: Original----------

# Cross-tabulate
crosstab("general_surgery")
crosstab("cardiac_surgery")
crosstab("hyperglycemia_final")
crosstab("anaemia_final")
crosstab("high_apache")
crosstab("fluid_balance_5L")
crosstab("ambulation")
crosstab("los_5")

# Score data - coefficients
coefficients_hammer <- log(0.005 / (1 - 0.005)) +
  ifelse(patients$sex == "M", 0.43, 0) +
  ifelse(patients$general_surgery, 0.64, 0) +
  ifelse(patients$cardiac_surgery, -1.01, 0) +
  ifelse(patients$hyperglycemia_final, 0.36, 0) +
  ifelse(patients$anaemia_final, 1.11, 0) +
  ifelse(patients$high_apache, 0.44, 0) +
  ifelse(patients$fluid_balance_5L, 0.52, 0) +
  ifelse(patients$ambulation == FALSE, 0.7, 0) +
  ifelse(patients$los_5, 0.88, 0)

# Convert coefficients to probs
probs_coefficients_hammer <- inverse_logit(coefficients_hammer)
probs_hammer_O <- probs_coefficients_hammer

# Martin: Original------------

patients %>% 
  group_by(readmission) %>% 
  summarise(mean = mean(respiratory_rate_final),
            sd = sd(respiratory_rate_final))

patients %>% 
  group_by(readmission) %>% 
  summarise(mean = mean(blood_urea_nitrogen_final),
            sd = sd(blood_urea_nitrogen_final))

patients %>% 
  group_by(readmission) %>% 
  summarise(mean = mean(serum_glucose_final),
            sd = sd(serum_glucose_final)) %>% 
  as.data.frame()

patients %>% 
  group_by(readmission) %>% 
  summarise(mean = mean(serum_chloride_final),
            sd = sd(serum_chloride_final))%>% 
  as.data.frame()

patients %>% 
  group_by(readmission) %>% 
  summarise(mean = mean(age),
            sd = sd(age))

crosstab("atrial_fibrillation")
crosstab("renal_insufficiency")

# Score data - coefficients
coefficients_martin <- -9.284491 +
  (0.04883 * patients$respiratory_rate_final) +
  (0.011588 * patients$age) +
  (0.036104 * patients$serum_chloride_final) +
  (0.004967 * patients$blood_urea_nitrogen_final) +
  (0.580153 * patients$atrial_fibrillation) +
  (0.458202 * patients$renal_insufficiency) +
  (0.003519 * patients$serum_glucose_final)

# Convert to probabilities
probs_martin_O <- inverse_logit(coefficients_martin)

# Frost: Original---------------

patients %>% 
  group_by(readmission) %>% 
  summarise(mean = mean(age),
            sd = sd(age))

crosstab("sex")
crosstab("admission_source")

patients %>% 
  group_by(readmission) %>% 
  summarise(mean = mean(apache_II),
            sd = sd(apache_II))

crosstab("los_7")
crosstab("after_hours_discharge")
crosstab("acute_renal_failure_total")

# Create scoring system from nomogram for age
age_score_system_input <- data.frame(
  input = seq(15, 90, 5),
  pix = seq(2, 287, length.out = 16)
)
age_score_system_output <- data.frame(
  output = seq(0, 8, 1),
  pix = seq(2, 290, length.out = 9)
)

# Create scoring system from nomogram for APACHE-II
apache_score_system_input <- data.frame(
  input = seq(0, 45, 5),
  pix = seq(12, 733, length.out = 10)
)
apache_score_system_output <- data.frame(
  output = seq(0, 20, 1),
  pix = seq(12, 733, length.out = 21)
)

# Score data - points
scores_frost <- nomogram_convert(patients$age, age_score_system_input,
                                 age_score_system_output) +
  ifelse(patients$sex == "M", 2, 0) +
  ifelse(patients$elective_admission == FALSE, 12.25, 0) +
  case_when(
    patients$admission_source == "OT" ~ 0,
    patients$admission_source == "Emergency" ~ 9.5,
    patients$admission_source == "OtherHosp" ~ 10.25,
    patients$admission_source == "Ward" ~ 14.75
  ) +
  nomogram_convert(patients$apache_II, apache_score_system_input,
                   apache_score_system_output) +
  ifelse(patients$los_7, 17, 0) +
  ifelse(patients$after_hours_discharge, 4, 0) +
  ifelse(patients$acute_renal_failure_total, 10, 0)

# Create nomogram system for total points
points_system_input <- data.frame(
  input = seq(0, 100, 5),
  pix = seq(7, 727, length.out = 21)
)
points_system_output <- data.frame(
  output = seq(0.05, 0.5, 0.05),
  pix = c(212, 330, 403, 460, 506,
          547, 583, 618, 652, 683)
)

# Convert using logarithmic transform
probs_frost_O <- nomogram_convert(scores_frost, points_system_input,
                                  points_system_output, log = TRUE)

# Discrimination----------

# Create prediction objects
prediction_hammer_A <- prediction(probs_hammer_A, patients$readmission)
prediction_martin_A <- prediction(probs_martin_A, patients$readmission)
prediction_frost_A <- prediction(probs_frost_A, patients$readmission)
prediction_hammer_O <- prediction(probs_hammer_O, patients$readmission)
prediction_martin_O <- prediction(probs_martin_O, patients$readmission)
prediction_frost_O <- prediction(probs_frost_O, patients$readmission)

# Create performance objects
performance_hammer_A <- performance(prediction_hammer_A, "tpr", "fpr")
performance_martin_A <- performance(prediction_martin_A, "tpr", "fpr")
performance_frost_A <- performance(prediction_frost_A, "tpr", "fpr")
performance_hammer_O <- performance(prediction_hammer_O, "tpr", "fpr")
performance_martin_O <- performance(prediction_martin_O, "tpr", "fpr")
performance_frost_O <- performance(prediction_frost_O, "tpr", "fpr")

# Create AUC objects
auc_hammer_A <- performance(prediction_hammer_A, measure = "auc")
auc_martin_A <- performance(prediction_martin_A, measure = "auc")
auc_frost_A <- performance(prediction_frost_A, measure = "auc")
auc_hammer_O <- performance(prediction_hammer_O, measure = "auc")
auc_martin_O <- performance(prediction_martin_O, measure = "auc")
auc_frost_O <- performance(prediction_frost_O, measure = "auc")


# Print AUC
auc_hammer_A@y.values[[1]]
auc_martin_A@y.values[[1]]
auc_frost_A@y.values[[1]]
auc_hammer_O@y.values[[1]]
auc_martin_O@y.values[[1]]
auc_frost_O@y.values[[1]]

# Plot AUC
data.frame(x = performance_hammer_A@x.values[[1]],
           y = performance_hammer_A@y.values[[1]],
           model = "Hammer_A") %>% 
  rbind(
    data.frame(x = performance_martin_A@x.values[[1]],
               y = performance_martin_A@y.values[[1]],
               model = "Martin_A"),
    data.frame(x = performance_frost_A@x.values[[1]],
               y = performance_frost_A@y.values[[1]],
               model = "Frost_A"),
    data.frame(x = performance_hammer_O@x.values[[1]],
               y = performance_hammer_O@y.values[[1]],
               model = "Hammer_O"),
    data.frame(x = performance_martin_O@x.values[[1]],
               y = performance_martin_O@y.values[[1]],
               model = "Martin_O"),
    data.frame(x = performance_frost_O@x.values[[1]],
               y = performance_frost_O@y.values[[1]],
               model = "Frost_O")
  ) %>% 
  ggplot(aes(x, y, colour = model))+
  geom_abline(slope = 1, intercept = 0,
              linetype = "dotted",
              size = 1)+
  geom_path(size = 1)+
  labs(x = "1 - Specificity",
       y = "Sensitivity")+
  theme_classic(20)+
  theme(legend.position = "top")+
  scale_color_brewer(palette = "Paired",
                     name = "")

# Calibration----------

# Split data into deciles
deciles_df <- tibble(
  patient_id = 1:nrow(patients),
  readmission = patients$readmission == "Readmitted to ICU",
  probs_hammer_A,
  decile_hammer_A = ntile(probs_hammer_A, 10),
  probs_martin_A,
  decile_martin_A = ntile(probs_martin_A, 10),
  probs_frost_A,
  decile_frost_A = ntile(probs_frost_A, 10),
  probs_hammer_O,
  decile_hammer_O = ntile(probs_hammer_O, 10),
  probs_martin_O,
  decile_martin_O = ntile(probs_martin_O, 10),
  probs_frost_O,
  decile_frost_O = ntile(probs_frost_O, 10)
)

# Calculate calibration
cal_hammer_A <- calibration(deciles_df, "probs_hammer_A", "decile_hammer_A") %>% 
  mutate(model = "hammer_A")
cal_martin_A <- calibration(deciles_df, "probs_martin_A", "decile_martin_A") %>% 
  mutate(model = "martin_A")
cal_frost_A <- calibration(deciles_df, "probs_frost_A", "decile_frost_A") %>% 
  mutate(model = "frost_A")
cal_hammer_O <- calibration(deciles_df, "probs_hammer_O", "decile_hammer_O") %>% 
  mutate(model = "hammer_O")
cal_martin_O <- calibration(deciles_df, "probs_martin_O", "decile_martin_O") %>% 
  mutate(model = "martin_O")
cal_frost_O <- calibration(deciles_df, "probs_frost_O", "decile_frost_O") %>% 
  mutate(model = "frost_O")

# Plot
rbind(cal_hammer_A %>% select(-decile_hammer_A), 
      cal_martin_A %>% select(-decile_martin_A),
      cal_frost_A %>% select(-decile_frost_A),
      cal_hammer_O %>% select(-decile_hammer_O), 
      cal_martin_O %>% select(-decile_martin_O),
      cal_frost_O %>% select(-decile_frost_O)) %>% 
  ggplot(aes(predicted, observed ,
             colour = model))+
  geom_abline(slope = 1, intercept = 0,
              size = 1)+
  geom_path(size = 1)+
  geom_pointrange(aes(ymin = observed - error ,
                      ymax = observed + error,
                      y = observed,
                      x = predicted,
                      colour = model)) +
  theme_classic(20)+
  theme(legend.position = "top")+
  labs(x = "Predicted readmission",
       y = "Observed readmission")+
  scale_colour_brewer(palette = "Paired",
                      name = "")+
  facet_wrap(~model, scales = "free_x")


# Calculate hosmer-lemeshow chi-squared
hoslem_hammer_A <- hoslem.test(patients$readmission == "Readmitted to ICU",
            probs_hammer_A, g = 10)
hoslem_martin_A <- hoslem.test(patients$readmission == "Readmitted to ICU",
            probs_martin_A, g = 10)
hoslem_frost_A <- hoslem.test(patients$readmission == "Readmitted to ICU",
            probs_frost_A, g = 10)
hoslem_hammer_O <- hoslem.test(patients$readmission == "Readmitted to ICU",
                               probs_hammer_O, g = 10)
hoslem_martin_O <- hoslem.test(patients$readmission == "Readmitted to ICU",
                               probs_martin_O, g = 10)
hoslem_frost_O <- hoslem.test(patients$readmission == "Readmitted to ICU",
                              probs_frost_O, g = 10)
hoslem_hammer_A$statistic
hoslem_martin_A$statistic
hoslem_frost_A$statistic
hoslem_hammer_O$statistic
hoslem_martin_O$statistic
hoslem_frost_O$statistic

# Write final patients file
write_csv(patients, "data/final_patients.csv")

# Combine and write----

list(model = "hammer",
     data = "MIMIC_A",
     discrimination = prediction_hammer_A,
     calibration = hoslem_hammer_A,
     deciles = cal_hammer_A) %>% 
  write_rds("scripts/readmission/shared/models/MIMIC_hammer_A.RDS")

list(model = "frost",
     data = "MIMIC_A",
     discrimination = prediction_frost_A,
     calibration = hoslem_frost_A,
     deciles = cal_frost_A) %>% 
  write_rds("scripts/readmission/shared/models/MIMIC_frost_A.RDS")

list(model = "martin",
     data = "MIMIC_A",
     discrimination = prediction_martin_A,
     calibration = hoslem_martin_A,
     deciles = cal_martin_A) %>% 
  write_rds("scripts/readmission/shared/models/MIMIC_martin_A.RDS")

list(model = "hammer",
     data = "MIMIC_O",
     discrimination = prediction_hammer_O,
     calibration = hoslem_hammer_O,
     deciles = cal_hammer_O) %>% 
  write_rds("scripts/readmission/shared/models/MIMIC_hammer_O.RDS")

list(model = "frost",
     data = "MIMIC_O",
     discrimination = prediction_frost_O,
     calibration = hoslem_frost_O,
     deciles = cal_frost_O) %>% 
  write_rds("scripts/readmission/shared/models/MIMIC_frost_O.RDS")

list(model = "martin",
     data = "MIMIC_O",
     discrimination = prediction_martin_O,
     calibration = hoslem_martin_O,
     deciles = cal_martin_O) %>% 
  write_rds("scripts/readmission/shared/models/MIMIC_martin_O.RDS")

# Bespoke model----

# Vector of acceptable features
feature_id <- c(7:13, 15:18, 20, 22, 24:33, 35:40, 44:58)

# Impute missing
results_mice <- read_csv("data/predictors.csv") %>% 
  select(any_of(feature_id)) %>% 
  mutate(age = ifelse(age == ">90", 90, age)) %>% 
  mutate(age = floor(as.numeric(age))) %>% 
  mice(m = 10, maxit = 10)

model_options <- foreach(i = 1:results_mice$m, .combine = "rbind") %do%
  {
    print(i)
    
    # Binarise readmission column
    results_imputed <- complete(results_mice, i) %>%
      as_tibble() %>% 
      mutate(readmission = ifelse(readmission == TRUE, 1, 0),
             high_apache  = as.logical(high_apache),
             hyperglycemia_initial = as.logical(hyperglycemia_initial),
             anaemia_initial = as.logical(anaemia_initial),
             ambulation = as.logical(ambulation),
             fluid_balance_5L = as.logical(fluid_balance_5L),
             glasgow_coma_below_15 = as.logical(glasgow_coma_below_15),
             discharge_delay = as.logical(discharge_delay),
             hyperglycemia_final = as.logical(hyperglycemia_final),
             anaemia_final = as.logical(anaemia_final),) %>% 
      na.omit()
    
    # Select predictor and outcome vectors
    x <- model.matrix(readmission~., results_imputed)[,-1]
    y <- results_imputed$readmission
    
    # Determine lambda
    cv <- cv.glmnet(x, y, alpha = 1, folds = nrow(results_imputed),
                    family = "binomial", type.measure = "auc")
    
    # Fit model
    initial_model <- glmnet(x, y, alpha = 1, lambda = quantile(cv$lambda, 0.8),
                            family = "binomial")
    
    # Output
    output <- as.matrix(initial_model$beta) %>% 
      as.data.frame()
    output$var <- rownames(output)
    output$iteration <- i
    output
  }

# Identify retained variables
model_options %>%
  group_by(var) %>% 
  summarise(prop = mean(s0 != 0)) %>% 
  filter(prop == 1)

# Loop through datasets
cl <- makeCluster(12)
registerDoParallel(cl)
ptm <- proc.time()
output <- foreach(i = 1:10000, .combine = "rbind",
                  .packages = c("dplyr", "ROCR",
                                "ResourceSelection")) %dopar%
  {
    # Seed
    print(i)
    set.seed(i)
    
    # Split data
    patients_train <- patients %>% 
      group_by(readmission) %>% 
      slice_sample(prop = 0.75)
    
    patients_validate <- patients %>% 
      filter(row_id %in% patients_train$row_id == FALSE)
    
    # Rebuild model
    final_model <- glm(readmission ~ acute_renal_failure_total +
                         age + atrial_fibrillation + days_before_ICU +
                         fluid_balance_5L + gcs + haemoglobin_final +
                         high_risk_speciality + serum_chloride_final,
                       data = patients_train,
                       family = "binomial")
    
    # Create predictions
    probs <- predict(final_model, newdata = patients_validate) %>% inverse_logit()
    
    # Assess discrimination
    pred <- prediction(probs[!is.na(probs)],
                       patients_validate$readmission[!is.na(probs)])
    AUC <- performance(pred, measure = "auc")@y.values[[1]]
    
    # Assess calibration
    cal <- hoslem.test(patients_validate$readmission[!is.na(probs)] == "Readmitted to ICU",
                       probs[!is.na(probs)], g = 10)
    
    # Output
    data.frame(i, disc = AUC,
               cal_chisq = cal$statistic,
               cal_p = cal$p.value)
  }
proc.time() - ptm 
stopCluster(cl)

# Summarise
output %>% 
  na.omit() %>% 
  summarise(AUC = mean(disc),
            AUC_error = sd(disc),
            cal = median(cal_chisq),
            cal_error_low = quantile(cal_chisq, 0.25),
            cal_error_high = quantile(cal_chisq, 0.75))

# AUC  AUC_error      cal cal_error_low cal_error_high
# 0.653 0.039          8.15     5.62      11.26

# Discrimination
set.seed(which.min(abs(output$disc - median(output$disc))))

# Split data
patients_train <- patients %>% 
  group_by(readmission) %>% 
  slice_sample(prop = 0.75)

patients_validate <- patients %>% 
  filter(subject_id %in% patients_train$subject_id == FALSE)

# Rebuild model
final_model <- glm(readmission ~ acute_renal_failure_total +
                     age + atrial_fibrillation + days_before_ICU +
                     fluid_balance_5L + gcs + haemoglobin_final +
                     high_risk_speciality + serum_chloride_final,
                    data = patients_train,
                    family = "binomial")
# Create predictions
probs <- predict(final_model, newdata = patients_validate) %>% inverse_logit()

# Assess discrimination
pred <- prediction(probs[!is.na(probs)],
                   patients_validate$readmission[!is.na(probs)])

# Calibration
set.seed(which.min(abs(output$cal_chisq - median(output$cal_chisq, na.rm = TRUE))))

# Split data
patients_train <- patients %>% 
  group_by(readmission) %>% 
  slice_sample(prop = 0.75)

patients_validate <- patients %>% 
  filter(subject_id %in% patients_train$subject_id == FALSE)

# Rebuild model
final_model <- glm(readmission ~ acute_renal_failure_total +
                     age + atrial_fibrillation + days_before_ICU +
                     fluid_balance_5L + gcs + haemoglobin_final +
                     high_risk_speciality + serum_chloride_final,
                    data = patients_train,
                    family = "binomial")

# Create predictions
probs <- predict(final_model, newdata = patients_validate) %>% inverse_logit()

# Split to deciles
decile_df <- tibble(
  patient_id = 1:nrow(patients_validate),
  readmission = patients_validate$readmission,
  probs,
  decile = ntile(probs, 10)) %>% 
  na.omit()

# Measure calibration
calibration_df <- decile_df %>% 
  select(readmission, probs, decile) %>% 
  group_by(decile) %>%
  summarise(N = length(readmission),
            observed = (sum(readmission == "Readmitted to ICU") / length(readmission)) * 100,
            predicted = mean(probs * 100),
            error = sd(probs * 100))

# Quantify calibration
cal_test <- hoslem.test(patients_validate$readmission[!is.na(probs)] == "Readmitted to ICU",
                        probs[!is.na(probs)], g = 10)

# Combine and write
list(model = "bespoke",
     data = "MIMIC",
     discrimination = pred,
     calibration = cal_test,
     deciles = calibration_df) %>% 
  write_rds("scripts/readmission/shared/models/MIMIC_bespoke.RDS")


