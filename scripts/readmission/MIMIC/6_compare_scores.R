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

# Load data
patients <- read_csv("data/predictors.csv")

# Load functions
source("functions/NA_count.R")
source("functions/inverse_logit.R")
source("functions/nomogram_convert.R")
source("functions/calibration.R")
source("functions/brier_extraction.R")

# Filter missing data----------

# Count and exclude missing chart data
sum(patients$chart_missing)
patients %<>% filter(chart_missing == FALSE)

# Count and exclude missing lab data
sum(patients$lab_missing)
patients %<>% filter(lab_missing == FALSE)

# Count and exclude missing events data
sum(is.na(patients$fluid_balance_5L))
patients %<>% filter(!is.na(patients$fluid_balance_5L))

# Filter patients missing respiratory rate assessments
sum(is.na(patients$respiratory_rate) | patients$respiratory_rate > 48)
patients %<>% filter(!is.na(patients$respiratory_rate),
                     patients$rr_time == 48)

# Filter patients missing functional mobility assessments
sum(is.na(patients$ambulation))
patients %<>% filter(!is.na(patients$ambulation))

# Filter patients missing information to calculate APACHE-II score
sum(is.na(patients$apache_II))
patients %<>% filter(!is.na(patients$apache_II))

# Filter patients with no blood labs
sum(is.na(patients$hyperglycemia) | is.na(patients$anaemia))
patients %<>% filter(!is.na(patients$hyperglycemia) & !is.na(patients$anaemia))

# Assess NA
sapply(patients, function(x) sum(is.na(x)))

# Filter missing physiology on discharge
patients %<>% 
  filter(!is.na(temperature) &
           !is.na(mean_arterial_pressure) &
           !is.na(respiratory_rate.1) &
           !is.na(creatinine) &
           !is.na(wbc) & !is.na(gcs))


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

# Hammer----------

# Score data - coefficients
coefficients_hammer <- log(0.005 / (1 - 0.005)) +
  ifelse(patients$sex == "M", 0.43, 0) +
  ifelse(patients$general_surgery, 0.64, 0) +
  ifelse(patients$cardiac_surgery, -1.01, 0) +
  ifelse(patients$hyperglycemia, 0.36, 0) +
  ifelse(patients$anaemia, 1.11, 0) +
  ifelse(patients$high_apache, 0.44, 0) +
  ifelse(patients$fluid_balance_5L, 0.52, 0) +
  ifelse(patients$ambulation == FALSE, 0.7, 0) +
  ifelse(patients$los_5, 0.88, 0)

# Convert coefficients to probs
probs_coefficients_hammer <- inverse_logit(coefficients_hammer)
probs_hammer <- probs_coefficients_hammer

# Martin------------

# Score data - coefficients
coefficients_martin <- -9.284491 +
  (0.04883 * patients$respiratory_rate) +
  (0.011588 * patients$age) +
  (0.036104 * patients$serum_chloride) +
  (0.004967 * patients$blood_urea_nitrogen) +
  (0.580153 * patients$atrial_fibrillation) +
  (0.458202 * patients$renal_insufficiency) +
  (0.003519 * patients$serum_glucose)
  
# Convert to probabilities
probs_martin <- inverse_logit(coefficients_martin)

# Frost---------------

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
  ifelse(patients$acute_renal_failure, 10, 0)

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
probs_frost <- nomogram_convert(scores_frost, points_system_input,
                 points_system_output, log = TRUE)

# Bespoke model----

# Vector of acceptable features
feature_id <- c(7:13, 15:18, 20:33, 35:40, 44:48)

# Binarise readmission column
patients_trim <-  patients %>% 
  mutate(readmission = readmission == "Readmitted to ICU") %>% 
  select(all_of(feature_id)) %>% 
  na.omit()

# Select predictor and outcome vectors
x <- model.matrix(readmission~., patients_trim)[,-1]
y <- patients_trim$readmission

# Determine lambda
cv <- cv.glmnet(x, y, alpha = 1, folds = nrow(patients_trim),
                family = "binomial", tytype.measure = "auc")

# Fit model
initial_model <- glmnet(x, y, alpha = 1, lambda = quantile(cv$lambda, 0.85),
                        family = "binomial")

# Identify retained variables
initial_model$beta

# Loop through datasets
ptm <- proc.time()
output <- foreach(i = 1:10000, .combine = "rbind") %do%
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
    final_model <- glm(readmission ~ cardiac_surgery +
                         acute_renal_failure + length_of_stay + 
                         days_before_ICU + blood_urea_nitrogen,
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
final_model <- glm(readmission ~ cardiac_surgery +
                     acute_renal_failure + length_of_stay + 
                     days_before_ICU, data = patients_train,
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
final_model <- glm(readmission ~ cardiac_surgery +
                     acute_renal_failure + length_of_stay + 
                     days_before_ICU, data = patients_train,
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

# Discrimination----------

# Create prediction objects
prediction_hammer <- prediction(probs_hammer, patients$readmission)
prediction_martin <- prediction(probs_martin, patients$readmission)
prediction_frost <- prediction(probs_frost, patients$readmission)

# Create performance objects
performance_hammer <- performance(prediction_hammer, "tpr", "fpr")
performance_martin <- performance(prediction_martin, "tpr", "fpr")
performance_frost <- performance(prediction_frost, "tpr", "fpr")

# Create AUC objects
auc_hammer <- performance(prediction_hammer, measure = "auc")
auc_martin <- performance(prediction_martin, measure = "auc")
auc_frost <- performance(prediction_frost, measure = "auc")

# Print AUC
auc_hammer@y.values[[1]]
auc_martin@y.values[[1]]
auc_frost@y.values[[1]]

# Plot AUC
data.frame(x = performance_hammer@x.values[[1]],
           y = performance_hammer@y.values[[1]],
           model = "Hammer") %>% 
  rbind(
    data.frame(x = performance_martin@x.values[[1]],
               y = performance_martin@y.values[[1]],
               model = "Martin"),
    data.frame(x = performance_frost@x.values[[1]],
               y = performance_frost@y.values[[1]],
               model = "Frost")
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
  scale_color_brewer(palette = "Set1",
                     name = "")

# Calibration----------

# Split data into deciles
deciles_df <- tibble(
  patient_id = 1:nrow(patients),
  readmission = patients$readmission == "Readmitted to ICU",
  probs_hammer,
  decile_hammer = ntile(probs_hammer, 10),
  probs_martin,
  decile_martin = ntile(probs_martin, 10),
  probs_frost,
  decile_frost = ntile(probs_frost, 10)
)

# Calculate calibration
cal_hammer <- calibration(deciles_df, "probs_hammer", "decile_hammer") %>% 
  mutate(model = "hammer")
cal_martin <- calibration(deciles_df, "probs_martin", "decile_martin") %>% 
  mutate(model = "martin")
cal_frost <- calibration(deciles_df, "probs_frost", "decile_frost") %>% 
  mutate(model = "frost")

# Plot
rbind(cal_hammer %>% select(-decile_hammer), 
      cal_martin %>% select(-decile_martin),
      cal_frost %>% select(-decile_frost)) %>% 
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
  scale_colour_brewer(palette = "Set1",
                      name = "")+
  facet_wrap(~model, scales = "free_x")


# Calculate hosmer-lemeshow chi-squared
hoslem_hammer <- hoslem.test(patients$readmission == "Readmitted to ICU",
            probs_hammer, g = 10)
hoslem_martin <- hoslem.test(patients$readmission == "Readmitted to ICU",
            probs_martin, g = 10)
hoslem_frost <- hoslem.test(patients$readmission == "Readmitted to ICU",
            probs_frost, g = 10)
hoslem_hammer
hoslem_martin
hoslem_frost

# Write final patients file
write_csv(patients, "data/final_patients.csv")

# Combine and write----

list(model = "hammer",
     data = "MIMIC",
     discrimination = prediction_hammer,
     calibration = hoslem_hammer,
     deciles = cal_hammer) %>% 
  write_rds("scripts/readmission/shared/models/MIMIC_hammer.RDS")

list(model = "frost",
     data = "MIMIC",
     discrimination = prediction_frost,
     calibration = hoslem_frost,
     deciles = cal_frost) %>% 
  write_rds("scripts/readmission/shared/models/MIMIC_frost.RDS")


list(model = "martin",
     data = "MIMIC",
     discrimination = prediction_martin,
     calibration = hoslem_martin,
     deciles = cal_martin) %>% 
  write_rds("scripts/readmission/shared/models/MIMIC_martin.RDS")

