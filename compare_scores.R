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
sum(is.na(patients$apache_II) | is.na(patients$apache_II_discharge))
patients %<>% filter(!is.na(patients$apache_II) & !is.na(patients$apache_II_discharge))

# Filter patients with no blood labs 24h before discharge
sum(is.na(patients$hyperglycemia) | is.na(patients$anaemia))
patients %<>% filter(!is.na(patients$hyperglycemia) & !is.na(patients$anaemia))

# Filter patients missing Fialho criteria
sum(is.na(patients$final_bp) | is.na(patients$final_lactate) |
      is.na(patients$final_platelets) | is.na(patients$final_pO2) |
      is.na(patients$final_pulse) | is.na(patients$final_temp))
patients %<>% filter(!is.na(final_bp) & !is.na(final_lactate) &
                       !is.na(final_platelets) & !is.na(final_pO2) &
                       !is.na(final_pulse) & !is.na(final_temp))

# Fix Fialho variables
patients %<>% 
  # Fix lactate
  mutate(final_lactate = final_lactate / 9) %>% 
  # Rename SpO2
  rename("final_SpO2" = final_pO2)

# Filter patients with physiologically implausible values
patients %<>% 
  mutate(physiologically_implausible = final_temp < 35 |
           final_temp > 40 | final_SpO2 > 100 | final_bp < 22 |
           final_bp > 136 | final_platelets > 1043)
sum(patients$physiologically_implausible)
patients %<>% filter(physiologically_implausible == FALSE)

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

# Table 1-------------

# Crosstabulation function
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
decile_cut <- function(in_var)
{
  quantile(in_var, probs = seq(0, 1, 0.1)) %>% 
    round()
}

### Hammer variables
crosstab("sex")
crosstab("general_surgery")
crosstab("cardiac_surgery")
crosstab("hyperglycemia")
crosstab("anaemia")
crosstab("high_apache")
crosstab("fluid_balance_5L")
crosstab("ambulation")
crosstab("los_5")

### Frost variables
patients %>%       
  mutate(age = cut(age, decile_cut(age))) %>% 
  crosstab("age", .)
crosstab("elective_admission")
crosstab("admission_source")
patients %>% 
  mutate(apache_II_discharge = cut(apache_II_discharge, decile_cut(apache_II_discharge),
                         right = FALSE)) %>% 
  crosstab("apache_II_discharge", .)
crosstab("los_7")
crosstab("after_hours_discharge")  
crosstab("acute_renal_failure")

### Martin variables
patients %>% 
  mutate(respiratory_rate = cut(respiratory_rate, decile_cut(respiratory_rate),
                                right = FALSE)) %>% 
  crosstab("respiratory_rate", .)
patients %>% 
  mutate(blood_urea_nitrogen = cut(blood_urea_nitrogen, 
                                   decile_cut(blood_urea_nitrogen),
                                right = FALSE)) %>% 
  crosstab("blood_urea_nitrogen", .)
patients %>% 
  mutate(serum_glucose = cut(serum_glucose, decile_cut(serum_glucose),
                                   right = FALSE)) %>% 
  crosstab("serum_glucose", .)
patients %>% 
  mutate(serum_choride = cut(serum_choride, decile_cut(serum_choride),
                             right = FALSE)) %>% 
  crosstab("serum_choride", .)
crosstab("atrial_fibrillation")
crosstab("renal_insufficiency")

### Fialho variables
patients %>% 
  mutate(final_pulse = cut(final_pulse, decile_cut(final_pulse),
                                right = FALSE)) %>% 
  crosstab("final_pulse", .)
patients %>% 
  mutate(final_temp = cut(final_temp, c(30,36, 36.5, 37,
                                        37.5, 38,40),
                           right = FALSE)) %>% 
  crosstab("final_temp", .)
patients %>% 
  mutate(final_SpO2 = cut(final_SpO2, c(75,90, 93:100),
                          right = TRUE)) %>% 
  crosstab("final_SpO2", .)
patients %>% 
  mutate(final_bp = cut(final_bp, c(20,60, 70, 80, 85, 90,
                                    100, 140),
                          right = TRUE)) %>% 
  crosstab("final_bp", .)
patients %>% 
  mutate(final_platelets = cut(final_platelets, decile_cut(final_platelets),
                           right = FALSE)) %>% 
  crosstab("final_platelets", .)
patients %>% 
  mutate(final_lactate = cut(final_lactate, c(0, 0.5, 1, 1.25, 1.5, 2, 
                                              2.5, 4, 8),
                           right = FALSE)) %>% 
  crosstab("final_lactate", .)

# Hammer----------

# Score data - points
points_hammer <- ifelse(patients$sex == "M", 1, 0) +
  ifelse(patients$general_surgery, 2, 0) +
  ifelse(patients$cardiac_surgery, -3, 0) +
  ifelse(patients$hyperglycemia, 1, 0) +
  ifelse(patients$anaemia, 3, 0) +
  ifelse(patients$high_apache, 1, 0) +
  ifelse(patients$fluid_balance_5L, 1, 0) +
  ifelse(patients$ambulation == FALSE, 2, 0) +
  ifelse(patients$los_5, 2, 0)
table(points_hammer)

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

# Convert scores to probs
scale_hammer <- data.frame(
  points = -3:12,
  probs = c(0.1, 0.2, 0.35, 0.5, 0.75, 1, 1.5, 2, 3.1,
            4.2, 6.4, 8.6, 12.55, 16.5, 23.05, 29.6)
)
probs_points_hammer <- scale_hammer$probs[match(points_hammer, scale_hammer$points)] / 100

# Correlate points with coefficients
#cor(probs_coefficients_hammer, probs_points_hammer)
probs_hammer <- probs_coefficients_hammer

# Martin------------

# Score data - coefficients
coefficients_martin <- -9.284491 +
  (0.04883 * patients$respiratory_rate) +
  (0.011588 * patients$age) +
  (0.036104 * patients$serum_choride) +
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
  nomogram_convert(patients$apache_II_discharge, apache_score_system_input,
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

# Fialho-----------

# Discrimination----------

# Create prediction objects
prediction_hammer <- prediction(probs_hammer, patients$readmission == "Readmitted to ICU")
prediction_martin <- prediction(probs_martin, patients$readmission == "Readmitted to ICU")
prediction_frost <- prediction(probs_frost, patients$readmission == "Readmitted to ICU")
prediction_apache <- prediction(patients$apache_II_discharge,
                                patients$readmission == "Readmitted to ICU")

# Create performance objects
performance_hammer <- performance(prediction_hammer, "tpr", "fpr")
performance_martin <- performance(prediction_martin, "tpr", "fpr")
performance_frost <- performance(prediction_frost, "tpr", "fpr")
performance_apache <- performance(prediction_apache, "tpr", "fpr")

# Create AUC objects
auc_hammer <- performance(prediction_hammer, measure = "auc")
auc_martin <- performance(prediction_martin, measure = "auc")
auc_frost <- performance(prediction_frost, measure = "auc")
auc_apache <- performance(prediction_apache, measure = "auc")

# Print AUC
auc_hammer@y.values[[1]]
auc_martin@y.values[[1]]
auc_frost@y.values[[1]]
auc_apache@y.values[[1]]

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
               model = "Frost"),
    data.frame(x = performance_apache@x.values[[1]],
               y = performance_apache@y.values[[1]],
               model = "APACHE-II")
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
  decile_frost = ntile(probs_frost, 10),
  probs_apache = patients$apache_II_discharge / 200,
  decile_apache = ntile(patients$apache_II_discharge, 10)
)

# Calculate calibration
cal_hammer <- calibration(deciles_df, "probs_hammer", "decile_hammer") %>% 
  mutate(model = "hammer")
cal_martin <- calibration(deciles_df, "probs_martin", "decile_martin") %>% 
  mutate(model = "martin")
cal_frost <- calibration(deciles_df, "probs_frost", "decile_frost") %>% 
  mutate(model = "frost")
cal_apache <- calibration(deciles_df, "probs_apache", "decile_apache") %>% 
  mutate(model = "apache")

# Plot
rbind(cal_hammer %>% select(-decile_hammer), 
      cal_martin %>% select(-decile_martin),
      cal_frost %>% select(-decile_frost), 
      cal_apache %>% select(-decile_apache)) %>% 
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
  facet_wrap(~model, scales = "fixed")


# Calculate hosmer-lemeshow chi-squared
hoslem.test(patients$readmission == "Readmitted to ICU",
            probs_hammer, g = 10)
hoslem.test(patients$readmission == "Readmitted to ICU",
            probs_martin, g = 10)
hoslem.test(patients$readmission == "Readmitted to ICU",
            probs_frost, g = 10)
hoslem.test(patients$readmission == "Readmitted to ICU",
            patients$apache_II_discharge / 200, g = 10)

# Calculate brier scores
brier_df <- rbind(
  brier_extraction(patients$readmission == "Readmitted to ICU", probs_hammer) %>% 
    mutate(model = "hammer"),
  brier_extraction(patients$readmission == "Readmitted to ICU", probs_martin) %>% 
    mutate(model = "martin"),  
  brier_extraction(patients$readmission == "Readmitted to ICU", probs_frost) %>% 
    mutate(model = "frost"),
  brier_extraction(patients$readmission == "Readmitted to ICU", 
                   patients$apache_II_discharge / 200) %>% 
    mutate(model = "apache")
)

# Plot brier scores
brier_df %>% 
  pivot_longer(1:3) %>% 
  ggplot(aes(x = model, y = value))+
  geom_bar(stat = "identity", colour = "black",
           aes(fill = model))+
  facet_wrap(~name, scales = "free_y")+
  theme_classic(20)+
  theme(legend.position = "none")+
  labs(x = "", y = "")+
  scale_fill_brewer(palette = "Set1",
                      name = "")

# Recalibration-----------

# Split data
set.seed(123)
patients_train <- patients %>% 
  group_by(readmission) %>% 
  slice_sample(prop = 0.5) %>% 
  mutate(type = "train")

patients_validate <- patients %>% 
  filter(row_id %in% patients_train$row_id == FALSE)

### Recalibrate models

# Hammer
hammer_recalibrate <- glm(readmission ~ sex + general_surgery + cardiac_surgery + hyperglycemia +
                            high_apache + fluid_balance_5L + ambulation + los_5,
                          data = patients_train, family = "binomial")

# Martin
martin_recalibrate <- glm(readmission ~ respiratory_rate + age + serum_choride + blood_urea_nitrogen +
                            atrial_fibrillation + renal_insufficiency + serum_glucose,
                          data = patients_train, family = "binomial")

# Frost
frost_recalibrate <- glm(readmission ~ age + sex + elective_admission + admission_source +
                            apache_II_discharge + los_7 + after_hours_discharge + acute_renal_failure,
                          data = patients_train, family = "binomial")

# APACHE
apache_recalibrate <- glm(readmission ~ apache_II_discharge,
                          data = patients_train, family = "binomial")

### Fit novel model automatically

# Build formula
formu <- paste(names(patients_train)[c(9:10, 12, 14,
                                       16:18, 20:35, 51)], collapse = " + ")

# Build model
dirty_model <- glm(paste("readmission", "~", formu, sep = " "),
                   data = patients_train, family = "binomial")


if(file.exists("data/cooper_model.RDS"))
{
  cooper_model <- readRDS("data/cooper_model.RDS")
}else
{
  # Automated model selection
  cooper_model <- step(dirty_model)
  write_rds(cooper_model, "data/cooper_model.RDS")
}

### Predict on new data
hammer_rc_probs <- predict(hammer_recalibrate, newdata = patients_validate) %>% inverse_logit()
martin_rc_probs <- predict(martin_recalibrate, newdata = patients_validate) %>% inverse_logit()
frost_rc_probs <- predict(frost_recalibrate, newdata = patients_validate) %>% inverse_logit()
apache_rc_probs <- predict(apache_recalibrate, newdata = patients_validate) %>% inverse_logit()
cooper_probs <- predict(cooper_model, newdata = patients_validate) %>% inverse_logit()

### Assess discrimination

# Create prediction objects
prediction_rc_hammer <- prediction(hammer_rc_probs, patients_validate$readmission)
prediction_rc_hammer <- prediction(hammer_rc_probs, patients_validate$readmission)

# Create performance objects
performance_rc_hammer <- performance(prediction_rc_hammer, "tpr", "fpr")

# Create AUC objects
auc_rc_hammer <- performance(prediction_rc_hammer, measure = "auc")

# Print AUC
auc_rc_hammer@y.values[[1]]


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
               model = "Frost"),
    data.frame(x = performance_apache@x.values[[1]],
               y = performance_apache@y.values[[1]],
               model = "APACHE-II")
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



