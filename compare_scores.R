# Packages
require(dplyr)
require(readr)
require(forcats)
require(magrittr)
require(lubridate)
require(ROCR)
require(ggplot2)

# Load data
patients <- read_csv("data/predictors.csv")

# Load functions
source("functions/NA_count.R")
source("functions/inverse_logit.R")
source("functions/nomogram_convert.R")

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

# Filter patients with no blood labs 24h before discharge
sum(is.na(patients$hyperglycemia) | is.na(patients$anaemia))
patients %<>% filter(!is.na(patients$hyperglycemia) & !is.na(patients$anaemia))

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

# Table 1-------------


# Crosstabulation function
crosstab <- function(varname, .df = patients)
{
  # Find
  .df  %>% 
    group_by_at(c("readmission", varname)) %>% 
    summarise(n = n()) %>% 
    mutate(`%` = (n/sum(n) * 100))
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
  mutate(age = cut(age, c(0, 25, 35, 45, 55,
                                   65, 75, 85, 89, 100),
                   c("18-25", "25-35", "35-45",
                     "45-55", "55-65", "65-75",
                     "75-85", "85-90", ">90"))) %>% 
  crosstab("age", .)
crosstab("sex")
crosstab("elective_admission")
crosstab("admission_source")
patients %>% 
  mutate(apache_II = cut(apache_II, seq(0,60,10),
                         right = FALSE)) %>% 
  crosstab("apache_II", .)
crosstab("los_7")
crosstab("after_hours_discharge")  
crosstab("acute_renal_failure")

### Martin variables
patients %>% 
  mutate(respiratory_rate = cut(respiratory_rate, seq(0, 90, 10),
                                right = FALSE)) %>% 
  crosstab("respiratory_rate", .)
patients %>% 
  mutate(blood_urea_nitrogen = cut(blood_urea_nitrogen, seq(0, 160, 40),
                                right = FALSE)) %>% 
  crosstab("blood_urea_nitrogen", .)
patients %>% 
  mutate(age = cut(age, seq(0, 100, 20), right = FALSE)) %>% 
  crosstab("age", .)
patients %>% 
  mutate(serum_glucose = cut(serum_glucose, seq(0, 650, 50),
                                   right = FALSE)) %>% 
  crosstab("serum_glucose", .)
patients %>% 
  mutate(serum_choride = cut(serum_choride, seq(80, 130, 10),
                             right = FALSE)) %>% 
  crosstab("serum_choride", .)
crosstab("atrial_fibrillation")
crosstab("renal_insufficiency")


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
coefficients_hammer <- log(0.005 / (1- 0.005)) +
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
cor(probs_coefficients_hammer, probs_points_hammer)
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

