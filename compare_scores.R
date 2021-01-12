# Packages
require(dplyr)
require(readr)
require(forcats)
require(magrittr)
require(lubridate)

# Load data
patients <- read_csv("data/predictors.csv")

# Load functions
source("functions/NA_count.R")
source("functions/inverse_logit.R")

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
coefficients_hammer <- -2.944439 +
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
plot(probs_coefficients_hammer, probs_points_hammer)
cor(probs_coefficients_hammer, probs_points_hammer)
