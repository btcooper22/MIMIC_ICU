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
