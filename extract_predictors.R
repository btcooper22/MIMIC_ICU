# Packages
require(tidyverse)
require(magrittr)
require(lubridate)

# Load MIMIC database
source("functions/mimic_load.R")

# Load preprocessed data
mimic_preproc <- read_rds("data/mimic_preprocessed.RDS")

# Load outcomes
outcomes <- read_csv("data/outcomes.csv")

# HAMMER------------
#     Sex
#     General Surgery
#     Cardiac Surgery
#     Hyperglycemia (>180 mg/dl in last 24h)
#     Anaemia (<7 mg/dl in last 24h)
#     APACHE II at admission (>20)
#     Positive fluid balance (>5L over ICU stay)
#     No ambulation in last 24h
#     ICU stay duration >5 days

# FROST------------
#---  Age
#     Sex
#     Elective admission
#     Source of admission (Operating theater/emergency/another hospital/ward)
#---  Apache II score at admission
#---  ICU stay duration >7 days
#     After-hours discharge (8:00-16:00)
#     Acute renal failure

# MARTIN------------
#     Respiratory rate (breaths/min)
#     Blood urea nitrogen (mg/dl)
#---  Age
#     Serum glucose (mg/dl)
#     Serum chloride (mmol/l)
#     History of atrial fibrillation
#     History of renal insufficiency

# Extraction loop----------

for(i in 1:100)
{
  # Identify subject
  # print(i)
  SUBJECT_ID <- mimic_preproc$patients[i]
  
  # Extract sex
  sex <- mimic$patients %>% 
    filter(subject_id == SUBJECT_ID) %>% 
    select(gender) %>% 
    collect %>% deframe()
}