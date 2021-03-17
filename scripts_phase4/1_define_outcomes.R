# Load packages
require(tidyverse)
require(magrittr)
require(lubridate)
require(doParallel)
require(tools)
options(dplyr.summarise.inform = FALSE)

# Load preprocessed data -> Make this into own script from scratch
mimic_preproc <- read_rds("data/mimic_preprocessed.RDS")

# Create outcome measures: Mortality-------------------

# Prepare parallel options
ptm <- proc.time()
psnice(value = 19)
registerDoParallel(ifelse(detectCores() <= 12,
                          detectCores() - 1,
                          12)
                   )

# Large loop for all patients
outcomes <- foreach(i = 1:length(mimic_preproc$patients),
                    .combine = "rbind",
                    .packages = c("dplyr", "magrittr")) %dopar%
{
  # Identify subject
  # print(i)
  SUBJECT_ID <- mimic_preproc$patients[i]
  
  # Extract stays (up to 124)
  stays <- mimic_preproc$stays %>% 
    filter(subject_id == SUBJECT_ID)
  
  # Identify surgical admission
  surgical_stay <- stays %>%
    group_by(hadm_id) %>%
    filter(any(surgical_hospitalisation == TRUE)) %>% 
    ungroup()
  
  # Flag readmission
  readmission_flag <- nrow(surgical_stay) > 1
  surgical_stay %<>% 
    mutate(readmission = readmission_flag)
  
  # Check for death
  if(!is.na(surgical_stay$dod[1]))
  {
    # Measure time between ICU discharge and death
    death_diff <- difftime(surgical_stay$dod[1], 
             surgical_stay$outtime[nrow(surgical_stay)],
             unit = "days")
    
    # Flag if <30 days
    surgical_stay$mortality_30d <- death_diff < 30
  }else
  {
    surgical_stay$mortality_30d <- FALSE
  }
  
  # # Write all information
  surgical_stay %>%
    select(subject_id, hadm_id, in_hospital_mortality,
           mortality_30d, readmission, in_unit_mortality) %>%
    mutate(subject_n = i) %>%
    group_by(hadm_id) %>%
    slice_head()
}
stopImplicitCluster()
proc.time() - ptm # 60 seconds @ 12 cores

# Clean and write results-------------

# REMOVE ORGAN DONORS
#WHERE lower(diagnosis) NOT LIKE '%organ donor%'

# All patients who died in hospital flagged as within 30 days?
any(outcomes$in_hospital_mortality & !outcomes$mortality_30d)

# All patients who died in ICU flagged as in hospital?
any(outcomes$in_unit_mortality & !outcomes$in_hospital_mortality)

# How many patients died in ICU?
table(outcomes$in_unit_mortality)

# Filter out in-unit mortality
# outcomes %<>% 
#   filter(in_unit_mortality == FALSE)

# How many patients died in hospital?
table(outcomes$in_hospital_mortality)

# How many died within 30 days
table(outcomes$mortality_30d)

# Write
outcomes %>% 
  ungroup() %>% 
  select(subject_id, hadm_id, readmission, in_unit_mortality,
         in_hospital_mortality, mortality_30d) %>% 
  write_csv("data/outcomes_mortality.csv")
