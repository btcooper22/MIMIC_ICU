# Load packages
require(tidyverse)
require(magrittr)
require(lubridate)
require(doParallel)
require(tools)
require(tidyr)
options(dplyr.summarise.inform = FALSE)

# Load preprocessed data
mimic_preproc <- read_rds("data/mimic_preprocessed_missing.RDS")

# Create outcome measures: Mortality-------------------

# Prepare parallel options
ptm <- proc.time()
psnice(value = 19)
registerDoParallel(ifelse(detectCores() <= 15,
                          detectCores() - 1,
                          15)
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
  
  # Check for death
  if(!is.na(stays$dod[1]))
  {
    # Measure time between ICU discharge and death
    death_diff <- difftime(stays$dod[1], 
                           stays$outtime[nrow(stays)],
             unit = "days")
    
    # Flag if <30 days
    stays$mortality_30d <- death_diff < 30
  }else
  {
    stays$mortality_30d <- FALSE
  }
  
  # Write all information
  stays %>%
    select(subject_id, hadm_id, icustay_id,
           in_hospital_mortality,
           mortality_30d, in_unit_mortality) %>%
    mutate(subject_n = i) %>%
    group_by(hadm_id) %>%
    slice_head()
}
stopImplicitCluster()
proc.time() - ptm # 60 seconds @ 15 cores

# Clean and write results-------------

outcomes %<>% na.omit(outcomes)

# All patients who died in-unit flagged as within 30 days?
any(outcomes$in_unit_mortality & outcomes$mortality_30d == FALSE)

# All patients who died in ICU flagged as in hospital?
any(outcomes$in_unit_mortality & !outcomes$in_hospital_mortality)

# In-unit mortality rate
mean(outcomes$in_unit_mortality) * 100

# In-hospital mortality rate
mean(subset(outcomes, in_unit_mortality == FALSE)$in_hospital_mortality) * 100

# 30-day mortality rate
mean(subset(outcomes, in_unit_mortality == FALSE)$mortality_30d) * 100

# Write
outcomes %>% 
  ungroup() %>% 
  select(subject_id, hadm_id, in_unit_mortality,
         in_hospital_mortality, mortality_30d) %>% 
  write_csv("data/outcomes_mortality.csv")
