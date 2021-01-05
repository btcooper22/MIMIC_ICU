# Load packages
require(tidyverse)
require(magrittr)
require(lubridate)
require(doParallel)
require(tools)
options(dplyr.summarise.inform = FALSE)

# Load preprocessed data
mimic_preproc <- read_rds("data/mimic_preprocessed.RDS")

# Create outcome measures: Readmission within same hospitalisation-------------------

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
  
  # Count stays (126-128)
  stays %<>% 
    group_by(hadm_id) %>% 
    summarise(counts = n()) %>% 
    left_join(stays, by  = "hadm_id") %>% 
    arrange(row_id.x)
  
  # Generate output points for each stay (130-132)
  stays %<>% 
    group_by(hadm_id) %>% 
    mutate(max_outtime = max(outtime) == outtime)
  
  # Detect back transfers to ICU (134)
  stays %<>% 
    mutate(readmission = (counts > 1) & (max_outtime == 0)) %>% 
    # Restrict to surgical hospitalisations
    filter(surgical_stay == TRUE)
    
  # Write all information
  stays %>% 
    select(subject_id, hadm_id, in_hospital_mortality, readmission) %>% 
    mutate(subject_n = i) %>% 
    group_by(hadm_id) %>% 
    slice_head()
}
stopImplicitCluster()
proc.time() - ptm # 70 seconds @ 12 cores

# Clean and write results-------------

# How many patients died in hospital?
table(outcomes$in_hospital_mortality)

# Remove those who died
outcomes %<>% filter(in_hospital_mortality == FALSE)

# How many readmissions?
table(outcomes$readmission)

# Write
outcomes %>% 
  ungroup() %>% 
  select(subject_id, hadm_id,
         readmission) %>% 
  write_csv("data/outcomes.csv")
