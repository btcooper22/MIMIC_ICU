# Load packages
require(tidyverse)
require(magrittr)
require(lubridate)
require(doParallel)
require(tools)
options(dplyr.summarise.inform = FALSE)

# Load preprocessed data
mimic_preproc <- read_rds("data/mimic_preprocessed.RDS")

# Create outcome measures: Any subsequent nonelective readmission within 30 days-------------------

# Prepare parallel options
ptm <- proc.time()
psnice(value = 19)
registerDoParallel(ifelse(detectCores() <= 12,
                          detectCores() - 1,
                          12)
                   )

# Large loop for all patients
outcomes <- foreach(i = 1:nrow(mimic_preproc$index_stays),
                    .combine = "rbind",
                    .packages = c("dplyr", "magrittr",
                                  "tibble")) %dopar%
{
  # Identify subject
  # print(i)
  SUBJECT_ID <- mimic_preproc$index_stays$subject_id[i]
  
  # Extract stays 
  stays <- mimic_preproc$stays %>% 
    filter(subject_id == SUBJECT_ID) %>% 
    arrange(intime)
  
  # Find admission ID of surgical stay
  surg_adm <- mimic_preproc$index_stays$hadm_id[mimic_preproc$index_stays$subject_id == SUBJECT_ID]
  
  # Filter out subsequent elective admissions
  stays %<>% 
    filter(admission_type != "ELECTIVE" |
             surgical_hospitalisation == TRUE)

  # Find all stays after surgical admission
  if(max(which(stays$hadm_id == surg_adm)) < nrow(stays))
  {
    subseq_stays <- stays[(max(which(stays$hadm_id == surg_adm))+1):nrow(stays),]
    
   surg_time <- stays %>% 
      filter(surgical_hospitalisation == TRUE) %>% 
      select(intime) %>% deframe() %>% min()
   
   readmission <- difftime(min(subseq_stays$intime),
            surg_time) < 30
  }else
  {
    readmission <- FALSE
  }

  # Allow readmission within same hospitalisation
  if(readmission == FALSE & (nrow(subset(stays, surgical_hospitalisation == TRUE)) > 1))
  {
    readmission <- TRUE
  }

  # Collect output
  stays %>% 
    filter(hadm_id == surg_adm) %>% 
    slice_min(intime) %>% 
    select(subject_id, hadm_id, in_hospital_mortality) %>% 
    mutate(subject_n = i,
           readmission)
}
stopImplicitCluster()
proc.time() - ptm # 8 seconds @ 12 cores

# Clean and write results-------------

# REMOVE ORGAN DONORS
#WHERE lower(diagnosis) NOT LIKE '%organ donor%'

# How many patients died in hospital?
table(outcomes$in_hospital_mortality)

# Remove those who died
outcomes %<>% filter(in_hospital_mortality == FALSE)

# How many readmissions?
sum(outcomes$readmission)
mean(outcomes$readmission) * 100

# Write
outcomes %>% 
  ungroup() %>% 
  select(subject_id, hadm_id,
         readmission) %>% 
  write_csv("data/outcomes.csv")
