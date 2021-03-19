# Adapted into R from https://github.com/Jeffreylin0925/MIMIC-III_ICU_Readmission_Analysis
# and https://github.com/YerevaNN/mimic3-benchmarks

# Load packages
require(tidyverse)
require(magrittr)
require(lubridate)
require(doParallel)
require(tools)

# Load MIMIC database
source("functions/mimic_load.R")

# Functions
source("functions/in_hospital_mortality.R")
source("functions/in_unit_mortality.R")

# Load valid patients
patient_df <- read_csv("data/valid_patients.csv") 
valid_patients <- patient_df %>% 
  select(subject_id) %>%  deframe()
surgical_hospitalisations <- patient_df %>% 
  select(hadm_id) %>%  deframe()


# Read and filter stay data-------------------------------------------------

# Merge and filter transfers (40-41)
transfers <- mimic$transfers %>% 
  inner_join(mimic$admissions,
             by = c("subject_id", "hadm_id")) %>% 
  inner_join(mimic$patients, by = "subject_id") %>% 
  filter(subject_id %in% valid_patients) %>% 
  collect()

# Merge and filter stays (43-44)
stays <- mimic$icustays %>% 
  inner_join(mimic$admissions,
             by = c("subject_id", "hadm_id")) %>% 
  inner_join(mimic$patients, by = "subject_id")%>% 
  filter(subject_id %in% valid_patients) %>% 
  collect()

# Add time events to transfers (46-49)
transfers %<>% 
  mutate(age = (intime - dob) / 364.25) %>% 
  mutate(age = ifelse(age > 90, 90, age)) %>% 
  mutate(in_unit_mortality = in_unit_mortality(.),
         in_hospital_mortality = in_hospital_mortality(.)) %>% 
  filter(age >= 18)

# Add time events to stays (54-57)
stays %<>% 
  mutate(age = (intime - dob) / 364.25) %>% 
  mutate(age = ifelse(age > 90, 90, age)) %>% 
  mutate(in_unit_mortality = in_unit_mortality(.),
         in_hospital_mortality = in_hospital_mortality(.),
         surgical_hospitalisation = ifelse(hadm_id %in% surgical_hospitalisations,
                                TRUE, FALSE)) %>% 
  filter(age >= 18)

# Find hospitalisations involving a surgical service
surgical_hosps_df <- stays %>% 
  filter(surgical_hospitalisation == TRUE)

# Flag index stays
index_stays <- foreach(i = 1:nrow(patient_df),
        .combine = "rbind",
        .packages = c("dplyr", "lubridate",
                      "tibble")) %do%
{
  # Extract all stays for patient
  patient_stays <- stays %>% 
    filter(subject_id == patient_df$subject_id[i]) %>% 
    arrange(intime)
  
  # Skip patients with no surgical stays
  if(any(patient_stays$surgical_hospitalisation))
  {
  # Take only elective surgical stays
  index_stay <- patient_stays %>% 
    filter(admission_type == "ELECTIVE" &
             surgical_hospitalisation == TRUE)
  
  # Filter multiple elective procedures in one hospitalisation
  if(nrow(index_stay) > 1 & length(unique(index_stay$hadm_id)) == 1)
  {
    # Take latest elective procedure
    index_stay <- index_stay %>% 
      filter(admission_type == "ELECTIVE") %>% 
      slice_tail(n = 1)
  }

  
  if(nrow(index_stay) > 1)
  {
    stop(i)
    # # Find time moved to surgical service
    # surgery_time <- patient_df %>%
    #   filter(hadm_id == unique_admissions[i]) %>% 
    #   select(transfertime) %>% deframe()
    # 
    # # Compare to ICU times
    # post_surgery <- difftime(admission_df$outtime, surgery_time,
    #                             units = "hours")
    # 
    # # Negative times == icu stays pre-surgery
    # if(any(post_surgery < 0))
    # {
    #   admission_df <- admission_df[post_surgery > 0,]
    # }
  }
  index_stay
  }else
  {
    output <- stays[i,]
    output[,1] <- NA
    output
  }
}

# Read and filter input data per stay----------------------------------

# Read and translate diagnoses (68)
diagnoses <- mimic$d_icd_diagnoses %>% 
  select(2:4) %>% 
  inner_join(mimic$diagnoses_icd, by = "icd9_code") %>% 
  collect()

# Filter diagnoses (69)
diagnoses %<>% 
  inner_join(stays %>% select(subject_id, hadm_id, icustay_id),
             by = c("subject_id", "hadm_id"))

# Count diagnosis codes (72)
diagnosis_counts <- diagnoses %>% 
  select(icd9_code, hadm_id) %>% 
  unique() %>% 
  group_by(icd9_code) %>% 
  summarise(count = n()) %>% 
  left_join(mimic$d_icd_diagnoses %>% 
              select(2:4) %>% 
              collect())

# Read and translate procedures (76)
procedures <- mimic$d_icd_procedures %>% 
  select(2:4) %>% 
  inner_join(mimic$procedures_icd, by = "icd9_code") %>% 
  collect()

# Filter procedures (77)
procedures %<>% 
  inner_join(stays %>% select(subject_id, hadm_id, icustay_id),
             by = c("subject_id", "hadm_id"))

# Count procedure codes (80)
procedure_counts <- procedures %>% 
  select(icd9_code, hadm_id) %>% 
  unique() %>% 
  group_by(icd9_code) %>% 
  summarise(count = n()) %>% 
  left_join(mimic$d_icd_procedures %>% 
              select(2:4) %>% 
              collect())

# Prepare data for writing-----------------------

# Create blank list
mimic_preproc <- list()

# Insert stay data
mimic_preproc[[1]] <- unique(stays$subject_id)
mimic_preproc[[2]] <- transfers
mimic_preproc[[3]] <- stays
names(mimic_preproc)[1:3] <- c("patients", "transfers", 
                               "stays")

# Insert input data
mimic_preproc[[4]] <- diagnoses
mimic_preproc[[5]] <- procedures
names(mimic_preproc)[4:5] <- c("diagnoses", "procedures")

# Insert counts
mimic_preproc[[6]] <- list("procedure" = procedure_counts,
                           "diagnosis" = diagnosis_counts)
names(mimic_preproc)[6] <- "counts"

# Insert index stays
mimic_preproc[[7]] <- index_stays %>% 
  filter(!is.na(row_id.x))
names(mimic_preproc)[7] <- "index_stays"

# Write data
write_rds(mimic_preproc, "data/mimic_preprocessed.RDS")
