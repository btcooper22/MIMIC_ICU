# Load packages
require(tidyverse)
require(magrittr)
require(doParallel)
require(tools)

# Functions
source("functions/in_hospital_mortality.R")
source("functions/in_unit_mortality.R")

# Load MIMIC dataset
source("functions/mimic_load.R")

# Find patients moved to surgical services
results <- mimic$services %>%
  filter(curr_service != "NB" &
           curr_service != "NBB") %>% 
  select(hadm_id, subject_id, curr_service,
         transfertime)

# Add admission type
results %<>% 
  left_join(mimic$admissions %>% 
              select(hadm_id, admittime,
                     admission_type),
            by = "hadm_id")

# Run query
show_query(results)
results %<>% collect()

# Select first admission for each patient
results %<>% 
  arrange(admittime) %>% 
  filter(!duplicated(subject_id))

# Load valid patients
valid_patients <- results %>% 
  select(subject_id) %>%  deframe()
valid_hospitalisations <- results %>% 
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
         in_hospital_mortality = in_hospital_mortality(.)) %>% 
  filter(age >= 18)

# Filter to first stay for each admission
stays_filtered <- stays %>% 
  arrange(intime) %>% 
  filter(!duplicated(subject_id)) %>% 
  filter(!grepl("organ donor", diagnosis, ignore.case = TRUE))

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
  inner_join(stays_filtered %>% select(subject_id, hadm_id, icustay_id),
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
mimic_preproc[[1]] <- unique(stays_filtered$subject_id)
mimic_preproc[[2]] <- transfers
mimic_preproc[[3]] <- stays_filtered
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

# Write data
write_rds(mimic_preproc, "data/mimic_preprocessed_missing.RDS")
