# Load packages
require(tidyverse)
require(magrittr)

# Initialise renv
# require(renv)
# renv::init()

# Load MIMIC dataset
source("functions/mimic_load.R")

# Find patients admitted under surgical services
surgical_services <- mimic$services %>%
  filter(curr_service %LIKE% "%SURG%" &
           is.na(prev_service)) %>% 
  select(hadm_id, subject_id, curr_service,
         transfertime)

# Add admission type
surgical_services %<>% 
  left_join(mimic$admissions %>% 
              select(hadm_id,
                     admission_type),
            by = "hadm_id") %>% 
  filter(admission_type == "ELECTIVE")

# Add surgery type
surgical_services %<>% 
  left_join(mimic$cptevents %>% 
              select(hadm_id, sectionheader,
                     subsectionheader, ticket_id_seq),
            by = "hadm_id") %>% 
  filter(sectionheader == "Surgery") %>% 
  select(-sectionheader)

# Run query
show_query(surgical_services)
surgical_services %<>% collect()

# Find number of surgeries
n_surg <- surgical_services %>% 
  group_by(subject_id) %>% 
  summarise(subseq = n() - 1)

# Clean and filter records
results <- surgical_services %>% 
  # Add additional surgeries
  left_join(n_surg, by = "subject_id") %>%
  # Select only first surgery
  group_by(subject_id) %>% 
  slice(which.min(ticket_id_seq)) %>% 
  ungroup()  %>% 
  rename(type = subsectionheader)

results %<>% 
  filter(type != "Operating microscope (deleted code)",
         type != "Laparoscopy, Surgical; Cholecystectomy (deleted code)")

# Write patient IDs
results %>% 
  write_csv("data/valid_patients.csv")
