# Load packages
require(tidyverse)
require(rms)
require(brms)
require(magrittr)
require(lubridate)

# Load MIMIC dataset
source("functions/mimic_load.R")

# Find patients moved to surgical services
surgical_services <- mimic$services %>%
  filter(curr_service %LIKE% "%SURG%") %>% 
  select(hadm_id, subject_id, curr_service,
         transfertime)

# Add admission type
surgical_services %<>% 
  left_join(mimic$admissions %>% 
              select(hadm_id,
                     admission_type),
            by = "hadm_id")

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

# Remove surgery types with n < 20
results %>% 
  group_by(type) %>% 
  summarise(n = n()) %>% 
  filter(n >= 20) %>% 
  select(type) %>% 
  deframe() -> surg_types

results %<>% 
  filter(type %in% surg_types) %>% 
  filter(type != "Operating microscope (deleted code)")

# Write patient IDs
results %>% 
  write_csv("data/valid_patients.csv")
