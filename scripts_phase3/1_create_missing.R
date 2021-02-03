# Packages
require(dplyr)
require(magrittr)
require(readr)
require(tibble)

# Load data and trim to APACHE-II
full_data <- read_csv("data/final_patients.csv") %>% 
  select(matches("_discharge$")) %>% 
  select(matches("^apache_")) %>% 
  cbind(read_csv("data/final_patients.csv") %>% 
          select(row_id, subject_id,
                 adm_id, readmission),.)
