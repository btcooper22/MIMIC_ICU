# Packages
require(bigrquery)
require(DBI)
require(dplyr)
require(magrittr)
require(tibble)
require(purrr)
require(readr)

# Find project name
bq_projects()

# Find datasets
bq_project_datasets("prd-proj-decovid-358375")

# Find tables
bq_dataset_tables("prd-proj-decovid-358375.icnarc_analytics")

# Set ICNARC table as connection
con <- dbConnect(
  bigrquery::bigquery(),
  project = "prd-proj-decovid-358375",
  dataset = "icnarc_analytics"
)
dbListTables(con)

# Read IDs of surgical patients discharged under normal conditions
surgical_ID <- tbl(con, "ICNARC") %>% 
  collect() %>% 
  filter(Source_ClassificationOfSurgery == "4. Elective" &
           AdmissionType == "04. Planned local surgical admission") %>%
  select(Identifiers_PatientPseudoId) %>% 
  deframe() %>% unique()

# Load dataset of surgical patients
icnarc <- tbl(con, "ICNARC") %>% 
  collect() %>% 
  filter(Identifiers_PatientPseudoId %in% surgical_ID)

# Remove cases where no physiological records
sum(icnarc$Physiology_AllDataMissing)
icnarc %<>%
  filter(Physiology_AllDataMissing == FALSE)

# Filter by past medical history available
icnarc %<>%
  filter(PastMedicalHistory_AssessmentEvidenceAvailable == TRUE)

# Filter empty demographics
icnarc %<>% 
  filter(Demographics_UnitAdmissionAgeBand != "")

# Screen merged patients (ensure all patient IDs have same core demographics, at least)
icnarc %<>%
  group_by(Identifiers_PatientPseudoId) %>%
  summarise(ethnicity = length(unique(Demographics_Ethnicity)),
            sex = length(unique(Demographics_Sex))) %>%
  filter(ethnicity == 1 & sex == 1) %>%
  ungroup() %>%
  select(Identifiers_PatientPseudoId) %>%
  left_join(icnarc)

# Count readmission rate
count_admissions <- icnarc %>% 
  group_by(Identifiers_PatientPseudoId) %>% 
  summarise(entries = n()) %>% 
  mutate(mult = entries > 1) %>% 
  ungroup()

count_admissions %>% 
  summarise(readmission_rate = mean(mult) * 100)

# Write
icnarc %>% 
  write_csv("data/icnarc_surgical_cohort.csv")

# Multiple in same month - which is elective?
# If can't tell, probably exclude

test <- icnarc %>% 
  select(PriorToAdmission_HospitalAdmissionDaysBefore,
         MonthYearOfAdmission,
         UnitDischarge_DischargedDiedOnDays,
         UnitDischarge_ClinicallyReadyForDischargeDays,
         UltimateHospitalDischarge_DischargedDays,
         HospitalDischarge_DateDischarged)

icnarc %>%
  group_by(Identifiers_PatientPseudoId) %>% 
  group_by(Identifiers_PatientPseudoId) %>%
  summarise(ethnicity = length(unique(Demographics_Ethnicity)),
            sex = length(unique(Demographics_Sex))) %>%
  filter(ethnicity != 1 | sex != 1) %>% # age bands non-consecutive
  ungroup() %>%
  select(Identifiers_PatientPseudoId) %>%
  left_join(icnarc)
  summarise(n_months = length(unique(MonthYearOfAdmission)),
            n_id = length(Identifiers_PatientPseudoId))
