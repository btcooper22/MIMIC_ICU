# Packages
require(readr)
require(dplyr)
require(tibble)
require(lubridate)
require(magrittr)
require(stringr)
require(foreach)
require(tidyr)

# Load cohort
icnarc <- read_csv("data/icnarc_surgical_cohort.csv") 
icnarc$row_id <- 1:nrow(icnarc)

# Fix admission month
admit_month <- c()
for(j in 1:nrow(icnarc))
{
  dlength <- nchar(icnarc$MonthYearOfAdmission[j])
  year <- substr(icnarc$MonthYearOfAdmission[j], dlength-3, dlength)
  month <- substr(icnarc$MonthYearOfAdmission[j], 1, dlength-4)
  admit_month[j] <- paste(year, str_pad(month, 2, "left", "0"),
                          sep = "") %>% as.numeric()
}
icnarc$admit_month <- as.Date(paste(admit_month, "01", sep = ""), format = "%Y%m%d")

# Count crude readmission rate
count_admissions <- icnarc %>% 
  group_by(Identifiers_PatientPseudoId) %>% 
  summarise(entries = n()) %>% 
  mutate(mult = entries > 1) %>% 
  ungroup()

count_admissions %>% 
  summarise(readmission_rate = mean(mult) * 100)

# Identify recurring patients
recurring_ID <- count_admissions %>% 
  filter(mult == TRUE) %>% 
  select(Identifiers_PatientPseudoId) %>% 
  deframe()

recurring_df <- icnarc %>% 
  filter(Identifiers_PatientPseudoId %in% recurring_ID)

isolated_df <- icnarc %>% 
  filter(Identifiers_PatientPseudoId %in% recurring_ID == FALSE)

# Confirm isolation
any(recurring_df$row_id %in% isolated_df$row_id)
any(isolated_df$row_id %in% recurring_df$row_id)


# Identify first stay for recurring patients
assessed_readmission <- foreach(i = 1:length(recurring_ID),
                               .combine = "rbind") %do%
{
  # Identify patient
  patient <- recurring_ID[i]
  patient_df <- recurring_df %>% 
    filter(Identifiers_PatientPseudoId == patient)
  
  # Find index admission
  index_admission <- patient_df %>% 
    filter(Source_ClassificationOfSurgery == "4. Elective" &
             AdmissionType == "04. Planned local surgical admission") %>% 
    slice_min(admit_month)
  
  if(nrow(index_admission) > 1)
  {
    patient_df$multiple_planned_elective <- TRUE
  }else
  {
    patient_df$multiple_planned_elective <- FALSE
  }
  
  patient_df$exclude_same_month <- FALSE
  
  # Mark index admission
  patient_df$index_admission <- patient_df$Identifiers_IcnarcPseudoId == index_admission$Identifiers_IcnarcPseudoId[1]
  
  # Extract index admission time
  index_time <- patient_df %>% 
    filter(index_admission == TRUE) %>% 
    select(admit_month) %>% deframe()
  
  # Remove cases before index admission
  patient_df %<>% filter(admit_month >= index_time)
  
  # Return
  patient_df$order_id <- i
  patient_df %>% 
    arrange(admit_month)
}

# Pull out records who were recurring but not "readmitted"
additional_isolated <- assessed_readmission %>% 
  group_by(Identifiers_PatientPseudoId) %>% 
  mutate(freq = n()) %>% 
  filter(freq == 1 & index_admission == TRUE)

# Add to isolated data frame
isolated_df <- additional_isolated %>% 
  select(-multiple_planned_elective,
         -exclude_same_month, -order_id,
         -freq, -index_admission) %>% 
  rbind(isolated_df) %>% 
  mutate(readmission = FALSE)

# Confirm no duplicates
any(duplicated(isolated_df$row_id))

# Identify all records of truly readmitted patients
readmitted_df <- assessed_readmission %>% 
  group_by(Identifiers_PatientPseudoId) %>% 
  mutate(freq = n()) %>% 
  filter(freq > 1 & multiple_planned_elective == FALSE
         & exclude_same_month == FALSE)

# Extract index admission of readmitted patients
results <- readmitted_df %>% 
  filter(index_admission == TRUE) %>% 
  select(-multiple_planned_elective,
         -exclude_same_month, -order_id,
         -freq, -index_admission) %>% 
  mutate(readmission = TRUE) %>% 
  # Add to total
  rbind(isolated_df)

# Remove deaths or discharges with intention of death
results %<>% 
  filter(UnitDischarge_ReasonDischarged != "F. Palliative care",
         UnitDischarge_ExpectedDependencyPostDischarge != "E. Discharged with the expectation of dying",
         UnitDischarge_UnitOutcome != "4. Died")

# Measure readmission timescales
readmitted_first <- results %>% 
  filter(readmission == TRUE)

readmission_delay <- c()
for(i in 1:nrow(readmitted_first))
{
  # Find readmission month
  readmit_time <- readmitted_df %>% 
    ungroup() %>% 
    filter(index_admission == FALSE &
             Identifiers_PatientPseudoId == readmitted_first$Identifiers_PatientPseudoId[i]) %>% 
    slice_min(admit_month) %>% 
    select(admit_month) %>% deframe() %>% unique()
  
  # Compare to admission month
  readmission_delay[i] <- difftime(readmitted_first$admit_month[i],
                                   readmit_time, unit = "days")
}

# Remove patients readmitted outside of a year
results$readmission[results$row_id %in% readmitted_first$row_id[abs(readmission_delay) > 32]] <- FALSE

# Confirm no duplicates
any(duplicated(results$row_id))

# Count unique patients
length(unique(results$Identifiers_PatientPseudoId)) == nrow(results)

# Count readmission rate
sum(results$readmission)
mean(results$readmission) * 100

# Write
results %>% 
  write_csv("data/icnarc_outcomes.csv")

#	Table on length of stay for single visit and readmission visit (needs rewrite)---- 
tabling_df <- results %>% 
  rbind(
    readmission_df %>% 
      select(-admit_month, -multiple_planned_elective,
             -exclude_same_month, -index_admission,
             -order_id) %>% 
      mutate(readmission = "more_visits")
  )

T1 <- tabling_df %>% 
  group_by(readmission) %>% 
  summarise(median = median(UnitDischarge_DischargedDiedOnDays),
            LQ = quantile(UnitDischarge_DischargedDiedOnDays, 0.25),
            UQ = quantile(UnitDischarge_DischargedDiedOnDays, 0.75),
            mean = mean(UnitDischarge_DischargedDiedOnDays),
            sd = sd(UnitDischarge_DischargedDiedOnDays)) %>% 
  mutate(name = "LengthOfStay")

# Total days support
T2 <- tabling_df %>% 
  mutate(total_support = OrganSupport_Ccmds_BasicRespiratoryDays +
                             OrganSupport_Ccmds_AdvancedRepiratoryDays +
                             OrganSupport_Ccmds_BasicCardiovascularDays +
                             OrganSupport_Ccmds_AdvancedCardiovascularDays +
                             OrganSupport_Ccmds_RenalDays +
                             OrganSupport_Ccmds_GastrointestinalDays +
                             OrganSupport_Ccmds_DermatologicalDays) %>% 
  group_by(readmission) %>% 
  summarise(median = median(total_support),
            LQ = quantile(total_support, 0.25),
            UQ = quantile(total_support, 0.75),
            mean = mean(total_support),
            sd = sd(total_support)) %>% 
  mutate(name = "TotalDaysSupport")


# Table on organ support for single visit and readmission visit
T3 <- tabling_df %>% 
  select(readmission, OrganSupport_Ccmds_BasicRespiratoryDays,
           OrganSupport_Ccmds_AdvancedRepiratoryDays,
           OrganSupport_Ccmds_BasicCardiovascularDays,
           OrganSupport_Ccmds_AdvancedCardiovascularDays,
           OrganSupport_Ccmds_RenalDays,
           OrganSupport_Ccmds_GastrointestinalDays,
           OrganSupport_Ccmds_DermatologicalDays) %>% 
  pivot_longer(2:8) %>% 
  group_by(readmission, name) %>% 
  summarise(median = median(value),
            LQ = quantile(value, 0.25),
            UQ = quantile(value, 0.75),
            mean = mean(value),
            sd = sd(value))

# Combine
rbind(T1, T2, T3) %>% 
  mutate(value = paste(mean %>% round(2), " +/-", 
                       sd %>% round(2), sep = "")) %>% 
  select(readmission, name, value) %>% 
  pivot_wider(names_from = "readmission",
              values_from = "value") %>% 
  write_csv("output/readmission_tables.csv")



