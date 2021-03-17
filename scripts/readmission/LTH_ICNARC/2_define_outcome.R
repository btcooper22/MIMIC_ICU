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

# Count crude readmission rate
count_admissions <- icnarc %>% 
  group_by(Identifiers_PatientPseudoId) %>% 
  summarise(entries = n()) %>% 
  mutate(mult = entries > 1) %>% 
  ungroup()

count_admissions %>% 
  summarise(readmission_rate = mean(mult) * 100)

# Identify readmitted patients
readmitted_ID <- count_admissions %>% 
  filter(mult == TRUE) %>% 
  select(Identifiers_PatientPseudoId) %>% 
  deframe()

readmitted_df <- icnarc %>% 
  filter(Identifiers_PatientPseudoId %in% readmitted_ID)


# Fix admission month
admit_month <- c()
for(j in 1:nrow(readmitted_df))
{
  dlength <- nchar(readmitted_df$MonthYearOfAdmission[j])
  year <- substr(readmitted_df$MonthYearOfAdmission[j], dlength-3, dlength)
  month <- substr(readmitted_df$MonthYearOfAdmission[j], 1, dlength-4)
  admit_month[j] <- paste(year, str_pad(month, 2, "left", "0"),
                                          sep = "") %>% as.numeric()
}
readmitted_df$admit_month <- admit_month

# Identify first stay
assessed_readmission <- foreach(i = 1:length(unique(readmitted_df$Identifiers_PatientPseudoId)),
                               .combine = "rbind") %do%
{
  # Identify patient
  patient <- unique(readmitted_df$Identifiers_PatientPseudoId)[i]
  patient_df <- readmitted_df %>% 
    filter(Identifiers_PatientPseudoId == patient)
  
  # Find index admission
  index_admission <- patient_df %>% 
    filter(Source_ClassificationOfSurgery == "4. Elective" &
             AdmissionType == "04. Planned local surgical admission")
  
  if(nrow(index_admission) > 1)
  {
    # Filter by first admission
    index_admission %<>% 
      slice_min(admit_month)
    
    patient_df$multiple_planned_elective <- TRUE
  }else
  {
    patient_df$multiple_planned_elective <- FALSE
  }
  
  # Check if multiple ICU stays in first admission
  if(nrow(index_admission) > 1)
  {
    patient_df$exclude_same_month <- TRUE
  }else
  {
    patient_df$exclude_same_month <- FALSE
  }
  
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

# Identify readmitted records
readmission_df <- assessed_readmission %>% 
  filter(index_admission == FALSE &
           exclude_same_month == FALSE &
           multiple_planned_elective == FALSE)

# Set aside single admissions
single_df <- icnarc %>% 
  filter(Identifiers_PatientPseudoId %in% readmitted_ID == FALSE) %>% 
  mutate(readmission = FALSE)

# Add first admissions of readmitted patients
first_readmitted <- assessed_readmission %>% 
  filter(index_admission == TRUE &
           exclude_same_month == FALSE &
           multiple_planned_elective == FALSE) %>% 
  mutate(readmission = TRUE)

results <- single_df %>% 
  rbind(first_readmitted %>% 
          select(-admit_month, -multiple_planned_elective,
                 -exclude_same_month, -index_admission,
                 -order_id))

# Count readmission rate
mean(results$readmission) * 100

# Count unique patients
length(unique(results$Identifiers_PatientPseudoId)) == nrow(results)

# Remove deaths or discharges with intention of death
results %<>% 
  filter(UnitDischarge_ReasonDischarged != "F. Palliative care",
         UnitDischarge_ExpectedDependencyPostDischarge != "E. Discharged with the expectation of dying",
         UnitDischarge_UnitOutcome != "4. Died")
# Write
results %>% 
  write_csv("data/icnarc_outcomes.csv")

#	Table on length of stay for single visit and readmission visit
tabling_df <- results %>% 
  rbind(
    readmission_df%>% 
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

