# Packages
require(tidyverse)
require(magrittr)
require(lubridate)
require(doParallel)
require(tools)

# Load MIMIC database
source("functions/mimic_load.R")
source("functions/event_filter_discharge.R")

# Load preprocessed data
mimic_preproc <- read_rds("data/mimic_preprocessed.RDS")

# Load outcomes
outcomes <- read_csv("data/outcomes.csv")

# Item definitions---------

ID_blood_glucose <- mimic$d_labitems %>% 
  filter(grepl("Glucose", label),
         fluid == "Blood") %>% 
  select(itemid) %>% 
  collect() %>% deframe()

ID_haemoglobin <- mimic$d_labitems %>% 
  filter(label == "Hemoglobin") %>% 
  select(itemid) %>% 
  collect() %>% deframe()

ID_acute_renal_failure <- mimic_preproc$diagnoses %>% 
  filter(grepl("Acute kidney failure", long_title)) %>% 
  select(icd9_code) %>% 
  deframe() %>% unique()

ID_respiratory_rate <- mimic$d_items %>% 
  filter(grepl("resp", label,
               ignore.case = TRUE)) %>% 
  filter(grepl("rate", label,
               ignore.case = TRUE)) %>% 
  select(itemid) %>% 
  collect() %>% deframe()

ID_blood_urea_nitrogen <-  mimic$d_labitems %>% 
  filter(grepl("urea nitrogen", label,
               ignore.case = TRUE),
         fluid == "Blood") %>% 
  select(itemid) %>% 
  collect() %>% deframe()

ID_serum_choride <-  mimic$d_labitems %>% 
  filter(grepl("chloride", label,
               ignore.case = TRUE),
               fluid == "Blood") %>% 
  select(itemid) %>% 
  collect() %>% deframe()

ID_atrial_fibrillation <- mimic$d_icd_diagnoses %>% 
  filter(grepl("atrial fibrillation", long_title,
               ignore.case = TRUE))%>% 
  select(icd9_code) %>% 
  collect() %>% deframe()

ID_renal_insufficiency <- mimic_preproc$diagnoses %>% 
  filter(grepl("renal", long_title,
               ignore.case = TRUE)) %>% 
  select(icd9_code, long_title) %>% 
  unique() %>% 
  mutate(row_id = 1:n())

ID_renal_insufficiency %<>% 
  filter(row_id %in% c(1:3, 10, 14, 19, 24:25,
                       29, 35)) %>% 
  select(icd9_code) %>% deframe()

ID_ambulation <- mimic$d_items %>% 
  filter(grepl("braden activity", label,
               ignore.case = TRUE)) %>% 
  collect() %>% 
  select(itemid) %>%  deframe()

# Extraction loop----------

# Load smaller database elements
mimic_icustays <- mimic$icustays %>% collect()
mimic_patients <- mimic$patients %>% collect()
mimic_services <- mimic$services %>% collect()
mimic_admissions <- mimic$admissions %>% collect()
mimic_diagnoses <- mimic$diagnoses_icd %>% collect()

# Prepare parallel options
ptm <- proc.time()
psnice(value = 19)
registerDoParallel(ifelse(detectCores() <= 15,
                          detectCores() - 1,
                          12)
)

# Run loop
predictors <- foreach(i = 1:nrow(outcomes), .combine = "rbind",
        .packages = c("magrittr", "readr", "dplyr",
                      "tibble", "lubridate")) %dopar%
{
  # Identify subject
  # print(i)
  subj <- outcomes$subject_id[i]
  adm <- outcomes$hadm_id[i]
  
  # Extract admission and discharge time
  in_out_times <- mimic_icustays %>% 
    filter(subject_id == subj) %>% 
    filter(hadm_id == adm) %>% 
    select(intime, outtime) %>% 
    slice_min(outtime) 
  
  admit_time <- in_out_times %>% 
    select(intime) %>% deframe()
  discharge_time <- in_out_times %>% 
    select(outtime) %>% deframe()
  
  # Detect after-hours discharge
  after_hours_discharge <-   hour(discharge_time) <= 7 | hour(discharge_time) >= 16 
  
  # Extract sex
  sex <- mimic_patients %>% 
    filter(subject_id == subj) %>% 
    select(gender) %>% 
    deframe()
  
  # Extract surgical type
  surg_type <- mimic_services %>% 
    filter(hadm_id == adm) %>%
    mutate(before_discharge = discharge_time - transfertime) %>% 
    filter(before_discharge > 0) %>%
    select(curr_service) %>% deframe() 
  
  # Define surgical types
  general_surgery <-  "SURG" %in% surg_type
  cardiac_surgery <-  "CSURG" %in% surg_type
  
  # Duration of stay
  los <- mimic_preproc$stays %>% 
    # Select patient and admission
    filter(subject_id == subj,
           hadm_id == adm) %>% 
    slice_head() %>% 
    select(los) %>% deframe()
  los_5 <- los > 5
  los_7 <- los > 7
  
  # Extract age
  age <- mimic_preproc$stays %>% 
    # Select patient and admission
    filter(subject_id == subj,
           hadm_id == adm) %>% 
    slice_head() %>%
    select(dob, intime) %>% 
    # Generate age
    mutate(age = (intime - dob) / 365) %>%
    select(age) %>% deframe() %>% 
    as.numeric()
  
  # Correct > 90
  age <- ifelse(age > 90, ">90", age)
  
  # Extract admission type
  admission_info <- mimic_preproc$stays %>% 
    # Select patient and admission
    filter(subject_id == subj,
           hadm_id == adm) %>% 
    slice_head()
  
  # Determine elective
  elective_admission <- admission_info %>% 
    select(admission_type) %>% 
    mutate(elective = case_when(
      admission_type == "ELECTIVE" ~ TRUE,
      admission_type != "ELECTIVE" ~ FALSE,
    )) %>%
    select(elective) %>% deframe()
  
  # Determine source
  admission_source <- admission_info %>%
    select(admission_location) %>% 
    mutate(frost_source = case_when(
      admission_location == "EMERGENCY ROOM ADMIT" ~ "Emergency",
      admission_location == "CLINIC REFERRAL/PREMATURE"|
        admission_location == "TRSF WITHIN THIS FACILITY" ~ "Ward",
      admission_location == "PHYS REFERRAL/NORMAL DELI" ~ "OT",
      TRUE ~ "OtherHosp"
    )) %>% 
    select(frost_source) %>%  deframe()
  
  # Extract lab events
  labfile <- paste("data/events/labevents_", adm, ".csv", sep = "")
  labs <- read_csv(labfile)
  
  if(nrow(labs) > 0)
  {
    # Extract serum glucose measurements
    serum_glucose <- labs %>% 
      filter(itemid %in% ID_blood_glucose) %>% 
      event_filter_discharge(discharge_time)
    
    # Hypoglycemia within 24h discharge
    serum_glucose_24h <- serum_glucose %>% 
      filter(before_discharge <= 24) %>% 
      select(valuenum) %>% deframe()
    hyperglycemia <- any(serum_glucose_24h > 180)
    
    # Final serum glucose measurement
    serum_glucose <- serum_glucose %>% 
      filter(before_discharge < 48) %>% 
      summarise(serum_glucose = mean(valuenum, na.rm = TRUE)) %>% 
      deframe()
    
    # Anaemia within 24h discharge
    anaemia <- labs %>% 
      filter(itemid %in% ID_haemoglobin) %>% 
      event_filter_discharge(discharge_time) %>% 
      filter(before_discharge <= 24) %>% 
      select(valuenum) %>%
      summarise(anaemia = any(valuenum < 7)) %>% 
      select(anaemia) %>%  deframe()
    
    # Extract blood urea nitrogen
    blood_urea_nitrogen <- labs %>% 
      filter(itemid %in% ID_blood_urea_nitrogen) %>% 
      event_filter_discharge(discharge_time) %>% 
      filter(before_discharge < 48) %>% 
      summarise(blood_urea_nitrogen = mean(valuenum, na.rm = TRUE)) %>% 
      deframe()
    
    # Extract serum chloride
    serum_choride <- labs %>% 
      filter(itemid %in% ID_serum_choride) %>% 
      event_filter_discharge(discharge_time) %>% 
      filter(before_discharge < 48) %>% 
      summarise(serum_choride = mean(valuenum, na.rm = TRUE)) %>% 
      deframe()
    
    lab_missing <- FALSE
  }else
  {
    hyperglycemia <- NA
    serum_glucose <- NA
    anaemia <- NA
    blood_urea_nitrogen <- NA
    serum_choride <- NA
    lab_missing <- TRUE
  }
  
  # Load chart data
  chartfile <- paste("data/events/chartevents_", adm, ".csv", sep = "")
  if(file.exists(chartfile))
  {
    charts <- read_csv(chartfile)
    
    # Determine respiratory rate 
    respiratory_rate <- charts %>% 
      filter(ITEMID %in% ID_respiratory_rate) %>% 
      # Calculate time difference from discharge
      mutate(before_discharge = discharge_time - CHARTTIME) %>% 
      # Select only events before discharge
      filter(before_discharge > 0) %>% 
      filter(before_discharge < 48) %>% 
      summarise(respiratory_rate = mean(VALUENUM, na.rm = TRUE)) %>% 
      deframe()
    
    ambulation <- charts %>% 
      filter(ITEMID %in% ID_ambulation) %>% 
      # Calculate time difference from discharge
      mutate(before_discharge = discharge_time - CHARTTIME) %>% 
      # Select only events 24h before discharge
      filter(before_discharge > 0) %>% 
      filter(before_discharge < 24) %>% 
      summarise(ambulation = any(grepl("walk", VALUE,
                                       ignore.case = TRUE) |
                                   VALUENUM >2, na.rm = TRUE)) %>% 
      deframe()
    
    chart_missing <- FALSE
  }else
  {
    respiratory_rate <- NA
    ambulation <- NA
    chart_missing <- TRUE
  }
  
  # Determine acute renal failure
  acute_renal_failure <- mimic_preproc$diagnoses %>% 
    # Select patient and admission
    filter(subject_id == subj,
           hadm_id == adm) %>% 
    # Search for diagnoses of acute renal failure
    summarise(ARF = any(icd9_code %in% ID_acute_renal_failure)) %>% 
    deframe()
  
  # Extract past diagnoses for patient
  history <- mimic_admissions %>% 
    filter(subject_id == subj,
           admittime < admit_time) %>% 
    select(hadm_id) %>% 
    left_join(mimic_diagnoses,
              by = "hadm_id")
  
  # Extract aFib
  atrial_fibrillation <- history %>% 
    summarise(afib = any(icd9_code == ID_atrial_fibrillation,
                         na.rm = TRUE)) %>% 
    deframe()
  
  # Extract renal insufficiency
  renal_insufficiency <- history %>% 
    summarise(renal = any(icd9_code %in% ID_renal_insufficiency,
                         na.rm = TRUE)) %>% 
    deframe()
  
  # HAMMER------------
  #     APACHE II at admission (>20)
  #     Positive fluid balance (>5L over ICU stay)
  #     No ambulation in last 24h
  
  # FROST------------
  #---  Apache II score at admission
  

  
  # Output
  data.frame(row_id = i,
             subject_id = subj,
             adm_id = adm,
             chart_missing,
             lab_missing,
             readmission = outcomes$readmission[i],
             # Hammer variables
             sex, general_surgery, cardiac_surgery,
             hyperglycemia, anaemia, los_5, age, ambulation,
             # Frost variables
             after_hours_discharge, los_7, elective_admission,
             admission_source, acute_renal_failure,
             # Martin variables
             serum_glucose, blood_urea_nitrogen, serum_choride,
             respiratory_rate, atrial_fibrillation, renal_insufficiency
             )
}
stopImplicitCluster()
proc.time() - ptm # 171s

# Quality control------------

# No duplicated patients or admissions
length(unique(predictors$subject_id)) == nrow(predictors)
length(unique(predictors$adm_id)) == nrow(predictors)

# Count and exclude missing chart data
sum(predictors$chart_missing)
predictors %<>% filter(chart_missing == FALSE)

# Count and exclude missing lab data
sum(predictors$lab_missing)
predictors %<>% filter(lab_missing == FALSE)

# Count cases of readmission
table(predictors$readmission)

# Tabulate sex
table(predictors$sex)

# Tabulate surgical types
table(predictors$general_surgery)
table(predictors$cardiac_surgery)

# Tabulate lengths of stay
table(predictors$los_5)
table(predictors$los_7)

# Tabulate admission and discharge information
table(predictors$after_hours_discharge)
table(predictors$elective_admission)
table(predictors$admission_source)

# Tabulate hammer binary variables
table(predictors$hyperglycemia)
table(predictors$anaemia)
table(predictors$ambulation)

# Tabulate remaining binaries
table(predictors$acute_renal_failure)
table(predictors$atrial_fibrillation)
table(predictors$renal_insufficiency)

# Clean and summarise age
predictors$age[predictors$age == ">90"] <- 90
predictors$age <- floor(as.numeric(predictors$age))
summary(predictors$age)

# Summarise continuous variables
summary(predictors$serum_glucose) # 115
summary(predictors$serum_choride) # 132
summary(predictors$blood_urea_nitrogen) # 120
summary(predictors$respiratory_rate) # 1084
