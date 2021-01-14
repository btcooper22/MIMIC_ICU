# Packages
require(tidyverse)
require(magrittr)
require(lubridate)
require(doParallel)
require(tools)

# Load MIMIC database and other functions
source("functions/mimic_load.R")
source("functions/event_filter_discharge.R")
source("functions/flag_within_discharge.R")
source("functions/apache_ii.R")
source("functions/fluid_balance.R")
source("functions/fialho_variables.R")
source("functions/NA_count.r")

# Load preprocessed data
mimic_preproc <- read_rds("data/mimic_preprocessed.RDS")

# Load outcomes
outcomes <- read_csv("data/outcomes.csv")

# Correct "both" datasources where inaccurate-----------

both_cases <- which(mimic_preproc$stays$dbsource == "both")
new_both <- foreach(i = 1:length(both_cases),.combine = "c") %do%
  {
    # Find admission
    hadm <- mimic_preproc$stays[both_cases[i],3]
    
    # Assess if admission in main data
    if(hadm %in% outcomes$hadm_id)
    {
      # Find files
      inputfile_mv <- paste("data/events/inputeventsmv_", hadm, ".csv", sep = "")
      inputfile_cv <- paste("data/events/inputeventscv_", hadm, ".csv", sep = "")
      
      # Create new dbsource
      case_when(
        file.exists(inputfile_mv) & file.exists(inputfile_cv) ~ "both",
        file.exists(inputfile_mv) ~ "metavision",
        file.exists(inputfile_cv) ~ "carevue"
      )
    }else
    {
      # Flag as not in data
      "irrelevant"
    }
  }
mimic_preproc$stays$dbsource[mimic_preproc$stays$dbsource == "both"] <- new_both
# test <- mimic_preproc$stays %>% filter(dbsource == "both")

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
                      "tibble", "lubridate", "purrr")) %dopar%
{
  # Identify subject
  # print(i)
  subj <- outcomes$subject_id[i]
  adm <- outcomes$hadm_id[i]
  
  # ICU stay and demographic predictors------------
  
  # Extract admission and discharge time
  in_out_times <- mimic_icustays %>% 
    filter(subject_id == subj) %>% 
    filter(hadm_id == adm) %>% 
    select(intime, outtime) %>% 
    slice_min(outtime) 
  
  admit_time <- in_out_times %>% 
    select(intime) %>% 
    deframe() %>% force_tz("UTC")
  discharge_time <- in_out_times %>% 
    select(outtime) %>% deframe() %>% force_tz("UTC")
  
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
    mutate(before_discharge = difftime(discharge_time, transfertime, unit = "hours")) %>% 
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
  
  # History predictors---------
  
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
  
  # Lab predictors----------
  
  # Extract lab events
  labfile <- paste("data/events/labevents_", adm, ".csv", sep = "")
  labs <- read_csv(labfile)
  
  if(nrow(labs) > 0)
  {
    # Extract serum glucose measurements
    serum_glucose_measurements <- labs %>% 
      filter(itemid %in% ID_blood_glucose) %>% 
      event_filter_discharge(discharge_time)
    
    # Extract haemoglobin measurements
    haemoglobin_measurements <- labs %>% 
      filter(itemid %in% ID_haemoglobin) %>% 
      event_filter_discharge(discharge_time) 
    
    # Extract blood urea nitrogen
    blood_urea_nitrogen_measurements <- labs %>% 
      filter(itemid %in% ID_blood_urea_nitrogen) %>% 
      event_filter_discharge(discharge_time)
    
    # Extract serum chloride
    serum_choride_measurements <- labs %>% 
      filter(itemid %in% ID_serum_choride) %>% 
      event_filter_discharge(discharge_time) 
    
    # If no serum glucose measurements
    if(nrow(serum_glucose_measurements) == 0)
    {
      hyperglycemia <- NA
      serum_glucose <- NA
    }else
    {
      # Hypoglycemia within 24h discharge
      serum_glucose_24h <- serum_glucose_measurements %>% 
        flag_within_discharge(24) %>% 
        select(valuenum) %>% deframe()
      hyperglycemia <- any(serum_glucose_24h > 180)
      
      # Final serum glucose measurement
      serum_glucose <- serum_glucose_measurements %>% 
        flag_within_discharge(48) %>% 
        summarise(serum_glucose = mean(valuenum, na.rm = TRUE)) %>% 
        deframe()
    }
    
    # If no haemoglobin measurements
    if(nrow(haemoglobin_measurements) == 0)
    {
      anaemia <- NA
    }else
    {
      # Anaemia within 24h discharge
      anaemia <- haemoglobin_measurements %>% 
        flag_within_discharge(24) %>% 
        select(valuenum) %>%
        summarise(anaemia = any(valuenum < 7)) %>% 
        select(anaemia) %>%  deframe()
    }
    
    # If no BUN measurements
    if(nrow(blood_urea_nitrogen_measurements) == 0)
    {
      blood_urea_nitrogen <- NA
    }else
    {
      # Final BUN measurement
      blood_urea_nitrogen <- blood_urea_nitrogen_measurements %>% 
        flag_within_discharge(48) %>% 
        summarise(blood_urea_nitrogen = mean(valuenum, na.rm = TRUE)) %>% 
        deframe()
    }
    
    # If no serum chloride measurements
    if(nrow(serum_choride_measurements) == 0)
    {
      serum_choride <- NA
    }else
    {
      # Final serum chloride measurement
      serum_choride <- serum_choride_measurements %>% 
        flag_within_discharge(48) %>% 
        summarise(serum_choride = mean(valuenum, na.rm = TRUE)) %>% 
        deframe()
    }
    
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
  
  # Charted predictors--------
  
  # Load chart data
  chartfile <- paste("data/events/chartevents_", adm, ".csv", sep = "")
  if(file.exists(chartfile))
  {
    charts <- read_csv(chartfile)
    
    # Determine respiratory rate 
    respiratory_rate <- charts %>% 
      filter(ITEMID %in% ID_respiratory_rate) %>% 
      # Calculate time difference from discharge
      mutate(before_discharge = difftime(discharge_time, CHARTTIME, unit = "hours")) %>% 
      # Select only events before discharge
      filter(before_discharge > 0) %>% 
      filter(before_discharge < 48)
    
    if(nrow(respiratory_rate) == 0)
    {
      # If no RR measurement within 48h of discharge
      respiratory_rate <- charts %>% 
        filter(ITEMID %in% ID_respiratory_rate) %>% 
        # Calculate time difference from discharge
        mutate(before_discharge = difftime(discharge_time, CHARTTIME, unit = "hours")) %>% 
        # Select only events before discharge
        filter(before_discharge > 0) %>% 
        # Find most recent
        slice_min(before_discharge)
      
      # If no RR measurements AT ALL before discharge
      if(nrow(respiratory_rate) == 0)
      {
        respiratory_rate <- NA
        rr_time <- NA
      }else
      {
        # Record time before discharge
        rr_time <- respiratory_rate$before_discharge %>% deframe()
        
        # Determine respiratory rate 
        respiratory_rate %<>% 
          summarise(respiratory_rate = mean(VALUENUM, na.rm = TRUE)) %>% 
          deframe()
      }

    }else
    {
      # Determine respiratory rate 
      respiratory_rate %<>% 
        summarise(respiratory_rate = mean(VALUENUM, na.rm = TRUE)) %>% 
        deframe()
      rr_time <- 48
    }
    
    # Extract ambulation data
    ambulation_df <- charts %>% 
      filter(ITEMID %in% ID_ambulation) %>% 
      # Calculate time difference from discharge
      mutate(before_discharge = difftime(discharge_time, CHARTTIME, unit = "hours")) %>% 
      filter(before_discharge > 0) 
    
    # Flag if no measurements 24h before discharge
    if(min(ambulation_df$before_discharge) > 24 | nrow(ambulation_df) == 0)
    {
      ambulation <- NA
    }else
    {
      # Select only events 24h before discharge
      ambulation <- ambulation_df %>% 
        filter(before_discharge < 24) %>% 
        summarise(ambulation = any(grepl("walk", VALUE,
                                         ignore.case = TRUE) |
                                     VALUENUM >2, na.rm = TRUE)) %>% 
        deframe()
    }

    chart_missing <- FALSE
  }else
  {
    respiratory_rate <- NA
    ambulation <- NA
    chart_missing <- TRUE
  }
  
  # APACHE II
  if(!chart_missing | !lab_missing | !is.na(respiratory_rate))
  {
    # Calculate apache score at admission
    apache_score_vector <- apacheII_score(labs_df = labs, chart_df = charts,
                                          patient_df = mimic_preproc$stays %>% 
                                            filter(subject_id == subj,
                                                   hadm_id == adm),
                                          admission_time = admit_time,
                                          elect_admit = elective_admission,
                                          prev_diagnoses = history,
                                          arf = acute_renal_failure)
    # Sum score
    apache_II <- sum(apache_score_vector)
    high_apache <- apache_II > 20
    
    # # Calculate apache score at discharge
    # apache_score_vector_discharge <- apacheII_score(labs_df = labs, chart_df = charts,
    #                                       patient_df = mimic_preproc$stays %>% 
    #                                         filter(subject_id == subj,
    #                                                hadm_id == adm),
    #                                       admission_time = discharge_time,
    #                                       elect_admit = elective_admission,
    #                                       prev_diagnoses = history,
    #                                       arf = acute_renal_failure)
    # # Sum score
    # apache_II_discharge <- sum(apache_score_vector_discharge)
  }else
  {
    apache_II <- NA
    high_apache <- NA
    apache_II_discharge <- NA
  }
  
  # Fialho's variables
  if(!chart_missing | !lab_missing | !is.na(respiratory_rate))
  {
    fialho_vars <- fialho_variables(discharge_time,
                                    charts,
                                    labs)
  }else
  {
    fialho_vars <- rep(NA, 6)
  }
  
  
  # I/O Predictors----------
  
  # Positive fluid balance (>5L)
  if(!chart_missing)
  {
    # Check if input and output files exist
    inputfile_mv <- paste("data/events/inputeventsmv_", adm, ".csv", sep = "")
    inputfile_cv <- paste("data/events/inputeventscv_", adm, ".csv", sep = "")
    outputfile <- paste("data/events/outputevents_", adm, ".csv", sep = "")
    
    if(file.exists(outputfile) & (file.exists(inputfile_mv) | file.exists(inputfile_cv)))
    {
      # Calculate fluid balance
      max_positive_fluid_balance <- fluid_balance(icustay_df = mimic_preproc$stays %>% 
                                                    filter(hadm_id == adm),
                                                  chart_df = charts,
                                                  readmission = outcomes$readmission[i])
      fluid_balance_5L <- max_positive_fluid_balance >= 5000
      
    } else
    {
      # Record as NA
      fluid_balance_5L <- NA
    }
  }else
  {
    fluid_balance_5L <- NA
  }
  
  # Final output----------
  
  # Output
  output <- data.frame(row_id = i,
             subject_id = subj,
             adm_id = adm,
             chart_missing,
             lab_missing,
             rr_time,
             readmission = outcomes$readmission[i],
             # Hammer variables
             sex, general_surgery, cardiac_surgery, high_apache,
             hyperglycemia, anaemia, los_5, age, ambulation,
             fluid_balance_5L,
             # Frost variables
             after_hours_discharge, los_7, elective_admission,
             admission_source, acute_renal_failure, apache_II,
             # Martin variables
             serum_glucose, blood_urea_nitrogen, serum_choride,
             respiratory_rate, atrial_fibrillation, renal_insufficiency,
             # Fialho variables
             final_pulse = fialho_vars[1], final_temp = fialho_vars[2],
             final_pO2 = fialho_vars[3], final_bp = fialho_vars[4],
             final_platelets = fialho_vars[5], final_lactate = fialho_vars[6],
             # Additional variables
             #apache_II_discharge,
             apache_temperature = apache_score_vector["temperature"],
             apache_map = apache_score_vector["map"],
             apache_pulse = apache_score_vector["pulse"],
             apache_respiratory = apache_score_vector["oxygen"],
             apache_artpH = apache_score_vector["arterialpH"],
             apache_sodium = apache_score_vector["sodium"],
             apache_potassium = apache_score_vector["potassium"],
             apache_creatinine = apache_score_vector["creatinine"],
             apache_hematocrit = apache_score_vector["haematocrit"],
             apache_wbc = apache_score_vector["whitebloodcount"],
             apache_gcs = apache_score_vector["glasgowcoma"],
             apache_age = apache_score_vector["age"],
             apache_oxygenation = apache_score_vector["oxygen"],
             apache_chronic = apache_score_vector["chronic"]
             )
  #row.names(output) <- i
  output
  #data.frame(i, cols = ncol(output))
}
stopImplicitCluster()
proc.time() - ptm # 1276s (~20 min)

# Write and quality control ------------

# Write to file
predictors %>% write_csv("data/predictors.csv")

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
table(predictors$readmission, useNA = "ifany")

# Summarise APACHE-II scores
summary(predictors$apache_II) # 7

# Tabulate fluid balance
table(predictors$fluid_balance_5L, useNA = "ifany") # 45

# Tabulate sex
table(predictors$sex, useNA = "ifany")

# Tabulate surgical types
table(predictors$general_surgery, useNA = "ifany")
table(predictors$cardiac_surgery, useNA = "ifany")

# Tabulate lengths of stay
table(predictors$los_5, useNA = "ifany")
table(predictors$los_7, useNA = "ifany")

# Tabulate admission and discharge information
table(predictors$after_hours_discharge, useNA = "ifany")
table(predictors$elective_admission, useNA = "ifany")
table(predictors$admission_source, useNA = "ifany")

# Tabulate hammer binary variables
table(predictors$hyperglycemia, useNA = "ifany")
table(predictors$anaemia, useNA = "ifany")
table(predictors$ambulation, useNA = "ifany")

# Tabulate remaining binaries
table(predictors$acute_renal_failure, useNA = "ifany")
table(predictors$atrial_fibrillation, useNA = "ifany")
table(predictors$renal_insufficiency, useNA = "ifany")

# Clean and summarise age
predictors$age[predictors$age == ">90"] <- 90
predictors$age <- floor(as.numeric(predictors$age))
summary(predictors$age)

# Summarise respiratory rate measurements
summary(predictors %>% filter(rr_time == 48) %>% 
          select(respiratory_rate))
summary(predictors %>% filter(rr_time > 48) %>% 
          select(respiratory_rate, rr_time))

# Summarise continuous variables
summary(predictors$serum_glucose)
summary(predictors$serum_choride)
summary(predictors$blood_urea_nitrogen)
summary(predictors$respiratory_rate)


