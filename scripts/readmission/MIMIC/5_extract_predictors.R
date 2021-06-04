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

filter_window_lab <- function(.df, time)
{
  # Filter within 24 of admission
    out_df <- .df %>% 
      mutate(window = difftime(charttime, time,
                               units = "hours"))  %>% 
      filter(window < 24 & 
               window >= 0)
  
  # Take nearest pre-admission if none available
  if(nrow(out_df) == 0)
  {
    out_df <- .df %>% 
      mutate(window = difftime(charttime, time,
                               units = "hours"))  %>% 
      filter(window < 0) %>% 
      slice_min(window)
  }
  return(out_df)
}

filter_window_lab_discharge <- function(.df, time)
{
  # Filter within 48h of discharge
  out_df <- .df %>% 
    mutate(window = difftime(charttime, time,
                             units = "hours"))  %>% 
    filter(window > -48 & 
             window <= 0)
  
  # Take nearest if none available
  if(nrow(out_df) == 0)
  {
    out_df <- .df %>% 
      mutate(window = difftime(charttime, time,
                               units = "hours"))  %>% 
      filter(window <= 0) %>% 
      slice_max(window)
  }
  
  # Return worse, if any
  if(nrow(out_df) > 0)
  {
    out_df %>% 
      select(valuenum) %>% 
      deframe() %>% max(na.rm = T) %>% 
      return()
  }else
  {
    return(NA)
  }
}

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
mimic_procedures <- mimic$cptevents %>% collect()
mimic_callouts <- mimic$callout %>% collect()

# Database of high-risk specialities
highrisk <- mimic_procedures %>% 
  filter(sectionheader == "Surgery",
         subsectionheader %in% c("Digestive system", "Mediastinum and diaphragm"))

# Prepare parallel options
ptm <- proc.time()
psnice(value = 19)
registerDoParallel(ifelse(detectCores() <= 15,
                          detectCores() - 1,
                          12)
)

# Run loop
predictors <- foreach(i = 1:nrow(outcomes),
                      .combine = "rbind",
        .packages = c("magrittr", "readr", "dplyr",
                      "tibble", "lubridate", "purrr",
                      "tidyr")) %dopar%
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
  
  # Identify callouts
  callout_time <- mimic_callouts %>% 
    filter(hadm_id == adm) %>% 
    filter(acknowledge_status == "Acknowledged") 
  
  # Take closest if multiple
  if(nrow(callout_time) > 1)
  {
    callout_time <- callout_time[which.min(abs(callout_time$outcometime - discharge_time)),]
  }
  callout_time <- callout_time %>% 
    select(createtime) %>% 
    deframe()
  
  # Measure ready for discharge days
  if(length(callout_time) > 0)
  {
    discharge_hours <- floor(difftime(discharge_time, callout_time,
                                     unit = "hours"))
  }else
  {
    discharge_hours <- NA
  }

  
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
  
  # Determine high-risk speciality
  high_risk_speciality <- nrow(highrisk %>% filter(hadm_id == adm)) > 0
  
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
  
  # Determine hospital days before ICU
  days_before_ICU <- admission_info %>% 
    mutate(days_before_ICU = difftime(intime, admittime,
                                      units = "days")) %>% 
    select(days_before_ICU) %>% 
    deframe() %>% ceiling()
  
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

  if(file.exists(labfile))
  {
    labs <- read_csv(labfile)
    
    # Extract serum glucose measurements
    serum_glucose_measurements <- labs %>% 
      filter(itemid %in% ID_blood_glucose) %>% 
      filter_window_lab(admit_time)
    
    # Extract haemoglobin measurements
    haemoglobin_measurements <- labs %>% 
      filter(itemid %in% ID_haemoglobin) %>% 
      filter_window_lab(admit_time) 
    
    # Extract blood urea nitrogen
    blood_urea_nitrogen_measurements <- labs %>% 
      filter(itemid %in% ID_blood_urea_nitrogen) %>% 
      filter_window_lab(admit_time)
    
    # Extract serum chloride
    serum_choride_measurements <- labs %>% 
      filter(itemid %in% ID_serum_choride) %>% 
      filter_window_lab(admit_time) 
    
    # If no serum glucose measurements
    if(nrow(serum_glucose_measurements) == 0)
    {
      hyperglycemia <- NA
      serum_glucose <- NA
    }else
    {
      # Hypoglycemia within 24h discharge
      serum_glucose <- serum_glucose_measurements %>% 
        summarise(valuenum = max(valuenum, na.rm = T)) %>% 
        select(valuenum) %>% deframe()
      hyperglycemia <- serum_glucose > 180
    }
    
    # If no haemoglobin measurements
    if(nrow(haemoglobin_measurements) == 0)
    {
      anaemia <- NA
    }else
    {
      # Anaemia within 24h discharge
      anaemia <- haemoglobin_measurements %>% 
        summarise(valuenum = min(valuenum, na.rm = T)) %>% 
        select(valuenum) %>%
        summarise(anaemia = valuenum < 9) %>% 
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
        summarise(valuenum = max(valuenum, na.rm = T)) %>% 
        select(valuenum) %>% deframe()
    }
    
    # If no serum chloride measurements
    if(nrow(serum_choride_measurements) == 0)
    {
      serum_choride <- NA
    }else
    {
      # Final serum chloride measurement
      serum_chloride <- serum_choride_measurements %>% 
        summarise(valuenum = max(valuenum, na.rm = T)) %>% 
        select(valuenum) %>% deframe()
    }
    
    # Final serum glucose
    serum_glucose_final <- labs %>% 
      filter(itemid %in% ID_blood_glucose) %>% 
      filter_window_lab_discharge(discharge_time)
    
    # Final haemoglobin
    haemoglobin_final <- labs %>% 
      filter(itemid %in% ID_haemoglobin) %>% 
      filter_window_lab_discharge(discharge_time)
    
    # Final BUN
    blood_urea_nitrogen_final <- labs %>% 
      filter(itemid %in% ID_blood_urea_nitrogen) %>% 
      filter_window_lab_discharge(discharge_time)
    
    # Final chloride
    serum_chloride_final <- labs %>% 
      filter(itemid %in% ID_serum_choride) %>% 
      filter_window_lab_discharge(discharge_time)
    
    
    lab_missing <- FALSE
  }else
  {
    hyperglycemia <- NA
    serum_glucose <- NA
    anaemia <- NA
    blood_urea_nitrogen <- NA
    serum_chloride <- NA
    lab_missing <- TRUE
    
    serum_glucose_final <- NA
    haemoglobin_final <- NA
    blood_urea_nitrogen_final <- NA
    serum_chloride_final <- NA
  }
  
  # Charted predictors--------
  
  # Load chart data
  chartfile <- paste("data/events/chartevents_", adm, ".csv", sep = "")
  if(file.exists(chartfile))
  {
    charts <- read_csv(chartfile)
    
    # Determine respiratory rate before discharge----
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
    
    # Determine respiratory rate upon admission----
    admission_resp_rate <- charts %>% 
      filter(ITEMID %in% ID_respiratory_rate) %>% 
      mutate(window = difftime(CHARTTIME, admit_time,
                               units = "hours"))  %>% 
      filter(window < 24 & 
               window >= 0)
    
    # Take nearest if none
    if(nrow(admission_resp_rate) == 0)
    {
      admission_resp_rate <- charts %>% 
        filter(ITEMID %in% ID_respiratory_rate) %>% 
        mutate(window = difftime(CHARTTIME, admit_time,
                                 units = "hours"))  %>% 
        filter(window < 24) %>% 
        slice_min(window)
    }
    
    if(nrow(admission_resp_rate) > 0)
    {
      admission_resp_rate %<>% 
        select(VALUENUM) %>% 
        deframe() %>% max(na.rm = T)
    }else
    {
      admission_resp_rate <- NA
    }
    
    
    # Extract ambulation data----
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
    
    # Extract GCS measurements----
    gcs_motor <- charts %>% 
      filter(ITEMID %in% c(454, 223901)) %>% 
      select(CHARTTIME, VALUENUM) %>% 
      mutate(type = "motor")
    
    gcs_verbal <- charts %>% 
      filter(ITEMID %in% c(723, 223900)) %>%
      filter(VALUE != "No Response-ETT") %>% 
      select(CHARTTIME, VALUENUM) %>% 
      mutate(type = "verbal")
    
    gcs_eye <- charts %>% 
      filter(ITEMID %in% c(184, 220739)) %>% 
      select(CHARTTIME, VALUENUM) %>% 
      mutate(type = "eye")
    
    # Combine and filter to first 24h
    gcs_df <- rbind(gcs_motor,
                         gcs_verbal,
                         gcs_eye) %>% 
      group_by(CHARTTIME) %>% 
      mutate(freq = n()) %>% 
      ungroup() %>% 
      filter(freq == 3) %>% 
      select(-freq) %>% 
      pivot_wider(names_from = "type",
                  values_from = "VALUENUM") %>% 
      filter(CHARTTIME > admit_time,
             CHARTTIME < (admit_time + 86400))
    
    # Measure if below 15
    if(nrow(gcs_df) > 0)
    {
      if(min(na.omit(rowSums(gcs_df[,2:4]))) < 15)
      {
        glasgow_coma_below_15 <- TRUE
      }else
      {
        glasgow_coma_below_15 <- FALSE
      }
    }else
    {
      glasgow_coma_below_15 <- NA
    }
    
    # Determine respiratory support----
    respiratory_support_df <- charts %>% 
      filter(ITEMID %in% c(3605, 225792,
                           225794,225303))
    
    respiratory_support <- nrow(respiratory_support_df) >0

    chart_missing <- FALSE
  }else
  {
    respiratory_rate <- NA
    rr_time <- NA
    ambulation <- NA
    glasgow_coma_below_15 <- NA
    respiratory_support <- NA
    admission_resp_rate <- NA
    chart_missing <- TRUE
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
  
  # APACHE II----
  if(!chart_missing & !lab_missing & !is.na(respiratory_rate))
  {
    # Estimate ARF during first 24h----
    if(file.exists(outputfile))
    {
      # Extract first 24h creatinine measurements
      creatinine <- labs %>% 
        filter(itemid %in% c(50912)) %>% 
        mutate(window = difftime(charttime, admit_time,
                                 units = "hours"))  %>% 
        filter(window < 24& 
                 window >= 0) %>% 
        select(valuenum) %>% deframe()
      
      if(length(creatinine) > 0)
      {
        highest_creatinine <- max(creatinine)
      }else
      {
        highest_creatinine <- NA
      }
      
      # Extract first 24h urine output
      urine <- read_csv(outputfile, col_types = c("ddddTddcTdlll")) %>% 
        # Carevue urine items
        filter(itemid %in% c(40055, 43175, 43175, 40069,
                             40094, 40715, 40473, 40085,
                             40057, 40056, 40405, 40428,
                             40086, 40096, 40651,
                             # Metavision urine items 
                             226559, 226560, 226561,
                             226584, 226563, 226564,
                             226565, 226567, 226557,
                             226558, 227488, 227489)) %>% 
        # Correct for input of GU irrigants
        mutate(value = ifelse(itemid == 227488, -value, value)) %>% 
        # Standardise units to ml and round
        mutate(amount = ifelse(valueuom == "L", 
                               value * 1000, value) %>% 
                 round()) %>% 
        # Filter to first 24h
        mutate(window = difftime(charttime, admit_time,
                                 units = "hours"))  %>% 
        filter(window < 24& 
                 window >= 0) %>% 
        # Take total
        select(value) %>% deframe() %>% sum(na.rm = T) 
      
      # Load prediction coefficients
      coef_arf <- read_rds("models/acute_renal_failure.RDS")
      
      # Predict ARF
      probs_arf <- coef_arf[1] + (coef_arf[2] * highest_creatinine) + (coef_arf[3] * urine) + 
        (coef_arf[4] * (highest_creatinine * urine))
      acute_renal_failure_24h <- probs_arf > coef_arf[5]
    }else
    {
      acute_renal_failure_24h <- acute_renal_failure
    }
    
    if(is.na(acute_renal_failure_24h))
    {
      acute_renal_failure_24h <- acute_renal_failure
    }
    
    # Calculate apache score at admission----
    apache_vector <- apacheII_score(labs_df = labs, chart_df = charts,
                                    patient_df = mimic_preproc$stays %>% 
                                      filter(subject_id == subj,
                                             hadm_id == adm),
                                    admission_time = admit_time,
                                    elect_admit = elective_admission,
                                    prev_diagnoses = history,
                                    arf = acute_renal_failure_24h)
    # Sum score
    if(is.na(apache_vector$full_scores["bicarbonate"]))
    {
      apache_II <- sum(apache_vector$full_scores[c("temperature", "map", "pulse", "respiratory",
                                                   "oxygen", "arterialpH", "sodium", "potassium",
                                                   "creatinine", "haematocrit", "whitebloodcount",
                                                   "glasgowcoma", "age", "chronic")], na.rm = T)
    }else
    {
      apache_II <- sum(apache_vector$full_scores[c("temperature", "map", "pulse", "respiratory",
                                                   "bicarbonate", "sodium", "potassium",
                                                   "creatinine", "haematocrit", "whitebloodcount",
                                                   "glasgowcoma", "age", "chronic")], na.rm = T)
    }
    
    high_apache <- apache_II > 20
  }else
  {
    acute_renal_failure_24h <- NA
    apache_II <- NA
    high_apache <- NA
    apache_vector <- list(full_values = c("temperature" = NA, "map" = NA, "pulse" = NA,
                                          "respiratory" = NA, "arterialpH" = NA, "sodium" = NA,
                                          "potassium" = NA, "creatinine" = NA, "haematocrit" = NA,
                                          "whitebloodcount" = NA, "glasgowcoma" = NA, "bicarbonate" = NA,
                                          "fractioninspiredoxygen" = NA, "arterialoxygen" = NA,
                                          "arterialcarbon" = NA))
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
             hyperglycemia_initial = hyperglycemia,
             anaemia_initial = anaemia, los_5, age, ambulation,
             fluid_balance_5L,
             # Frost variables
             after_hours_discharge, los_7, elective_admission,
             admission_source, apache_II,
             acute_renal_failure_initial = acute_renal_failure_24h,
             # Martin variables
             serum_glucose_initial = serum_glucose, 
             blood_urea_nitrogen_initial = blood_urea_nitrogen,
             serum_chloride_initial = serum_chloride,
             respiratory_rate_initial = admission_resp_rate,
             atrial_fibrillation, renal_insufficiency,
             # APACHE-II variables
             temperature = apache_vector$full_values["temperature"],
             mean_arterial_pressure = apache_vector$full_values["map"],
             pulse = apache_vector$full_values["pulse"],
             respiratory_rate = apache_vector$full_values["respiratory"],
             arterial_ph = apache_vector$full_values["arterialpH"],
             sodium = apache_vector$full_values["sodium"],
             potassium = apache_vector$full_values["potassium"],
             creatinine = apache_vector$full_values["creatinine"],
             hematocrit = apache_vector$full_values["haematocrit"],
             wbc = apache_vector$full_values["whitebloodcount"],
             gcs = apache_vector$full_values["glasgowcoma"],
             arterial_o2 = apache_vector$full_values["arterialoxygen"],
             arterial_co2 = apache_vector$full_values["arterialcarbon"],
             bicarbonate = apache_vector$full_values["bicarbonate"],
             # Additional variables
             glasgow_coma_below_15, days_before_ICU, respiratory_support,
             high_risk_speciality, length_of_stay = ceiling(los),
             discharge_hours = ifelse(discharge_hours >= 0, discharge_hours,
                                     NA),
             discharge_delay = discharge_hours > 24,
             acute_renal_failure_total = acute_renal_failure,
             serum_glucose_final, haemoglobin_final,
             blood_urea_nitrogen_final, serum_chloride_final,
             respiratory_rate_final = respiratory_rate,
             hyperglycemia_final = serum_glucose_final > 180,
             anaemia_final = haemoglobin_final < 9
             )
  #row.names(output) <- i
  output
  #data.frame(i, cols = ncol(output))
}
row.names(predictors) <- c()
stopImplicitCluster()
proc.time() - ptm # 320 s

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


