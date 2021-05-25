# Packages
require(readr)
require(dplyr)
require(magrittr)
require(lubridate)
require(doParallel)
require(tools)
require(tibble)

# Load MIMIC database and other functions
source("functions/mimic_load.R")
source("functions/event_filter_discharge.R")
source("functions/flag_within_discharge.R")
source("functions/apache_ii_mort.R")
source("functions/NA_count.r")

# Load preprocessed data
mimic_preproc <- read_rds("data/mimic_preprocessed_missing.RDS")

# Load outcomes
outcomes <- read_csv("data/outcomes_mortality.csv")

# Extraction loop----------

ID_acute_renal_failure <- mimic_preproc$diagnoses %>% 
  filter(grepl("Acute kidney failure", long_title)) %>% 
  select(icd9_code) %>% 
  deframe() %>% unique()

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
                          15)
)

# EXAMPLES
# i <- 1
# i <- which(outcomes$in_unit_mortality == TRUE)[1]
# i <- which(outcomes$in_hospital_mortality == TRUE & 
#              outcomes$in_unit_mortality == FALSE)[1]
# i <- which(outcomes$mortality_30d == TRUE &
#              outcomes$in_hospital_mortality == FALSE)[1]

# Run loop
predictors <- foreach(i = 1:nrow(outcomes), .combine = "rbind",
        .packages = c("magrittr", "readr", "dplyr",
                      "tibble", "lubridate", "purrr"),
        .errorhandling = "remove") %dopar%
{
  # Identify subject and file
  subj <- outcomes$subject_id[i]
  adm <- outcomes$hadm_id[i]
  filename <- paste("data/missing_full/", adm,
                    ".csv", sep = "")
  
  if(!file.exists(filename))
  {
  # Demographics and stay----
  
  # Extract admission and discharge times
  in_out_times <- mimic_icustays %>% 
    filter(subject_id == subj) %>% 
    filter(hadm_id == adm) %>% 
    select(intime, outtime) %>% 
    slice_min(outtime) 
  
  admit_time <- in_out_times %>% 
    select(intime) %>% 
    deframe() %>% force_tz("UTC")
  
  # If outtime is NA?

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
  
  # Medical history---------
  
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

  # Lab predictors----------
  
  # Extract lab events
  labfile <- paste("data/events/labevents_", adm, ".csv", sep = "")
  labs <- read_csv(labfile)
  lab_missing <- !nrow(labs) > 0
  
  # Charted predictors--------
  
  # Load chart data
  chartfile <- paste("data/events/chartevents_", adm, ".csv", sep = "")
  charts <- read_csv(chartfile)
  chart_missing <- !file.exists(chartfile)
  
  # APACHE II
  if(!chart_missing & !lab_missing)
  {
    # Calculate apache score at admission
    apache_admission <- apacheII_score(labs_df = labs, chart_df = charts,
                                       patient_df = mimic_preproc$stays %>% 
                                         filter(subject_id == subj,
                                                hadm_id == adm),
                                       event_time = admit_time,
                                       elect_admit = elective_admission,
                                       prev_diagnoses = history,
                                       arf = acute_renal_failure,
                                       event_type = "admission")
    
    # Check if survived ICU stay
    if(outcomes$in_unit_mortality[i] == FALSE)
    {
      # Extract discharge time
      discharge_time <- in_out_times %>% 
        select(outtime) %>% deframe() %>% force_tz("UTC")
      
      # Calculate apache score at discharge
      apache_discharge <- apacheII_score(labs_df = labs, chart_df = charts,
                                         patient_df = mimic_preproc$stays %>%
                                           filter(subject_id == subj,
                                                  hadm_id == adm),
                                         event_time = discharge_time,
                                         elect_admit = elective_admission,
                                         prev_diagnoses = history,
                                         arf = acute_renal_failure,
                                         event_type = "discharge")
    }else
    {
      apache_discharge <- NA
    }
    
  }else
  {
    apache_admission <- NA
    apache_discharge <- NA
  }
  
  # Final output----------
  
  # Output
  output <- data.frame(row_id = i,
             subject_id = subj,
             adm_id = adm,
             chart_missing,
             lab_missing,
             acute_renal_failure,
             mort_inunit = outcomes$in_unit_mortality[i],
             mort_inhosp = outcomes$in_hospital_mortality[i],
             mort_30 = outcomes$mortality_30d[i],
             a_temperature = apache_admission["temperature"],
             a_systolicbp = apache_admission["systolicbp"],
             a_diastolicbp = apache_admission["diastolicbp"],
             a_pulse = apache_admission["pulse"],
             a_respiratory = apache_admission["respiratory"],
             a_arterialpH = apache_admission["arterialpH"],
             a_sodium = apache_admission["sodium"],
             a_potassium = apache_admission["potassium"],
             a_creatinine = apache_admission["creatinine"],
             a_haematocrit = apache_admission["haematocrit"],
             a_whitebloodcount = apache_admission["whitebloodcount"],
             a_glasgowcomaeye = apache_admission["glasgowcomaeye"],
             a_glasgowcomamotor = apache_admission["glasgowcomamotor"],
             a_glasgowcomaverbal = apache_admission["glasgowcomaverbal"],
             a_age = apache_admission["age"],
             a_chronic = apache_admission["chronic"],
             a_bicarbonate = apache_admission["bicarbonate"],
             a_fractioninspiredoxygen = apache_admission["fractioninspiredoxygen"],
             a_arterialoxygen = apache_admission["arterialoxygen"],
             a_arterialcarbon = apache_admission["arterialcarbon"],
             d_temperature = apache_discharge["temperature"],
             d_systolicbp = apache_discharge["systolicbp"],
             d_diastolicbp = apache_discharge["diastolicbp"],
             d_pulse = apache_discharge["pulse"],
             d_respiratory = apache_discharge["respiratory"],
             d_arterialpH = apache_discharge["arterialpH"],
             d_sodium = apache_discharge["sodium"],
             d_potassium = apache_discharge["potassium"],
             d_creatinine = apache_discharge["creatinine"],
             d_haematocrit = apache_discharge["haematocrit"],
             d_whitebloodcount = apache_discharge["whitebloodcount"],
             d_glasgowcomaeye = apache_discharge["glasgowcomaeye"],
             d_glasgowcomamotor = apache_discharge["glasgowcomamotor"],
             d_glasgowcomaverbal = apache_discharge["glasgowcomaverbal"],
             d_age = apache_discharge["age"],
             d_chronic = apache_discharge["chronic"],
             d_bicarbonate = apache_discharge["bicarbonate"],
             d_fractioninspiredoxygen = apache_discharge["fractioninspiredoxygen"],
             d_arterialoxygen = apache_discharge["arterialoxygen"],
             d_arterialcarbon = apache_discharge["arterialcarbon"])
  #row.names(output) <- i
  write_csv(output, filename)
  rm(charts,labs)
  gc()
  }
  
  output
  #data.frame(i, cols = ncol(output))
}
row.names(predictors) <- c()
stopImplicitCluster()
proc.time() - ptm

#
predictor_files <- dir("data/missing_full", full.names = TRUE)
cl <- makeCluster(15)
predictors <- do.call(rbind.data.frame, parLapply(cl, predictor_files, read_csv))
stopCluster(cl)


# Write and quality control ------------

# Write to file
predictors %>% write_csv("data/predictors_mortality.csv")

# No duplicated patients or admissions
length(unique(predictors$subject_id)) == nrow(predictors)
length(unique(predictors$adm_id)) == nrow(predictors)

# Count and exclude missing chart data
sum(predictors$chart_missing)
predictors %<>% filter(chart_missing == FALSE)

# Count and exclude missing lab data
sum(predictors$lab_missing)
predictors %<>% filter(lab_missing == FALSE)

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


