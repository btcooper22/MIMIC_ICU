# Modified from https://github.com/bgulbis/icuriskr/blob/master/R/apache2.R
# MIMIC extraction code modified from https://github.com/MIT-LCP/mimic-code/tree/master/concepts/cookbook
require(purrr)

apacheII_score <- function(labs_df, chart_df, patient_df, event_time,
                           elect_admit, prev_diagnoses, arf, event_type)
{
  # # For debug
  # labs_df <- labs
  # chart_df <- charts
  # patient_df <- mimic_preproc$stays %>%
  #   filter(subject_id == subj,
  #          hadm_id == adm)
  # event_time <- admit_time
  # elect_admit <- elective_admission
  # prev_diagnoses <- history
  # arf <- acute_renal_failure
  # event_type <-  "admission"
  
  # event_time <- discharge_time
  # event_type <-  "discharge"
  
  # Helper function: Find closest chart measurement to admission
  filter_window <- function(.df, time, type = "chart", event = event_type)
  {
    # Filter within 24 of admission
    if(event == "admission")
    {
      if(type == "chart")
      {
        .df %>% 
          mutate(window = difftime(CHARTTIME, time,
                                            units = "hours"))  %>% 
          filter(window < 24)
      } else
      {
        .df %>% 
          mutate(window = difftime(charttime, time,
                                            units = "hours"))  %>% 
          filter(window < 24)
      }
    }else
    {
      if(type == "chart")
      {
        .df %>% 
          mutate(window = difftime(CHARTTIME, time,
                                            units = "hours"))  %>% 
          filter(window < 0, window > -24)
      } else
      {
        .df %>% 
          mutate(window = difftime(charttime, time,
                                            units = "hours"))  %>% 
          filter(window < 0, window > -24)
      }
    }

  }
  
  # Take the worst value inside the window
  worst_value <- function(val, bad_direction)
  {
    if(all(is.na(val)))
    {
      return(NA)
    }else 
    {
      if(bad_direction == "max")
      {
        return(max(val, na.rm = TRUE))
      }else if(bad_direction == "min")
      {
        return(min(val, na.rm = TRUE))
      } else
      {
        return(val[which.max(abs(val - bad_direction))])
      }
    }
  }
  
  # NA value if missing score
  missing_to_na <- function(val)
  {
    if(is_empty(val))
    {
      return(NA)
    }else
    {
      return(val)
    }
  }
  
  # Filter out implausible values
  implausability_filter <- function(val, good_range)
  {
    val[val < good_range[1]] <- NA
    val[val > good_range[2]] <- NA
    return(val)
  }
  
  # Calculate temperature---------
  
  # Gather all temperature measurements
  temp_measurements <- chart_df %>% 
    filter(ITEMID %in% c(676, 677, 678, 679, 223761, 223762)) %>% 
    # Convert any F to C
    mutate(value_degC = ifelse(ITEMID %in% c(678, 679, 223761),
                               (VALUENUM - 32) * 5/9,
                               VALUENUM)
           )

  # Find closest to ICU admission 
  temp_value <- temp_measurements %>% 
    filter_window(event_time) %>% 
    select(value_degC) %>% deframe()
  
  temp_value %<>%
    implausability_filter(c(30,40)) %>% 
    worst_value(37.2)
  
  # Calculate mean arterial pressure----------
  
  # Find first systolic BP measurement
  systolic_bp_value <- chart_df %>% 
    filter(ITEMID %in% c(6, 51, 455, 6701, 220050, 220179)) %>% 
    # Find closest
    filter_window(event_time) %>% 
    # Extract sBP
    select(VALUENUM) %>% deframe()
  
  # Find first diastolic BP measurement
  diastolic_bp_value <- chart_df %>% 
    filter(ITEMID %in% c(8364, 8368, 8441, 
                         8555, 220051, 220180)) %>% 
    # Find closest
    filter_window(event_time) %>% 
    # Extract sBP
    select(VALUENUM) %>% deframe()
  
  # Find worst
  systolic_bp_value %<>% 
    implausability_filter(c(70,220)) %>% 
    worst_value("max")
  
  diastolic_bp_value %<>%
    implausability_filter(c(40,140)) %>% 
    worst_value("max")
  
  # Calculate heart rate---------------
  
  # Gather all pulse measurements
  pulse_measurements <- chart_df %>% 
    filter(ITEMID %in% c(211, 220045)) 
  
  # Find closest to ICU admission 
  pulse_value <- pulse_measurements %>% 
    filter_window(event_time) %>% 
    select(VALUENUM) %>% deframe()
  
  # Find worst
  pulse_value %<>%
    implausability_filter(c(20,200)) %>% 
    worst_value(90)
  
  # Calculate respiratory rate-----------
  
  # Gather all RR measurements
  respiratory_measurements <- chart_df %>% 
    filter(ITEMID %in% c(614, 615, 618, 219, 619, 653, 1884,
                         8113, 1635, 3603, 224688, 224689,
                         224690, 220210))
  
  # Find closest to ICU admission 
  respiratory_value <- respiratory_measurements %>% 
    filter_window(event_time) %>% 
    select(VALUENUM) %>% deframe()
  
  # Find worst
  respiratory_value %<>%
    implausability_filter(c(0,60)) %>% 
    worst_value("max")
  
  # Calculate arterial pH---------

  # Gather all arterial pH measurements
  artpH_measurements <- chart_df %>% 
    filter(ITEMID %in% c(1126, 780,4753, 223830))
  
  # Find closest to ICU admission 
  artpH_value <- artpH_measurements %>% 
    filter_window(event_time) %>% 
    select(VALUENUM) %>% deframe()
  
  # Find worst
  artpH_value %<>% 
    implausability_filter(c(7,8)) %>% 
    worst_value(7.55)
  
  # Calculate serum sodium (mMol/L)------------

  # Gather all serum sodium measurements
  sodium_measurements <- labs_df %>% 
    filter(itemid %in% c(50824, 50983))
  
  # Find closest to ICU admission 
  sodium_value <- sodium_measurements %>% 
    filter_window(event_time, "labs") %>% 
    select(valuenum) %>% deframe()
  
  # Find worst
  sodium_value %<>% 
    implausability_filter(c(100,160)) %>% 
    worst_value(140)
  
  # Calculate serum potassium (mMol/L)---------
  
  # Gather all serum potassium measurements
  potassium_measurements <- labs_df %>% 
    filter(itemid %in% c(50822, 50971))
  
  # Find closest to ICU admission 
  potassium_value <- potassium_measurements %>% 
    filter_window(event_time, "labs") %>% 
    select(valuenum) %>% deframe()
  
  # Find worst
  potassium_value %<>%
    implausability_filter(c(0,8.5)) %>% 
    worst_value(4.5)
  
  # Calculate serum creatinine (mg/dL)------------
  
  # Gather all serum creatinine measurements
  creatinine_measurements <- labs_df %>% 
    filter(itemid %in% c(50912))
  
  # Find closest to ICU admission 
  creatinine_value <- creatinine_measurements %>% 
    filter_window(event_time, "labs") %>% 
    select(valuenum) %>% deframe()
  
  # Find worst
  creatinine_value %<>% 
    implausability_filter(c(0,6)) %>% 
    worst_value("max")
  
  # Calculate haematocrit (%)-----------------
  
  # Gather all haematocrit measurements
  hematocrit_measurements <- labs_df %>% 
    filter(itemid %in% c(51221, 50810))
  
  # Find closest to ICU admission 
  hematocrit_value <- hematocrit_measurements %>% 
    filter_window(event_time, "labs") %>% 
    select(valuenum) %>% deframe()
  
  # Find worst
  hematocrit_value %<>%
    implausability_filter(c(15,60)) %>% 
    worst_value(44)
  
  # Calculate white blood count (1000s)---------
  
  # Gather all WBC measurements
  wbc_measurements <- labs_df %>% 
    filter(itemid %in% c(51300, 51301))
  
  # Find closest to ICU admission 
  wbc_value <- wbc_measurements %>% 
    filter_window(event_time, "labs") %>% 
    select(valuenum) %>% deframe()
  
  # Find worst
  wbc_value %<>%
    implausability_filter(c(0,50)) %>% 
    worst_value("max")
  
  # Calculate Glasgow coma score---------
  
  # Gather all GCS motor measurements
  gcs_motor_measurements <- chart_df %>% 
    filter(ITEMID %in% c(454, 223901))
  
  # Gather all GCS verbal measurements
  gcs_verbal_measurements <- chart_df %>% 
    filter(ITEMID %in% c(723, 223900))
  
  # Gather all GCS eye measurements
  gcs_eye_measurements <- chart_df %>% 
    filter(ITEMID %in% c(184, 220739))
  
  # Find closest motor measurement
  gcs_motor_value <- gcs_motor_measurements %>% 
    filter_window(event_time) %>% 
    select(VALUENUM) %>%  deframe()
  gcs_motor_value %<>% worst_value("min")
  
  # Find closest verbal measurement
  gcs_verbal_value <- gcs_verbal_measurements %>% 
    filter_window(event_time) %>% 
    select(VALUENUM) %>%  deframe()
  gcs_verbal_value %<>% worst_value("min")
  
  # Find closest eye measurement
  gcs_eye_value <- gcs_eye_measurements %>% 
    filter_window(event_time) %>% 
    select(VALUENUM) %>%  deframe()
  gcs_eye_value %<>% worst_value("min")
  
  # Calculate age-----------------
  
  # Find patient data
  age_value <- patient_df %>% 
    slice_head() %>%
    select(dob, intime) %>% 
    # Generate age
    mutate(age = (intime - dob) / 365) %>%
    select(age) %>% deframe() %>% 
    as.numeric()
  
  # Score
  age_score <- case_when(
    age_value >= 75 ~ 6L,
    age_value >= 65 ~ 5L,
    age_value >= 55 ~ 3L,
    age_value >= 45 ~ 2L,
    is.numeric(age_value) ~ 0L
  )
  
  # Calculate oxygenation score--------
  
  # Gather all FiO2 measurements
  fio2_measurements <- chart_df %>% 
    filter(ITEMID %in% c(223835, 2981, 3420))
  
  # Find closest to ICU admission 
  fio2_value <- fio2_measurements %>% 
    filter_window(event_time) %>% 
    select(VALUENUM) %>% deframe()
  
  # Find worst
  fio2_value %<>% 
    implausability_filter(c(21,100)) %>% 
    worst_value("max")
  
  # Calculate arterial o2
  pao2_measurements <- chart_df %>% 
    filter(ITEMID %in%  c(4203, 3785, 3837, 3828, 220224, 779))
  
  # Find closest to ICU admission 
  pao2_value <- pao2_measurements %>% 
    filter_window(event_time) %>% 
    select(VALUENUM) %>% deframe()
  
  # Find worst
  pao2_value %<>% 
    implausability_filter(c(20,250)) %>% 
    worst_value("min")
  
  # Calculate Pco2
  paco2_measurements <- chart_df %>% 
    filter(ITEMID %in%  c(779, 220235, 777, 778, 3784, 3835, 3836))
  
  # Find closest to ICU admission 
  paco2_value <- paco2_measurements %>% 
    filter_window(event_time) %>% 
    select(VALUENUM) %>% deframe()
  
  # Find worst
  paco2_value %<>%
    implausability_filter(c(8,200)) %>% 
    worst_value(40)
  
  # Gather all bicarbonate measurements
  bicarbonate_measurements <- labs_df %>% 
    filter(itemid %in% c(50803, 50882))
  
  # Find closest to ICU admission 
  bicarbonate_value <- bicarbonate_measurements %>% 
    filter_window(event_time, "labs") %>% 
    select(valuenum) %>% deframe()
  
  # Find worst
  bicarbonate_value %<>%
    implausability_filter(c(0,60)) %>% 
    worst_value(27)

  # Calculate chronic health score----------
  
  # Load comborbitity icd9 codes (from icuriskr)
  load("functions/comorbidity_map_icd9.rda")
  
  # Unlist codes
  comorbidity_icd9_codes <- unlist(comorbidity_map_icd9)
  
  # Search for any codes in patient history
  comorbs <- any(comorbidity_icd9_codes %in% prev_diagnoses$icd9_code)
  
  # Score
  if (comorbs == FALSE) {
    chronic_score <- 0L
  } else {
    if (elect_admit == TRUE) {
      chronic_score <- 2L
    } else {
      chronic_score <- 5L
    }
  }
  
  # Return full vector of values
  full_values <- c("temperature" = missing_to_na(temp_value),
                   "systolicbp" = missing_to_na(systolic_bp_value),
                   "diastolicbp" = missing_to_na(diastolic_bp_value),
                   "pulse" = missing_to_na(pulse_value),
                   "respiratory" = missing_to_na(respiratory_value),
                   "arterialpH" = missing_to_na(artpH_value),
                   "sodium" = missing_to_na(sodium_value),
                   "potassium" = missing_to_na(potassium_value),
                   "creatinine" = missing_to_na(creatinine_value),
                   "haematocrit" = missing_to_na(hematocrit_value),
                   "whitebloodcount" = missing_to_na(wbc_value),
                   "glasgowcomaeye" = missing_to_na(gcs_eye_value),
                   "glasgowcomamotor" = missing_to_na(gcs_motor_value),
                   "glasgowcomaverbal" = missing_to_na(gcs_verbal_value),
                   "age" = missing_to_na(age_score),
                   "chronic" = missing_to_na(chronic_score),
                   "bicarbonate" = missing_to_na(bicarbonate_value),
                   "fractioninspiredoxygen" = missing_to_na(fio2_value),
                   "arterialoxygen" = missing_to_na(pao2_value),
                   "arterialcarbon" = missing_to_na(paco2_value))
  return(full_values)
}


