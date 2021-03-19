# Modified from https://github.com/bgulbis/icuriskr/blob/master/R/apache2.R
# MIMIC extraction code modified from https://github.com/MIT-LCP/mimic-code/tree/master/concepts/cookbook
require(purrr)

apacheII_score <- function(labs_df, chart_df, patient_df, admission_time,
                           elect_admit, prev_diagnoses, arf)
{
  # For debug
  # labs_df <- labs
  # chart_df <- charts
  # patient_df <- mimic_preproc$stays %>%
  #   filter(subject_id == subj,
  #          hadm_id == adm)
  # admission_time <- admit_time
  # elect_admit <- elective_admission
  # prev_diagnoses <- history
  # arf <- acute_renal_failure
  
  # Helper function: Find measurements within 24h of admission (or before, if not availble)
  filter_window <- function(.df, time, type = "chart")
  {
    # Filter within 24 of admission
      if(type == "chart")
      {
        out_df <- .df %>% 
          mutate(window = difftime(CHARTTIME, time,
                                   units = "hours"))  %>% 
          filter(window < 24 & 
                   window >= 0)
      } else
      {
        out_df <- .df %>% 
          mutate(window = difftime(charttime, time,
                                   units = "hours"))  %>% 
          filter(window < 24& 
                   window >= 0)
      }
    
    # Take nearest pre-admission if none available
    if(nrow(out_df) == 0)
    {
      if(type == "chart")
      {
        out_df <- .df %>% 
          mutate(window = difftime(CHARTTIME, time,
                                   units = "hours"))  %>% 
          filter(window < 0) %>% 
          slice_min(window)
      }else
      {
        out_df <- .df %>% 
          mutate(window = difftime(charttime, time,
                                   units = "hours"))  %>% 
          filter(window < 0) %>% 
          slice_min(window)
      }
    }
    
    # Return NA if none available at all
    if(nrow(out_df) == 0)
    {
      ifelse(type == "chart",
             return(data.frame(VALUENUM = NA)),
             return(data.frame(valuenum = NA)))
      
    }
    return(out_df)
  }
  
  # Filter out implausible values
  implausability_filter <- function(val, good_range)
  {
    val[val < good_range[1]] <- NA
    val[val > good_range[2]] <- NA
    return(val)
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
  
  # Calculate temperature---------
  
  # Gather all temperature measurements
  temp_measurements <- chart_df %>% 
    filter(ITEMID %in% c(676, 677, 678, 679, 223761, 223762)) %>% 
    # Convert any F to C
    mutate(VALUENUM = ifelse(ITEMID %in% c(678, 679, 223761),
                               (VALUENUM - 32) * 5/9,
                               VALUENUM)
           )
  
  # Find closest to ICU admission 
  temp_value <- temp_measurements %>% 
    filter_window(admission_time) %>% 
    select(VALUENUM) %>% deframe()
  
  temp_value %<>%
    implausability_filter(c(30,40)) %>% 
    worst_value(37.2)
  
  # Score
  temp_score <- case_when(
    temp_value >= 41 | temp_value <= 29.9 ~ 4L,
    temp_value >= 39 | temp_value <= 31.9 ~ 3L,
    temp_value <= 33.9 ~ 2L,
    temp_value >= 38.5 | temp_value <= 35.9 ~ 1L,
    is.numeric(temp_value) ~ 0L
  )
  
  # Calculate mean arterial pressure----------
  
  # Find first systolic BP measurement
  systolic_bp <- chart_df %>% 
    filter(ITEMID %in% c(6, 51, 455, 6701, 220050, 220179)) %>% 
    # Find closest
    filter_window(admission_time) %>% 
    # Extract sBP
    select(VALUENUM) %>% deframe()
  
  # Find first diastolic BP measurement
  diastolic_bp <- chart_df %>% 
    filter(ITEMID %in% c(8364, 8368, 8441, 
                         8555, 220051, 220180)) %>% 
    # Find closest
    filter_window(admission_time) %>% 
    # Extract sBP
    select(VALUENUM) %>% deframe()
  
  # Find worst
  systolic_bp %<>% 
    implausability_filter(c(70,220)) %>% 
    worst_value("max")
  
  diastolic_bp %<>%
    implausability_filter(c(40,140)) %>% 
    worst_value("max")
  
  # Calculate MAP and score
  map_value <- (systolic_bp + (2 * diastolic_bp) )/ 3
  map_score <- case_when(
    map_value >= 160 | map_value <= 49 ~ 4L,
    map_value >= 130 ~ 3L,
    map_value >= 110 | map_value <= 69 ~ 2L,
    is.numeric(map_value) ~ 0L
  )
  
  # Calculate heart rate---------------
  
  # Gather all pulse measurements
  pulse_measurements <- chart_df %>% 
    filter(ITEMID %in% c(211, 220045)) 
  
  # Find closest to ICU admission 
  pulse_value <- pulse_measurements %>% 
    filter_window(admission_time) %>% 
    select(VALUENUM) %>% deframe()
  
  # Find worst
  pulse_value %<>%
    implausability_filter(c(20,200)) %>% 
    worst_value(90)
  
  # Score
  pulse_score <- case_when(
    pulse_value >= 180 | pulse_value <= 39 ~ 4L,
    pulse_value >= 140 | pulse_value <= 54 ~ 3L,
    pulse_value >= 110 | pulse_value <= 69 ~ 2L,
    is.numeric(pulse_value) ~ 0L
  )
  
  # Calculate respiratory rate-----------
  
  # Gather all RR measurements
  respiratory_measurements <- chart_df %>% 
    filter(ITEMID %in% c(614, 615, 618, 219, 619, 653, 1884,
                         8113, 1635, 3603, 224688, 224689,
                         224690, 220210))
  
  # Find closest to ICU admission 
  respiratory_value <- respiratory_measurements %>% 
    filter_window(admission_time) %>% 
    select(VALUENUM) %>% deframe()
  
  # Find worst
  respiratory_value %<>%
    implausability_filter(c(0,60)) %>% 
    worst_value("max")
  
  # Score
  respiratory_score <-  case_when(
    respiratory_value >= 50 | respiratory_value <= 5 ~ 4L,
    respiratory_value >= 35 ~ 3L,
    respiratory_value <= 9 ~ 2L,
    respiratory_value >= 25 | respiratory_value <= 11 ~ 1L,
    is.numeric(respiratory_value) ~ 0L
  )
  
  # Calculate arterial pH---------

  # Gather all arterial pH measurements
  artpH_measurements <- chart_df %>% 
    filter(ITEMID %in% c(1126, 780,4753, 223830))
  
  # Find closest to ICU admission 
  artpH_value <- artpH_measurements %>% 
    filter_window(admission_time) %>% 
    select(VALUENUM) %>% deframe()
  
  # Find worst
  artpH_value %<>% 
    implausability_filter(c(7,8)) %>% 
    worst_value(7.55)
  
  # Score
  artpH_score <-  case_when(
    artpH_value >= 7.7 | artpH_value < 7.15 ~ 4L,
    artpH_value >= 7.6 | artpH_value <= 7.24 ~ 3L,
    artpH_value <= 7.32 ~ 2L,
    artpH_value >= 7.5 ~ 1L,
    is.numeric(artpH_value) ~ 0L
  )
  
  # Calculate serum sodium (mMol/L)------------

  # Gather all serum sodium measurements
  sodium_measurements <- labs_df %>% 
    filter(itemid %in% c(50824, 50983))
  
  # Find closest to ICU admission 
  sodium_value <- sodium_measurements %>% 
    filter_window(admission_time, "labs") %>% 
    select(valuenum) %>% deframe()
  
  # Find worst
  sodium_value %<>% 
    implausability_filter(c(100,160)) %>% 
    worst_value(140)
  
  # Score
  sodium_score <-  case_when(
    sodium_value >= 180 | sodium_value <= 110 ~ 4L,
    sodium_value >= 160 | sodium_value <= 119 ~ 3L,
    sodium_value >= 155 | sodium_value <= 129 ~ 2L,
    sodium_value >= 150 ~ 1L,
    is.numeric(sodium_value) ~ 0L
  )
  
  # Calculate serum potassium (mMol/L)---------
  
  # Gather all serum potassium measurements
  potassium_measurements <- labs_df %>% 
    filter(itemid %in% c(50822, 50971))
  
  # Find closest to ICU admission 
  potassium_value <- potassium_measurements %>% 
    filter_window(admission_time, "labs") %>% 
    select(valuenum) %>% deframe()
  
  # Find worst
  potassium_value %<>%
    implausability_filter(c(0,8.5)) %>% 
    worst_value(4.5)
  
  # Score
  potassium_score <- case_when(
    potassium_value >= 7 | potassium_value < 2.5 ~ 4L,
    potassium_value >= 6 ~ 3L,
    potassium_value <= 2.9 ~ 2L,
    potassium_value >= 5.5 | potassium_value <= 3.4 ~ 1L,
    is.numeric(potassium_value) ~ 0L
  )

  # Calculate serum creatinine (mg/dL)------------
  
  # Gather all serum creatinine measurements
  creatinine_measurements <- labs_df %>% 
    filter(itemid %in% c(50912))
  
  # Find closest to ICU admission 
  creatinine_value <- creatinine_measurements %>% 
    filter_window(admission_time, "labs") %>% 
    select(valuenum) %>% deframe()
  
  # Find worst
  creatinine_value %<>% 
    implausability_filter(c(0,6)) %>% 
    worst_value("max")
  
  # Score
  creatinine_score <- case_when(
    creatinine_value >= 3.5 ~ 4L,
    creatinine_value >= 2 ~ 3L,
    creatinine_value >= 1.5 | creatinine_value < 0.6 ~ 2L,
    is.numeric(creatinine_value) ~ 0L
  )
  
  # Double for acute renal failure
  if(arf)
  {
    creatinine_score <- creatinine_score * 2
  }
  
  # Calculate haematocrit (%)-----------------
  
  # Gather all haematocrit measurements
  hematocrit_measurements <- labs_df %>% 
    filter(itemid %in% c(51221, 50810))
  
  # Find closest to ICU admission 
  hematocrit_value <- hematocrit_measurements %>% 
    filter_window(admission_time, "labs") %>% 
    select(valuenum) %>% deframe()
  
  # Find worst
  hematocrit_value %<>%
    implausability_filter(c(15,60)) %>% 
    worst_value(44)
  
  # Score
  hematocrit_score <- case_when(
    hematocrit_value >= 60 | hematocrit_value < 20 ~ 4L,
    hematocrit_value >= 50 | hematocrit_value <= 29.9 ~ 2L,
    hematocrit_value >= 46 ~ 1L,
    is.numeric(hematocrit_value) ~ 0L
  )
  
  # Calculate white blood count (1000s)---------
  
  # Gather all WBC measurements
  wbc_measurements <- labs_df %>% 
    filter(itemid %in% c(51300, 51301))
  
  # Find closest to ICU admission 
  wbc_value <- wbc_measurements %>% 
    filter_window(admission_time, "labs") %>% 
    select(valuenum) %>% deframe()

  # Find worst
  wbc_value %<>%
    implausability_filter(c(0,50)) %>% 
    worst_value("max")
  
  # Score
  wbc_score <- case_when(
    wbc_value >= 40 | wbc_value < 1 ~ 4L,
    wbc_value >= 20 | wbc_value <= 2.9 ~ 2L,
    wbc_value >= 15 ~ 1L,
    is.numeric(wbc_value) ~ 0L
  )
  
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
    filter_window(admission_time) %>% 
    select(VALUENUM) %>%  deframe()
  gcs_motor_value %<>% worst_value("min")
  
  # Find closest verbal measurement
  gcs_verbal_value <- gcs_verbal_measurements %>% 
    filter_window(admission_time) %>% 
    select(VALUENUM) %>%  deframe()
  gcs_verbal_value %<>% worst_value("min")
  
  # Find closest eye measurement
  gcs_eye_value <- gcs_eye_measurements %>% 
    filter_window(admission_time) %>% 
    select(VALUENUM) %>%  deframe()
  gcs_eye_value %<>% worst_value("min")
  
  # Sum and score
  gcs_value <- gcs_motor_value + gcs_verbal_value + gcs_eye_value
  gcs_score <- 15 - gcs_value
  
  
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
  
  # Set to 21 if measurements missing
  if(nrow(fio2_measurements) == 0)
  {
    fio2_value <- 21
  }else
  {
    # Find closest to ICU admission 
    fio2_value <- fio2_measurements %>% 
      filter_window(admission_time) %>% 
      select(VALUENUM) %>% deframe()
    
    # Find worst
    fio2_value %<>% 
      implausability_filter(c(21,100)) %>% 
      worst_value("max")
  }
  
  # Calculate arterial o2
  pao2_measurements <- chart_df %>% 
    filter(ITEMID %in%  c(4203, 3785, 3837, 3828, 220224, 779))
  
  # Find closest to ICU admission 
  pao2_value <- pao2_measurements %>% 
    filter_window(admission_time) %>% 
    select(VALUENUM) %>% deframe()

  # Find worst
  pao2_value %<>% 
    implausability_filter(c(20,250)) %>% 
    worst_value("min")
  
  # Calculate Pco2
  paco2_measurements <- chart_df %>% 
    filter(ITEMID %in%  c(220235, 777, 778, 3784, 3835, 3836))
  
  # Find closest to ICU admission 
  paco2_value <- paco2_measurements %>% 
    filter_window(admission_time) %>% 
    select(VALUENUM) %>% deframe()
  
  # Find worst
  paco2_value %<>%
    implausability_filter(c(8,200)) %>% 
    worst_value(40)
  
  # If fio2 > 0.5 use aa_grad
  if(fio2_value > 50 & !is.na(fio2_value))
  {
    # Calculate AA gradient
    elev <-  5
    patm <- 760 * exp(elev / -7000)
    ph2o <- 47 * exp((temp_value - 37) / 18.4)
    fio2 <- purrr::map_dbl(fio2_value, ~ dplyr::if_else(.x > 1, .x / 100, .x))
    aa_grad <- (fio2 * (patm - ph2o) - (paco2_value / 0.8)) - pao2_value
    
    # Score AA gradient
    oxygenation_score <- case_when(
      aa_grad >= 500 ~ 4L,
      aa_grad >= 350 ~ 3L,
      aa_grad >= 200 ~ 2L,
      is.numeric(aa_grad) ~ 0L
    )
    
  }else # Else use Pao2
  {
    oxygenation_score <- case_when(
      pao2_value < 55 ~ 4L,
      pao2_value <= 60 ~ 3L,
      pao2_value <= 70 ~ 1L,
      is.numeric(pao2_value) ~ 0L
    )
  }
  
  # Serum bicarbonate (venous mMol/L, only use if no arterial blood gases)-----------
  if(is_empty(oxygenation_score) & is_empty(artpH_score))
  {
    # Gather all bicarbonate measurements
    bicarbonate_measurements <- labs_df %>% 
      filter(itemid %in% c(50803, 50882))
    
    # Find closest to ICU admission 
    bicarbonate_value <- bicarbonate_measurements %>% 
      filter_window(admission_time, "labs") %>% 
      select(valuenum) %>% deframe()
    
    # Find worst
    bicarbonate_value %<>%
      implausability_filter(c(0,60)) %>% 
      worst_value(27)
    
    # Score
    bicarbonate_score <- case_when(
            bicarbonate_value >= 52 | bicarbonate_value < 15 ~ 4L,
            bicarbonate_value >= 41 | bicarbonate_value <= 17.9 ~ 3L,
            bicarbonate_value <= 21.9 ~ 2L,
            bicarbonate_value >= 32 ~ 1L,
            is.numeric(bicarbonate_value) ~ 0L
          )
  }else
  {
    bicarbonate_value <- NA
    bicarbonate_score <- NA
  }
  
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
  
  # Create score vector
  full_scores <- c("temperature" = missing_to_na(temp_score),
                   "map" = missing_to_na(map_score),
                   "pulse" = missing_to_na(pulse_score),
                   "respiratory" = missing_to_na(respiratory_score),
                   "oxygen" = missing_to_na(oxygenation_score),
                   "arterialpH" = missing_to_na(artpH_score),
                   "sodium" = missing_to_na(sodium_score),
                   "potassium" = missing_to_na(potassium_score),
                   "creatinine" = missing_to_na(creatinine_score),
                   "haematocrit" = missing_to_na(hematocrit_score),
                   "whitebloodcount" = missing_to_na(wbc_score),
                   "glasgowcoma" = missing_to_na(gcs_score),
                   "age" = missing_to_na(age_score),
                   "chronic" = missing_to_na(chronic_score),
                   "bicarbonate" = missing_to_na(bicarbonate_score))
  
  # Create value vector
  full_values <- c("temperature" = missing_to_na(temp_value),
                   "map" = missing_to_na(map_value),
                   "pulse" = missing_to_na(pulse_value),
                   "respiratory" = missing_to_na(respiratory_value),
                   "arterialpH" = missing_to_na(artpH_value),
                   "sodium" = missing_to_na(sodium_value),
                   "potassium" = missing_to_na(potassium_value),
                   "creatinine" = missing_to_na(creatinine_value),
                   "haematocrit" = missing_to_na(hematocrit_value),
                   "whitebloodcount" = missing_to_na(wbc_value),
                   "glasgowcoma" = missing_to_na(gcs_value),
                   "bicarbonate" = missing_to_na(bicarbonate_value),
                   "fractioninspiredoxygen" = missing_to_na(fio2_value),
                   "arterialoxygen" = missing_to_na(pao2_value),
                   "arterialcarbon" = missing_to_na(paco2_value))
  return(list(full_scores = full_scores, 
              full_values = full_values))
}


