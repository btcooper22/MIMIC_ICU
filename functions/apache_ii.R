# Modified from https://github.com/bgulbis/icuriskr/blob/master/R/apache2.R
# MIMIC extraction code modified from https://github.com/MIT-LCP/mimic-code/tree/master/concepts/cookbook
require(purrr)

apacheII_score <- function(labs_df, chart_df, admission_time)
{
  # For debug
  labs_df <- labs
  chart_df <- charts
  admission_time <- admit_time
  
  # Helper function: Find closest chart measurement to admission (or up to 12h before)
  filter_closest <- function(.df, time, type = "chart")
  {
    if(type == "chart")
    {
      .df %>% 
        mutate(after_admission = difftime(CHARTTIME, time,
                                          units = "hours"))  %>% 
        filter(after_admission > -12) %>% 
        slice_min(abs(after_admission))
    } else
    {
      .df %>% 
        mutate(after_admission = difftime(charttime, time,
                                          units = "hours")) %>% 
        filter(after_admission > -12) %>% 
        slice_min(abs(after_admission))
    }
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
    filter_closest(admission_time) %>% 
    select(value_degC) %>% deframe()
  
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
    filter_closest(admission_time) %>% 
    # Extract sBP
    select(VALUENUM) %>% deframe()
  
  # Find first diastolic BP measurement
  diastolic_bp <- chart_df %>% 
    filter(ITEMID %in% c(8364, 8368, 8441, 
                         8555, 220051, 220180)) %>% 
    # Find closest
    filter_closest(admission_time) %>% 
    # Extract sBP
    select(VALUENUM) %>% deframe()
  
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
    filter_closest(admission_time) %>% 
    select(VALUENUM) %>% deframe()
  
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
    filter_closest(admission_time) %>% 
    select(VALUENUM) %>% deframe()
  
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
    filter_closest(admission_time) %>% 
    select(VALUENUM) %>% deframe()
  
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
    filter_closest(admission_time, "labs") %>% 
    select(valuenum) %>% deframe()
  
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
    filter_closest(admission_time, "labs") %>% 
    select(valuenum) %>% deframe()
  
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
    filter_closest(admission_time, "labs") %>% 
    select(valuenum) %>% deframe()
  
  # Score
  creatinine_score <- case_when(
    creatinine_value >= 3.5 ~ 4L,
    creatinine_value >= 2 ~ 3L,
    creatinine_value >= 1.5 | creatinine_value < 0.6 ~ 2L,
    is.numeric(creatinine_value) ~ 0L
  )
  
  # Calculate haematocrit (%)-----------------
  
  # Gather all haematocrit measurements
  hematocrit_measurements <- labs_df %>% 
    filter(itemid %in% c(51221, 50810))
  
  # Find closest to ICU admission 
  hematocrit_value <- hematocrit_measurements %>% 
    filter_closest(admission_time, "labs") %>% 
    select(valuenum) %>% deframe()
  
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
    filter_closest(admission_time, "labs") %>% 
    select(valuenum) %>% deframe()
  
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
    filter_closest(admission_time) %>% 
    select(VALUENUM) %>%  deframe()
  
  # Find closest verbal measurement
  gcs_verbal_value <- gcs_verbal_measurements %>% 
    filter_closest(admission_time) %>% 
    select(VALUENUM) %>%  deframe()
  
  # Find closest eye measurement
  gcs_eye_value <- gcs_eye_measurements %>% 
    filter_closest(admission_time) %>% 
    select(VALUENUM) %>%  deframe()
  
  # Sum and score
  gcs_value <- gcs_motor_value + gcs_verbal_value + gcs_eye_value
  gcs_score <- 15 - gcs_value
  
  
  # Calculate age-----------------
  apache2_score.age <- function(x, ...) {
    score <- function(y) {
      dplyr::case_when(
        y >= 75 ~ 6L,
        y >= 65 ~ 5L,
        y >= 55 ~ 3L,
        y >= 45 ~ 2L,
        is.numeric(y) ~ 0L
      )
    }
    
    purrr::map_int(x, score)
  }
  
  # Serum bicarbonate (venous mMol/L, only use if no aaDo2 or Pa02)
  # apache2_score.hco3 <- function(x, ...) {
  #   score <- function(y) {
  #     dplyr::case_when(
  #       y >= 52 | y < 15 ~ 4L,
  #       y >= 41 | y <= 17.9 ~ 3L,
  #       y <= 21.9 ~ 2L,
  #       y >= 32 ~ 1L,
  #       is.numeric(y) ~ 0L
  #     )
  #   }
  #   
  #   purrr::map_int(x, score)
  # }
  
  # Arterial O2 (mmHg)
  apache2_score.pao2 <- function(x, ...) {
    score <- function(y) {
      dplyr::case_when(
        y < 55 ~ 4L,
        y <= 60 ~ 3L,
        y <= 70 ~ 1L,
        is.numeric(y) ~ 0L
      )
    }
    
    purrr::map_int(x, score)
  }
  
  # Alveolar–arterial gradient (mmHg)
  apache2_score.aa_grad <- function(x, ...) {
    
    score <- function(y) {
      dplyr::case_when(
        y >= 500 ~ 4L,
        y >= 350 ~ 3L,
        y >= 200 ~ 2L,
        is.numeric(y) ~ 0L
      )
    }
    
    purrr::map_int(x, score)
  }
  
  # Calculate AA gradient
  aa_gradient <- function(pco2, pao2, fio2 = 21, temp = 37, elev = 0) {
    # pco2: partial pressure of carbon dioxide in arterial blood
    # pao2: partial pressure of oxygen in arterial blood
    # fio2:  Fraction of Inspired Oxygen
    # temp: patient temperature in degrees Celsius
    # elev  elevation above sea level in meters
    
    # Aa DO2 = (FIO2 * (Patm - PH2O) - (PaCO2 / 0.8)) - PaO2
    # Patm = 760 * exp(Elevation / -7000); Houston elevation = 43 feet (13.106 m)
    # PH2O = 47 * exp((Temp - 37) / 18.4)
    
    patm <- 760 * exp(elev / -7000)
    ph2o <- 47 * exp((temp - 37) / 18.4)
    
    fio2 <- purrr::map_dbl(fio2, ~ dplyr::if_else(.x > 1, .x / 100, .x))
    
    (fio2 * (patm - ph2o) - (pco2 / 0.8)) - pao2
  }
  
  
  # Need elective admission + comorbs
  apache2_score.admit <- function(x, ..., comorbidity) {
    score <- function(y, z) {
      if (is.na(z) | z == FALSE) {
        0L
      } else {
        if (is.na(y) | y == "elective") {
          2L
        } else {
          5L
        }
      }
    }
    
    purrr::map2_int(x, comorbidity, score)
  }
}


