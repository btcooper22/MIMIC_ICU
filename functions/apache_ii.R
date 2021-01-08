# Modified from https://github.com/bgulbis/icuriskr/blob/master/R/apache2.R
# MIMIC extraction code modified from https://github.com/MIT-LCP/mimic-code/tree/master/concepts/cookbook
require(purrr)

apacheII_score <- function(labs_df, chart_df, admission_time)
{
  # For debug
  labs_df <- labs
  chart_df <- charts
  admission_time <- admit_time
  
  # Helper function: Find closest chart measurement to admission
  filter_closest <- function(.df, time)
  {
    .df %>% 
      mutate(after_admission = difftime(CHARTTIME, time,
                                        units = "hours")) %>% 
      filter(after_admission > 0) %>% 
      slice_min(after_admission)
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
  respiratory_measurements <- charts %>% 
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
  
  # Arterial pH
  apache2_score.ph <- function(x, ...) {
    score <- function(y) {
      dplyr::case_when(
        y >= 7.7 | y < 7.15 ~ 4L,
        y >= 7.6 | y <= 7.24 ~ 3L,
        y <= 7.32 ~ 2L,
        y >= 7.5 ~ 1L,
        is.numeric(y) ~ 0L
      )
    }
    
    purrr::map_int(x, score)
  }
  
  # Serum sodium (mMol/L)
  apache2_score.sodium <- function(x, ...) {
    score <- function(y) {
      dplyr::case_when(
        y >= 180 | y <= 110 ~ 4L,
        y >= 160 | y <= 119 ~ 3L,
        y >= 155 | y <= 129 ~ 2L,
        y >= 150 ~ 1L,
        is.numeric(y) ~ 0L
      )
    }
    
    purrr::map_int(x, score)
  }
  
  # Serum potassium (mMol/L)
  apache2_score.potassium <- function(x, ...) {
    score <- function(y) {
      dplyr::case_when(
        y >= 7 | y < 2.5 ~ 4L,
        y >= 6 ~ 3L,
        y <= 2.9 ~ 2L,
        y >= 5.5 | y <= 3.4 ~ 1L,
        is.numeric(y) ~ 0L
      )
    }
    
    purrr::map_int(x, score)
  }
  
  # Serum Creatinine (mg/dL)
  apache2_score.scr <- function(x, ...) {
    score <- function(y) {
      dplyr::case_when(
        y >= 3.5 ~ 4L,
        y >= 2 ~ 3L,
        y >= 1.5 | y < 0.6 ~ 2L,
        is.numeric(y) ~ 0L
      )
    }
    
    purrr::map_int(x, score)
  }
  
  # Haematocrit (%)
  apache2_score.hct <- function(x, ...) {
    score <- function(y) {
      dplyr::case_when(
        y >= 60 | y < 20 ~ 4L,
        y >= 50 | y <= 29.9 ~ 2L,
        y >= 46 ~ 1L,
        is.numeric(y) ~ 0L
      )
    }
    
    purrr::map_int(x, score)
  }
  
  # White blood count (1000s)
  apache2_score.wbc <- function(x, ...) {
    score <- function(y) {
      dplyr::case_when(
        y >= 40 | y < 1 ~ 4L,
        y >= 20 | y <= 2.9 ~ 2L,
        y >= 15 ~ 1L,
        is.numeric(y) ~ 0L
      )
    }
    
    purrr::map_int(x, score)
  }
  
  # Glasgow coma score (units)
  apache2_score.gcs <- function(x, ...) {
    purrr::map_int(x, ~ 15L - as.integer(.x))
  }
  
  # Serum bicarbonate (venous mMol/L, only use if no aaDo2 or Pa02)
  apache2_score.hco3 <- function(x, ...) {
    score <- function(y) {
      dplyr::case_when(
        y >= 52 | y < 15 ~ 4L,
        y >= 41 | y <= 17.9 ~ 3L,
        y <= 21.9 ~ 2L,
        y >= 32 ~ 1L,
        is.numeric(y) ~ 0L
      )
    }
    
    purrr::map_int(x, score)
  }
  
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
  
  # Age
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


