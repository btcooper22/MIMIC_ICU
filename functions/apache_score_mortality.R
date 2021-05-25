require(purrr)

apache_score <- function(data_in, mort_type = "inunit")
{
  # Temperature
  temperature_scores <- case_when(
    data_in$temperature >= 41 | data_in$temperature <= 29.9 ~ 4,
    data_in$temperature >= 39 | data_in$temperature <= 31.9 ~ 3,
    data_in$temperature <= 33.9 ~ 2,
    data_in$temperature >= 38.5 | data_in$temperature <= 35.9 ~ 1,
    is.na(data_in$temperature) ~ NaN,
    is.numeric(data_in$temperature) ~ 0,
  )
  
  # Mean arterial pressure
  map_values <- (data_in$systolicbp + (2 * data_in$diastolicbp) )/ 3
  map_scores <- case_when(
    map_values >= 160 | map_values <= 49 ~ 4,
    map_values >= 130 ~ 3,
    map_values >= 110 | map_values <= 69 ~ 2,
    is.na(map_values) ~ NaN,
    is.numeric(map_values) ~ 0
  )
  
  # Heart rate
  pulse_score <- case_when(
    data_in$pulse >= 180 | data_in$pulse <= 39 ~ 4,
    data_in$pulse >= 140 | data_in$pulse <= 54 ~ 3,
    data_in$pulse >= 110 | data_in$pulse <= 69 ~ 2,
    is.na(data_in$pulse) ~ NaN,
    is.numeric(data_in$pulse) ~ 0
  )
  
  # Respiratory rate
  respiratory_score <-  case_when(
    data_in$respiratory >= 50 | data_in$respiratory <= 5 ~ 4,
    data_in$respiratory >= 35 ~ 3,
    data_in$respiratory <= 9 ~ 2,
    data_in$respiratory >= 25 | data_in$respiratory <= 11 ~ 1,
    is.na(data_in$respiratory) ~ NaN,
    is.numeric(data_in$respiratory) ~ 0
  )
  
  # Arterial pH
  artpH_score <-  case_when(
    data_in$arterialpH >= 7.7 | data_in$arterialpH < 7.15 ~ 4,
    data_in$arterialpH >= 7.6 | data_in$arterialpH <= 7.24 ~ 3,
    data_in$arterialpH <= 7.32 ~ 2,
    data_in$arterialpH >= 7.5 ~ 1,
    is.na(data_in$arterialpH) ~ NaN,
    is.numeric(data_in$arterialpH) ~ 0
  )
  
  # Sodium
  sodium_score <-  case_when(
    data_in$sodium >= 180 | data_in$sodium <= 110 ~ 4,
    data_in$sodium >= 160 | data_in$sodium <= 119 ~ 3,
    data_in$sodium >= 155 | data_in$sodium <= 129 ~ 2,
    data_in$sodium >= 150 ~ 1,
    is.na(data_in$sodium) ~ NaN,
    is.numeric(data_in$sodium) ~ 0
  )
  
  # Potassium
  potassium_score <- case_when(
    data_in$potassium >= 7 | data_in$potassium < 2.5 ~ 4,
    data_in$potassium >= 6 ~ 3,
    data_in$potassium <= 2.9 ~ 2,
    data_in$potassium >= 5.5 | data_in$potassium <= 3.4 ~ 1,
    is.na(data_in$potassium) ~ NaN,
    is.numeric(data_in$potassium) ~ 0
  )
  
  # Creatinine
  creatinine_score <- case_when(
    data_in$creatinine >= 3.5 ~ 4,
    data_in$creatinine >= 2 ~ 3,
    data_in$creatinine >= 1.5 | data_in$creatinine < 0.6 ~ 2,
    is.na(data_in$creatinine) ~ NaN,
    is.numeric(data_in$creatinine) ~ 0
  )
  
  creatinine_score <- ifelse(data_in$acute_renal_failure,
                             creatinine_score * 2, 
                             creatinine_score)
  
  # Haematocrit
  haematocrit_score <- case_when(
    data_in$haematocrit >= 60 | data_in$haematocrit < 20 ~ 4,
    data_in$haematocrit >= 50 | data_in$haematocrit <= 29.9 ~ 2,
    data_in$haematocrit >= 46 ~ 1,
    is.na(data_in$haematocrit) ~ NaN,
    is.numeric(data_in$haematocrit) ~ 0
  )
  
  # White blood cell count
  wbc_score <- case_when(
    data_in$whitebloodcount >= 40 | data_in$whitebloodcount < 1 ~ 4,
    data_in$whitebloodcount >= 20 | data_in$whitebloodcount <= 2.9 ~ 2,
    data_in$whitebloodcount >= 15 ~ 1,
    is.na(data_in$whitebloodcount) ~ NaN,
    is.numeric(data_in$whitebloodcount) ~ 0
  )
  
  # PaO2
  pao2_score <- case_when(
    data_in$arterialoxygen < 55 ~ 4,
    data_in$arterialoxygen <= 60 ~ 3,
    data_in$arterialoxygen <= 70 ~ 1,
    is.na(data_in$arterialoxygen) ~ NaN,
    is.numeric(data_in$arterialoxygen) ~ 0
  )
  
  # Extract FiO2
  fiO2 <- ifelse(is.na(data_in$fractioninspiredoxygen), 21, 
                 data_in$fractioninspiredoxygen)
  
  # Alveolar-arterial gradient
  elev <-  5
  patm <- 760 * exp(elev / -7000)
  ph2o <- 47 * exp((data_in$temperature - 37) / 18.4)
  fio2_remap <- purrr::map_dbl(fiO2, ~ dplyr::if_else(.x > 1, .x / 100, .x))
  aa_grad <- (fio2_remap * (patm - ph2o) - (data_in$arterialcarbon / 0.8)) - data_in$arterialoxygen
  
  aagrad_score <- case_when(
    aa_grad >= 500 ~ 4,
    aa_grad >= 350 ~ 3,
    aa_grad >= 200 ~ 2,
    is.na(aa_grad) ~ NaN,
    is.numeric(aa_grad) ~ 0
  )
  
  # If FiO2 > 0.5, use aa_grad for oxygenation score, else use PaO2
  oxygenation_score <- ifelse(fiO2 > 50, aagrad_score, pao2_score)
  
  # Bicarbonate
  bicarbonate_score <- case_when(
    data_in$bicarbonate >= 52 | data_in$bicarbonate < 15 ~ 4,
    data_in$bicarbonate >= 41 | data_in$bicarbonate <= 17.9 ~ 3,
    data_in$bicarbonate <= 21.9 ~ 2,
    data_in$bicarbonate >= 32 ~ 1,
    is.na(data_in$bicarbonate) ~ NaN,
    is.numeric(data_in$bicarbonate) ~ 0
  )
  
  # Glasgow coma scale
  gcs_score <- 15 - (data_in$glasgowcomaeye + data_in$glasgowcomamotor +
                       data_in$glasgowcomaverbal) %>% round()
  
  # Age
  age_score <- data_in$age
  
  # Chronic health
  chronic_score <- data_in$chronic
  
  # If arterial blood gasses missing, use bicarbonate instead
  abg_score <- 
    ifelse(data_in$missing_abg, bicarbonate_score,
           oxygenation_score + artpH_score)
  
  # Sum scores
  apache_scores <- temperature_scores + map_scores + pulse_score +
    respiratory_score + abg_score + sodium_score + potassium_score +
    creatinine_score + haematocrit_score + wbc_score + gcs_score +
    age_score + chronic_score
  
  # Extract mortality
  mort_varname <- grep("mort", names(data_in), value = TRUE)
  
  return(data.frame(apache_scores,
                    mortality = data_in %>% select(any_of(mort_varname))))
}