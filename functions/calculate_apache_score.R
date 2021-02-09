calculate_apache_scores <- function(apache_df)
{
  # Debug
  # apache_df <- full_data
  
  # Value extraction
  temp_value <- apache_df$apache_temperature_discharge
  map_value <- apache_df$apache_map_discharge
  pulse_value <- apache_df$apache_pulse_discharge
  respiratory_value <- apache_df$apache_respiratory_discharge
  artpH_value <- apache_df$apache_artpH_discharge
  sodium_value <- apache_df$apache_sodium_discharge
  potassium_value <- apache_df$apache_potassium_discharge
  creatinine_value <- apache_df$apache_creatinine_discharge
  hematocrit_value <- apache_df$apache_hematocrit_discharge
  wbc_value <- apache_df$apache_wbc_discharge
  gcs_value <- apache_df$apache_gcs_discharge
  
  # Temperature-----
  temp_score <- case_when(
    temp_value >= 41 | temp_value <= 29.9 ~ 4L,
    temp_value >= 39 | temp_value <= 31.9 ~ 3L,
    temp_value <= 33.9 ~ 2L,
    temp_value >= 38.5 | temp_value <= 35.9 ~ 1L,
    is.numeric(temp_value) ~ 0L
  )
  
  # Mean arterial pressure----
  map_score <- case_when(
    map_value >= 160 | map_value <= 49 ~ 4L,
    map_value >= 130 ~ 3L,
    map_value >= 110 | map_value <= 69 ~ 2L,
    is.numeric(map_value) ~ 0L
  )
  
  # Heart rate----
  pulse_score <- case_when(
    pulse_value >= 180 | pulse_value <= 39 ~ 4L,
    pulse_value >= 140 | pulse_value <= 54 ~ 3L,
    pulse_value >= 110 | pulse_value <= 69 ~ 2L,
    is.numeric(pulse_value) ~ 0L
  )
  
  # Respiratory rate----
  respiratory_score <-  case_when(
    respiratory_value >= 50 | respiratory_value <= 5 ~ 4L,
    respiratory_value >= 35 ~ 3L,
    respiratory_value <= 9 ~ 2L,
    respiratory_value >= 25 | respiratory_value <= 11 ~ 1L,
    is.numeric(respiratory_value) ~ 0L
  )
  
  # Arterial pH----
  artpH_score <-  case_when(
    artpH_value >= 7.7 | artpH_value < 7.15 ~ 4L,
    artpH_value >= 7.6 | artpH_value <= 7.24 ~ 3L,
    artpH_value <= 7.32 ~ 2L,
    artpH_value >= 7.5 ~ 1L,
    is.numeric(artpH_value) ~ 0L
  )
  
  # Sodium----
  sodium_score <-  case_when(
    sodium_value >= 180 | sodium_value <= 110 ~ 4L,
    sodium_value >= 160 | sodium_value <= 119 ~ 3L,
    sodium_value >= 155 | sodium_value <= 129 ~ 2L,
    sodium_value >= 150 ~ 1L,
    is.numeric(sodium_value) ~ 0L
  )
  
  # Potassium----
  potassium_score <- case_when(
    potassium_value >= 7 | potassium_value < 2.5 ~ 4L,
    potassium_value >= 6 ~ 3L,
    potassium_value <= 2.9 ~ 2L,
    potassium_value >= 5.5 | potassium_value <= 3.4 ~ 1L,
    is.numeric(potassium_value) ~ 0L
  )
  
  # Creatinine----
  creatinine_score <- case_when(
    creatinine_value >= 3.5 ~ 4L,
    creatinine_value >= 2 ~ 3L,
    creatinine_value >= 1.5 | creatinine_value < 0.6 ~ 2L,
    is.numeric(creatinine_value) ~ 0L
  )
  # Double if ARF
  creatinine_score <- ifelse(apache_df$acute_renal_failure, 
         creatinine_score * 2,
         creatinine_score)
  
  # Hematocrit----
  hematocrit_score <- case_when(
    hematocrit_value >= 60 | hematocrit_value < 20 ~ 4L,
    hematocrit_value >= 50 | hematocrit_value <= 29.9 ~ 2L,
    hematocrit_value >= 46 ~ 1L,
    is.numeric(hematocrit_value) ~ 0L
  )
  
  # White bloob cells----
  wbc_score <- case_when(
    wbc_value >= 40 | wbc_value < 1 ~ 4L,
    wbc_value >= 20 | wbc_value <= 2.9 ~ 2L,
    wbc_value >= 15 ~ 1L,
    is.numeric(wbc_value) ~ 0L
  )
  
  # Glasgow Coma Scale
  gcs_score <- 15 - gcs_value
  
  # Sum and output----
  apache_II_score <- temp_score + map_score + pulse_score +
    respiratory_score + sodium_score + potassium_score + creatinine_score +
    hematocrit_score + wbc_score + gcs_score + apache_df$apache_age_discharge +
    apache_df$apache_oxygenation_discharge + apache_df$apache_chronic_discharge
  
  return(apache_II_score)
}