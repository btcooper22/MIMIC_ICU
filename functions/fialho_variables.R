fialho_variables <- function(time_of_discharge, chart_df,
                             lab_df)
{
  # Debug
  # time_of_discharge <- discharge_time
  # labs_df <- labs
  # chart_df <- charts
  
  # Function to filter to 24h pre-discharge
  discharge_filter <- function(.df, time = time_of_discharge,
                               type = "chart")
  {
    if(type == "chart")
    {
      .df %>% 
        mutate(before_discharge = difftime(CHARTTIME, time,
                                          units = "hours"))  %>% 
        filter(before_discharge > -24,
               before_discharge <= 0)
    } else
    {
      .df %>% 
        mutate(before_discharge = difftime(charttime, time,
                                          units = "hours")) %>% 
        filter(before_discharge > -24,
               before_discharge <= 0)
    }
  }
  
  # Heart rate----------
  
  # Gather all pulse measurements
  pulse_measurements <- chart_df %>% 
    filter(ITEMID %in% c(211, 220045)) 
  
  # Average all within 24h of discharge
  pulse_value <- pulse_measurements %>% 
    discharge_filter(time_of_discharge) %>% 
    select(VALUENUM) %>%
    deframe() %>% mean(na.rm = TRUE)
  
  # Temperature--------
  
  # Gather all temperature measurements
  temp_measurements <- chart_df %>% 
    filter(ITEMID %in% c(676, 677, 678, 679, 223761, 223762)) %>% 
    # Convert any F to C
    mutate(value_degC = ifelse(ITEMID %in% c(678, 679, 223761),
                               (VALUENUM - 32) * 5/9,
                               VALUENUM)
    )
  
  # Average all within 24h of discharge
  temp_value <- temp_measurements %>% 
    discharge_filter(time_of_discharge) %>% 
    select(value_degC) %>%
    deframe() %>% mean(na.rm = TRUE)
  
  # SpO2----------------
  
  # Gather all SpO2 measurements
  spo2_measurements <- chart_df %>% 
    filter(ITEMID %in% c(646, 220277)) 
  
  # Average all within 24h of discharge
  spo2_value <- spo2_measurements %>% 
    discharge_filter(time_of_discharge) %>% 
    select(VALUENUM) %>%
    deframe() %>% mean(na.rm = TRUE)
  
  # Non-invasive arterial blood pressure------------
  
  # Gather all art. BP measurements
  BP_measurements <- chart_df %>% 
    filter(ITEMID %in% c(52, 6702, 6702, 6927,
                         220052)) 
  
  # Average all within 24h of discharge
  BP_value <- BP_measurements %>% 
    discharge_filter(time_of_discharge) %>% 
    select(VALUENUM) %>%
    deframe() %>% mean(na.rm = TRUE)
  
  # Platelets-----------
  
  # Gather all platelet measures
  platelet_measurements <- labs_df %>% 
    filter(itemid %in% c(51265)) 
  
  # Average all within 24h of discharge
  platelet_value <- platelet_measurements %>% 
    discharge_filter(time_of_discharge,
                     type = "labs") %>% 
    select(valuenum) %>%
    deframe() %>% mean(na.rm = TRUE)
  
  # Lactic acid------------
  
  # Gather all lactate measures
  lactate_measurements <- labs_df %>% 
    filter(itemid %in% c(50813)) 
  
  # Average all within 24h of discharge
  lactate_value <- lactate_measurements %>% 
    discharge_filter(time_of_discharge,
                     type = "labs") %>% 
    select(valuenum) %>%
    deframe() %>% mean(na.rm = TRUE)
  
  # Convert to mg/dl
  lactate_value <- lactate_value * 0.111
  
  # Output------------
  return(
    c(pulse_value,
      temp_value,
      spo2_value,
      BP_value,
      platelet_value,
      lactate_value)
  )
}