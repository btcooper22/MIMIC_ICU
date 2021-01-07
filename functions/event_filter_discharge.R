event_filter_discharge <- function(.df, discharge_time)
{
  return(
    .df %>% 
      # Calculate time difference from discharge
      mutate(before_discharge = difftime(discharge_time, charttime,
                                         units = "hours")) %>% 
      # Select only events before discharge
      filter(before_discharge > 0)
  )
}