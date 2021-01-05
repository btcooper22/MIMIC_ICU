event_filter_discharge <- function(.df, discharge_time)
{
  return(
    .df %>% 
      # Calculate time difference from discharge
      mutate(before_discharge = discharge_time - charttime) %>% 
      # Select only events before discharge
      filter(before_discharge > 0)
  )
}