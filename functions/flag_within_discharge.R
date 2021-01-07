flag_within_discharge <- function(.df, window = 48)
{
  # Attempt full window filter
  filt_df <- .df %>% 
    filter(before_discharge <= window)
  
  # If filter worked
  if(nrow(filt_df) > 0)
  {
    # Return filtered df
    return(filt_df)
  }else
  {
    return(
    .df %>% 
      slice_head() %>% 
      mutate(valuenum = NA)
    )
  }
  
}