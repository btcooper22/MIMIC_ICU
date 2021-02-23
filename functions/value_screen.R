value_screen <- function(.df, variable_name,
                         safe_range, return_safe = FALSE)
{
  # Select variable of interest
  pivoted_df <- .df %>% 
    select(any_of(c(paste(c("a", "d"), variable_name,
                          sep = "_"),
                    "row_id"))) %>% 
    pivot_longer(1:2) %>% 
    na.omit()
  
  
  # Perform filtering
  if(return_safe)
  {
    filt_df <- pivoted_df %>% 
      filter(value > safe_range[1] &
               value < safe_range[2])
  }else
  {
    filt_df <- pivoted_df %>% 
      filter(value < safe_range[1] |
               value > safe_range[2])
  }
  return(filt_df)
}