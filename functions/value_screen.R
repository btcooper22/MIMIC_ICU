value_screen <- function(.df, variable_name,
                         safe_range)
{
  .df %>% 
    select(any_of(c(paste(c("a", "d"), variable_name,
                          sep = "_"),
                    "row_id"))) %>% 
    na.omit() %>% 
    pivot_longer(1:2) %>% 
    filter(value < safe_range[1] |
             value > safe_range[2]) %>% 
    return()
}