crosstab <- function(varname, .df = patients)
{
  options(dplyr.summarise.inform=F)
  
  # Find
  .df  %>% 
    group_by_at(c("readmission", varname)) %>% 
    summarise(n = n()) %>% 
    mutate(`%` = (n/sum(n) * 100)) %>% 
    mutate(`n (%)` = paste(n, " (", round(`%`,1),
                           "%)", sep = "")) %>% 
    select(any_of(c("readmission", varname, "n (%)"))) %>% 
    pivot_wider(names_from = "readmission",
                values_from = "n (%)") %>% 
    as.data.frame()
}