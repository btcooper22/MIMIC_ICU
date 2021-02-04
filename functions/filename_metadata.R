filename_metadata <- function(filename_in)
{
  metadata <- filename_in %>% 
    str_split("/", simplify = TRUE) %>% .[,3] %>% 
    str_split("_", simplify = TRUE)
  
  split <- metadata[,1] %>% 
    substr(2,4) %>% 
    as.numeric()
  
  n <- metadata[,2] %>% 
    substr(2,4) %>% 
    as.numeric()
}