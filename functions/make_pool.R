make_pool <- function(pooling_variables, target, pool_df)
{
  # Debug
  # pool_df <- present_variable
  # target <- NA_variable[i,]
  # pooling_variables <- pooling_vars[[1]]
  
  # Extract pooling values
  pooling_values <- target %>% 
    select(any_of(pooling_variables)) %>% 
    slice_head()
  
  # Create distance to pool
  dist_matrix <- foreach(p = 1:length(pooling_variables),
          .combine = "cbind") %do%
    {
      if(pooling_variables[p] == "diagnosis")
      {
        pool_dist <- adist(pooling_values[p], pool_df$diagnosis)[1,]
      }else if(pooling_variables[p] == "admission_type")
      {
        admission_class <- pool_df %>% 
          mutate(admission_type = fct_relevel(admission_type, "ELECTIVE",
                                              "URGENT", "EMERGENCY")) %>% 
          select(admission_type) %>% deframe()
        
        pool_dist <- abs(which(c("ELECTIVE", "URGENT", "EMERGENCY") == 
                                  as.character(pooling_values[p])) -
                          as.numeric(admission_class))
      }else
      {
        pool_dist <- abs(pool_df %>% 
          select(any_of(pooling_variables[p])) %>% 
          deframe() - pooling_values[p] %>% deframe())
      }
      
      pool_dist
    }
  
  # Select closest
  closest_rows <- dist_matrix %>% 
    as_tibble() %>% 
    mutate(dist = rowSums(dist_matrix),
           row_id = 1:nrow(dist_matrix)) %>% 
    slice_min(dist) %>% 
    select(row_id) %>% 
    deframe()
  
  # Output pool
  pool_df[closest_rows,]
}