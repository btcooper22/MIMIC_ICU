# Packages
require(dplyr)
require(magrittr)
require(readr)
require(foreach)
require(mice)
require(stringr)
require(randomForest)
require(tidyr)

# Functions

# Main loop through files-------

# Load data
source("functions/data_loader_splitter.R")

# Set methods
mice_methods <- c("pmm", "norm.boot", "norm.predict", "midastouch",
                  "sample", "norm", "norm.nob",
                  "quadratic", "ri")


# Loop through methods
imputed_dfs <- foreach(j = 1:length(mice_methods)) %do%
  {
    print(j)
    
    # Perform imputation
    MICE_results <- apache_scores %>% 
      mice(method = mice_methods[j])
    
    # Return imputed data
    list(method = mice_methods[j],
         df1 = complete(MICE_results, 1),
         df2 = complete(MICE_results, 2),
         df3 = complete(MICE_results, 3),
         df4 = complete(MICE_results, 4),
         df5 = complete(MICE_results, 5))
  }

# Disentangle
disentangled_df <- foreach(k = 1:length(imputed_dfs)) %do%
  {
    # Loop through datasets within a method
    value_matrix <- foreach(d = 1:5,
                           .combine = "cbind") %do%
      {
        pivoted_df <- imputed_dfs[[k]][[d+1]] %>% 
          as.data.frame() %>% 
          mutate(id = 1:nrow(imputed_dfs[[k]][[d+1]])) %>% 
          pivot_longer(1:17)
        
        if(d == 1)
        {
          pivoted_df
        }else
        {
          pivoted_df[,3]
        }
      }
    
    # Calculate mean
    means_df <- data.frame(id = value_matrix[,1],
                           var = value_matrix[,2],
               value = rowMeans(value_matrix[,3:7]))
    
    # Re-pivot
    means_df %>% 
      pivot_wider(names_from = var,
                  values_from = value) %>% 
      select(-id)
  }
names(disentangled_df) <- paste("MICE", mice_methods, sep = "_")

write_rds(disentangled_df, "data/impute_discharge/MICE.RDS")
