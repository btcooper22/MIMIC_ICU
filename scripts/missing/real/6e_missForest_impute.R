# Packages
require(dplyr)
require(magrittr)
require(readr)
require(foreach)
require(missForest)
require(stringr)
require(tidyr)

# Load data-------
source("functions/data_loader_splitter.R")

# Set mtry
mtry_range <- 2:15

# Set up parallel
source("functions/parallel_setup.R")
parallel_setup(14)

# Loop through forest configuration
imputed_dfs <- foreach(j = 1:length(mtry_range)) %do%
  {
    print(j)
    
    # Perform random forest
    ptm <- proc.time()
    forest_results <- apache_scores %>%
      as.data.frame() %>% 
      missForest(mtry = mtry_range[j],
                 parallelize = "variables")
    proc.time() - ptm
    
    # Return imputed data
    forest_results$ximp
  }
names(imputed_dfs) <- paste("RF", mtry_range, sep = "_")

write_rds(imputed_dfs, "data/impute_discharge/forest.RDS",
          compress = "gz")


# Load data-------
source("functions/data_loader_splitter_admission.R")

# Set mtry
mtry_range <- 2:15

# Set up parallel
source("functions/parallel_setup.R")
parallel_setup(14)

# Loop through forest configuration
imputed_dfs <- foreach(j = 1:length(mtry_range)) %do%
  {
    print(j)
    
    # Perform random forest
    ptm <- proc.time()
    forest_results <- apache_scores %>%
      as.data.frame() %>% 
      missForest(mtry = mtry_range[j],
                 parallelize = "variables")
    proc.time() - ptm
    
    # Return imputed data
    forest_results$ximp
  }
names(imputed_dfs) <- paste("RF", mtry_range, sep = "_")

write_rds(imputed_dfs, "data/impute_mortality/forest.RDS",
          compress = "gz")

