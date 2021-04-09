# Packages
require(dplyr)
require(magrittr)
require(readr)
require(foreach)
require(mice)
require(stringr)
require(randomForest)
require(tidyr)
require(doParallel)
require(tools)

# Functions

# Main loop through files-------

# Load data
source("functions/data_loader_splitter.R")

# Set methods
mice_methods <- c("pmm", "norm.boot", "norm.predict", "midastouch",
                  "sample", "norm", "norm.nob",
                  "quadratic", "ri")

# Prepare parallel options
psnice(value = 19)
registerDoParallel(length(mice_methods))

# Loop through methods
imputed_dfs <- foreach(j = 1:length(mice_methods),
                       .packages = c("mice", "magrittr")) %dopar%
  {
    print(j)
    
    # Perform imputation
    MICE_results <- apache_scores %>% 
      mice(method = mice_methods[j],
           m = 5)
    
    # Return imputed data
    MICE_results
  }
names(imputed_dfs) <- paste("MICE", mice_methods, sep = "_")

write_rds(imputed_dfs, "data/impute_discharge/MICE.RDS")


# Main loop through files: Admission-------

# Load data
source("functions/data_loader_splitter_admission.R")

# Set methods
mice_methods <- c("pmm", "norm.boot", "norm.predict", "midastouch",
                  "sample", "norm", "norm.nob",
                  "quadratic", "ri")


# Prepare parallel options
psnice(value = 19)
registerDoParallel(length(mice_methods))

# Loop through methods
imputed_dfs <- foreach(j = 1:length(mice_methods),
                       .packages = c("mice", "magrittr")) %dopar%
  {
    print(j)
    
    # Perform imputation
    MICE_results <- apache_scores %>% 
      mice(method = mice_methods[j],
           m = 5)
    
    # Return imputed data
    MICE_results
  }
names(imputed_dfs) <- paste("MICE", mice_methods, sep = "_")

write_rds(imputed_dfs, "data/impute_mortality/MICE.RDS")

stopImplicitCluster()