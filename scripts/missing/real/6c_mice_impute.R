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
require(tibble)

# Functions

# Set methods
mice_methods <- c("pmm", "norm.boot", "norm.predict",
                  "sample", "norm", "norm.nob",
                  "quadratic", "ri")

m_profile <- c(5,10,20)

mice_profile <- expand_grid(mice_methods, m_profile) %>% 
  mutate(name = paste(mice_methods, m_profile,
                      sep = "-")) %>% 
  select(name) %>% 
  deframe()

# Main loop through files-------

# Load data
source("functions/data_loader_splitter.R")


# Prepare parallel options
psnice(value = 19)
registerDoParallel(12)

# Loop through methods
imputed_dfs <- foreach(j = 1:length(mice_profile),
                       .packages = c("mice", "magrittr")) %dopar%
  { 
    print(j)
    
    # Extract names
    hyperparameters <- strsplit(mice_profile[j], "-")[[1]]
    
    # Perform imputation
    MICE_results <- apache_scores %>% 
      mice(method = hyperparameters[1],
           m = as.numeric(hyperparameters[2]))
    
    # Return imputed data
    MICE_results
  }
names(imputed_dfs) <- paste("MICE", mice_profile, sep = "_")
stopImplicitCluster()

write_rds(imputed_dfs, "data/impute_discharge/MICE.RDS",
          compress = "gz")
rm(imputed_dfs, apache_additional,
   apache_scores)

# Main loop through files: Admission-------

# Load data
source("functions/data_loader_splitter_admission.R")

# Prepare parallel options
psnice(value = 19)
registerDoParallel(12)

# Loop through methods
imputed_dfs <- foreach(j = 1:length(mice_profile),
                       .packages = c("mice", "magrittr")) %dopar%
  { 
    print(j)
    
    # Extract names
    hyperparameters <- strsplit(mice_profile[j], "-")[[1]]
    
    # Perform imputation
    MICE_results <- apache_scores %>% 
      mice(method = hyperparameters[1],
           m = as.numeric(hyperparameters[2]))
    
    # Return imputed data
    MICE_results
  }
names(imputed_dfs) <- paste("MICE", mice_profile, sep = "_")
stopImplicitCluster()


write_rds(imputed_dfs, "data/impute_mortality/MICE.RDS",
          compress = "gz")

stopImplicitCluster()