# Packages
require(dplyr)
require(magrittr)
require(readr)
require(foreach)
require(missMDA)
require(stringr)

# Load data-------
source("functions/data_loader_splitter.R")

# Set number of components to profile
ncomps <- 2:16

# Loop through ncomps
imputed_dfs <- foreach(j = 1:length(ncomps)) %do%
  {
    print(j)
    
    # Perform PCA
    PCA_results <- apache_scores %>% 
      as.data.frame() %>% 
      imputePCA(ncp = ncomps[j])
    
    # Return imputed data
    PCA_results$completeObs
  }
    
names(imputed_dfs) <- paste("PCA", ncomps, sep = "_")

write_rds(imputed_dfs, "data/impute_discharge/PCA.RDS",
          compress = "gz")
