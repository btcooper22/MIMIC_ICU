# Packages
require(dplyr)
require(magrittr)
require(readr)
require(foreach)
require(missMDA)
require(stringr)

# Functions
source("functions/calculate_apache_score.R")
source("functions/filename_metadata.R")
source("functions/parallel_setup.R")


# Load data-------
full_data <- read_csv("data/impute/complete_cases.csv") 

# Load model for prediction
apache_model <- read_rds("models/apache_model.RDS")

# Main loop through files-------

# Locate files
files <- dir("data/MCAR/", full.names = TRUE)

# Set number of components to profile
ncomps <- 2:11

# Set up parallel
ptm <- proc.time()
parallel_setup(15)

output <- foreach(f = 1:length(files),
                  .packages = c("foreach", "dplyr",
                               "missMDA", "stringr")) %dopar%
  {
    # Load file
    mcar_df <- read.csv(files[f])
    
    # Loop through ncomps
    imputed_dfs <- foreach(j = 1:length(ncomps)) %do%
      {
        # Select apache variables
        PCA_results <- mcar_df %>% 
          select(1:12) %>% 
          # Perform PCA
          imputePCA(ncp = ncomps[j])
        
        # Return imputed data
        list(ncomps = j,
          PCA_results$completeObs)
      }
    
    # Extract metadata
    metadata <- filename_metadata(files[f])

    # Return
    list(n = metadata[1],
         split = metadata[2],
         data = imputed_dfs)
  }
stopImplicitCluster()
proc.time() - ptm # 

write_rds(output, "data/impute/PCA.RDS",
          compress = "gz")