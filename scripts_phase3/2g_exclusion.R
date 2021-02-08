# Packages
require(dplyr)
require(magrittr)
require(readr)
require(tibble)
require(foreach)
require(stringr)
require(doParallel)
require(tools)
require(tidyr)

# Functions
source("functions/parallel_setup.R")
source("functions/filename_metadata.R")

# Load data-------
full_data <- read_csv("data/impute/complete_cases.csv") 

# Main loop through files-------

# Load model for prediction
apache_model <- read_rds("models/apache_model.RDS")

# Locate files
files <- dir("data/MCAR/", full.names = TRUE)

# Set up parallel
ptm <- proc.time()
parallel_setup(14)

output <- foreach(f = 1:length(files), .packages = c("dplyr", "tidyr",
                                                     "stringr", "tibble")) %dopar%
  {
    # Load file
    mcar_df <- read.csv(files[f])
    
    # Remove NA and predict
    probs <- mcar_df %>% 
      mutate(apache_II_discharge = rowSums(mcar_df)) %>% 
      predict(apache_model, newdata = .)
    
    # Extract metadata
    metadata <- filename_metadata(files[f])
    
    # Output data frame for list
    data.frame(split = metadata[2], 
               n = metadata[1], probs)
  }
stopImplicitCluster()
proc.time() - ptm # 9s 

write_rds(output, "data/impute/exclusion.RDS",
          compress = "gz")
