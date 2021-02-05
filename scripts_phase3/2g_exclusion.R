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

# Load data-------
full_data <- read_csv("data/impute/complete_cases.csv") 

# Main loop through files-------

# Load model for prediction
apache_model <- read_rds("models/apache_model.RDS")

# Locate files
files <- dir("data/MCAR/", full.names = TRUE)

# Set up parallel
ptm <- proc.time()
psnice(value = 19)
registerDoParallel(ifelse(detectCores() <= 16,
                          detectCores() - 1,
                          16)
)

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
    metadata <- files[f] %>% 
      str_split("/", simplify = TRUE) %>% .[,3] %>% 
      str_split("_", simplify = TRUE)
    
    split <- metadata[,1] %>% 
      substr(2,5) %>% 
      as.numeric()
    
    n <- metadata[,2] %>% 
      substr(2,4) %>% 
      as.numeric()
    
    # Output data frame for list
    data.frame(split, n, probs)
  }
stopImplicitCluster()
proc.time() - ptm # 9s 

write_rds(output, "data/impute/exclusion.RDS",
          compress = "gz")
