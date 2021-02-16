# Packages
require(dplyr)
require(magrittr)
require(readr)
require(foreach)
require(missForest)
require(stringr)
require(tidyr)

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
files <- dir("data/MCAR/", full.names = TRUE)#[seq(1,900,200)]

# Set mtry
mtry_range <- 3

# Set up parallel
ptm <- proc.time()
parallel_setup(15)

all_forest <- foreach(f = 1:length(files),
                   .packages = c("foreach", "dplyr",
                                 "missForest", "stringr",
                                 "readr")) %dopar%
  {
    # Extract metadata
    metadata <- filename_metadata(files[f])
    
    # Create filename
    filename <- paste("data/impute/forest/forest_S",
                      str_pad(metadata[2], 4, "right", 0),"_N",
                      str_pad(metadata[1],3,"left",0), ".RDS", sep = "")
    
    if(!file.exists(filename))
    {
      # Load file
      mcar_df <- read.csv(files[f])
      
      # Loop through forest configuration
      imputed_dfs <- foreach(j = 1:length(mtry_range)) %do%
        {
          print(j)

          # Select apache variables
          forest_results <- mcar_df %>% 
            select(1:12) %>% 
            # Perform random forest
            missForest(mtry = mtry_range[j])

          # Return imputed data
          list(mtry = mtry_range[j],
               data = forest_results$ximp,
               error = forest_results$OOBerror)
        }
      
      # Perform predictions and disentangle
      probs <- foreach(k = 1:length(imputed_dfs),
                       .combine = "cbind") %do%
        {
          imputed_dfs[[k]][[2]] %>% 
            as.data.frame() %>% 
            cbind(mcar_df[,13:15]) %>% 
            mutate(apache_II_discharge = calculate_apache_scores(.)) %>% 
            predict(apache_model, newdata = .)
        }
      
      # Combine and write
      output <- list(n = metadata[1],
                     split = metadata[2],
                     mtry_range,
                     probs)
      
      write_rds(output, filename)
    } else
    {
      output <- read_rds(filename)
    }
    
    # Return
    output
  }
stopImplicitCluster()
proc.time() - ptm # 15 parallel iterations of sequential forests =  68 min
# 15 sequential iterations of parallel forests = 

write_rds(all_forest, "data/impute/forest.RDS",
          compress = "gz")