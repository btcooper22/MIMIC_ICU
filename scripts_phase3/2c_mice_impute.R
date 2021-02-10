# Packages
require(dplyr)
require(magrittr)
require(readr)
require(foreach)
require(mice)
require(stringr)
require(randomForest)

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

# Set methods
mice_methods <- c("pmm","norm.boot",
                  "norm.predict",
                  "midastouch")

# Set up parallel
ptm <- proc.time()
parallel_setup(14)

all_mice <- foreach(f = 1:1200,
                   .packages = c("foreach", "dplyr",
                                 "mice", "stringr",
                                 "readr")) %dopar%
  {
    # Extract metadata
    metadata <- filename_metadata(files[f])
    
    # Create filename
    filename <- paste("data/impute/MICE/MICE_S",
                      str_pad(metadata[2], 4, "right", 0),"_N",
                      str_pad(metadata[1],3,"left",0), ".RDS", sep = "")
    
    if(!file.exists(filename))
    {
      # Load file
      mcar_df <- read.csv(files[f])
      
      # Loop through methods
      imputed_dfs <- foreach(j = 1:length(mice_methods)) %do%
        {
          print(j)
          
          # Select apache variables
          MICE_results <- mcar_df %>% 
            select(1:12) %>% 
            # Perform imputation
            mice(method = mice_methods[j])

          # Return imputed data
          list(method = mice_methods[j],
               df1 = complete(MICE_results, 1),
               df2 = complete(MICE_results, 2),
               df3 = complete(MICE_results, 3),
               df4 = complete(MICE_results, 4),
               df5 = complete(MICE_results, 5))
        }
      
      # Perform predictions and disentangle
      probs <- foreach(k = 1:length(imputed_dfs),
                       .combine = "cbind") %do%
        {
          # Loop through datasets within a method
          prob_matrix <- foreach(d = 1:5,
                  .combine = "cbind") %do%
          {
            imputed_dfs[[k]][[d+1]] %>% 
              as.data.frame() %>% 
              cbind(mcar_df[,13:15]) %>% 
              mutate(apache_II_discharge = calculate_apache_scores(.)) %>% 
              predict(apache_model, newdata = .)
          }
          # Return mean
          rowMeans(prob_matrix)
        }
      
      # Combine and write
      output <- list(n = metadata[1],
                     split = metadata[2],
                     mice_methods,
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
proc.time() - ptm #

write_rds(all_mice, "data/impute/MICE.RDS",
          compress = "gz")