# Packages
require(dplyr)
require(magrittr)
require(readr)
require(foreach)
require(caret)
require(stringr)
require(RANN)

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

# Set number of K to profile
k_list <- 2:8

# Set up parallel
ptm <- proc.time()
parallel_setup(15)

all_KNN <- foreach(f = 1:length(files),
                   .packages = c("foreach", "dplyr",
                                 "caret", "stringr",
                                 "readr", "RANN")) %dopar%
  {
    # Extract metadata
    metadata <- filename_metadata(files[f])
    
    # Create filename
    filename <- paste("data/impute/KNN/KNN_S",
                      str_pad(metadata[2], 4, "right", 0),"_N",
                      str_pad(metadata[1],3,"left",0), ".RDS", sep = "")
    
    if(!file.exists(filename))
    {
      # Load file
      mcar_df <- read.csv(files[f])
      
      # Loop through K
      imputed_dfs <- foreach(k = 1:length(k_list),
                             .errorhandling = "remove") %do%
        {
          print(k)
          
          # Prepare KNN imputation
          KNN_preproc <- mcar_df %>% 
            select(1:12) %>% 
            preProcess(method = "knnImpute",
                       k = k_list[k],
                       knnSummary = mean)
          
          # Perform imputation
          KNN_results <- KNN_preproc %>% 
            predict(mcar_df %>% 
                      select(1:12),
                    na.action = na.pass)
        
          # De-normalise data
          procNames <- data.frame(
            col = names(KNN_preproc$mean),
            mean = KNN_preproc$mean,
            sd = KNN_preproc$std
          )
          
          for (i in procNames$col) {
            KNN_results[i] <-
              KNN_results[i] * KNN_preproc$std[i] + KNN_preproc$mean[i]
          }
          
          # Return imputed data
          list(k = k_list[k],
               KNN_results)
        }
      
      if(length(imputed_dfs) > 0)
      {
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
        
        # Extract successful K
        k_out <- foreach(k = 1:length(imputed_dfs),
                         .combine = "c") %do%
          {
            imputed_dfs[[k]][[1]]
          }
        
        # Combine and write
        output <- list(n = metadata[1],
                       split = metadata[2],
                       k_out,
                       probs)
      }else
      {
        # Combine and write
        output <- list(n = metadata[1],
                       split = metadata[2],
                       k_out = NA,
                       probs = NA)
      }
        
        write_rds(output, filename)

    } else
    {
      output <- read_rds(filename)
    }
    
    # Return
    output
  }
stopImplicitCluster()
proc.time() - ptm # 115 min

write_rds(all_KNN, "data/impute/KNN.RDS",
          compress = "gz")