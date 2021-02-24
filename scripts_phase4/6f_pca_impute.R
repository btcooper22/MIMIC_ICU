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

all_PCA <- foreach(f = 1:length(files),
                   .packages = c("foreach", "dplyr",
                                 "missMDA", "stringr",
                                 "readr")) %dopar%
  {
    # Extract metadata
    metadata <- filename_metadata(files[f])
    
    # Create filename
    filename <- paste("data/impute/PCA/PCA_S",
                      str_pad(metadata[2], 4, "right", 0),"_N",
                      str_pad(metadata[1],3,"left",0), ".RDS", sep = "")
    
    if(!file.exists(filename))
    {
      # Load file
      mcar_df <- read.csv(files[f])
      
      # Loop through ncomps
      imputed_dfs <- foreach(j = 1:length(ncomps)) %do%
        {
          print(j)
          # Select apache variables
          PCA_results <- mcar_df %>% 
            select(1:12) %>% 
            # Perform PCA
            imputePCA(ncp = ncomps[j])
          
          # Return imputed data
          list(ncomps = j,
               PCA_results$completeObs)
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
                     ncomps,
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
proc.time() - ptm # 115 min

write_rds(all_PCA, "data/impute/PCA.RDS",
          compress = "gz")