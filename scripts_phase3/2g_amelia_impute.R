# Packages
require(dplyr)
require(magrittr)
require(readr)
require(foreach)
require(Amelia)
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
files <- dir("data/MCAR/", full.names = TRUE)


# Set up parallel
ptm <- proc.time()
parallel_setup(15)

all_amelia <- foreach(f = 1:length(files),
                      .packages = c("foreach", "dplyr",
                                    "Amelia", "stringr",
                                    "readr")) %dopar%
  {
    # Extract metadata
    metadata <- filename_metadata(files[f])
    
    # Create filename
    filename <- paste("data/impute/amelia/amelia_S",
                      str_pad(metadata[2], 4, "right", 0),"_N",
                      str_pad(metadata[1],3,"left",0), ".RDS", sep = "")
    
    if(!file.exists(filename))
    {
      # Load file
      mcar_df <- read.csv(files[f])
      
          
          # Select apache variables
          amelia_results <- mcar_df %>% 
            select(1:12) %>% 
            # Perform amelia
            amelia(parallel = "no",
                   m = 100)
          
          # Extract imputed datasets
          amelia_array <- array(unlist(amelia_results$imputations), 
                                dim = c(dim(amelia_results$imputations[[1]]), 100))
          rm(amelia_results)
          amelia_output <- apply(amelia_array, c(1,2), mean) %>% 
            as.data.frame()
          rm(amelia_array)
          names(amelia_output) <- names(mcar_df)[1:12]

      
      # Perform predictions
         probs <-  amelia_output %>% 
            cbind(mcar_df[,13:15]) %>% 
            mutate(apache_II_discharge = calculate_apache_scores(.)) %>% 
            predict(apache_model, newdata = .)
         rm(amelia_output)
      
      # Combine and write
      output <- list(n = metadata[1],
                     split = metadata[2],
                     probs)
      gc()
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

write_rds(all_amelia, "data/impute/amelia.RDS",
          compress = "gz")

