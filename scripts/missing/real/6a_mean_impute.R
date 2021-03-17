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

# Main loop through files-------

# Load files
source("functions/data_loader_splitter.R")

# Create zero values for each variable
zero_values <- c("temperature" = 37, "systolicbp" = 120, "diastolicbp" = 80,
                 "pulse" = 80, "respiratory" = 15, "arterialpH" = 7.35,
                 "sodium" = 135, "potassium" = 3.6, "creatinine" = 0.8,
                 "haematocrit" = 35, "whitebloodcount" = 10,
                 "glasgowcomaeye" = 4, "glasgowcomamotor" = 6,
                 "glasgowcomaverbal" = 5, "bicarbonate" = 25,
                 "arterialoxygen" = 90, "arterialcarbon" = 40)

# Loop through columns
imputed_cols <- foreach(i = 1:ncol(apache_scores)) %do%
  {
    # Replace NA with zero values
    values <- apache_scores[,i] %>% deframe()
    values[is.na(values)] <- zero_values[i]
    values_zero <- values
    
    # Replace NA with mean
    values <- apache_scores[,i] %>% deframe()
    values[is.na(values)] <- mean(values, na.rm = T)
    values_mean <- values
    
    # Replace NA with median
    values <- apache_scores[,i] %>% deframe()
    values[is.na(values)] <- median(values, na.rm = T)
    values_median <- values
    
    # Output column list
    data.frame(zero = values_zero,
         mean = values_mean,
         median = values_median)
  }

# Disentangle
output <- foreach(i = 1:3) %:%
  foreach(j = 1:17, .combine = "cbind") %do%
  {
    df <- data.frame(imputed_cols[[j]][,i])
    names(df) <- names(zero_values)[j]
    df
  }

names(output) <- paste("average_", c("zero", "mean", "median"), sep = "")

write_rds(output, "data/impute_discharge/average.RDS")
