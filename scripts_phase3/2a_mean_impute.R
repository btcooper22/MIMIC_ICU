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

# Generate external mean and median
external_averages <- full_data %>% 
  select(7:20) %>% 
  pivot_longer(1:14) %>% 
  group_by(name) %>% 
  summarise(mean = mean(value),
            median = median(value))

# Create global replacement function
global_replace <- function(key)
{
  # Key = 2col: name, value
  return(
    list(apache_artpH_discharge = key %>% filter(name == "apache_artpH_discharge") %>% select(value) %>% deframe(),
         apache_creatinine_discharge = key %>% filter(name == "apache_creatinine_discharge") %>% select(value) %>% deframe(),
         apache_gcs_discharge = key %>% filter(name == "apache_gcs_discharge") %>% select(value) %>% deframe(),
         apache_hematocrit_discharge = key %>% filter(name == "apache_hematocrit_discharge") %>% select(value) %>% deframe(),
         apache_map_discharge = key %>% filter(name == "apache_map_discharge") %>% select(value) %>% deframe(),
         apache_oxygenation_discharge = key %>% filter(name == "apache_oxygenation_discharge") %>% select(value) %>% deframe(),
         apache_potassium_discharge = key %>% filter(name == "apache_potassium_discharge") %>% select(value) %>% deframe(),
         apache_pulse_discharge = key %>% filter(name == "apache_pulse_discharge") %>% select(value) %>% deframe(),
         apache_respiratory_discharge = key %>% filter(name == "apache_respiratory_discharge") %>% select(value) %>% deframe(),
         apache_sodium_discharge = key %>% filter(name == "apache_sodium_discharge") %>% select(value) %>% deframe(),
         apache_temperature_discharge = key %>% filter(name == "apache_temperature_discharge") %>% select(value) %>% deframe(),
         apache_wbc_discharge = key %>% filter(name == "apache_wbc_discharge") %>% select(value) %>% deframe())
  )
}

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
    
    # Generate internal mean and median
    internal_averages <- mcar_df %>% 
      pivot_longer(1:14) %>% 
      group_by(name) %>% 
      summarise(mean = mean(value, na.rm = T),
                median = median(value, na.rm = T))
    
    # Assume 0
    zero_df <- mcar_df
    zero_df[is.na(zero_df)] <- 0
    
    # Internal mean
    int_mean_df <- mcar_df %>% 
      replace_na(global_replace(internal_averages %>% 
                                  select(name, mean) %>% 
                                  rename(value = "mean")))
    
    # External mean
    ext_mean_df <- mcar_df %>% 
      replace_na(global_replace(external_averages %>% 
                                  select(name, mean) %>% 
                                  rename(value = "mean")))
    
    # Internal median
    int_median_df <- mcar_df %>% 
      replace_na(global_replace(internal_averages %>% 
                                  select(name, median) %>% 
                                  rename(value = "median")))
    
    # External median
    ext_median_df <- mcar_df %>% 
      replace_na(global_replace(external_averages %>% 
                                  select(name, median) %>% 
                                  rename(value = "median")))
    
    # Sum and predict
    zero_df %>% mutate(apache_II_discharge = rowSums(zero_df)) %>% 
      predict(apache_model, newdata = .) -> zero_probs
    
    int_mean_df %>% mutate(apache_II_discharge = rowSums(int_mean_df)) %>% 
      predict(apache_model, newdata = .) -> int_mean_probs
    
    ext_mean_df %>% mutate(apache_II_discharge = rowSums(ext_mean_df)) %>% 
      predict(apache_model, newdata = .) -> ext_mean_probs
    
    int_median_df %>% mutate(apache_II_discharge = rowSums(int_median_df)) %>% 
      predict(apache_model, newdata = .) -> int_median_probs
    
    ext_median_df %>% mutate(apache_II_discharge = rowSums(ext_median_df)) %>% 
      predict(apache_model, newdata = .) -> ext_median_probs
    
    # Extract metadata
    metadata <- files[f] %>% 
      str_split("/", simplify = TRUE) %>% .[,3] %>% 
      str_split("_", simplify = TRUE)
    
    split <- metadata[,1] %>% 
      substr(2,4) %>% 
      as.numeric()
    
    n <- metadata[,2] %>% 
      substr(2,4) %>% 
      as.numeric()
    
    # Output data frame for list
    data.frame(split, n, zero_probs, int_mean_probs,
               ext_mean_probs, int_median_probs,
               ext_median_probs)
  }
stopImplicitCluster()
proc.time() - ptm # 32s 

write_rds(output, "data/impute/average.RDS")