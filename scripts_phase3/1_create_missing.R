# Packages
require(dplyr)
require(magrittr)
require(readr)
require(tibble)
require(ggplot2)
require(ROCR)
require(foreach)
require(stringr)
require(doParallel)
require(tools)

source("functions/inverse_logit.R")
source("functions/calibration.R")

# Load data----- 
# Trim to APACHE-II
full_data <- read_csv("data/final_patients.csv") %>% 
  select(matches("_discharge$")) %>% 
  select(matches("^apache_")) %>% 
  cbind(read_csv("data/final_patients.csv") %>% 
          select(row_id, subject_id,
                 adm_id, readmission,
                 age),.)

# Binarise readmission
full_data %<>% 
  mutate(readmission = readmission == "Readmitted to ICU") %>% 
  # Map bicarbonate onto oxygenation
  mutate(apache_oxygenation_discharge = case_when(
    is.na(apache_oxygenation_discharge) ~ apache_bicarbonate_discharge,
    !is.na(apache_oxygenation_discharge) ~ apache_oxygenation_discharge)) %>% 
  select(-apache_bicarbonate_discharge) %>% 
  # Remove NA
  na.omit()
  
# Confirm no missing values
table(full_data$readmission, useNA = "ifany")
summary(full_data)

# Remove implausible values
full_data %<>% 
  filter(apache_temperature_discharge < 40, 
         apache_temperature_discharge > 30,
         apache_map_discharge > 0,
         apache_map_discharge < 180,
         apache_pulse_discharge > 0,
         apache_pulse_discharge < 160,
         apache_respiratory_discharge > 0,
         apache_respiratory_discharge <= 60)

# Write
write_csv(full_data,"data/impute/complete_cases.csv")

# Simulate MCAR (missing completely at random)---------
# Set n
N <- 100

# Set splits
splits <- c(0.01, 0.05, seq(0.1, 0.8, 0.1))

# Pull APACHE matrix
apache_matrix <- full_data[,c(7:17,19)]
apache_additional <- full_data[,c(18,20)]

# Set up parallel
ptm <- proc.time()
psnice(value = 19)
registerDoParallel(ifelse(detectCores() <= 16,
                          detectCores() - 1,
                          16)
)

# Outer loop over N
foreach(n = 1:N, .packages = c("foreach", "stringr")) %dopar%
  {
    # Inner loop for splits
    foreach(s = 1:length(splits)) %do%
      {
        # Define shape
        shape <- nrow(apache_matrix) * ncol(apache_matrix)
        
        # Draw randomisation array
        rand_array <- rbinom(shape, 1, splits[s]) %>% 
          as.logical()
        
        # Reshape array
        NA_array <- matrix(rand_array, 
               nrow(apache_matrix),
               ncol(apache_matrix))
        
        # Remove values
        data_out <- apache_matrix
        data_out[NA_array] <- NA
        
        # Reattach non-NA'd values
        data_out <- cbind(data_out, apache_additional)
        
        # Write to file
        filename <- paste("data/MCAR/S", str_pad(splits[s], 4, "right", 0),
                          "_N", str_pad(n,3,"left",0), ".csv", sep = "")
        write.csv(data_out, filename,
                  row.names = FALSE)
      }
  }
stopImplicitCluster()
proc.time() - ptm # 17s par, 112 seq


# Basic methods -> Assume 0, Mean, hot-deck

# Intermediate methods -> MICE, KNN, MissForest

# Complex methods from competiton