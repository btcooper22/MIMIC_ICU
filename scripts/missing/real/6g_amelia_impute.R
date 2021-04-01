# Packages
require(dplyr)
require(magrittr)
require(readr)
require(foreach)
require(Amelia)
require(stringr)
require(tidyr)

# Load data-------
source("functions/data_loader_splitter.R")

# Perform for 5
amelia_results_5 <- apache_scores %>% 
  as.data.frame() %>% 
  amelia(parallel = "no",
         m = 5)

# Perform for 10
amelia_results_10 <- apache_scores %>% 
  as.data.frame() %>% 
  amelia(parallel = "no",
         m = 10)

# Perform for 20
amelia_results_20 <- apache_scores %>% 
  as.data.frame() %>% 
  amelia(parallel = "no",
         m = 20)

# Save to file
all_amelia <- list(amelia_results_5,
                   amelia_results_10,
                   amelia_results_20)
names(all_amelia) <- c("amelia_5", "amelia_10", "amelia_20")
write_rds(all_amelia, "data/impute_discharge/amelia.RDS",
          compress = "gz")


# Load data-------
source("functions/data_loader_splitter_admission.R")

# Perform for 5
amelia_results_5 <- apache_scores %>% 
  as.data.frame() %>% 
  amelia(parallel = "no",
         m = 5)

# Perform for 10
amelia_results_10 <- apache_scores %>% 
  as.data.frame() %>% 
  amelia(parallel = "no",
         m = 10)

# Perform for 20
amelia_results_20 <- apache_scores %>% 
  as.data.frame() %>% 
  amelia(parallel = "no",
         m = 20)

# Save to file
all_amelia <- list(amelia_results_5,
                   amelia_results_10,
                   amelia_results_20)
names(all_amelia) <- c("amelia_5", "amelia_10", "amelia_20")
write_rds(all_amelia, "data/impute_mortality/amelia.RDS",
          compress = "gz")

