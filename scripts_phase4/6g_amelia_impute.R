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

# Select apache variables
amelia_results <- apache_scores %>% 
  as.data.frame() %>% 
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
names(amelia_output) <- names(apache_scores)

# Save to file
all_amelia <- list(amelia_output)
names(all_amelia) <- c("amelia_100")
write_rds(all_amelia, "data/impute_mortality/amelia.RDS",
          compress = "gz")

