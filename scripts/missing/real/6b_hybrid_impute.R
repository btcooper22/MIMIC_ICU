# Packages
require(dplyr)
require(magrittr)
require(readr)
require(foreach)
require(stringr)
require(tidyr)
require(tibble)
require(lubridate)
require(doParallel)
require(tools)
require(mice)
require(Amelia)

# Load recent imputed data
results <- read_rds("data/impute_discharge/recent.RDS")[["recent_worsttotal"]]

# Perform MICE imputation
MICE_results <- results %>% 
  mice(method = "quadratic",
       m = 10, maxit = 10) # 10

# Perform Amelia imputation
amelia_results <- results %>% 
  amelia(parallel = "no",
         m = 20) # 5

# Write results
write_rds(MICE_results, "data/impute_discharge/hybrid_mice.RDS",
          compress = "gz")

write_rds(amelia_results, "data/impute_discharge/hybrid_amelia.RDS",
          compress = "gz")