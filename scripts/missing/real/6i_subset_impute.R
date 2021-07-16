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

source("functions/apache_score_mortality.R")

# Load data
apache_df <- read_csv("data/apache_discharge_missing.csv")
apache_scores_30d <- apache_df %>% 
  select(-row_id, -adm_id,
         -mort_30, -age,
         -chronic, -fractioninspiredoxygen,
         -apache_II, -missing_abg, -acute_renal_failure)
apache_additional_30d <- apache_df %>% 
  select(row_id, adm_id,
         mort_30, age,
         chronic, fractioninspiredoxygen,
         apache_II, missing_abg, acute_renal_failure)
rm(apache_df)

# Load recent imputed data
results <- read_rds("data/impute_discharge/recent.RDS")[["recent_worsttotal"]]

# Score recent
new_missing_abg <- is.na(results[["arterialpH"]]) &
  (is.na(results[["arterialoxygen"]]) |
     is.na(results[["arterialcarbon"]]))

apache_additional_long <- apache_additional_30d %>% 
  mutate(missing_abg = new_missing_abg)

scores <- apache_score(cbind(results, apache_additional_long))

# Highlight complete cases
complete_cases <- !is.na(apache_additional_long$apache_II)

# Highlight still missing recent data and filter
still_missing <- is.na(scores$apache_scores) & !complete_cases
recent_out <- scores %>% 
  filter(!still_missing & !complete_cases)

# Perform MICE imputation
MICE_results <- apache_scores_30d %>% 
  filter(!still_missing) %>% 
  mice(method = "quadratic",
       m = 10, maxit = 10) # 10

# Perform Amelia imputation
amelia_results <- apache_scores_30d %>% 
  filter(!still_missing) %>% 
  as.data.frame() %>% 
  amelia(parallel = "no",
         m = 20) # 20

# Write results
write_rds(recent_out, "data/impute_discharge/subset_recent.RDS")

write_rds(list(results = MICE_results,
               additional = apache_additional_long %>% 
                 mutate(complete_cases) %>% 
                 filter(!still_missing)), 
          "data/impute_discharge/subset_mice.RDS",
          compress = "gz")

write_rds(list(results = amelia_results,
               additional = apache_additional_long %>% 
                 mutate(complete_cases) %>% 
                 filter(!still_missing)), "data/impute_discharge/subset_amelia.RDS",
          compress = "gz")
