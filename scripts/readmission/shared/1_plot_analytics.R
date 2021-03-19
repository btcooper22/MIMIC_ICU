# Packages
require(readr)
require(ROCR)
require(foreach)
require(dplyr)
require(tidyr)

# Load models
ICNARC_bespoke <- read_rds("scripts/readmission/shared/models/LTH_ICNARC_bespoke.RDS")
ICNARC_hammer <- read_rds("scripts/readmission/shared/models/LTH_ICNARC_hammer.RDS")
ICNARC_martin <- read_rds("scripts/readmission/shared/models/LTH_ICNARC_martin.RDS")
ICNARC_frost <- read_rds("scripts/readmission/shared/models/LTH_ICNARC_frost.RDS")

MIMIC_bespoke <- read_rds("scripts/readmission/shared/models/MIMIC_bespoke.RDS")
MIMIC_hammer <- read_rds("scripts/readmission/shared/models/MIMIC_hammer.RDS")
MIMIC_martin <- read_rds("scripts/readmission/shared/models/MIMIC_martin.RDS")
MIMIC_frost <- read_rds("scripts/readmission/shared/models/MIMIC_frost.RDS")

model_list <- list(ICNARC_bespoke, ICNARC_hammer, ICNARC_martin, ICNARC_frost,
                   MIMIC_bespoke, MIMIC_hammer, MIMIC_martin, MIMIC_frost)

# Table of values----
performance_df <- foreach(i = 1:length(model_list), .combine = "rbind") %do%
  {
    # Extract prediction object
    pred <- model_list[[i]][3]$discrimination
    
    # Extract AUC
    AUC <- performance(pred, "auc")@y.values[[1]]
    
    # Extract calibration test
    hoslem <- model_list[[i]][4]$calibration
    
    # Output
    data.frame(model = model_list[[i]][1]$model,
               dataset = model_list[[i]][2]$data,
               AUC, chisq = hoslem$statistic %>% unname(),
               p = hoslem$p.value > 0.05)
  }
performance_df %>% 
  pivot_wider(names_from = "dataset",
              values_from = c("AUC", "chisq", "p")) 

# Plot discrimination----

