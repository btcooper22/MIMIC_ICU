# Packages
require(readr)
require(dplyr)
require(magrittr)
require(ROCR)
require(ResourceSelection)
require(foreach)
require(tidyr)

source("functions/inverse_logit.R")
source("functions/crosstab.R")

# Load data
results <- read_csv("data/icnarc_predictors.csv") %>% 
  filter(surgical_type != "221. Paediatric Cardiac Surgery") %>% 
  select(id, readmission,age, sex, admission_source, apache_II,
         length_of_stay, out_of_hours_discharge, 
         acute_renal_failure, respiratory_rate,
         age, sodium, potassium, urea, glucose, sex, general_surgery,
         cardiac_surgery, hyperglycaemia,
         anaemia, apache_II, posthospital_dependency,
         length_of_stay) %>% 
  na.omit()

crosstab("general_surgery", results)
crosstab("cardiac_surgery", results)
crosstab("hyperglycaemia", results)
crosstab("anaemia", results)
crosstab("high_apache", results %>% 
           mutate(high_apache = apache_II >= 20))
crosstab("ambulation", results %>% 
           mutate(ambulation = posthospital_dependency == FALSE))
crosstab("los_5", results %>% 
           mutate(los_5 = length_of_stay >= 5))

# Generate hammer scores
coefficients_hammer <- log(0.005 / (1 - 0.005)) +
  ifelse(results$sex == "M", 0.43, 0) +
  ifelse(results$general_surgery, 0.64, 0) +
  ifelse(results$cardiac_surgery, -1.01, 0) +
  ifelse(results$hyperglycaemia, 0.36, 0) +
  ifelse(results$anaemia, 1.11, 0) +
  ifelse(results$apache_II >= 20, 0.44, 0) +
  ifelse(results$posthospital_dependency == TRUE, 0.7, 0) +
  ifelse(results$length_of_stay >= 5, 0.88, 0) + 0.26

# Convert coefficients to probs
probs_hammer <- inverse_logit(coefficients_hammer)

# Assess discrimination
pred <- prediction(probs_hammer, results$readmission)
performance(pred, measure = "auc")@y.values[[1]]

# Split to deciles
decile_df <- tibble(
  patient_id = 1:nrow(results),
  readmission = results$readmission,
  probs = probs_hammer,
  decile = ntile(probs_hammer, 10)) %>% 
  na.omit()

# Measure calibration
calibration_df <- decile_df %>% 
  select(readmission, probs, decile) %>% 
  group_by(decile) %>%
  summarise(N = length(readmission),
            observed = (sum(readmission) / length(readmission)) * 100,
            predicted = mean(probs * 100),
            error = sd(probs * 100))

# Assess calibration
hoslem_hammer <- hoslem.test(results$readmission,
            probs_hammer, g = 10)

list(model = "hammer",
     data = "LTH_ICNARC",
     discrimination = pred,
     calibration = hoslem_hammer,
     deciles = calibration_df) %>% 
  write_rds("scripts/readmission/shared/models/LTH_ICNARC_hammer.RDS")

