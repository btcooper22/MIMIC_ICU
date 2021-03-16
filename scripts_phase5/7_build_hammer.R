# Packages
require(readr)
require(dplyr)
require(magrittr)
require(ROCR)
require(ResourceSelection)
require(foreach)

source("functions/inverse_logit.R")

# Load data
results <- read_csv("data/icnarc_predictors.csv") %>% 
  filter(surgical_type != "221. Paediatric Cardiac Surgery") %>% 
  select(id, readmission, sex, general_surgery,
         cardiac_surgery, hyperglycaemia,
         anaemia, apache_II, posthospital_dependency,
         length_of_stay) %>% 
  na.omit()

# Generate hammer scores
coefficients_hammer <- log(0.005 / (1 - 0.005)) +
  ifelse(results$sex == "M", 0.43, 0) +
  ifelse(results$general_surgery, 0.64, 0) +
  ifelse(results$cardiac_surgery, -1.01, 0) +
  ifelse(results$hyperglycaemia, 0.36, 0) +
  ifelse(results$anaemia, 1.11, 0) +
  ifelse(results$apache_II >- 20, 0.44, 0) +
  ifelse(results$posthospital_dependency == TRUE, 0.7, 0) +
  ifelse(results$length_of_stay >= 5, 0.88, 0)

# Convert coefficients to probs
probs_hammer <- inverse_logit(coefficients_hammer)

# Assess discrimination
pred <- prediction(probs_hammer, results$readmission)
performance(pred, measure = "auc")@y.values[[1]]

# Assess calibration
hoslem.test(results$readmission,
            probs_hammer, g = 10)
