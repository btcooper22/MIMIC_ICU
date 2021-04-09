# Packages
require(readr)
require(dplyr)
require(magrittr)
require(ROCR)
require(ResourceSelection)
require(foreach)
require(tidyr)

source("functions/inverse_logit.R")
source("functions/nomogram_convert.R")
source("functions/crosstab.R")

# Load data (2083)
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

# Create scoring system from nomogram for APACHE-II
apache_score_system_input <- data.frame(
  input = seq(0, 45, 5),
  pix = seq(12, 733, length.out = 10)
)
apache_score_system_output <- data.frame(
  output = seq(0, 20, 1),
  pix = seq(12, 733, length.out = 21)
)

# Score for age
age_convert <- data.frame(age_band = unique(results$age),
                          points = c(4.5, 4, 6, 3, 2,
                                     5, 5.5, 7, 6.5,
                                     3.5, 2.5, 1.5,
                                     7.5, 0, 8),
                          med_age = c(73, 58, 78, 83, 68,
                                      53, 38, 48, 90, 63, 
                                      28, 43, 88, 18, 33))
results$age_est <- NA
for(i in 1:nrow(results))
{
  results$age_est[i] <- age_convert$med_age[which(age_convert$age_band == results$age[i])] 
  results$age[i] <- age_convert$points[which(age_convert$age_band == results$age[i])] 
}


# Cross tabulate
table(results$readmission)
(table(results$readmission) / nrow(results) * 100) %>% round(1)

results %>% 
  group_by(readmission) %>% 
  summarise(mean = mean(age_est),
            sd = sd(age_est))

crosstab("sex", results)
crosstab("admission_source", results)

results %>% 
  group_by(readmission) %>% 
  summarise(mean = mean(apache_II),
            sd = sd(apache_II))

crosstab("los_7", results %>% 
           mutate(los_7 = length_of_stay >=7 ))
crosstab("out_of_hours_discharge", results)
crosstab("acute_renal_failure", results)


# Build score
scores_frost <-  as.numeric(results$age) + 
  ifelse(results$sex == "M", 2, 0) +
  ifelse(results$admission_source == "Ward", 14.75, 1) +
  nomogram_convert(results$apache_II, apache_score_system_input,
                   apache_score_system_output) +
  ifelse(results$length_of_stay >= 7, 17, 0) +
  ifelse(results$out_of_hours_discharge, 4, 0) +
  ifelse(results$acute_renal_failure, 10, 0)

# Create nomogram system for total points
points_system_input <- data.frame(
  input = seq(0, 100, 5),
  pix = seq(7, 727, length.out = 21)
)
points_system_output <- data.frame(
  output = seq(0.05, 0.5, 0.05),
  pix = c(212, 330, 403, 460, 506,
          547, 583, 618, 652, 683)
)

# Convert using logarithmic transform
probs_frost <- nomogram_convert(scores_frost, points_system_input,
                                points_system_output, log = TRUE)

# Assess discrimination
pred <- prediction(probs_frost, results$readmission)
perf <- performance(pred, "tpr", "fpr")
performance(pred, measure = "auc")@y.values[[1]]

# Split to deciles
decile_df <- tibble(
  patient_id = 1:nrow(results),
  readmission = results$readmission,
  probs = probs_frost,
  decile = ntile(probs_frost, 10)) %>% 
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
hoslem_frost <- hoslem.test(results$readmission,
                            probs_frost, g = 10)

list(model = "frost",
     data = "LTH_ICNARC",
     discrimination = pred,
     calibration = hoslem_frost,
     deciles = calibration_df) %>% 
  write_rds("scripts/readmission/shared/models/LTH_ICNARC_frost.RDS")



