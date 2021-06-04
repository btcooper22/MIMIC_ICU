# Packages
require(readr)
require(dplyr)
require(magrittr)
require(ROCR)
require(ResourceSelection)
require(foreach)
require(tidyr)
require(ggplot2)
require(doParallel)

source("functions/inverse_logit.R")
source("functions/nomogram_convert.R")
source("functions/crosstab.R")

# Load data----
results <- read_csv("data/icnarc_predictors.csv") %>% 
  filter(surgical_type != "221. Paediatric Cardiac Surgery") %>% 
  select(id, readmission,age, sex, admission_source, apache_II,
         length_of_stay, out_of_hours_discharge, 
         acute_renal_failure, respiratory_rate,
         age, sodium, potassium, urea, glucose, sex, general_surgery,
         cardiac_surgery, hyperglycaemia,
         anaemia, apache_II, posthospital_dependency,
         length_of_stay,
         clinically_ready_discharge_days) %>% 
  na.omit()

# Build frost scores----

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

# Bootstrap assessment----

# Prepare data
results %<>% 
  mutate(probs_frost = log(probs_frost / (1 - probs_frost)))

# Prepare parallel
ptm <- proc.time()
cl <- makeCluster(4)
registerDoParallel(cl)
nboot <- 1000

model_results <- foreach(i = 1:nboot, .combine = "rbind",
                         .packages = c("dplyr", "ROCR",
                                       "ResourceSelection")) %dopar%
  {
    set.seed(i)
    print(i)
    
    # Split data
    patients_train <- results %>% 
      group_by(readmission) %>% 
      slice_sample(prop = 0.75)
    
    patients_validate <- results %>% 
      filter(id %in% patients_train$id == FALSE)
    
    # Fit frost+
    frost_plus_model <- glm(readmission ~ probs_frost + clinically_ready_discharge_days,
                            data = patients_train,
                       family = "binomial")
    
    # Predict
    probs_frost_plus <- predict(frost_plus_model, newdata = patients_validate) %>% inverse_logit()
    
    # Make prediction object
    pred_frost <- prediction(patients_validate$probs_frost %>% inverse_logit(),
                             patients_validate$readmission)
    pred_frost_plus <- prediction(probs_frost_plus, patients_validate$readmission)
    
    # Extract AUC
    AUC_frost <- performance(pred_frost, measure = "auc")@y.values[[1]]
    AUC_frost_plus <- performance(pred_frost_plus, measure = "auc")@y.values[[1]]
    
    # Calibration
    cal_frost <- hoslem.test(patients_validate$readmission,
                             patients_validate$probs_frost %>% inverse_logit(),
                             g = 10)
    cal_frost_plus <- hoslem.test(patients_validate$readmission,
                             probs_frost_plus, g = 10)
    
    # Output
    data.frame(i, AUC_frost,
               AUC_frost_plus,
               cal_frost = cal_frost$statistic,
               cal_frost_plus = cal_frost_plus$statistic)
  }

proc.time() - ptm 
stopCluster(cl)

# Plot----

# Discrimination
model_results %>% 
  mutate(frost_AUC_diff = AUC_frost_plus - AUC_frost) %>% 
  ggplot(aes(x = frost_AUC_diff))+
  geom_density(fill = "#e41a1c", alpha = 0.8)+
  geom_vline(linetype = "dashed",
             xintercept = 0)+
  theme_classic(20)

# Calibration
model_results %>% 
  mutate(frost_cal_diff = cal_frost_plus - cal_frost) %>% 
  ggplot(aes(x = frost_cal_diff))+
  geom_density(fill = "#377eb8", alpha = 0.8)+
  geom_vline(linetype = "dashed",
             xintercept = 0)+
  theme_classic(20)
