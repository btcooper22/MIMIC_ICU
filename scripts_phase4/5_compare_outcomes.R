# Setup----
require(readr)
require(dplyr)
require(tidyr)
require(ROCR)
require(ggplot2)
require(ResourceSelection)

source("functions/inverse_logit.R")

# Load data
results <- read_csv("data/apache_data_full.csv") %>% 
  filter(!is.na(apache_II))

# Split data
admission_results <- results %>% 
  filter(type == "admission")

discharge_results <- results %>% 
  filter(type == "discharge")

# Describe data
admission_results %>% 
  filter(mort_inhosp == FALSE &
           mort_inunit == FALSE) %>% 
  summarise(readmission = mean(readmission) * 100)

admission_results %>% 
  filter(mort_inunit == FALSE) %>% 
  summarise(`30-day mortality` = mean(mort_30) * 100)

admission_results %>% 
  filter(mort_inunit == FALSE) %>% 
  summarise(`In-hospital mortality` = mean(mort_inhosp) * 100)

admission_results %>% 
  summarise(`In-unit mortality` = mean(mort_inunit) * 100)

# Build models----

# Readmission, admission score
model_admission_readmission <- glm(readmission ~ apache_II,
                                   data = admission_results,
                                   subset = mort_inhosp == FALSE &
                                     mort_inunit == FALSE,
                                   family = "binomial")

# Readmission, discharge score
model_discharge_readmission <- glm(readmission ~ apache_II,
                                   data = discharge_results,
                                   subset = mort_inhosp == FALSE &
                                     mort_inunit == FALSE,
                                   family = "binomial")

# 30-day mortality, admission score
model_admission_mort30 <- glm(mort_30 ~ apache_II,
                              data = admission_results,
                              subset = mort_inunit == FALSE,
                              family = "binomial")

# 30-day mortality, discharge score
model_discharge_mort30 <- glm(mort_30 ~ apache_II,
                              data = discharge_results,
                              subset = mort_inunit == FALSE,
                              family = "binomial")

# In-hospital mortality, admission score
model_admission_mortInHospital <- glm(mort_inhosp ~ apache_II,
                              data = admission_results,
                              subset = mort_inunit == FALSE,
                              family = "binomial")

# In-hospital mortality, discharge score
model_discharge_mortInHospital <- glm(mort_inhosp ~ apache_II,
                                      data = discharge_results,
                                      subset = mort_inunit == FALSE,
                                      family = "binomial")

# In-unit mortality, admission score
model_admission_mortInUnit <- glm(mort_inunit ~ apache_II,
                                      data = admission_results,
                                      family = "binomial")

# Generate predictions----

# admission_readmission
data.frame(probs = predict(model_admission_readmission, 
                           newdata = model_admission_readmission$data),
           outcome = model_admission_readmission$data$readmission) %>% 
  mutate(probs = inverse_logit(probs)) -> probs_admission_readmission

# discharge_readmission
data.frame(probs = predict(model_discharge_readmission, 
                           newdata = model_discharge_readmission$data),
           outcome = model_discharge_readmission$data$readmission) %>% 
  mutate(probs = inverse_logit(probs)) -> probs_discharge_readmission

# admission_30d
data.frame(probs = predict(model_admission_mort30, 
                           newdata = model_admission_mort30$data),
           outcome = model_admission_mort30$data$mort_30) %>% 
  mutate(probs = inverse_logit(probs)) -> probs_admission_mort30

# discharge_30d
data.frame(probs = predict(model_discharge_mort30, 
                           newdata = model_discharge_mort30$data),
           outcome = model_discharge_mort30$data$mort_30) %>% 
  mutate(probs = inverse_logit(probs)) -> probs_discharge_mort30

# admission_inhosp
data.frame(probs = predict(model_admission_mortInHospital, 
                           newdata = model_admission_mortInHospital$data),
           outcome = model_admission_mortInHospital$data$mort_inhosp) %>% 
  mutate(probs = inverse_logit(probs)) -> probs_admission_mortInHospital

# discharge_inhosp
data.frame(probs = predict(model_discharge_mortInHospital, 
                           newdata = model_discharge_mortInHospital$data),
           outcome = model_discharge_mortInHospital$data$mort_inhosp) %>% 
  mutate(probs = inverse_logit(probs)) -> probs_discharge_mortInHospital

# admission_inunit
data.frame(probs = predict(model_admission_mortInUnit, 
                           newdata = model_admission_mortInUnit$data),
           outcome = model_admission_mortInUnit$data$mort_inunit) %>% 
  mutate(probs = inverse_logit(probs)) -> probs_admission_mortInUnit

# Assess discrimination----

# Generate prediction objects
prediction_admission_readmission <- prediction(probs_admission_readmission$probs,
                                               probs_admission_readmission$outcome)

prediction_discharge_readmission <- prediction(probs_discharge_readmission$probs,
                                               probs_discharge_readmission$outcome)

prediction_admission_mort30 <- prediction(probs_admission_mort30$probs,
                                               probs_admission_mort30$outcome)

prediction_discharge_mort30 <- prediction(probs_discharge_mort30$probs,
                                               probs_discharge_mort30$outcome)

prediction_admission_mortInHospital <- prediction(probs_admission_mortInHospital$probs,
                                               probs_admission_mortInHospital$outcome)

prediction_discharge_mortInHospital <- prediction(probs_discharge_mortInHospital$probs,
                                               probs_discharge_mortInHospital$outcome)

prediction_admission_mortInUnit <- prediction(probs_admission_mortInUnit$probs,
                                               probs_admission_mortInUnit$outcome)


# Generate AUC
performance(prediction_admission_readmission, measure = "auc")@y.values[[1]]
performance(prediction_discharge_readmission, measure = "auc")@y.values[[1]]
performance(prediction_admission_mort30, measure = "auc")@y.values[[1]]
performance(prediction_discharge_mort30, measure = "auc")@y.values[[1]]
performance(prediction_admission_mortInHospital, measure = "auc")@y.values[[1]]
performance(prediction_discharge_mortInHospital, measure = "auc")@y.values[[1]]
performance(prediction_admission_mortInUnit, measure = "auc")@y.values[[1]]

# Generate performance objects
performance_admission_mort30 <- performance(prediction_admission_mort30, "tpr", "fpr")
performance_discharge_mort30 <- performance(prediction_discharge_mort30, "tpr", "fpr")
performance_admission_mortInHospital <- performance(prediction_admission_mortInHospital, "tpr", "fpr")
performance_discharge_mortInHospital <- performance(prediction_discharge_mortInHospital, "tpr", "fpr")
performance_admission_mortInUnit <- performance(prediction_admission_mortInUnit, "tpr", "fpr")

# Plot AUC
data.frame(x = performance_admission_mort30@x.values[[1]],
           y = performance_admission_mort30@y.values[[1]],
           model = "admission_30") %>% 
  rbind(
    data.frame(x = performance_discharge_mort30@x.values[[1]],
               y = performance_discharge_mort30@y.values[[1]],
               model = "discharge_30"),
    data.frame(x = performance_admission_mortInHospital@x.values[[1]],
               y = performance_admission_mortInHospital@y.values[[1]],
               model = "admission_inhospital"),
    data.frame(x = performance_discharge_mortInHospital@x.values[[1]],
               y = performance_discharge_mortInHospital@y.values[[1]],
               model = "discharge_inhospital"),
    data.frame(x = performance_admission_mortInUnit@x.values[[1]],
               y = performance_admission_mortInUnit@y.values[[1]],
               model = "admission_inunit")
  ) %>% 
  ggplot(aes(x, y, colour = model))+
  geom_abline(slope = 1, intercept = 0,
              linetype = "dotted",
              size = 1)+
  geom_path(size = 1)+
  labs(x = "1 - Specificity",
       y = "Sensitivity")+
  theme_classic(20)+
  scale_color_brewer(palette = "Set1",
                     name = "")

# Assess calibration-----

# Split data into deciles
deciles_admission_mort30 <- tibble(
  mortality = probs_admission_mort30$outcome,
  probs = probs_admission_mort30$probs,
  decile = ntile(probs_admission_mort30$probs, 10)
)

deciles_discharge_mort30 <- tibble(
  mortality = probs_discharge_mort30$outcome,
  probs = probs_discharge_mort30$probs,
  decile = ntile(probs_discharge_mort30$probs, 10)
)


deciles_admission_mortInHospital <- tibble(
  mortality = probs_admission_mortInHospital$outcome,
  probs = probs_admission_mortInHospital$probs,
  decile = ntile(probs_admission_mortInHospital$probs, 10)
)


deciles_discharge_mortInHospital <- tibble(
  mortality = probs_discharge_mortInHospital$outcome,
  probs = probs_discharge_mortInHospital$probs,
  decile = ntile(probs_discharge_mortInHospital$probs, 10)
)

deciles_admission_mortInUnit <- tibble(
  mortality = probs_admission_mortInUnit$outcome,
  probs = probs_admission_mortInUnit$probs,
  decile = ntile(probs_admission_mortInUnit$probs, 10)
)

# Create quick calibration function
calibration <- function(.df)
{
  .df %>% 
    group_by(decile) %>%
    summarise(N = length(mortality),
              observed = (sum(mortality) / length(mortality)) * 100,
              predicted = mean(probs * 100),
              error = sd(probs * 100))
}

# Calculate calibration
calib_admission_mort30 <- calibration(deciles_admission_mort30)
calib_discharge_mort30 <- calibration(deciles_discharge_mort30)
calib_admission_mortInHospital <- calibration(deciles_admission_mortInHospital)
calib_discharge_mortInHospital <- calibration(deciles_discharge_mortInHospital)
calib_admission_mortInUnit <- calibration(deciles_admission_mortInUnit)

# Plot calibration
rbind(calib_admission_mort30 %>% mutate(model = "admission_mort30"), 
      calib_discharge_mort30 %>% mutate(model = "discharge_mort30"),
      calib_admission_mortInHospital %>% mutate(model = "admission_mortInHospital"),
      calib_discharge_mortInHospital %>% mutate(model = "discharge_mortInHospital"),
      calib_admission_mortInUnit %>% mutate(model = "admission_mortInUnit")) %>% 
  ggplot(aes(predicted, observed ,
             colour = model))+
  geom_abline(slope = 1, intercept = 0,
              size = 1)+
  geom_path(size = 1)+
  geom_pointrange(aes(ymin = observed - error ,
                      ymax = observed + error,
                      y = observed,
                      x = predicted,
                      colour = model)) +
  theme_classic(20)+
  theme(legend.position = "top")+
  labs(x = "Predicted readmission",
       y = "Observed readmission")+
  scale_colour_brewer(palette = "Set1",
                      name = "")+
  facet_wrap(~model, scales = "free")

# Quantify calibration
hoslem.test(probs_admission_mort30$outcome,
            probs_admission_mort30$probs, g = 10)
hoslem.test(probs_discharge_mort30$outcome,
            probs_discharge_mort30$probs, g = 10)
hoslem.test(probs_admission_mortInHospital$outcome,
            probs_admission_mortInHospital$probs, g = 10)
hoslem.test(probs_discharge_mortInHospital$outcome,
            probs_discharge_mortInHospital$probs, g = 10)
hoslem.test(probs_admission_mortInUnit$outcome,
            probs_admission_mortInUnit$probs, g = 10)

# Write best data----

# Load data
results_out <- 
  read_csv("data/apache_data_full.csv") %>% 
  filter(type == "admission") %>% 
  mutate(fractioninspiredoxygen = ifelse(is.na(fractioninspiredoxygen), 21, 
                fractioninspiredoxygen)) %>% 
  select(-subject_id, -type, -readmission, 
         -mort_inhosp, -mort_30)

# Count NA in each column
results_out %>% 
  select(-row_id, -adm_id,
         -mort_inunit, -age,
         -chronic, -fractioninspiredoxygen) %>% 
  sapply(function(x) sum(is.na(x)))

# Table of missingness
missing_hist <- results_out %>% 
  select(-adm_id, -mort_inunit, -age, -apache_II, 
         -chronic, -fractioninspiredoxygen) %>% 
  pivot_longer(2:18) %>% 
  group_by(row_id) %>% 
  summarise(n = sum(is.na(value)))

table(missing_hist$n)
  
# Write
results_out %>% 
  write_csv("data/apache_real_missing")


