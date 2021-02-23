# Setup----
require(readr)
require(dplyr)
require(tidyr)
require(ROCR)

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

# Generate AUC values----

# Generate prediction objects
prediction_admission_readmission <- prediction(probs_admission_readmission$probs,
                                               probs_admission_readmission$outcome)

# Generate AUC objects
