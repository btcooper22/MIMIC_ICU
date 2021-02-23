# Setup----
require(readr)
require(dplyr)
require(tidyr)

# Load data
results <- read_csv("data/apache_data_full.csv") %>% 
  filter(!is.na(apache_II))

# Split data
admission_results <- results %>% 
  filter(type == "admission")

discharge_results <- results %>% 
  filter(type == "discharge")

# Describe data
mean(admission_results$readmission) * 100
mean(admission_results$mort_30) * 100
mean(admission_results$mort_inhosp) * 100
mean(admission_results$mort_inunit) * 100

# Build models----

# Readmission, admission score
model_admission_readmission <- glm(readmission ~ apache_II,
                                   data = admission_results,
                                   family = "binomial")

# Readmission, discharge score
model_discharge_readmission <- glm(readmission ~ apache_II,
                                   data = discharge_results,
                                   family = "binomial")

# 30-day mortality, admission score
model_admission_mort30 <- glm(mort_30 ~ apache_II,
                                   data = admission_results,
                                   family = "binomial")

# 30-day mortality, discharge score
model_discharge_mort30 <- glm(mort_30 ~ apache_II,
                              data = discharge_results,
                              family = "binomial")

# In-hospital mortality, admission score
model_admission_mortInHospital <- glm(mort_inhosp ~ apache_II,
                              data = admission_results,
                              family = "binomial")

# In-hospital mortality, discharge score
model_discharge_mortInHospital <- glm(mort_inhosp ~ apache_II,
                                      data = discharge_results,
                                      family = "binomial")

# In-unit mortality, admission score
model_admission_mortInUnit <- glm(mort_inunit ~ apache_II,
                                      data = admission_results,
                                      family = "binomial")
