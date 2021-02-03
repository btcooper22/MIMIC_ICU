# Packages
require(dplyr)
require(magrittr)
require(readr)
require(tibble)
require(ggplot2)
require(ROCR)

source("functions/inverse_logit.R")
source("functions/calibration.R")

# Load data----- 
# Trim to APACHE-II
full_data <- read_csv("data/final_patients.csv") %>% 
  select(matches("_discharge$")) %>% 
  select(matches("^apache_")) %>% 
  cbind(read_csv("data/final_patients.csv") %>% 
          select(row_id, subject_id,
                 adm_id, readmission,
                 age),.)

# Binarise readmission
full_data %<>% 
  mutate(readmission = readmission == "Readmitted to ICU")

# Confirm no missing values
table(full_data$readmission, useNA = "ifany")
summary(full_data$apache_II_discharge)

# Create and assess APACHE model------
# Create model and predict
apache_model <- glm(readmission ~ apache_II_discharge,
                    data = full_data, family = "binomial")
probs <- predict(apache_model, newdata = full_data) %>% inverse_logit()

# Assess ROC
pred <- prediction(probs, full_data$readmission)
perf <- performance(pred, "tpr", "fpr")
auc <- performance(pred, measure = "auc")
auc@y.values[[1]]

# Plot ROC
data.frame(x = perf@x.values[[1]],
           y = perf@y.values[[1]]) %>% 
  ggplot(aes(x, y))+
  geom_abline(slope = 1, intercept = 0,
              linetype = "dotted",
              size = 1)+
  geom_path(size = 1)+
  labs(x = "1 - Specificity",
       y = "Sensitivity")+
  theme_classic(20)

# Split into deciles
deciles_df <- tibble(
  patient_id = 1:nrow(full_data),
  readmission = full_data$readmission,
  probs,
  decile = ntile(probs, 10))

# Calculate calibration
cal <- calibration(deciles_df, "probs", "decile")

# Plot calibration
cal %>% 
  ggplot(aes(predicted, observed))+
  geom_abline(slope = 1, intercept = 0,
              size = 1)+
  geom_path(size = 1)+
  geom_pointrange(aes(ymin = observed - error ,
                      ymax = observed + error,
                      y = observed,
                      x = predicted)) +
  theme_classic(20)+
  theme(legend.position = "top")+
  labs(x = "Predicted readmission",
       y = "Observed readmission")

# Simulate MCAR (missing completely at random)---------

# Basic methods -> Mean, MICE, KNN MissForest

# Complex methods from competiton