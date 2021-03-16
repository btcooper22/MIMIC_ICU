# Packages
require(readr)
require(dplyr)
require(caret)
require(glmnet)
require(magrittr)

source("functions/inverse_logit.R")

# Load and filter data
results <- read_csv("data/icnarc_predictors.csv") %>% 
  filter(surgical_type != "221. Paediatric Cardiac Surgery") %>% 
  select(-height, -weight, -surgical_type)

# Assess NA
sapply(results, function(x) sum(is.na(x)))


# Take complete cases only
results_complete <- results %>% 
  select(-apache_II, -glucose,
         -hyperglycaemia,
         -haematocrit,
         -lactate, -follow_up_days,
         -CV_support_prop,
         -total_support_prop) %>% 
  na.omit()

sum(results_complete$readmission)

# Build model----

# Vector of acceptable features
feature_id <- 2:ncol(results_complete)

# Binarise readmission column
results_complete %<>%
  mutate(readmission = ifelse(readmission == TRUE, 1, 0))

# Select predictor and outcome vectors
x <- model.matrix(readmission~., results_complete[,feature_id])[,-1]
y <- results_complete$readmission

# Find optimal lambda
cv <- cv.glmnet(x, y, alpha = 1)

# Fit model
initial_model <- glmnet(x, y, alpha = 1, lambda = cv$lambda.min)
initial_model$beta

# Rebuild model on complete cases for relevant variables----

# Split data for training
set.seed(123)
patients_train <- results %>% 
  group_by(readmission) %>% 
  slice_sample(prop = 0.75) %>% 
  mutate(type = "train")

patients_validate <- results %>% 
  filter(id %in% patients_train$id == FALSE)

# Rebuild model
final_model <- glm(readmission ~ glasgow_coma_below_15 + 
                     respiratory_support +
                     days_before_ICU + high_risk_speciality +
                     out_of_hours_discharge, data = patients_train,
                   family = "binomial")

# Predict and assess----

# Create predictions
probs <- predict(final_model, newdata = patients_validate) %>% inverse_logit()

# Assess discrimination

# Assess calibration