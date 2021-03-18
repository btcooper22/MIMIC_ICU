# Packages
require(readr)
require(dplyr)
require(caret)
require(glmnet)
require(magrittr)
require(ROCR)
require(ResourceSelection)
require(tidyr)
require(foreach)
require(doParallel)

source("functions/inverse_logit.R")

# Load and filter data
results <- read_csv("data/icnarc_predictors.csv") %>% 
  filter(surgical_type != "221. Paediatric Cardiac Surgery") %>% 
  select(-height, -weight, -surgical_type)

# Assess NA
sapply(results, function(x) sum(is.na(x)))

# Take complete cases only
results_complete <- results %>% 
  select(-glucose,
         -hyperglycaemia,
         -haematocrit,
         -lactate, -follow_up_days,
         -CV_support_prop,
         -total_support_prop,
         -admission_source,
         -age) %>% 
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

# Determine lambda
cv <- cv.glmnet(x, y, alpha = 1, folds = nrow(results_complete),
                family = "binomial")

# Fit model
initial_model <- glmnet(x, y, alpha = 1, lambda = quantile(cv$lambda, 0.9),
                        family = "binomial")

# Identify retained variables
data.frame(varname = initial_model$beta@Dimnames[[1]],
           included = 1:length(initial_model$beta@Dimnames[[1]]) %in%
             (initial_model$beta@i + 1)) %>% 
  pivot_wider(names_from = "varname",
              values_from = "included") %>% 
  pivot_longer(1:34) %>% 
  group_by(name) %>% 
  summarise(prop = mean(value)) %>% 
  filter(prop != 0)

# Rebuild model on complete cases for relevant variables----

# Split data for training
set.seed(123)
patients_train <- results %>% 
  group_by(readmission) %>% 
  slice_sample(prop = 0.75) 
patients_validate <- results %>% 
  filter(id %in% patients_train$id == FALSE)

# Rebuild model
final_model <- glm(readmission ~ 
                     anaemia + apache_II + 
                     glasgow_coma_below_15 +
                     high_risk_speciality + length_of_stay +
                     out_of_hours_discharge + pulse_rate,
                     data = patients_train,
                   family = "binomial")

# Predict and assess----

# Create predictions
probs <- predict(final_model, newdata = patients_validate) %>% inverse_logit()

# Assess discrimination
pred <- prediction(probs[!is.na(probs)],
                   patients_validate$readmission[!is.na(probs)])
perf <- performance(pred, "tpr", "fpr")
performance(pred, measure = "auc")@y.values[[1]]

# Assess calibration
hoslem.test(patients_validate$readmission[!is.na(probs)],
            probs[!is.na(probs)], g = 10)
