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
require(mice)

source("functions/inverse_logit.R")

# Load and filter data
results <- read_csv("data/icnarc_predictors.csv") %>% 
  filter(surgical_type != "221. Paediatric Cardiac Surgery") %>% 
  select(-height, -weight, -surgical_type) %>% 
  select(-CV_support_prop,
         -total_support_prop,
         -admission_source,
         -age)

# Assess NA
sapply(results, function(x) sum(is.na(x))) %>% 
  as.data.frame() %>% 
  filter(. != 0)

# Impute missing
results_mice <- results %>% 
  mice(m = 10)

sum(complete(results_mice, 1)$readmission)

# Build model----

# Vector of acceptable features
feature_id <- 2:ncol(results)

model_options <- foreach(i = 1:results_mice$m, .combine = "rbind") %do%
  {
    print(i)
    
    # Binarise readmission column
    results_imputed <- complete(results_mice, i) %>%
      as_tibble() %>% 
      mutate(readmission = ifelse(readmission == TRUE, 1, 0),
             acute_renal_failure = as.logical(acute_renal_failure),
             glasgow_coma_below_15 = as.logical(glasgow_coma_below_15),
             hyperglycaemia = as.logical(hyperglycaemia),
             anaemia = as.logical(anaemia),
             sedation = as.logical(sedation),
             discharge_delay = as.logical(discharge_delay))
    
    # Select predictor and outcome vectors
    x <- model.matrix(readmission~., results_imputed[,feature_id])[,-1]
    y <- results_imputed$readmission
    
    # Determine lambda
    cv <- cv.glmnet(x, y, alpha = 1, folds = nrow(results_imputed),
                    family = "binomial", type.measure = "auc")
    
    # Fit model
    initial_model <- glmnet(x, y, alpha = 1, lambda = quantile(cv$lambda, 0.8),
                            family = "binomial")
    
    # Output
    output <- as.matrix(initial_model$beta) %>% 
      as.data.frame()
    output$var <- rownames(output)
    output$iteration <- i
    output
  }

# Identify retained variables
model_options %>%
  group_by(var) %>% 
  summarise(prop = mean(s0 != 0)) %>% 
  filter(prop == 1)

# Rebuild model on complete cases for relevant variables----

# Split data for training
set.seed(123)
patients_train <- results %>% 
  group_by(readmission) %>% 
  slice_sample(prop = 0.75) 
patients_validate <- results %>% 
  filter(id %in% patients_train$id == FALSE)

# Rebuild model
final_model <- glm(readmission ~ apache_II +
                     general_surgery + 
                     glasgow_coma_below_15 +
                     respiratory_support +
                     anaemia + lactate +
                     out_of_hours_discharge + 
                     total_support,
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
