# Packages
require(readr)
require(dplyr)
require(magrittr)
require(ROCR)
require(ResourceSelection)
require(foreach)

source("functions/inverse_logit.R")

# Load and filter data
results <- read_csv("data/icnarc_predictors.csv") %>% 
  filter(surgical_type != "221. Paediatric Cardiac Surgery") %>% 
  select(-height, -weight, -surgical_type) %>% 
  select(id, readmission, glasgow_coma_below_15,
         respiratory_support, high_risk_speciality,
         out_of_hours_discharge, days_before_ICU) %>% 
  na.omit()


# Loop through datasets
ptm <- proc.time()
output <- foreach(i = 1:10000, .combine = "rbind") %do%
  {
    # Seed
    set.seed(i)
    
    # Split data
    patients_train <- results %>% 
      group_by(readmission) %>% 
      slice_sample(prop = 0.75)
    
    patients_validate <- results %>% 
      filter(id %in% patients_train$id == FALSE)
    
    # Rebuild model
    final_model <- glm(readmission ~ glasgow_coma_below_15 + 
                         respiratory_support +
                         days_before_ICU + high_risk_speciality +
                         out_of_hours_discharge, data = patients_train,
                       family = "binomial")
    
    
    # Create predictions
    probs <- predict(final_model, newdata = patients_validate) %>% inverse_logit()
    
    # Assess discrimination
    pred <- prediction(probs[!is.na(probs)],
                       patients_validate$readmission[!is.na(probs)])
    AUC <- performance(pred, measure = "auc")@y.values[[1]]
    
    # Assess calibration
    cal <- hoslem.test(patients_validate$readmission[!is.na(probs)],
                probs[!is.na(probs)], g = 10)
    
    # Output
    data.frame(i, disc = AUC,
               cal_chisq = cal$statistic,
               cal_p = cal$p.value)
  }
proc.time() - ptm # 4 minutes (220 seconds)

# Summarise
output %>% 
   na.omit() %>% 
  summarise(AUC = mean(disc),
            AUC_error = sd(disc),
            cal = mean(cal_chisq),
            cal_error = sd(cal_chisq))
# AUC 0.658 (0.044)
# Cal 8.79 (4.68)
