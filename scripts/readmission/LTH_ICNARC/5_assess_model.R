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
  select(id, readmission, anaemia,  
         glasgow_coma_below_15,
         high_risk_speciality,
         out_of_hours_discharge,
         respiratory_support,
         clinically_ready_discharge_days) %>% 
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
                         anaemia + respiratory_support +
                         high_risk_speciality + 
                         clinically_ready_discharge_days +
                         out_of_hours_discharge,
                       data = patients_train,
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
            cal = median(cal_chisq),
            cal_error_low = quantile(cal_chisq, 0.25),
            cal_error_high = quantile(cal_chisq, 0.75))
# AUC 0.629 pm 0.053
# Cal 9.59 [6.81, 13.88]

# AUC  AUC_error      cal cal_error_low cal_error_high
# 1 0.6182097 0.05428703 9.877884      6.820818       14.20504

# Extract median discrimination and calibration data----

# Discrimination
set.seed(which.min(abs(output$disc - median(output$disc))))

# Split data
patients_train <- results %>% 
  group_by(readmission) %>% 
  slice_sample(prop = 0.75)

patients_validate <- results %>% 
  filter(id %in% patients_train$id == FALSE)

# Rebuild model
final_model <- glm(readmission ~ glasgow_coma_below_15 +
                     anaemia + respiratory_support +
                     high_risk_speciality + 
                     clinically_ready_discharge_days +
                     out_of_hours_discharge,
                   data = patients_train,
                   family = "binomial")

# Create predictions
probs <- predict(final_model, newdata = patients_validate) %>% inverse_logit()

# Assess discrimination
pred <- prediction(probs[!is.na(probs)],
                   patients_validate$readmission[!is.na(probs)])

# Calibration
set.seed(which.min(abs(output$cal_chisq - median(output$cal_chisq, na.rm = T))))

# Split data
patients_train <- results %>% 
  group_by(readmission) %>% 
  slice_sample(prop = 0.75)

patients_validate <- results %>% 
  filter(id %in% patients_train$id == FALSE)

# Rebuild model
final_model <- glm(readmission ~ glasgow_coma_below_15 +
                     anaemia + respiratory_support +
                     high_risk_speciality + 
                     clinically_ready_discharge_days +
                     out_of_hours_discharge,
                   data = patients_train,
                   family = "binomial")

# Create predictions
probs <- predict(final_model, newdata = patients_validate) %>% inverse_logit()

# Split to deciles
decile_df <- tibble(
  patient_id = 1:nrow(patients_validate),
  readmission = patients_validate$readmission,
  probs,
  decile = ntile(probs, 10)) %>% 
  na.omit()

# Measure calibration
calibration_df <- decile_df %>% 
  select(readmission, probs, decile) %>% 
  group_by(decile) %>%
  summarise(N = length(readmission),
            observed = (sum(readmission) / length(readmission)) * 100,
            predicted = mean(probs * 100),
            error = sd(probs * 100))

# Quantify calibration
cal_test <- hoslem.test(patients_validate$readmission[!is.na(probs)],
            probs[!is.na(probs)], g = 10)

# Combine and write
list(model = "bespoke",
     data = "LTH_ICNARC",
     discrimination = pred,
     calibration = cal_test,
     deciles = calibration_df) %>% 
  write_rds("scripts/readmission/shared/models/LTH_ICNARC_bespoke.RDS")
