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
require(stringr)

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

# Generate hammer scores----
coefficients_hammer <- log(0.005 / (1 - 0.005)) +
  ifelse(results$sex == "M", 0.43, 0) +
  ifelse(results$general_surgery, 0.64, 0) +
  ifelse(results$cardiac_surgery, -1.01, 0) +
  ifelse(results$hyperglycaemia, 0.36, 0) +
  ifelse(results$anaemia, 1.11, 0) +
  ifelse(results$apache_II >= 20, 0.44, 0) +
  ifelse(results$posthospital_dependency == TRUE, 0.7, 0) +
  ifelse(results$length_of_stay >= 5, 0.88, 0) + 0.26

# Generate Martin scores----
# Rebuild age
age_band <- results$age
results$age <- case_when(
  age_band == "under 25" ~ 21 %>% as.character(),
  age_band == "over 90" ~ 90 %>% as.character(),
  TRUE ~ substr(age_band, 1, 2)
) %>% as.numeric()


# Estimate serum chloride from sodium and potassium
chloride_coef <- read_rds("models/serum_chloride.RDS")

results$serum_chloride <- chloride_coef[1] + (chloride_coef[2] * results$sodium) + 
  (chloride_coef[3] * results$potassium) + 
  (chloride_coef[4] * (results$sodium * results$potassium))

# Convert urea to BUN
results$blood_urea_nitrogen <- results$urea / 0.357

# Extract medical history
medical_history <- read_csv("data/icnarc_outcomes.csv") %>% 
  select(contains("PastMedicalHistory")) %>% 
  select("PastMedicalHistory_VerySevereCardiovascularDisease",
         "PastMedicalHistory_ChronicRenalReplacement",
         "PastMedicalHistory_Condition") %>% 
  cbind(read_csv("data/icnarc_outcomes.csv") %>% 
          select(Identifiers_PatientPseudoId) %>% 
          rename(id = "Identifiers_PatientPseudoId")) %>% 
  mutate(renal_insufficiency = PastMedicalHistory_ChronicRenalReplacement == TRUE |
           PastMedicalHistory_Condition == "Chronic renal failure",
         atrial_fibrillation = PastMedicalHistory_VerySevereCardiovascularDisease == TRUE |
           PastMedicalHistory_Condition == "Supra-ventricular tachycardia, atrial fibrillation or flutter") %>% 
  select(id, renal_insufficiency, atrial_fibrillation)

# Attach
results %<>% 
  left_join(medical_history)
results$renal_insufficiency[is.na(results$renal_insufficiency)] <- FALSE
results$atrial_fibrillation[is.na(results$atrial_fibrillation)] <- FALSE

# Score data - coefficients
coefficients_martin <- -9.284491 +
  (0.04883 * results$respiratory_rate) +
  (0.011588 * results$age) +
  (0.036104 * results$serum_chloride) +
  (0.004967 * results$blood_urea_nitrogen) +
  (0.580153 * results$atrial_fibrillation) +
  (0.458202 * results$renal_insufficiency) +
  (0.003519 * results$glucose)

# Bootstrap assessment----

# Prepare data
results %<>% 
  mutate(probs_frost = log(probs_frost / (1 - probs_frost)),
         probs_hammer = coefficients_hammer,
         probs_martin = coefficients_martin)

# Prepare parallel
ptm <- proc.time()
cl <- makeCluster(4)
registerDoParallel(cl)
nboot <- 100000

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
    frost_plus_model <- glm(readmission ~ probs_frost + clinically_ready_discharge_days + 0,
                            data = patients_train,
                       family = "binomial")
    
    # Fit hammer+
    hammer_plus_model <- glm(readmission ~ probs_hammer + clinically_ready_discharge_days + 0,
                            data = patients_train,
                            family = "binomial")
    
    # Fit martin+
    martin_plus_model <- glm(readmission ~ probs_martin + clinically_ready_discharge_days + 0,
                             data = patients_train,
                             family = "binomial")
    
    
    # Predict
    probs_frost_plus <- predict(frost_plus_model, newdata = patients_validate) %>% inverse_logit()
    probs_hammer_plus <- predict(hammer_plus_model, newdata = patients_validate) %>% inverse_logit()
    probs_martin_plus <- predict(martin_plus_model, newdata = patients_validate) %>% inverse_logit()
    
    # Make prediction object
    pred_frost <- prediction(patients_validate$probs_frost %>% inverse_logit(),
                             patients_validate$readmission)
    pred_frost_plus <- prediction(probs_frost_plus, patients_validate$readmission)
    pred_hammer <- prediction(patients_validate$probs_hammer %>% inverse_logit(),
                             patients_validate$readmission)
    pred_hammer_plus <- prediction(probs_hammer_plus, patients_validate$readmission)
    pred_martin <- prediction(patients_validate$probs_martin %>% inverse_logit(),
                              patients_validate$readmission)
    pred_martin_plus <- prediction(probs_martin_plus, patients_validate$readmission)
    
    # Extract AUC
    AUC_frost <- performance(pred_frost, measure = "auc")@y.values[[1]]
    AUC_frost_plus <- performance(pred_frost_plus, measure = "auc")@y.values[[1]]
    AUC_hammer <- performance(pred_hammer, measure = "auc")@y.values[[1]]
    AUC_hammer_plus <- performance(pred_hammer_plus, measure = "auc")@y.values[[1]]
    AUC_martin <- performance(pred_martin, measure = "auc")@y.values[[1]]
    AUC_martin_plus <- performance(pred_martin_plus, measure = "auc")@y.values[[1]]
    
    # Calibration
    cal_frost <- hoslem.test(patients_validate$readmission,
                             patients_validate$probs_frost %>% inverse_logit(),
                             g = 10)
    cal_frost_plus <- hoslem.test(patients_validate$readmission,
                             probs_frost_plus, g = 10)
    cal_hammer <- hoslem.test(patients_validate$readmission,
                             patients_validate$probs_hammer %>% inverse_logit(),
                             g = 10)
    cal_hammer_plus <- hoslem.test(patients_validate$readmission,
                                  probs_hammer_plus, g = 10)
    cal_martin <- hoslem.test(patients_validate$readmission,
                              patients_validate$probs_martin %>% inverse_logit(),
                              g = 10)
    cal_martin_plus <- hoslem.test(patients_validate$readmission,
                                   probs_martin_plus, g = 10)
    
    
    # Output
    data.frame(i, AUC_frost,
               AUC_frost_plus,
               cal_frost = cal_frost$statistic,
               cal_frost_plus = cal_frost_plus$statistic,
               AUC_hammer,
               AUC_hammer_plus,
               cal_hammer = cal_hammer$statistic,
               cal_hammer_plus = cal_hammer_plus$statistic,
               AUC_martin,
               AUC_martin_plus,
               cal_martin = cal_martin$statistic,
               cal_martin_plus = cal_martin_plus$statistic)
  }

write_rds(model_results, "models/context_bootstrap.RDS",
          compress = "gz")
proc.time() - ptm 
stopCluster(cl)


# Plot----
model_results <- read_rds("models/context_bootstrap.RDS")

# Discrimination
model_results %>% 
  mutate(frost_AUC_diff = AUC_frost_plus - AUC_frost,
         hammer_AUC_diff = AUC_hammer_plus - AUC_hammer,
         martin_AUC_diff = AUC_martin_plus - AUC_martin) %>%
  select(frost_AUC_diff, hammer_AUC_diff, martin_AUC_diff) %>% 
  pivot_longer(1:3) %>% 
  ggplot(aes(x = value))+
  geom_density(fill = "#e41a1c", alpha = 0.8)+
  geom_vline(linetype = "dashed",
             xintercept = 0)+
  theme_classic(20)+
    facet_wrap(~name, scales = "free")

# Calibration
model_results %>% 
  mutate(frost_cal_diff = cal_frost_plus - cal_frost,
         hammer_cal_diff = cal_hammer_plus - cal_hammer,
         martin_cal_diff = cal_martin_plus - cal_martin) %>% 
  select(frost_cal_diff, hammer_cal_diff, martin_cal_diff) %>% 
  pivot_longer(1:3) %>% 
  ggplot(aes(x = value))+
  geom_density(fill = "#377eb8", alpha = 0.8)+
  geom_vline(linetype = "dashed",
             xintercept = 0)+
  theme_classic(20)+
  facet_wrap(~name, scales = "free")

# Summarise
model_results %>% 
  mutate(cal_frost_diff = cal_frost_plus - cal_frost,
         cal_hammer_diff = cal_hammer_plus - cal_hammer,
         cal_martin_diff = cal_martin_plus - cal_martin,
         AUC_frost_diff = AUC_frost_plus - AUC_frost,
         AUC_hammer_diff = AUC_hammer_plus - AUC_hammer,
         AUC_martin_diff = AUC_martin_plus - AUC_martin) %>% 
  pivot_longer(2:19) %>% 
  group_by(name) %>% 
  summarise(median = median(value, na.rm = TRUE),
            Q1 = quantile(value, 0.25, na.rm = TRUE),
            Q3 = quantile(value, 0.75, na.rm = TRUE)) %>% 
  mutate(model = str_split(name, "_", simplify = TRUE)[,2],
         measure = str_split(name, "_", simplify = TRUE)[,1],
         type = case_when(
           grepl("plus", name) ~ "plus",
           grepl("diff", name) ~ "diff",
           is.character(name) ~ "standard")
         ) %>% 
  mutate(value = paste(signif(median, 3), " [",
                       signif(Q1, 3), ", ",
                       signif(Q3, 3), "]", sep = "")) %>% 
  select(-median, -Q1, -Q3, -name) %>% 
  pivot_wider(names_from = c("measure", "type"),
              values_from = "value") %>% 
  relocate("model", "AUC_standard", "AUC_plus", "AUC_diff",
           "cal_standard", "cal_plus", "cal_diff")
