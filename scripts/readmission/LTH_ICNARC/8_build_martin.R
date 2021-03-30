# Packages
require(readr)
require(dplyr)
require(magrittr)
require(ROCR)
require(ResourceSelection)
require(foreach)
require(stringr)

source("functions/inverse_logit.R")

# Load data
results <- read_csv("data/icnarc_predictors.csv") %>% 
  filter(surgical_type != "221. Paediatric Cardiac Surgery") %>% 
  select(id, readmission,age, sex, admission_source, apache_II,
         length_of_stay, out_of_hours_discharge, 
         acute_renal_failure, respiratory_rate,
         age, sodium, potassium, urea, glucose, sex, general_surgery,
         cardiac_surgery, hyperglycaemia,
         anaemia, apache_II, posthospital_dependency,
         length_of_stay) %>% 
  na.omit()
  
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

# Convert to probabilities
probs_martin <- inverse_logit(coefficients_martin)

# Assess discrimination
pred <- prediction(probs_martin, results$readmission)
performance(pred, measure = "auc")@y.values[[1]]

# Split to deciles
decile_df <- tibble(
  patient_id = 1:nrow(results),
  readmission = results$readmission,
  probs = probs_martin,
  decile = ntile(probs_martin, 10)) %>% 
  na.omit()

# Measure calibration
calibration_df <- decile_df %>% 
  select(readmission, probs, decile) %>% 
  group_by(decile) %>%
  summarise(N = length(readmission),
            observed = (sum(readmission) / length(readmission)) * 100,
            predicted = mean(probs * 100),
            error = sd(probs * 100))

# Assess calibration
hoslem_martin <- hoslem.test(results$readmission,
            probs_martin, g = 10)

list(model = "martin",
     data = "LTH_ICNARC",
     discrimination = pred,
     calibration = hoslem_martin,
     deciles = calibration_df) %>% 
  write_rds("scripts/readmission/shared/models/LTH_ICNARC_martin.RDS")
