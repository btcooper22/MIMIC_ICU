# Packages
require(readr)
require(magrittr)
require(dplyr)
require(tibble)
require(stringr)
require(lubridate)

# Load data
df <- read_csv("data/icnarc_outcomes.csv")

# Define functions
quicktab <- function(varname)
{
  results %>% 
    group_by_at(all_of(varname)) %>% 
      summarise(n = n(), percent = mean(readmission) * 100)
}

unfold <- function(x, repeats = 2)
{
    x[is.na(x)] <- -1
    score <- data.frame(x, id = rep(df$Identifiers_PatientPseudoId, repeats)) %>%
      group_by(id) %>%
      summarise(score = max(x, na.rm = T)) %>%
      select(score) %>% deframe()
    
    score[score == -1] <- NA
    return(score)
}

# Prepare results
results <- data.frame(id = df$Identifiers_PatientPseudoId,
                      readmission = df$readmission)

# Acute renal failure----
# Extract urine and creatinine
creatinine <- ifelse(is.na(df$SerumCreatinine_High),
                     df$SerumCreatinine_Low, df$SerumCreatinine_High) / 88.42

creatinine <- ifelse(is.na(creatinine), df$SerumCreatinine_PreAdmit / 88.42,
       creatinine)

urine <- df$UrineOutput_Total1st24Hrs

# Load prediction coefficients
coef_arf <- read_rds("models/acute_renal_failure.RDS")

# Predict ARF
probs_arf <- coef_arf[1] + (coef_arf[2] * creatinine) + (coef_arf[3] * urine) + 
  (coef_arf[4] * (creatinine * urine))
results$acute_renal_failure <- probs_arf > coef_arf[5]
with(results, table(readmission, acute_renal_failure))

# Medical history vars----
medical_history <- df %>% 
  select(contains("PastMedicalHistory")) %>% 
  select(2:19)

results$poor_chronic_health <- rowSums(medical_history) >0
with(results, table(readmission, poor_chronic_health))

# APACHE-II----

# Temperature
temp_low <- df$Temperature_LowCentral
temp_high <- df$Temperature_HighCentral

temp_all <- c(temp_low, temp_high)

temperature_scores <- case_when(
  temp_all >= 41 | temp_all <= 29.9 ~ 4,
  temp_all >= 39 | temp_all <= 31.9 ~ 3,
  temp_all <= 33.9 ~ 2,
  temp_all >= 38.5 | temp_all <= 35.9 ~ 1,
  is.na(temp_all) ~ NaN,
  is.numeric(temp_all) ~ 0,
) %>% unfold(2)

# Mean arterial pressure
map_low <- (df$BloodPressure_LowSystolic + (2 * df$BloodPressure_DiastolicPairedLowSystolic) )/ 3
map_high <- (df$BloodPressure_HighSystolic + (2 * df$BloodPressure_DiastolicPairedHighSystolic) )/ 3

map_full <- c(map_low, map_high)

map_scores <- case_when(
  map_full >= 160 | map_full <= 49 ~ 4,
  map_full >= 130 ~ 3,
  map_full >= 110 | map_full <= 69 ~ 2,
  is.na(map_full) ~ NaN,
  is.numeric(map_full) ~ 0
) %>% unfold()

# Heart rate
pulse_low <- df$VentricularRate_Low
pulse_high <- df$VentricularRate_High

pulse_full <- c(pulse_low, pulse_high)

pulse_score <- case_when(
  pulse_full >= 180 | pulse_full <= 39 ~ 4,
  pulse_full >= 140 | pulse_full <= 54 ~ 3,
  pulse_full >= 110 | pulse_full <= 69 ~ 2,
  is.na(pulse_full) ~ NaN,
  is.numeric(pulse_full) ~ 0
) %>% unfold()

# Respiratory rate
respiratory_full <- c(
  df$RespiratoryRate_HighNonventilated,
  df$RespiratoryRate_LowNonventilated,
  df$RespiratoryRate_LowVentilated,
  df$RespiratoryRate_HighVentilated
)

respiratory_score <-  case_when(
  respiratory_full >= 50 | respiratory_full <= 5 ~ 4,
  respiratory_full >= 35 ~ 3,
  respiratory_full <= 9 ~ 2,
  respiratory_full >= 25 | respiratory_full <= 11 ~ 1,
  is.na(respiratory_full) ~ NaN,
  is.numeric(respiratory_full) ~ 0
) %>% unfold(4)

# Arterial pH
artpH_full <- df$ArterialBloodGases_LowestPaO2_HPaH

artpH_score <-  case_when(
  artpH_full >= 7.7 | artpH_full < 7.15 ~ 4,
  artpH_full >= 7.6 | artpH_full <= 7.24 ~ 3,
  artpH_full <= 7.32 ~ 2,
  artpH_full >= 7.5 ~ 1,
  is.na(artpH_full) ~ NaN,
  is.numeric(artpH_full) ~ 0
)

# Sodium
sodium_full <- c(df$SerumSodium_High, df$SerumSodium_Low,
                 df$SerumSodium_PreAdmit)

sodium_score <-  case_when(
  sodium_full >= 180 | sodium_full <= 110 ~ 4,
  sodium_full >= 160 | sodium_full <= 119 ~ 3,
  sodium_full >= 155 | sodium_full <= 129 ~ 2,
  sodium_full >= 150 ~ 1,
  is.na(sodium_full) ~ NaN,
  is.numeric(sodium_full) ~ 0
) %>% unfold(3)


# Potassium
potassium_full <- c(df$SerumPotassium_Low, df$SerumPotassium_High,
                    df$SerumPotassium_PreAdmit)

potassium_score <- case_when(
  potassium_full >= 7 | potassium_full < 2.5 ~ 4,
  potassium_full >= 6 ~ 3,
  potassium_full <= 2.9 ~ 2,
  potassium_full >= 5.5 | potassium_full <= 3.4 ~ 1,
  is.na(potassium_full) ~ NaN,
  is.numeric(potassium_full) ~ 0
) %>% unfold(3)

# Creatinine
creatinine_full <- c(df$SerumCreatinine_High, df$SerumCreatinine_Low,
                     df$SerumCreatinine_PreAdmit) / 88.42

creatinine_score <- case_when(
  creatinine_full >= 3.5 ~ 4,
  creatinine_full >= 2 ~ 3,
  creatinine_full >= 1.5 | creatinine_full < 0.6 ~ 2,
  is.na(creatinine_full) ~ NaN,
  is.numeric(creatinine_full) ~ 0
) %>% unfold(3)

creatinine_score <- ifelse(results$acute_renal_failure,
                           creatinine_score * 2, 
                           creatinine_score)

# Haematocrit
haematocrit_full <- c(df$Haemoglobin_High, df$Haemoglobin_Low,
                        df$Haemoglobin_PreAdmit) * 3

haematocrit_score <- case_when(
  haematocrit_full >= 60 | haematocrit_full < 20 ~ 4,
  haematocrit_full >= 50 | haematocrit_full <= 29.9 ~ 2,
  haematocrit_full >= 46 ~ 1,
  is.na(haematocrit_full) ~ NaN,
  is.numeric(haematocrit_full) ~ 0
) %>% unfold(3)

# White blood cell count
wbc_full <- c(df$WhiteCellCount_High, df$WhiteCellCount_Low,
              df$WhiteCellCount_PreAdmit)

wbc_score <- case_when(
  wbc_full >= 40 | wbc_full < 1 ~ 4,
  wbc_full >= 20 | wbc_full <= 2.9 ~ 2,
  wbc_full >= 15 ~ 1,
  is.na(wbc_full) ~ NaN,
  is.numeric(wbc_full) ~ 0
) %>% unfold(3)

# PaO2
pao2_full <- df$ArterialBloodGases_LowestPaO2_Any / 0.133333

pao2_score <- case_when(
  pao2_full < 55 ~ 4,
  pao2_full <= 60 ~ 3,
  pao2_full <= 70 ~ 1,
  is.na(pao2_full) ~ NaN,
  is.numeric(pao2_full) ~ 0
)

# Alveolar-arterial gradient
elev <-  63
patm <- 760 * exp(elev / -7000)
ph2o <- 47 * exp((median(temp_all, na.rm = T)- 37) / 18.4)
fio2 <- df$ArterialBloodGases_LowestPaO2_O2PerCent /100
aa_grad <- (fio2 * (patm - ph2o) - ((df$ArterialBloodGases_LowestPaO2_LowestPaCO2 / 0.13333) / 0.8)) -
  pao2_full

aagrad_score <- case_when(
  aa_grad >= 500 ~ 4,
  aa_grad >= 350 ~ 3,
  aa_grad >= 200 ~ 2,
  is.na(aa_grad) ~ NaN,
  is.numeric(aa_grad) ~ 0
)

# If FiO2 > 0.5, use aa_grad for oxygenation score, else use PaO2
oxygenation_score <- ifelse(fio2 > 0.5, aagrad_score, pao2_score)

# Bicarbonate
bicarb_full <- c(df$SerumBicarbonate_High, df$SerumBicarbonate_Low,
                 df$SerumBicarbonate_PreAdmit)

bicarbonate_score <- case_when(
  bicarb_full >= 52 | bicarb_full < 15 ~ 4,
  bicarb_full >= 41 | bicarb_full <= 17.9 ~ 3,
  bicarb_full <= 21.9 ~ 2,
  bicarb_full >= 32 ~ 1,
  is.na(bicarb_full) ~ NaN,
  is.numeric(bicarb_full) ~ 0
) %>% unfold(3)

# Glasgow coma scale
gcs_extract <- str_split(df$Sedation_GcsLowestElements, "",
                         simplify = TRUE) %>% 
  as.data.frame() %>% 
  mutate(sum = as.numeric(V1) + 
           as.numeric(V2) + 
           as.numeric(V3))

gcs_score <- 15 - gcs_extract$sum

# Age
age_band <- df$Demographics_UnitAdmissionAgeBand

age_score <- case_when(
  age_band %in% c("66 to 70", "71 to 75",
                  "76 to 80", "81 to 85",
                  "86 to 90", "over 90") ~ 5,
  age_band %in% c("56 to 60", "61 to 65") ~ 3,
  age_band %in% c("46 to 50", "51 to 55") ~ 2,
  is.character(age_band) ~ 0
)

# Chronic health
chronic_score <- ifelse(results$poor_chronic_health, 2, 0)

# If arterial blood gasses missing, use bicarbonate instead
abg_score <- 
  ifelse(df$ArterialBloodGases_NotAvailable, bicarbonate_score,
         oxygenation_score + artpH_score)

# Set NA scores to 0
abg_score[is.na(abg_score)] <- 0
sodium_score[is.na(sodium_score)] <- 0
potassium_score[is.na(potassium_score)] <- 0
creatinine_score[is.na(creatinine_score)] <- 0
wbc_score[is.na(wbc_score)] <- 0
gcs_score[is.na(gcs_score)] <- 0
haematocrit_score[is.na(haematocrit_score)] <- 0

# Sum scores
apache_scores <- temperature_scores + map_scores + pulse_score +
  respiratory_score + abg_score + sodium_score + potassium_score +
  creatinine_score + wbc_score + haematocrit_score + gcs_score +
  age_score + chronic_score

# Add to results
results %<>% 
  mutate(apache_II = apache_scores)


# Apache vars individually----
worst_value <- function(val, bad_direction)
{
  if(all(is.na(val)))
  {
    return(NA)
  }else 
  {
    if(bad_direction == "max")
    {
      return(max(val, na.rm = TRUE))
    }else if(bad_direction == "min")
    {
      return(min(val, na.rm = TRUE))
    } else
    {
      return(val[which.max(abs(val - bad_direction))])
    }
  }
}

unfold_value <- function(x, repeats = 2, bd)
{
  score <- data.frame(x, id = rep(df$Identifiers_PatientPseudoId, repeats)) %>%
    group_by(id) %>%
    summarise(score = worst_value(x, bd)) %>%
    select(score) %>% deframe()
  return(score)
}

# Temperature
results$temperature <- unfold_value(temp_all, bd = 37.2)

# Mean arterial pressure
results$mean_arterial_pressure <- unfold_value(map_full, bd = "max")

# Pulse
results$pulse_rate <- unfold_value(pulse_full, bd = 90)

# Respiratory
results$respiratory_rate <- unfold_value(respiratory_full, 4, "max")

# Sodium
results$sodium <- unfold_value(sodium_full, 3, 140)

# Potassium
results$potassium <- unfold_value(potassium_full, 3, 4.5)

# Creatinine
results$creatinine <- creatinine

# Haematocrit
results$haematocrit <- unfold_value(haematocrit_full, 3, 44)

# White cell count
results$white_cell_count <- unfold_value(wbc_full, 3, "max")

# Glasgow coma scale
results$glasgow_coma_below_15 <- gcs_extract$sum < 15

# Age
results$age_over_65 <- age_band %in% c("66 to 70", "71 to 75",
                                       "76 to 80", "81 to 85",
                                       "86 to 90", "over 90")

# Other variables----

# Glucose
results$glucose <-  case_when(
  is.na(df$SerumGlucose_High) & 
    is.na(df$SerumGlucose_Low) ~ df$SerumGlucose_PreAdmit,
  is.na(df$SerumGlucose_High) ~ df$SerumGlucose_Low,
  is.numeric(df$SerumGlucose_High) ~ df$SerumGlucose_High
) * 18

# Hyperglycaemia
results$hyperglycaemia <- results$glucose > 180 

# Anaemia
results$anaemia <- case_when(
  is.na(df$Haemoglobin_Low) & 
    is.na(df$Haemoglobin_High) ~ df$Haemoglobin_PreAdmit,
  is.na(df$Haemoglobin_Low) ~ df$Haemoglobin_High,
  is.numeric(df$Haemoglobin_Low) ~ df$Haemoglobin_Low
) < 9

# Lactate
results$lactate <- ifelse(is.na(df$SerumLactate_High), df$SerumLactate_PreAdmit,
       df$SerumLactate_High)

# Urea
results$urea <- ifelse(is.na(df$SerumUrea_High), df$SerumUrea_PreAdmit,
                          df$SerumUrea_High)
                             
# Platelets  
results$platelets <- ifelse(is.na(df$Platelets_Low), df$Platelets_PreAdmit,
                       df$Platelets_Low)

# Organ support----

# > 0 days on resp. support
results$respiratory_support <-  (df$OrganSupport_Ccmds_BasicRespiratoryDays + 
  df$OrganSupport_Ccmds_AdvancedRepiratoryDays) > 0

# > 0  days on advanced CV support
results$advanced_cv_support <- df$OrganSupport_Ccmds_AdvancedCardiovascularDays > 0

# Days on CV support
results$CV_support <- df$OrganSupport_Ccmds_BasicCardiovascularDays

# Basic CV support for  % of stay
results$CV_support_prop <- df %>% 
  mutate(CV_support_prop = OrganSupport_Ccmds_BasicCardiovascularDays /
           UnitDischarge_DischargedDiedOnDays) %>% 
  mutate(CV_support_prop = case_when(is.infinite(CV_support_prop) ~ 1,
                                     is.numeric(CV_support_prop) ~ CV_support_prop)) %>% 
  select(CV_support_prop) %>% deframe()

# Total support
results$total_support <- df %>% 
  mutate(total_support = OrganSupport_Ccmds_BasicRespiratoryDays +
           OrganSupport_Ccmds_AdvancedRepiratoryDays +
           OrganSupport_Ccmds_BasicCardiovascularDays +
           OrganSupport_Ccmds_AdvancedCardiovascularDays +
           OrganSupport_Ccmds_RenalDays +
           OrganSupport_Ccmds_GastrointestinalDays +
           OrganSupport_Ccmds_DermatologicalDays) %>% 
  select(total_support) %>% deframe()

# Total support as % of stay
results$total_support_prop <- df %>% 
  mutate(total_support = OrganSupport_Ccmds_BasicRespiratoryDays +
           OrganSupport_Ccmds_AdvancedRepiratoryDays +
           OrganSupport_Ccmds_BasicCardiovascularDays +
           OrganSupport_Ccmds_AdvancedCardiovascularDays +
           OrganSupport_Ccmds_RenalDays +
           OrganSupport_Ccmds_GastrointestinalDays +
           OrganSupport_Ccmds_DermatologicalDays) %>% 
  mutate(total_support_prop = total_support /
           UnitDischarge_DischargedDiedOnDays) %>% 
  mutate(total_support_prop = case_when(is.infinite(total_support_prop) ~ total_support,
                                     is.numeric(total_support_prop) ~ total_support_prop)) %>% 
  select(total_support_prop) %>% deframe()

# Any organ support required
results$organ_support <- results$total_support > 0

# Simple variables----

# Sex
results$sex <- df$Demographics_Sex

# BMI
results$height <- df$Demographics_UnitAdmissionHeightCm
results$weight <- df$Demographics_UnitAdmissionWeightKg
results$BMI <-  df$Demographics_UnitAdmissionWeightKg / 
  ((df$Demographics_UnitAdmissionHeightCm / 100) ^ 2)

# LOS
results$length_of_stay <- df$UnitDischarge_DischargedDiedOnDays

# Time in hospital before ICU
results$days_before_ICU <- abs(df$PriorToAdmission_HospitalAdmissionDaysBefore)

# Surgery types
results$surgical_type <- df$PriorToAdmission_Specialty
results$cardiac_surgery <- results$surgical_type %in% c("172. Cardiac Surgery", 
                             "170. Cardiothoracic Surgery", "320. Cardiology")
results$high_risk_speciality <- results$surgical_type %in%
  c("106. Upper Gastrointestinal Surgery", "104. Colorectal Surgery",
    "105. Hepatobiliary & Pancreatic Surgery", "306. Hepatology",
    "173. Thoracic Surgery", "107. Vascular Surgery",
    "301. Gastroenterology")

# Diabetes
results$diabetes <- grepl("diabetes", df$PastMedicalHistory_Condition,
      ignore.case = TRUE)

# Prehospital dependency
results$prehospital_dependency <-  df$PastMedicalHistory_PrehospitalDependency != 
  "A. Able to live without assistance in daily activities"

# Sedation
results$sedation <- df$Sedation != "5. No sedation or paralysis during scored period"

# Clinically ready for discharge days == 0
results$clinically_ready_discharge_days <-  df$UnitDischarge_ClinicallyReadyForDischargeDays

# Out-of-hours discharge
discharge_time <- substr(df$UnitDischarge_DischargedDiedAtTime, 1, 2) %>% 
  as.numeric()
results$out_of_hours_discharge <-  discharge_time < 8 | discharge_time > 17

# Some expected dependency post-discharge
results$posthospital_dependency <- df$UnitDischarge_ExpectedDependencyPostDischarge != 
  "A. Able to live without assistance in daily activities"

# Follow up care days >5
results$follow_up_days <- df$FollowUpCare_FirstNumberOfDays

# Admission source
results$admission_source <- ifelse(df$Source_AdmitFrom == "B. Recovery/theatre (following surgery/anaesthetic procedure)",
                                   "OT", "Ward")

# Write results
results %>% 
  write_csv("data/icnarc_predictors.csv")