# Packages and functions----
require(readr)
require(dplyr)
require(tidyr)

source("functions/value_screen.R")

# Load and check implausible values----
predictors <- read_csv("data/predictors_mortality.csv")

# Temperature
hist(value_screen(predictors, "temperature",
                  c(30,40), TRUE)$value)

value_screen(predictors, "temperature",
             c(30,40))

# Systolic blood pressure
hist(value_screen(predictors, "systolicbp",
                  c(70,220), TRUE)$value)

value_screen(predictors, "systolicbp",
             c(70,220))


# Diastolic blood pressure
hist(value_screen(predictors, "diastolicbp",
                  c(40,140), TRUE)$value)

value_screen(predictors, "diastolicbp",
             c(40,140))

# Pulse
hist(value_screen(predictors, "pulse",
                  c(20,200), TRUE)$value)

value_screen(predictors, "pulse",
             c(20,200))

# Respiratory rate
hist(value_screen(predictors, "respiratory",
                  c(0,60), TRUE)$value)

value_screen(predictors, "respiratory",
             c(0,60))

# Arterial pH
hist(value_screen(predictors, "arterialpH",
                  c(7,8), TRUE)$value)

value_screen(predictors, "arterialpH",
             c(7,8))

# Sodium
hist(value_screen(predictors, "sodium",
                  c(100,160), TRUE)$value)

value_screen(predictors, "sodium",
             c(100, 160))

# Potassium
hist(value_screen(predictors, "potassium",
                  c(0,8.5), TRUE)$value)

value_screen(predictors, "potassium",
             c(0, 8.5))

# Creatinine
hist(value_screen(predictors, "creatinine",
                  c(0,6), TRUE)$value)

value_screen(predictors, "creatinine",
             c(0,6))

# Haematocrit
hist(value_screen(predictors, "haematocrit",
                  c(15,60), TRUE)$value)

value_screen(predictors, "haematocrit",
             c(15,60))

# WBC
hist(value_screen(predictors, "whitebloodcount",
                  c(0,50), TRUE)$value)

value_screen(predictors, "whitebloodcount",
             c(0,50))

# Bicarbonate
hist(value_screen(predictors, "bicarbonate",
                  c(0,60), TRUE)$value)

value_screen(predictors, "bicarbonate",
             c(0,60))

# FiO2
hist(value_screen(predictors, "fractioninspiredoxygen",
                  c(21,100), TRUE)$value)

value_screen(predictors, "fractioninspiredoxygen",
             c(21,100))

# PaO2
hist(value_screen(predictors, "arterialoxygen",
                  c(20,250), TRUE)$value)

value_screen(predictors, "arterialoxygen",
             c(20,250))

# PaCO2
hist(value_screen(predictors, "arterialcarbon",
                  c(8,200), TRUE)$value)

value_screen(predictors, "arterialcarbon",
             c(8,200))

# Calculate apache scores


# Pivot and recombine data----

# Seperate admission and discharge apache scores
results_admission <- predictors %>% 
  select(1:29) %>% 
  mutate(type = "admission")

results_discharge <- predictors %>% 
  select(c(1:9, 30:49)) %>% 
  mutate(type = "discharge")

# Rename variables
names(results_admission)[10:29] <- c("temperature", "systolicbp", "diastolicbp",
                                     "pulse", "respiratory", "arterialpH",
                                     "sodium", "potassium", "creatinine",
                                     "haematocrit", "whitebloodcount",
                                     "glasgowcomaeye", "glasgowcomamotor",
                                     "glasgowcomaverbal", "age", "chronic",
                                     "bicarbonate", "fractioninspiredoxygen",
                                     "arterialoxygen", "arterialcarbon")

names(results_discharge)[10:29] <- c("temperature", "systolicbp", "diastolicbp",
                                     "pulse", "respiratory", "arterialpH",
                                     "sodium", "potassium", "creatinine",
                                     "haematocrit", "whitebloodcount",
                                     "glasgowcomaeye", "glasgowcomamotor",
                                     "glasgowcomaverbal", "age", "chronic",
                                     "bicarbonate", "fractioninspiredoxygen",
                                     "arterialoxygen", "arterialcarbon")

# Recombine
results <- rbind(results_admission, results_discharge)

# Calculate apache score----

# Temperature
temperature_scores <- case_when(
  results$temperature >= 41 | results$temperature <= 29.9 ~ 4,
  results$temperature >= 39 | results$temperature <= 31.9 ~ 3,
  results$temperature <= 33.9 ~ 2,
  results$temperature >= 38.5 | results$temperature <= 35.9 ~ 1,
  is.na(results$temperature) ~ NaN,
  is.numeric(results$temperature) ~ 0,
)

# Mean arterial pressure
map_values <- (results$systolicbp + (2 * results$diastolicbp) )/ 3
map_scores <- case_when(
  map_values >= 160 | map_values <= 49 ~ 4,
  map_values >= 130 ~ 3,
  map_values >= 110 | map_values <= 69 ~ 2,
  is.na(map_values) ~ NaN,
  is.numeric(map_values) ~ 0
)

# Heart rate
pulse_score <- case_when(
  results$pulse >= 180 | results$pulse <= 39 ~ 4,
  results$pulse >= 140 | results$pulse <= 54 ~ 3,
  results$pulse >= 110 | results$pulse <= 69 ~ 2,
  is.na(map_values) ~ NaN,
  is.numeric(results$pulse) ~ 0
)

# Respiratory rate
respiratory_score <-  case_when(
  results$respiratory >= 50 | results$respiratory <= 5 ~ 4,
  results$respiratory >= 35 ~ 3,
  results$respiratory <= 9 ~ 2,
  results$respiratory >= 25 | results$respiratory <= 11 ~ 1,
  is.na(results$respiratory) ~ NaN,
  is.numeric(results$respiratory) ~ 0
)

# Arterial pH
artpH_score <-  case_when(
  results$arterialpH >= 7.7 | results$arterialpH < 7.15 ~ 4,
  results$arterialpH >= 7.6 | results$arterialpH <= 7.24 ~ 3,
  results$arterialpH <= 7.32 ~ 2,
  results$arterialpH >= 7.5 ~ 1,
  is.na(results$arterialpH) ~ NaN,
  is.numeric(results$arterialpH) ~ 0
)

# Sodium
sodium_score <-  case_when(
  results$sodium >= 180 | results$sodium <= 110 ~ 4,
  results$sodium >= 160 | results$sodium <= 119 ~ 3,
  results$sodium >= 155 | results$sodium <= 129 ~ 2,
  results$sodium >= 150 ~ 1,
  is.na(results$sodium) ~ NaN,
  is.numeric(results$sodium) ~ 0
)

# Potassium
potassium_score <- case_when(
  results$potassium >= 7 | results$potassium < 2.5 ~ 4,
  results$potassium >= 6 ~ 3,
  results$potassium <= 2.9 ~ 2,
  results$potassium >= 5.5 | results$potassium <= 3.4 ~ 1,
  is.na(results$potassium) ~ NaN,
  is.numeric(results$potassium) ~ 0
)

# Creatinine
creatinine_score <- case_when(
  results$creatinine >= 3.5 ~ 4,
  results$creatinine >= 2 ~ 3,
  results$creatinine >= 1.5 | results$creatinine < 0.6 ~ 2,
  is.na(results$creatinine) ~ NaN,
  is.numeric(results$creatinine) ~ 0
)

creatinine_score <- ifelse(results$acute_renal_failure,
       creatinine_score * 2, 
       creatinine_score)

# Haematocrit
haematocrit_score <- case_when(
  results$haematocrit >= 60 | results$haematocrit < 20 ~ 4,
  results$haematocrit >= 50 | results$haematocrit <= 29.9 ~ 2,
  results$haematocrit >= 46 ~ 1,
  is.na(results$haematocrit) ~ NaN,
  is.numeric(results$haematocrit) ~ 0
)

# White blood cell count
wbc_score <- case_when(
  results$whitebloodcount >= 40 | results$whitebloodcount < 1 ~ 4,
  results$whitebloodcount >= 20 | results$whitebloodcount <= 2.9 ~ 2,
  results$whitebloodcount >= 15 ~ 1,
  is.na(results$whitebloodcount) ~ NaN,
  is.numeric(results$whitebloodcount) ~ 0
)

# PaO2
pao2_score <- case_when(
  results$arterialoxygen < 55 ~ 4,
  results$arterialoxygen <= 60 ~ 3,
  results$arterialoxygen <= 70 ~ 1,
  is.na(results$arterialoxygen) ~ NaN,
  is.numeric(results$arterialoxygen) ~ 0
)

# Extract FiO2
fiO2 <- ifelse(is.na(results$fractioninspiredoxygen), 21, 
               results$fractioninspiredoxygen)

# Alveolar-arterial gradient
elev <-  5
patm <- 760 * exp(elev / -7000)
ph2o <- 47 * exp((results$temperature - 37) / 18.4)
fio2_remap <- purrr::map_dbl(fiO2, ~ dplyr::if_else(.x > 1, .x / 100, .x))
aa_grad <- (fio2_remap * (patm - ph2o) - (results$arterialcarbon / 0.8)) - results$arterialoxygen

aagrad_score <- case_when(
  aa_grad >= 500 ~ 4,
  aa_grad >= 350 ~ 3,
  aa_grad >= 200 ~ 2,
  is.na(aa_grad) ~ NaN,
  is.numeric(aa_grad) ~ 0
)

# If FiO2 > 0.5, use aa_grad for oxygenation score, else use PaO2
oxygenation_score <- ifelse(fiO2 > 50, aagrad_score, pao2_score)

# Bicarbonate
bicarbonate_score <- case_when(
  results$bicarbonate >= 52 | results$bicarbonate < 15 ~ 4,
  results$bicarbonate >= 41 | results$bicarbonate <= 17.9 ~ 3,
  results$bicarbonate <= 21.9 ~ 2,
  results$bicarbonate >= 32 ~ 1,
  is.na(results$bicarbonate) ~ NaN,
  is.numeric(results$bicarbonate) ~ 0
)

# Glasgow coma scale
gcs_score <- 15 - (results$glasgowcomaeye + results$glasgowcomamotor +
        results$glasgowcomaverbal)

# Age
age_score <- results$age

# Chronic health
chronic_score <- results$chronic

# If arterial blood gasses missing, use bicarbonate instead
abg_score <- 
ifelse(is.na(oxygenation_score) & is.na(artpH_score),
       bicarbonate_score, oxygenation_score + artpH_score)

# Sum scores
apache_scores <- temperature_scores + map_scores + pulse_score +
  respiratory_score + abg_score + sodium_score + potassium_score +
  creatinine_score + wbc_score + gcs_score + age_score + chronic_score

# Attach and write
results %>% 
  mutate(apache_II = apache_scores,
         missing_abg = is.na(oxygenation_score) & is.na(artpH_score)) %>% 
  filter(lab_missing == FALSE) %>% 
  select(-chart_missing, -lab_missing) %>% 
  relocate(row_id, subject_id,
           adm_id, type, apache_II) %>% 
  write_csv("data/apache_data_full.csv")
