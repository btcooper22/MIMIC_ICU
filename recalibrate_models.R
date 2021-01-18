# Packages and functions
require(dplyr)
require(readr)
require(forcats)
require(magrittr)
require(lubridate)
require(ROCR)
require(ggplot2)
require(tidyr)
require(ResourceSelection)
require(s2dverification)
require(scales)

source("functions/inverse_logit.R")
source("functions/calibration.R")
source("functions/brier_extraction.R")

# Process data----------

# Load data
patients <- read_csv("data/final_patients.csv")

# Fix variables
patients %<>% 
  # Readmission
  mutate(readmission = readmission == "Readmitted to ICU") %>% 
  # Admission source
  mutate(admission_source = fct_relevel(admission_source, "OT"))



# Create categorised versions of continous variables----------

# Factorise variables
patients %<>% 
  # Serum choride (98 - 107)
  mutate(serum_chloride = case_when(
    serum_choride > 107 ~ "High",
    serum_choride > 98 ~ "Normal",
    is.numeric(serum_choride) ~ "Low"
  ),
  # BUN (7 - 20)
  blood_urea_nitrogen = case_when(
    blood_urea_nitrogen > 20 ~ "High",
    blood_urea_nitrogen > 7 ~ "Normal",
    is.numeric(blood_urea_nitrogen) ~ "Low"
  ),
  # Serum glucose (< 140)
  serum_glucose = case_when(
    serum_glucose > 140 ~ "High",
    is.numeric(serum_glucose) ~ "Normal"
  ),
  # final_pulse (60-100)
  final_pulse = case_when(
    final_pulse > 100 ~ "High",
    final_pulse > 60 ~ "Normal",
    is.numeric(final_pulse) ~ "Low"
  ),
  # final_temp (>36.5)
  final_temp = case_when(
    final_temp > 36.49 ~ "Normal",
    is.numeric(final_temp) ~ "Low"
  ),
  # SpO2 (< 95)
  final_SpO2 = case_when(
    final_SpO2 < 95 ~ "Low",
    is.numeric(final_SpO2) ~ "Normal"
  ),
  # bp *(> 100)
  final_bp = case_when(
    final_bp > 100 ~ "High",
    is.numeric(final_bp) ~ "Normal"
  ),
  # platelets (150-450)
  final_platelets = case_when(
    final_platelets > 450 ~ "High",
    final_platelets > 150 ~ "Normal",
    is.numeric(final_platelets) ~ "Low"
  ),
  # Lactate (> 2)
  final_lactate = case_when(
    final_lactate > 2 ~ "High",
    is.numeric(final_lactate) ~ "Normal"
  )
  )

# Reorder factor levels
patients %<>% 
  mutate(serum_chloride = fct_relevel(serum_chloride, "Normal"),
         blood_urea_nitrogen = fct_relevel(blood_urea_nitrogen, "Normal"),
         serum_glucose = fct_relevel(serum_glucose, "Normal"),
         final_pulse = fct_relevel(final_pulse, "Normal"),
         final_temp = fct_relevel(final_temp, "Normal"),
         final_SpO2 = fct_relevel(final_SpO2, "Normal"),
         final_bp = fct_relevel(final_bp, "Normal"),
         final_platelets = fct_relevel(final_platelets, "Normal"),
         final_lactate = fct_relevel(final_lactate, "Normal"))
  

# Recalibrate models------------

# Split data
set.seed(123)
patients_train <- patients %>% 
  group_by(readmission) %>% 
  slice_sample(prop = 0.5) %>% 
  mutate(type = "train")

patients_validate <- patients %>% 
  filter(row_id %in% patients_train$row_id == FALSE)

# Hammer
hammer_recalibrate <- glm(readmission ~ sex + general_surgery + cardiac_surgery + hyperglycemia +
                            high_apache + fluid_balance_5L + ambulation + los_5,
                          data = patients_train, family = "binomial")

# Martin
martin_recalibrate <- glm(readmission ~ respiratory_rate + age + serum_chloride + blood_urea_nitrogen +
                            atrial_fibrillation + renal_insufficiency + serum_glucose,
                          data = patients_train, family = "binomial")

# Frost
frost_recalibrate <- glm(readmission ~ age + sex + elective_admission + admission_source +
                           apache_II_discharge + los_7 + after_hours_discharge + acute_renal_failure,
                         data = patients_train, family = "binomial")

# APACHE
apache_recalibrate <- glm(readmission ~ apache_II_discharge,
                          data = patients_train, family = "binomial")

# Fialho
fialho_recalibrate <- glm(readmission ~ final_pulse + final_temp + final_SpO2 +
                            final_bp + final_platelets + final_lactate,
                          data = patients_train, family = "binomial")

# Fit novel model-----------

# Build formula
formu <- paste(names(patients_train)[c(9:10, 12, 14, 16:18, 
                                       20:25, 27:35, 51, 74)],
               collapse = " + ")

# Build model
dirty_model <- glm(paste("readmission", "~", formu, sep = " "),
                   data = patients_train, family = "binomial")


if(file.exists("data/cooper_model.RDS"))
{
  cooper_model <- readRDS("data/cooper_model.RDS")
}else
{
  # Automated model selection
  cooper_model <- step(glm(readmission ~ 1, data = patients_train, family = "binomial"),
                       scope = paste("readmission", "~", formu, sep = " "))
  write_rds(cooper_model, "data/cooper_model.RDS")
}

# Compare coefficients and predict---------

### Predict on new data
hammer_rc_probs <- predict(hammer_recalibrate, newdata = patients_validate) %>% inverse_logit()
martin_rc_probs <- predict(martin_recalibrate, newdata = patients_validate) %>% inverse_logit()
frost_rc_probs <- predict(frost_recalibrate, newdata = patients_validate) %>% inverse_logit()
apache_rc_probs <- predict(apache_recalibrate, newdata = patients_validate) %>% inverse_logit()
cooper_rc_probs <- predict(cooper_model, newdata = patients_validate) %>% inverse_logit()
fialho_rc_probs <- predict(fialho_recalibrate, newdata = patients_validate) %>% inverse_logit()

# Assess discrimination----------

# Create prediction objects
prediction_rc_hammer <- prediction(hammer_rc_probs, patients_validate$readmission)
prediction_rc_martin <- prediction(martin_rc_probs, patients_validate$readmission)
prediction_rc_frost <- prediction(frost_rc_probs, patients_validate$readmission)
prediction_rc_apache <- prediction(apache_rc_probs, patients_validate$readmission)
prediction_rc_cooper <- prediction(cooper_rc_probs, patients_validate$readmission)

# Create performance objects
performance_rc_hammer <- performance(prediction_rc_hammer, "tpr", "fpr")
performance_rc_martin <- performance(prediction_rc_martin, "tpr", "fpr")
performance_rc_frost <- performance(prediction_rc_frost, "tpr", "fpr")
performance_rc_apache <- performance(prediction_rc_apache, "tpr", "fpr")
performance_rc_cooper <- performance(prediction_rc_cooper, "tpr", "fpr")

# Create AUC objects
auc_rc_hammer <- performance(prediction_rc_hammer, measure = "auc")
auc_rc_martin <- performance(prediction_rc_martin, measure = "auc")
auc_rc_frost <- performance(prediction_rc_frost, measure = "auc")
auc_rc_apache <- performance(prediction_rc_apache, measure = "auc")
auc_rc_cooper <- performance(prediction_rc_cooper, measure = "auc")

# Print AUC
auc_rc_hammer@y.values[[1]]
auc_rc_martin@y.values[[1]]
auc_rc_frost@y.values[[1]]
auc_rc_apache@y.values[[1]]
auc_rc_cooper@y.values[[1]]


# Plot AUC
data.frame(x = performance_rc_hammer@x.values[[1]],
           y = performance_rc_hammer@y.values[[1]],
           model = "Hammer") %>% 
  rbind(
    data.frame(x = performance_rc_martin@x.values[[1]],
               y = performance_rc_martin@y.values[[1]],
               model = "Martin"),
    data.frame(x = performance_rc_frost@x.values[[1]],
               y = performance_rc_frost@y.values[[1]],
               model = "Frost"),
    data.frame(x = performance_rc_apache@x.values[[1]],
               y = performance_rc_apache@y.values[[1]],
               model = "APACHE-II"),
    data.frame(x = performance_rc_cooper@x.values[[1]],
               y = performance_rc_cooper@y.values[[1]],
               model = "Cooper")
  ) %>% 
  ggplot(aes(x, y, colour = model))+
  geom_abline(slope = 1, intercept = 0,
              linetype = "dotted",
              size = 1)+
  geom_path(size = 1)+
  labs(x = "1 - Specificity",
       y = "Sensitivity")+
  theme_classic(20)+
  theme(legend.position = "top")+
  scale_color_brewer(palette = "Set1",
                     name = "")

# Assess calibration-----------

# Split data into deciles
deciles_rc <- tibble(
  patient_id = 1:nrow(patients_validate),
  readmission = patients_validate$readmission,
  hammer_rc_probs,
  decile_hammer = ntile(hammer_rc_probs, 10),
  martin_rc_probs,
  decile_martin = ntile(martin_rc_probs, 10),
  frost_rc_probs,
  decile_frost = ntile(frost_rc_probs, 10),
  apache_rc_probs,
  decile_apache = ntile(apache_rc_probs, 10),
  cooper_rc_probs,
  decile_cooper = ntile(cooper_rc_probs, 10),
)

# Calculate calibration
cal_hammer <- calibration(deciles_rc, "hammer_rc_probs", "decile_hammer") %>% 
  mutate(model = "hammer")
cal_martin <- calibration(deciles_rc, "martin_rc_probs", "decile_martin") %>% 
  mutate(model = "martin")
cal_frost <- calibration(deciles_rc, "frost_rc_probs", "decile_frost") %>% 
  mutate(model = "frost")
cal_apache <- calibration(deciles_rc, "apache_rc_probs", "decile_apache") %>% 
  mutate(model = "apache")
cal_cooper <- calibration(deciles_rc, "cooper_rc_probs", "decile_cooper") %>% 
  mutate(model = "cooper")

# Plot
rbind(cal_hammer %>% select(-decile_hammer), 
      cal_martin %>% select(-decile_martin),
      cal_frost %>% select(-decile_frost), 
      cal_apache %>% select(-decile_apache),
      cal_cooper %>% select(-decile_cooper)) %>% 
  ggplot(aes(predicted, observed ,
             colour = model))+
  geom_abline(slope = 1, intercept = 0,
              size = 1)+
  geom_path(size = 1)+
  geom_pointrange(aes(ymin = observed - error ,
                      ymax = observed + error,
                      y = observed,
                      x = predicted,
                      colour = model)) +
  theme_classic(20)+
  theme(legend.position = "top")+
  labs(x = "Predicted readmission",
       y = "Observed readmission")+
  scale_colour_brewer(palette = "Set1",
                      name = "")+
  facet_wrap(~model, scales = "fixed")


# Calculate hosmer-lemeshow chi-squared
hoslem_hammer <- hoslem.test(patients_validate$readmission,
                             hammer_rc_probs, g = 10)
hoslem_martin <- hoslem.test(patients_validate$readmission,
                             martin_rc_probs, g = 10)
hoslem_frost <- hoslem.test(patients_validate$readmission,
                            frost_rc_probs, g = 10)
hoslem_apache <- hoslem.test(patients_validate$readmission,
                             apache_rc_probs, g = 10)
hoslem_cooper <- hoslem.test(patients_validate$readmission,
                             cooper_rc_probs, g = 10)
hoslem_hammer
hoslem_martin
hoslem_frost
hoslem_apache
hoslem_cooper

# Calculate brier scores
brier_df <- rbind(
  brier_extraction(patients_validate$readmission, hammer_rc_probs) %>% 
    mutate(model = "hammer"),
  brier_extraction(patients_validate$readmission, martin_rc_probs) %>% 
    mutate(model = "martin"),  
  brier_extraction(patients_validate$readmission, frost_rc_probs) %>% 
    mutate(model = "frost"),
  brier_extraction(patients_validate$readmission, apache_rc_probs) %>% 
    mutate(model = "apache"),
  brier_extraction(patients_validate$readmission, cooper_rc_probs) %>% 
    mutate(model = "cooper")
)

# Plot brier scores
brier_df %>% 
  pivot_longer(1:3) %>% 
  ggplot(aes(x = model, y = value))+
  geom_bar(stat = "identity", colour = "black",
           aes(fill = model))+
  facet_wrap(~name, scales = "free_y")+
  theme_classic(20)+
  theme(legend.position = "none")+
  labs(x = "", y = "")+
  scale_fill_brewer(palette = "Set1",
                    name = "")

# Plot combined
data.frame(
  model = c("hammer", "martin", "frost", "apache", "cooper"),
  discrimination = c(auc_rc_hammer@y.values[[1]],
                     auc_rc_martin@y.values[[1]],
                     auc_rc_frost@y.values[[1]],
                     auc_rc_apache@y.values[[1]],
                     auc_rc_cooper@y.values[[1]]),
  calibration = c(hoslem_hammer$statistic,
                  hoslem_martin$statistic,
                  hoslem_frost$statistic,
                  hoslem_apache$statistic,
                  hoslem_cooper$statistic)
) %>% 
  ggplot(aes(x = discrimination,
             y = calibration))+
  geom_point(aes(fill = model),
             size = 6,
             shape = 21)+
  theme_classic(20)+
  scale_fill_brewer(palette = "Set1",
                    name = "")+
  theme(legend.position = "top")+
  scale_y_reverse()+
  labs(x = "Discrimination (AUROC)",
       y = expression(paste("Calibration (", chi^2, ")")))

