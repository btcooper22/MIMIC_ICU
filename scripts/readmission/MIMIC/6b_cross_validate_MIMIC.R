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
require(caret)
require(glmnet)
require(mice)
require(doParallel)
require(tools)
require(snakecase)

source("functions/NA_count.R")
source("functions/inverse_logit.R")
source("functions/nomogram_convert.R")
source("functions/calibration.R")

# Load and clean data-----
patients <- read_csv("data/predictors.csv")

# Count and exclude missing chart data
patients %<>% filter(chart_missing == FALSE)

# Count and exclude missing events data
patients %<>% filter(!is.na(patients$fluid_balance_5L))

# Filter patients missing respiratory rate assessments
patients$respiratory_rate_initial[patients$respiratory_rate_initial > 60] <- NA
patients %<>% filter(!is.na(respiratory_rate_initial))

# Filter patients missing functional mobility assessments
patients %<>% filter(!is.na(patients$ambulation))

# Filter patients with no initial or final blood labs
patients %<>% filter(!is.na(patients$hyperglycemia_initial) & !is.na(patients$anaemia_initial) &
                       !is.na(patients$hyperglycemia_final) & !is.na(patients$anaemia_final) &
                       !is.na(patients$blood_urea_nitrogen_initial) & !is.na(patients$blood_urea_nitrogen_final) &
                       !is.na(patients$serum_chloride_initial) & !is.na(patients$serum_chloride_final))

# Clear up outcome
patients %<>% 
  mutate(readmission = fct_recode(as.character(readmission),
                                  `No Readmission` = "FALSE",
                                  `Readmitted to ICU` = "TRUE"))

# Fix age
patients %<>% 
  mutate(age = ifelse(age == ">90", 90, age)) %>% 
  mutate(age = floor(as.numeric(age)))

# Fix admission source
patients %<>% 
  mutate(admission_source = fct_relevel(admission_source, "OT"))

# Create splits
splits_df <- data.frame(patients,
                        split = floor(seq(0,10,length.out = nrow(patients))) + 1)
splits_df$split[splits_df$split == 11] <- 10

# Frost----
frost_xval <- foreach(i = 1:10, .combine = "rbind") %do%
  {
    # Extract validation values
    validation_df <- splits_df %>% 
      filter(split == i)
    
    # Create scoring system from nomogram for age
    age_score_system_input <- data.frame(
      input = seq(15, 90, 5),
      pix = seq(2, 287, length.out = 16)
    )
    age_score_system_output <- data.frame(
      output = seq(0, 8, 1),
      pix = seq(2, 290, length.out = 9)
    )
    
    # Create scoring system from nomogram for APACHE-II
    apache_score_system_input <- data.frame(
      input = seq(0, 45, 5),
      pix = seq(12, 733, length.out = 10)
    )
    apache_score_system_output <- data.frame(
      output = seq(0, 20, 1),
      pix = seq(12, 733, length.out = 21)
    )
    
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
    
    # Score data - original
    scores_original <- nomogram_convert(validation_df$age, age_score_system_input,
                                     age_score_system_output) +
      ifelse(validation_df$sex == "M", 2, 0) +
      ifelse(validation_df$elective_admission == FALSE, 12.25, 0) +
      case_when(
        validation_df$admission_source == "OT" ~ 0,
        validation_df$admission_source == "Emergency" ~ 9.5,
        validation_df$admission_source == "OtherHosp" ~ 10.25,
        validation_df$admission_source == "Ward" ~ 14.75
      ) +
      nomogram_convert(validation_df$apache_II, apache_score_system_input,
                       apache_score_system_output) +
      ifelse(validation_df$los_7, 17, 0) +
      ifelse(validation_df$after_hours_discharge, 4, 0) +
      ifelse(validation_df$acute_renal_failure_total, 10, 0)
    
    # Score data - aligned
    scores_aligned <- nomogram_convert(validation_df$age, age_score_system_input,
                                     age_score_system_output) +
      ifelse(validation_df$sex == "M", 2, 0) +
      ifelse(validation_df$elective_admission == FALSE, 12.25, 0) +
      case_when(
        validation_df$admission_source == "OT" ~ 0,
        validation_df$admission_source == "Emergency" ~ 9.5,
        validation_df$admission_source == "OtherHosp" ~ 10.25,
        validation_df$admission_source == "Ward" ~ 14.75
      ) +
      nomogram_convert(validation_df$apache_II, apache_score_system_input,
                       apache_score_system_output) +
      ifelse(validation_df$los_7, 17, 0) +
      ifelse(validation_df$after_hours_discharge, 4, 0) +
      ifelse(validation_df$acute_renal_failure_initial, 10, 0)
    
    # Convert to probs
    probs_frost_O <- nomogram_convert(scores_original, points_system_input,
                                      points_system_output, log = TRUE)
    probs_frost_A <- nomogram_convert(scores_aligned, points_system_input,
                                      points_system_output, log = TRUE)
    
    # Calculate AUC
    pred_O <- prediction(probs_frost_O, validation_df$readmission == "Readmitted to ICU")
    pred_A <- prediction(probs_frost_A, validation_df$readmission == "Readmitted to ICU")
    
    AUC_O <- performance(pred_O, measure = "auc")@y.values[[1]]
    AUC_A <- performance(pred_A, measure = "auc")@y.values[[1]]
    
    # Calculate chi-squared
    hoslem_A <- hoslem.test(validation_df$readmission == "Readmitted to ICU",
                             probs_frost_A, g = 10)
    hoslem_O <- hoslem.test(validation_df$readmission == "Readmitted to ICU",
                            probs_frost_O, g = 10)
    
    # Output
    data.frame(model = "frost", AUC_O, AUC_A,
               chi_O = hoslem_O$statistic,
               chi_A = hoslem_A$statistic)
  }

# Martin----
martin_xval <- foreach(i = 1:10, .combine = "rbind") %do%
  {
    # Extract validation values
    validation_df <- splits_df %>% 
      filter(split == i)
    
    # Score data - Original
    coefficients_martin_O <- -9.284491 +
      (0.04883 * validation_df$respiratory_rate_final) +
      (0.011588 * validation_df$age) +
      (0.036104 * validation_df$serum_chloride_final) +
      (0.004967 * validation_df$blood_urea_nitrogen_final) +
      (0.580153 * validation_df$atrial_fibrillation) +
      (0.458202 * validation_df$renal_insufficiency) +
      (0.003519 * validation_df$serum_glucose_final)
    
    # Score data - Aligned
    coefficients_martin_A <- -9.284491 +
      (0.04883 * validation_df$respiratory_rate_initial) +
      (0.011588 * validation_df$age) +
      (0.036104 * validation_df$serum_chloride_initial) +
      (0.004967 * validation_df$blood_urea_nitrogen_initial) +
      (0.580153 * validation_df$atrial_fibrillation) +
      (0.458202 * validation_df$renal_insufficiency) +
      (0.003519 * validation_df$serum_glucose_initial)
    
    # Convert to probs
    probs_martin_O <- inverse_logit(coefficients_martin_O)
    probs_martin_A <- inverse_logit(coefficients_martin_A)
    
    # Calculate AUC
    pred_O <- prediction(probs_martin_O, validation_df$readmission == "Readmitted to ICU")
    pred_A <- prediction(probs_martin_A, validation_df$readmission == "Readmitted to ICU")
    
    AUC_O <- performance(pred_O, measure = "auc")@y.values[[1]]
    AUC_A <- performance(pred_A, measure = "auc")@y.values[[1]]
    
    # Calculate chi-squared
    hoslem_A <- hoslem.test(validation_df$readmission == "Readmitted to ICU",
                            probs_martin_A, g = 10)
    hoslem_O <- hoslem.test(validation_df$readmission == "Readmitted to ICU",
                            probs_martin_O, g = 10)
    
    # Output
    data.frame(model = "martin", AUC_O, AUC_A,
               chi_O = hoslem_O$statistic,
               chi_A = hoslem_A$statistic)
  }

# Hammer----
hammer_xval <- foreach(i = 1:10, .combine = "rbind") %do%
  {
    # Extract validation values
    validation_df <- splits_df %>% 
      filter(split == i)
    
    # Score data - Original
    coefficients_hammer_O <- log(0.005 / (1 - 0.005)) +
      ifelse(validation_df$sex == "M", 0.43, 0) +
      ifelse(validation_df$general_surgery, 0.64, 0) +
      ifelse(validation_df$cardiac_surgery, -1.01, 0) +
      ifelse(validation_df$hyperglycemia_final, 0.36, 0) +
      ifelse(validation_df$anaemia_final, 1.11, 0) +
      ifelse(validation_df$high_apache, 0.44, 0) +
      ifelse(validation_df$fluid_balance_5L, 0.52, 0) +
      ifelse(validation_df$ambulation == FALSE, 0.7, 0) +
      ifelse(validation_df$los_5, 0.88, 0)
    
    # Score data - Aligned
    coefficients_hammer_A <- log(0.005 / (1 - 0.005)) +
      ifelse(validation_df$sex == "M", 0.43, 0) +
      ifelse(validation_df$general_surgery, 0.64, 0) +
      ifelse(validation_df$cardiac_surgery, -1.01, 0) +
      ifelse(validation_df$hyperglycemia_initial, 0.36, 0) +
      ifelse(validation_df$anaemia_initial, 1.11, 0) +
      ifelse(validation_df$high_apache, 0.44, 0) +
      ifelse(validation_df$fluid_balance_5L, 0.52, 0) +
      ifelse(validation_df$ambulation == FALSE, 0.7, 0) +
      ifelse(validation_df$los_5, 0.88, 0)
    
    # Convert to probs
    probs_hammer_O <- inverse_logit(coefficients_hammer_O)
    probs_hammer_A <- inverse_logit(coefficients_hammer_A)
    
    # Calculate AUC
    pred_O <- prediction(probs_hammer_O, validation_df$readmission == "Readmitted to ICU")
    pred_A <- prediction(probs_hammer_A, validation_df$readmission == "Readmitted to ICU")
    
    AUC_O <- performance(pred_O, measure = "auc")@y.values[[1]]
    AUC_A <- performance(pred_A, measure = "auc")@y.values[[1]]
    
    # Calculate chi-squared
    hoslem_A <- hoslem.test(validation_df$readmission == "Readmitted to ICU",
                            probs_hammer_A, g = 10)
    hoslem_O <- hoslem.test(validation_df$readmission == "Readmitted to ICU",
                            probs_hammer_O, g = 10)
    
    # Output
    data.frame(model = "hammer", AUC_O, AUC_A,
               chi_O = hoslem_O$statistic,
               chi_A = hoslem_A$statistic)
  }

# Output and plot-----

# Combine and tidy
xval_df <- rbind(frost_xval, martin_xval,
                 hammer_xval) %>% 
  pivot_longer(2:5) %>% 
  mutate(type = ifelse(grepl("AUC", name), "Discrimination", "Calibration"),
         dataset = ifelse(grepl("O", name), "MIMIC-O", "MIMIC-A"),
         model = to_any_case(model, "title"))
  
# Output means and SD
xval_summarised <- xval_df %>% 
  group_by(model, type, dataset) %>% 
  summarise(mean = mean(value, na.rm = T),
            sd = sd(value, na.rm = T))
xval_summarised

# Plot
xval_summarised %>% 
  ggplot(aes(model, mean, fill = dataset)) + 
  geom_bar(stat = "identity", color = "black", 
           position=position_dodge(), alpha = 0.8)+
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
                width = 0.2, position = position_dodge(0.9))+
  facet_wrap(~type, nrow = 1,
             scales = "free")+
  scale_fill_manual(name = "Dataset",
                    values = c("#984ea3", "#ff7f00"))+
  theme_classic(20)+
  theme(legend.position = "top")+
  labs(x = "", y = "")

ggsave("figures/cross_validation.png",
       height = 8, width = 15)  
