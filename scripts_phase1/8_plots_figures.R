# Source line reader
source("functions/run_lines.R")

# Prepare data
run_lines("scripts_phase1/6_compare_scores.R", 1:124)

# Cross-tabulate model data---------

# N
crosstab("readmission")

# Hammer
crosstab("sex")
crosstab("general_surgery")
crosstab("cardiac_surgery")
crosstab("hyperglycemia")
crosstab("anaemia")
crosstab("high_apache")
crosstab("fluid_balance_5L")
crosstab("ambulation")
crosstab("los_5")

# Frost
patients %>%       
  group_by(readmission) %>% 
  summarise(mean = mean(age),
            sd = sd(age))
crosstab("elective_admission")
crosstab("admission_source")
patients %>%       
  group_by(readmission) %>% 
  summarise(mean = mean(apache_II_discharge),
            sd = sd(apache_II_discharge))
crosstab("los_7")
crosstab("after_hours_discharge")  
crosstab("acute_renal_failure")

# Martin
patients %>%       
  group_by(readmission) %>% 
  summarise(mean = mean(respiratory_rate),
            sd = sd(respiratory_rate))
patients %>%       
  group_by(readmission) %>% 
  summarise(mean = mean(blood_urea_nitrogen),
            sd = sd(blood_urea_nitrogen))
patients %>%       
  group_by(readmission) %>% 
  summarise(mean = mean(serum_glucose),
            sd = sd(serum_glucose))
patients %>%       
  group_by(readmission) %>% 
  summarise(mean = mean(serum_choride),
            sd = sd(serum_choride))
crosstab("atrial_fibrillation")
crosstab("renal_insufficiency")

# Fialho
patients %>%       
  group_by(readmission) %>% 
  summarise(mean = mean(final_pulse),
            sd = sd(final_pulse))
patients %>%       
  group_by(readmission) %>% 
  summarise(mean = mean(final_temp),
            sd = sd(final_temp))
patients %>%       
  group_by(readmission) %>% 
  summarise(mean = mean(final_SpO2),
            sd = sd(final_SpO2))
patients %>%       
  group_by(readmission) %>% 
  summarise(mean = mean(final_bp),
            sd = sd(final_bp))
patients %>%       
  group_by(readmission) %>% 
  summarise(mean = mean(final_platelets),
            sd = sd(final_platelets))
patients %>%       
  group_by(readmission) %>% 
  summarise(mean = mean(final_lactate),
            sd = sd(final_lactate))

# Discrimination plots--------

# Generate standard predictions
run_lines("scripts_phase1/6_compare_scores.R", 202:478)
discrimination_standard <- data.frame(x = performance_hammer@x.values[[1]],
                                      y = performance_hammer@y.values[[1]],
                                      model = "Hammer") %>% 
  rbind(
    data.frame(x = performance_martin@x.values[[1]],
               y = performance_martin@y.values[[1]],
               model = "Martin"),
    data.frame(x = performance_frost@x.values[[1]],
               y = performance_frost@y.values[[1]],
               model = "Frost"),
    data.frame(x = performance_apache@x.values[[1]],
               y = performance_apache@y.values[[1]],
               model = "APACHE-II"),
    data.frame(x = performance_fialho@x.values[[1]],
               y = performance_fialho@y.values[[1]],
               model = "Fialho")
  ) %>% 
  mutate(type = "Standard")

# Generate recalibrated predictions
run_lines("scripts_phase1/7_recalibrate_models.R", 1:184)
discrimination_recalibrate <- data.frame(x = performance_rc_hammer@x.values[[1]],
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
               model = "Cooper"),
    data.frame(x = performance_rc_fialho@x.values[[1]],
               y = performance_rc_fialho@y.values[[1]],
               model = "Fialho")
  ) %>% 
  mutate(type = "Recalibrated")

# Combine and plot
discrimination_standard %>% 
  rbind(discrimination_recalibrate) %>% 
  mutate(type = fct_relevel(type, "Standard")) %>% 
  ggplot(aes(x, y, colour = model))+
  geom_abline(slope = 1, intercept = 0,
              linetype = "dotted",
              size = 1)+
  geom_path(size = 1)+
  labs(x = "1 - Specificity",
       y = "Sensitivity")+
  theme_classic(20)+
  theme(legend.position = "top")+
  facet_wrap(~type, nrow = 1)+
  scale_colour_manual(values = c("#1f78b4", "#b15928", "#e31a1c",
                                 "#ff7f00", "#6a3d9a", "#33a02c"),
                      name = "")

# Save
ggsave("writeup/figures/discrimination.png", 
       width = 14, height = 8)

# Calibration plots-------

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

# Generate standard calibration
run_lines("scripts_phase1/6_compare_scores.R", 1:539)
calibration_standard <- rbind(cal_hammer %>% select(-decile_hammer), 
                              cal_martin %>% select(-decile_martin),
                              cal_frost %>% select(-decile_frost), 
                              cal_apache %>% select(-decile_apache),
                              cal_fialho %>% select(-decile_fialho)) %>% 
  mutate(type = "Standard")

# Generate recalibration
run_lines("scripts_phase1/7_recalibrate_models.R", 1:250)
calibration_recalibrated <- rbind(cal_hammer %>% select(-decile_hammer), 
                              cal_martin %>% select(-decile_martin),
                              cal_frost %>% select(-decile_frost), 
                              cal_apache %>% select(-decile_apache),
                              cal_fialho %>% select(-decile_fialho),
                              cal_cooper %>% select(-decile_cooper)) %>% 
  mutate(type = "Recalibrated")

# Combine and clean
calibration_plot <- calibration_standard %>% 
  rbind(calibration_recalibrated) %>% 
  # Convert model names to upper case
  mutate(model = firstup(model)) %>% 
  mutate(model = ifelse(model == "Apache", "APACHE-II", model)) %>% 
  # Create combined model-type variable
  mutate(model_type = paste(model, type, sep = ": ")) %>% 
  # Reorder variable
  mutate(model_type = ifelse(model_type == "Cooper: Recalibrated",
                "Cooper", model_type)) %>% 
  mutate(model_type = fct_relevel(model_type, "APACHE-II: Standard",
                                  "APACHE-II: Recalibrated",
                                  "Cooper",
                                  "Fialho: Standard",
                                  "Fialho: Recalibrated",
                                  "Frost: Standard",
                                  "Frost: Recalibrated",
                                  "Hammer: Standard",
                                  "Hammer: Recalibrated",
                                  "Martin: Standard",
                                  "Martin: Recalibrated",
                                  ))
  
# Plot
calibration_plot %>% 
  ggplot(aes(predicted, observed ,
             colour = model_type))+
  geom_abline(slope = 1, intercept = 0,
              size = 1)+
  geom_path(size = 1)+
  geom_pointrange(aes(ymin = observed - error ,
                      ymax = observed + error,
                      y = observed,
                      x = predicted,
                      colour = model_type)) +
  theme_classic(20)+
  theme(legend.position = "top")+
  labs(x = "Predicted readmission",
       y = "Observed readmission")+
  facet_wrap(~model, scales = "free")+
  scale_colour_manual(values = c("#1f78b4", "#a6cee3", "#b15928",
                                 "#e31a1c", "#fb9a99",
                                 "#ff7f00", "#fdbf6f",
                                 "#6a3d9a", "#cab2d6",
                                 "#33a02c", "#b2df8a"),
                      name = "")

# Save
ggsave("writeup/figures/calibration.png", 
       width = 14, height = 8)

# Model performance table---------
