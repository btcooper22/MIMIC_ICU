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
ggsave("writeup/figures/calibration.png")

