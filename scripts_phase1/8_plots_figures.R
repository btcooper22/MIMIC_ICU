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
