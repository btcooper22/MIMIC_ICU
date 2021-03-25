# Packages
require(readr)
require(ROCR)
require(foreach)
require(dplyr)
require(tidyr)
require(ggplot2)
require(patchwork)

# Load models
ICNARC_bespoke <- read_rds("scripts/readmission/shared/models/LTH_ICNARC_bespoke.RDS")
ICNARC_hammer <- read_rds("scripts/readmission/shared/models/LTH_ICNARC_hammer.RDS")
ICNARC_martin <- read_rds("scripts/readmission/shared/models/LTH_ICNARC_martin.RDS")
ICNARC_frost <- read_rds("scripts/readmission/shared/models/LTH_ICNARC_frost.RDS")

MIMIC_bespoke <- read_rds("scripts/readmission/shared/models/MIMIC_bespoke.RDS")
MIMIC_hammer <- read_rds("scripts/readmission/shared/models/MIMIC_hammer.RDS")
MIMIC_martin <- read_rds("scripts/readmission/shared/models/MIMIC_martin.RDS")
MIMIC_frost <- read_rds("scripts/readmission/shared/models/MIMIC_frost.RDS")

model_list <- list(ICNARC_bespoke, ICNARC_hammer, ICNARC_martin, ICNARC_frost,
                   MIMIC_bespoke, MIMIC_hammer, MIMIC_martin, MIMIC_frost)

# Table of values----
performance_df <- foreach(i = 1:length(model_list), .combine = "rbind") %do%
  {
    # Extract prediction object
    pred <- model_list[[i]][3]$discrimination
    
    # Extract AUC
    AUC <- performance(pred, "auc")@y.values[[1]]
    
    # Extract calibration test
    hoslem <- model_list[[i]][4]$calibration
    
    # Output
    data.frame(model = model_list[[i]][1]$model,
               dataset = model_list[[i]][2]$data,
               AUC, chisq = hoslem$statistic %>% unname(),
               p = hoslem$p.value > 0.05)
  }
performance_df %>% 
  pivot_wider(names_from = "dataset",
              values_from = c("AUC", "chisq", "p")) 

ggplot(performance_df,
       aes(AUC, chisq))+
  geom_point(size = 8,
             aes(fill = paste(model, dataset)),
             shape = 21)+
  theme_classic(20)+
  theme(legend.position = "top")+
  scale_fill_manual(values = c("#a6cee3", "#1f78b4",
                               "#b2df8a", "#33a02c",
                               "#fb9a99", "#e31a1c",
                               "#cab2d6", "#6a3d9a"),
                    name = "")+
  scale_y_reverse()

# Plot discrimination----
discrimination_df <- foreach(i = 1:length(model_list), .combine = "rbind") %do%
  {
    # Extract prediction object
    pred <- model_list[[i]][3]$discrimination
    
    # Calculate performance
    perf <- performance(pred, "tpr", "fpr")
    
    # Output
    data.frame(x = perf@x.values[[1]],
               y = perf@y.values[[1]],
               model = model_list[[i]][1]$model,
               dataset = model_list[[i]][2]$data)
  } %>% 
  mutate(identifier = paste(dataset, model, sep = "-"))

# Plot combined
discrimination_df %>% 
  ggplot(aes(x, y, colour = model))+
  geom_abline(slope = 1, intercept = 0,
              linetype = "dotted",
              size = 1)+
  geom_path(size = 1)+
  labs(x = "1 - Specificity",
       y = "Sensitivity")+
  theme_classic(20)+
  theme(legend.position = "top")+
  facet_wrap(~dataset)

# Bespoke
dplot_bespoke <- discrimination_df %>% 
  filter(model == "bespoke") %>% 
  ggplot(aes(x, y, colour = dataset))+
  geom_abline(slope = 1, intercept = 0,
              linetype = "dotted",
              size = 1)+
  geom_path(size = 2)+
  labs(x = "1 - Specificity",
       y = "Sensitivity")+
  theme_classic(20)+
  theme(legend.position = "top")+
  scale_colour_manual(values = c("#a6cee3", "#1f78b4"),
                      name = "Bespoke")

# Hammer
dplot_hammer <- discrimination_df %>% 
  filter(model == "hammer") %>% 
  ggplot(aes(x, y, colour = dataset))+
  geom_abline(slope = 1, intercept = 0,
              linetype = "dotted",
              size = 1)+
  geom_path(size = 2)+
  labs(x = "1 - Specificity",
       y = "Sensitivity")+
  theme_classic(20)+
  theme(legend.position = "top")+
  scale_colour_manual(values = c("#b2df8a", "#33a02c"),
                      name = "Hammer")

# Martin
dplot_martin <- discrimination_df %>% 
  filter(model == "martin") %>% 
  ggplot(aes(x, y, colour = dataset))+
  geom_abline(slope = 1, intercept = 0,
              linetype = "dotted",
              size = 1)+
  geom_path(size = 2)+
  labs(x = "1 - Specificity",
       y = "Sensitivity")+
  theme_classic(20)+
  theme(legend.position = "top")+
  scale_colour_manual(values = c("#fb9a99", "#e31a1c"),
                      name = "Martin")

# Frost
dplot_frost <- discrimination_df %>% 
  filter(model == "frost") %>% 
  ggplot(aes(x, y, colour = dataset))+
  geom_abline(slope = 1, intercept = 0,
              linetype = "dotted",
              size = 1)+
  geom_path(size = 2)+
  labs(x = "1 - Specificity",
       y = "Sensitivity")+
  theme_classic(20)+
  theme(legend.position = "top")+
  scale_colour_manual(values = c("#cab2d6", "#6a3d9a"),
                      name = "Frost")

# Write
dplot_colour <- (dplot_bespoke + dplot_frost) / (dplot_hammer + dplot_martin)
ggsave("figures/Discrimination_colour.png", dplot_colour,
       height = 15, width = 15)

# B/W Version
bw_plot <- discrimination_df %>% 
  ggplot(aes(x, y, colour = identifier))+
  geom_abline(slope = 1, intercept = 0,
              linetype = "dotted",
              size = 1)+
  geom_path(size = 2)+
  labs(x = "1 - Specificity",
       y = "Sensitivity")+
  theme_classic(20)+
  theme(legend.position = "top")+
  scale_colour_manual(values = rep("#f0f0f0",8),
                      name = "BW")

# Write
dplot_bw <- (bw_plot + bw_plot) / (bw_plot + bw_plot)
ggsave("figures/Discrimination_BW.png", dplot_bw,
       height = 15, width = 15)

# Plot calibration----
calibration_df <- foreach(i = 1:length(model_list), .combine = "rbind") %do%
  {
    # Extract prediction object
    deciles <- model_list[[i]][5]$deciles
    
    
    # Output
    data.frame(observed = deciles$observed,
               predicted = deciles$predicted,
               error = deciles$error,
               model = model_list[[i]][1]$model,
               dataset = model_list[[i]][2]$data)
  } %>% 
  mutate(identifier = paste(dataset, model, sep = "-"))

# Plot combined
calibration_df %>% 
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
  facet_wrap(~dataset, scales = "free")+
  scale_colour_brewer(palette = "Set1")

# Bespoke
cplot_bespoke <- calibration_df %>% 
  filter(model == "bespoke") %>% 
  ggplot(aes(predicted, observed ,
             colour = dataset))+
  geom_abline(slope = 1, intercept = 0,
              size = 1)+
  geom_path(size = 2)+
  geom_pointrange(aes(ymin = observed - error ,
                      ymax = observed + error,
                      y = observed,
                      x = predicted,
                      colour = dataset), size = 1) +
  theme_classic(20)+
  theme(legend.position = "top")+
  labs(x = "Predicted readmission",
       y = "Observed readmission")+
  scale_colour_manual(values = c("#a6cee3", "#1f78b4"),
                      name = "Bespoke")

# Hammer
cplot_hammer <- calibration_df %>% 
  filter(model == "hammer") %>% 
  ggplot(aes(predicted, observed ,
             colour = dataset))+
  geom_abline(slope = 1, intercept = 0,
              size = 1)+
  geom_path(size = 2)+
  geom_pointrange(aes(ymin = observed - error ,
                      ymax = observed + error,
                      y = observed,
                      x = predicted,
                      colour = dataset), size = 1) +
  theme_classic(20)+
  theme(legend.position = "top")+
  labs(x = "Predicted readmission",
       y = "Observed readmission")+
  scale_colour_manual(values = c("#b2df8a", "#33a02c"),
                      name = "Hammer")

# Martin
cplot_martin <- calibration_df %>% 
  filter(model == "martin") %>% 
  ggplot(aes(predicted, observed ,
             colour = dataset))+
  geom_abline(slope = 1, intercept = 0,
              size = 1)+
  geom_path(size = 2)+
  geom_pointrange(aes(ymin = observed - error ,
                      ymax = observed + error,
                      y = observed,
                      x = predicted,
                      colour = dataset), size = 1) +
  theme_classic(20)+
  theme(legend.position = "top")+
  labs(x = "Predicted readmission",
       y = "Observed readmission")+
  scale_colour_manual(values = c("#fb9a99", "#e31a1c"),
                      name = "Martin")

# Frost
cplot_frost <- calibration_df %>% 
  filter(model == "frost") %>% 
  ggplot(aes(predicted, observed ,
             colour = dataset))+
  geom_abline(slope = 1, intercept = 0,
              size = 1)+
  geom_path(size = 2)+
  geom_pointrange(aes(ymin = observed - error ,
                      ymax = observed + error,
                      y = observed,
                      x = predicted,
                      colour = dataset), size = 1) +
  theme_classic(20)+
  theme(legend.position = "top")+
  labs(x = "Predicted readmission",
       y = "Observed readmission")+
  scale_colour_manual(values = c("#cab2d6", "#6a3d9a"),
                      name = "Frost")

# Write
cplot_colour <- (cplot_bespoke + cplot_frost) / (cplot_hammer + cplot_martin)
ggsave("figures/Calibration_colour.png", cplot_colour,
       height = 15, width = 15)
