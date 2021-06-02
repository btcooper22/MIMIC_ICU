# Packages
require(readr)
require(ROCR)
require(foreach)
require(dplyr)
require(tidyr)
require(ggplot2)
require(cowplot)
require(forcats)

# Load models
ICNARC_bespoke <- read_rds("scripts/readmission/shared/models/LTH_ICNARC_bespoke.RDS")
ICNARC_hammer <- read_rds("scripts/readmission/shared/models/LTH_ICNARC_hammer.RDS")
ICNARC_martin <- read_rds("scripts/readmission/shared/models/LTH_ICNARC_martin.RDS")
ICNARC_frost <- read_rds("scripts/readmission/shared/models/LTH_ICNARC_frost.RDS")

MIMIC_bespoke <- read_rds("scripts/readmission/shared/models/MIMIC_bespoke.RDS")
MIMIC_hammer_A <- read_rds("scripts/readmission/shared/models/MIMIC_hammer_A.RDS")
MIMIC_martin_A <- read_rds("scripts/readmission/shared/models/MIMIC_martin_A.RDS")
MIMIC_frost_A <- read_rds("scripts/readmission/shared/models/MIMIC_frost_A.RDS")

MIMIC_hammer_O <- read_rds("scripts/readmission/shared/models/MIMIC_hammer_O.RDS")
MIMIC_martin_O <- read_rds("scripts/readmission/shared/models/MIMIC_martin_O.RDS")
MIMIC_frost_O <- read_rds("scripts/readmission/shared/models/MIMIC_frost_O.RDS")


model_list <- list(ICNARC_bespoke, ICNARC_hammer, ICNARC_martin, ICNARC_frost,
                   MIMIC_bespoke, MIMIC_hammer_A, MIMIC_martin_A, MIMIC_frost_A,
                   MIMIC_hammer_O, MIMIC_martin_O, MIMIC_frost_O)

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
               p = hoslem$p.value %>% round(3))
  } %>% 
  mutate(dataset = ifelse(dataset == "MIMIC","MIMIC_O",
                          dataset))
performance_df %>% 
  pivot_wider(names_from = "dataset",
              values_from = c("AUC", "chisq", "p")) 

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
  mutate(identifier = paste(dataset, model, sep = "-"),
         dataset = ifelse(dataset == "LTH_ICNARC", "ICNARC", dataset))


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
  scale_colour_manual(values = c("#bdc9e1", "#0570b0"),
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
  scale_colour_manual(values = c("#c2e699", "#78c679",
                                 "#238443"),
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
  scale_colour_manual(values = c("#fecc5c","#fb9a99",
                                 "#e31a1c"),
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
  scale_colour_manual(values = c("#fbb4b9", "#f768a1",
                                 "#ae017e"),
                      name = "Frost")

# Write
dplot_colour <- plot_grid(dplot_frost, dplot_martin,
                          dplot_hammer, dplot_bespoke,
                          labels = "AUTO", label_size = 32)
ggsave("figures/Discrimination_colour.png", dplot_colour,
       height = 10, width = 12.5)

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

# Bespoke
cplot_bespoke <- calibration_df %>% 
  filter(model == "bespoke") %>% 
  ggplot(aes(predicted, observed ,
             colour = dataset))+
  geom_abline(slope = 1, intercept = 0,
              size = 1)+
  geom_path(size = 1, alpha = 0.8)+
  geom_pointrange(aes(ymin = observed - error ,
                      ymax = observed + error,
                      y = observed,
                      x = predicted,
                      colour = dataset), size = 0.25,
                  alpha = 0.8) +
  theme_classic(20)+
  theme(legend.position = "top")+
  labs(x = "Predicted readmission",
       y = "Observed readmission")+
  scale_colour_manual(values = c("#bdc9e1", "#0570b0"),
                      name = "Bespoke")+
  coord_cartesian(xlim = c(0,25),
                  ylim = c(0,30))


# Hammer
cplot_hammer <- calibration_df %>% 
  filter(model == "hammer") %>% 
  ggplot(aes(predicted, observed ,
             colour = dataset))+
  geom_abline(slope = 1, intercept = 0,
              size = 1)+
  geom_path(size = 1, alpha = 0.8)+
  geom_pointrange(aes(ymin = observed - error ,
                      ymax = observed + error,
                      y = observed,
                      x = predicted,
                      colour = dataset), size = 0.25,
                  alpha = 0.8) +
  theme_classic(20)+
  theme(legend.position = "top")+
  labs(x = "Predicted readmission",
       y = "Observed readmission")+
  scale_colour_manual(values = c("#c2e699", "#78c679",
                                 "#238443"),
                      name = "Hammer")+
  coord_cartesian(xlim = c(0,25),
                  ylim = c(0,30))

# Martin
cplot_martin <- calibration_df %>% 
  filter(model == "martin") %>% 
  ggplot(aes(predicted, observed ,
             colour = dataset))+
  geom_abline(slope = 1, intercept = 0,
              size = 1)+
  geom_path(size = 1)+
  geom_pointrange(aes(ymin = observed - error ,
                      ymax = observed + error,
                      y = observed,
                      x = predicted,
                      colour = dataset), size = 0.25,
                  alpha = 0.8) +
  theme_classic(20)+
  theme(legend.position = "top")+
  labs(x = "Predicted readmission",
       y = "Observed readmission")+
  scale_colour_manual(values = c("#fecc5c","#fb9a99",
                                 "#e31a1c"),
                      name = "Martin")+
  coord_cartesian(xlim = c(0,25),
                  ylim = c(0,30))


# Frost
cplot_frost <- calibration_df %>% 
  filter(model == "frost") %>% 
  ggplot(aes(predicted, observed ,
             colour = dataset))+
  geom_abline(slope = 1, intercept = 0,
              size = 1)+
  geom_path(size = 1)+
  geom_pointrange(aes(ymin = observed - error ,
                      ymax = observed + error,
                      y = observed,
                      x = predicted,
                      colour = dataset), size = 0.25,
                  alpha = 0.8) +
  theme_classic(20)+
  theme(legend.position = "top")+
  labs(x = "Predicted readmission",
       y = "Observed readmission")+
  scale_colour_manual(values = c("#fbb4b9", "#f768a1",
                                 "#ae017e"),
                      name = "Frost")+
  coord_cartesian(xlim = c(0,25),
                  ylim = c(0,30))

# Write
cplot_colour <- plot_grid(cplot_frost, cplot_martin,
                          cplot_hammer, cplot_bespoke,
                          labels = "AUTO", label_size = 32)
ggsave("figures/Calibration_colour.png", cplot_colour,
       height = 8, width = 15)
