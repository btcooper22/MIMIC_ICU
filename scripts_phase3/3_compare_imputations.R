# Packages
require(dplyr)
require(magrittr)
require(readr)
require(tibble)
require(foreach)
require(stringr)
require(doParallel)
require(tools)
require(tidyr)
require(ROCR)
require(ggplot2)
require(ResourceSelection)
require(cowplot)
require(RColorBrewer)

# Functions
source("functions/inverse_logit.R")
source("functions/calibration.R")

# Load data
full_data <- read_csv("data/impute/complete_cases.csv") 

# Complete cases------

# Create model and predict
apache_model <- glm(readmission ~ apache_II_discharge,
                    data = full_data, family = "binomial")
write_rds(apache_model, "models/apache_model.RDS")
probs <- predict(apache_model, newdata = full_data) %>% inverse_logit()

# Assess ROC
pred <- prediction(probs, full_data$readmission)
perf <- performance(pred, "tpr", "fpr")
auc <- performance(pred, measure = "auc")
auc@y.values[[1]]

# Plot ROC
data.frame(x = perf@x.values[[1]],
           y = perf@y.values[[1]]) %>% 
  ggplot(aes(x, y))+
  geom_abline(slope = 1, intercept = 0,
              linetype = "dotted",
              size = 1)+
  geom_path(size = 1)+
  labs(x = "1 - Specificity",
       y = "Sensitivity")+
  theme_classic(20)

# Split into deciles
deciles_df <- tibble(
  patient_id = 1:nrow(full_data),
  readmission = full_data$readmission,
  probs,
  decile = ntile(probs, 10))

# Calculate calibration
cal <- calibration(deciles_df, "probs", "decile")

# Plot calibration
cal %>% 
  ggplot(aes(predicted, observed))+
  geom_abline(slope = 1, intercept = 0,
              size = 1)+
  geom_path(size = 1)+
  geom_pointrange(aes(ymin = observed - error ,
                      ymax = observed + error,
                      y = observed,
                      x = predicted)) +
  theme_classic(20)+
  theme(legend.position = "top")+
  labs(x = "Predicted readmission",
       y = "Observed readmission")

# Assess calibration
hoslem_full <- hoslem.test(full_data$readmission,
            probs, 10)


# Tabulate
results <- data.frame(method = "full_dataset",
                      split = 0,
                      discrimination = auc@y.values[[1]],
                      discrimination_error = 0,
                      calibration = hoslem_full$statistic,
                      calibration_error = 0)

# Mean imputation-----
results_mean <- read_rds("data/impute/average.RDS")

# Extraction function
#x <- results_mean[[1]]
extraction_mean <- function(x)
{
  # Load packages
  require(magrittr)
  require(ROCR)
  require(ResourceSelection)
  
  # Calculate discrimination
  auc_zero <- prediction(x["zero_probs"] %>% inverse_logit(),
                     full_data$readmission) %>% 
    performance(measure = "auc")
  auc_mean <- prediction(x["int_mean_probs"] %>% inverse_logit(),
                         full_data$readmission) %>% 
    performance(measure = "auc")
  auc_median <- prediction(x["int_median_probs"] %>% inverse_logit(),
                         full_data$readmission) %>% 
    performance(measure = "auc")
  
  # Calculate calibration
  cal_zero <- hoslem.test(full_data$readmission,
                          x["zero_probs"][,1] %>% inverse_logit(), 10)
  cal_mean <- hoslem.test(full_data$readmission,
                          x["int_mean_probs"][,1] %>% inverse_logit(), 10)
  cal_median <- hoslem.test(full_data$readmission,
                            x["int_median_probs"][,1] %>% inverse_logit(), 10)
  
  # Output
  output <- data.frame(method = c("assume_zero", "mean", "median"),
                       split = x["split"][1,1], n = x["n"][1,1],
                       AUC = c(auc_zero@y.values[[1]],
                               auc_mean@y.values[[1]],
                               auc_median@y.values[[1]]),
                       hoslem = c(cal_zero$statistic,
                                 cal_mean$statistic,
                                 cal_median$statistic))
  return(output)
}

# Run extraction
psnice(value = 19)
cl <- makeCluster(ifelse(detectCores() <= 12,
                         detectCores() - 1,
                         12))
clusterExport(cl, varlist = c("full_data", "inverse_logit"))
metrics_mean <- parLapply(cl, 
                          results_mean,
                          extraction_mean)
stopCluster(cl)

# Summarise
summary_mean <- do.call(rbind.data.frame, metrics_mean) %>% 
  na.omit() %>% 
  group_by(method, split) %>% 
  summarise(discrimination = mean(AUC),
            discrimination_error = sd(AUC),
            calibration = mean(hoslem),
            calibration_error = sd(hoslem))

# Plot discrimination
summary_mean %>% 
  ggplot(aes(x = split,
             y = discrimination))+
  geom_ribbon(aes(ymin = discrimination - discrimination_error,
                  ymax = discrimination + discrimination_error),
              fill = "Gray90")+
  geom_path()+
  geom_hline(linetype = "dashed",
             colour =  "red",
             yintercept = results$discrimination)+
  facet_wrap(~method)+
  theme_classic(20)

# Plot calibration
summary_mean %>% 
  ggplot(aes(x = split,
             y = calibration))+
  geom_ribbon(aes(ymin = calibration - calibration_error,
                  ymax = calibration + calibration_error),
              fill = "Gray90")+
  geom_path()+
  geom_hline(linetype = "dashed",
             colour =  "red",
             yintercept = results$calibration)+
  facet_wrap(~method)+
  theme_classic(20)

# Exclusion-------
results_exclusion <- read_rds("data/impute/exclusion.RDS")

# Extraction function
x <- results_exclusion[[1]]
extraction_exclusion <- function(x)
{
  # Load packages
  require(magrittr)
  require(ROCR)
  require(ResourceSelection)
  
  # Combine and trim data
  exclude_df <- data.frame(readmission = full_data$readmission,
                           x["probs"]) %>% 
    na.omit()
  
  # Filter high-exclusion cases
  if(length(unique(exclude_df$readmission)) == 2)
  {
    # Calculate discrimination
    auc_exclude <- prediction(exclude_df$probs %>% inverse_logit(),
                              exclude_df$readmission) %>% 
      performance(measure = "auc")
    
    # Calculate calibration
    cal_exclude <- hoslem.test(exclude_df$readmission,
                               exclude_df$probs %>% inverse_logit(), 10)
    
    # Output
    output <- data.frame(method = "exclude",
                         split = x["split"][1,1], n = x["n"][1,1],
                         AUC = auc_exclude@y.values[[1]],
                         hoslem = cal_exclude$statistic)
  }else
  {
    # Output
    output <- data.frame(method = "exclude",
                         split = x["split"][1,1], n = x["n"][1,1],
                         AUC = NA,
                         hoslem = NA)
  }
  

  return(output)
}

# Run extraction
psnice(value = 19)
cl <- makeCluster(ifelse(detectCores() <= 12,
                         detectCores() - 1,
                         12))
clusterExport(cl, varlist = c("full_data", "inverse_logit"))
metrics_exclusion <- parLapply(cl, 
                          results_exclusion,
                          extraction_exclusion)
stopCluster(cl)

# Summarise
summary_exclusion <- do.call(rbind.data.frame, metrics_exclusion) %>% 
  na.omit() %>% 
  group_by(method, split) %>% 
  summarise(discrimination = mean(AUC),
            discrimination_error = sd(AUC),
            calibration = mean(hoslem),
            calibration_error = sd(hoslem))

# Plot discrimination
summary_exclusion %>% 
  ggplot(aes(x = split,
             y = discrimination))+
  geom_ribbon(aes(ymin = discrimination - discrimination_error,
                  ymax = discrimination + discrimination_error),
              fill = "Gray90")+
  geom_path()+
  geom_hline(linetype = "dashed",
             colour =  "red",
             yintercept = results$discrimination)+
  facet_wrap(~method)+
  theme_classic(20)

# Plot calibration
summary_exclusion %>% 
  ggplot(aes(x = split,
             y = calibration))+
  geom_ribbon(aes(ymin = calibration - calibration_error,
                  ymax = calibration + calibration_error),
              fill = "Gray90")+
  geom_path()+
  geom_hline(linetype = "dashed",
             colour =  "red",
             yintercept = results$calibration)+
  facet_wrap(~method)+
  theme_classic(20)

# Hot-deck imputation--------
results_hotdeck <- read_rds("data/impute/hotdeck.RDS")

# Extraction function
#x <- results_hotdeck[[1]]
extraction_hotdeck <- function(x)
{
  # Load packages
  require(magrittr)
  require(ROCR)
  require(ResourceSelection)
  
  # List pooling methods
  pooling_methods <- c("admit.chronic.age",
                       "admit.age", "admit.chronic",
                       "age.chronic", "admit",
                       "age", "chronic")
  
  # Calculate discrimination
  auc_admit.chronic.age <- prediction(x$probs[,1] %>% inverse_logit(),
                         full_data$readmission) %>% 
    performance(measure = "auc") #admit.chronic.age
  
  auc_admit.age <- prediction(x$probs[,2] %>% inverse_logit(),
                                      full_data$readmission) %>% 
    performance(measure = "auc") #admit.age
  
  auc_admit.chronic <- prediction(x$probs[,3] %>% inverse_logit(),
                                      full_data$readmission) %>% 
    performance(measure = "auc") #admit.chronic
  
  auc_age.chronic <- prediction(x$probs[,4] %>% inverse_logit(),
                                      full_data$readmission) %>% 
    performance(measure = "auc") #age.chronic
  
  auc_admit <- prediction(x$probs[,5] %>% inverse_logit(),
                                      full_data$readmission) %>% 
    performance(measure = "auc") #admit
  
  auc_age <- prediction(x$probs[,6] %>% inverse_logit(),
                                      full_data$readmission) %>% 
    performance(measure = "auc") #age
  
  auc_chronic <- prediction(x$probs[,7] %>% inverse_logit(),
                                      full_data$readmission) %>% 
    performance(measure = "auc") #chronic


  # Calculate calibration
  cal_admit.chronic.age <- hoslem.test(full_data$readmission,
                          x$probs[,1] %>% inverse_logit(), 10)
  
  cal_admit.age <- hoslem.test(full_data$readmission,
                                       x$probs[,2] %>% inverse_logit(), 10)
  
  cal_admit.chronic <- hoslem.test(full_data$readmission,
                                       x$probs[,3] %>% inverse_logit(), 10)
  
  cal_age.chronic <- hoslem.test(full_data$readmission,
                                       x$probs[,4] %>% inverse_logit(), 10)
  
  cal_admit <- hoslem.test(full_data$readmission,
                                       x$probs[,5] %>% inverse_logit(), 10)
  
  cal_chronic <- hoslem.test(full_data$readmission,
                                       x$probs[,6] %>% inverse_logit(), 10)
  
  cal_age <- hoslem.test(full_data$readmission,
                                       x$probs[,7] %>% inverse_logit(), 10)
  
  # Output
  output <- data.frame(method = pooling_methods,
                       split = x$split[1], n = x$N[1],
                       AUC = c(auc_admit.chronic.age@y.values[[1]],
                               auc_admit.age@y.values[[1]],
                               auc_admit.chronic@y.values[[1]],
                               auc_age.chronic@y.values[[1]],
                               auc_admit@y.values[[1]],
                               auc_age@y.values[[1]],
                               auc_chronic@y.values[[1]]),
                       hoslem = c(cal_admit.chronic.age$statistic,
                                  cal_admit.age$statistic,
                                  cal_admit.chronic$statistic,
                                  cal_age.chronic$statistic,
                                  cal_admit$statistic,
                                  cal_chronic$statistic,
                                  cal_age$statistic))
  return(output)
}

# Run extraction
psnice(value = 19)
cl <- makeCluster(ifelse(detectCores() <= 12,
                         detectCores() - 1,
                         12))
clusterExport(cl, varlist = c("full_data", "inverse_logit"))
metrics_hotdeck <- parLapply(cl, 
                          results_hotdeck,
                          extraction_hotdeck)
stopCluster(cl)

# Summarise
summary_hotdeck <- do.call(rbind.data.frame, metrics_hotdeck) %>% 
  na.omit() %>% 
  group_by(method, split) %>% 
  summarise(discrimination = mean(AUC),
            discrimination_error = sd(AUC),
            calibration = mean(hoslem),
            calibration_error = sd(hoslem))

# Plot discrimination -> All pooling methods
summary_hotdeck %>% 
  ggplot(aes(x = split,
             y = discrimination))+
  geom_path(aes(colour = method), size = 1)+
  geom_hline(linetype = "dashed",
             colour =  "red",
             yintercept = results$AUROC)+
  theme_classic(20)+
  theme(legend.position = "top")+
  scale_colour_manual(values = brewer.pal(8,"Set1")[-6],
                      name = "")

# Plot calibration -> All pooling methods
summary_hotdeck %>% 
  ggplot(aes(x = split,
             y = calibration))+
  geom_path(aes(colour = method), size = 1)+
  geom_hline(linetype = "dashed",
             colour =  "red",
             yintercept = results$chisq)+
  theme_classic(20)+
  theme(legend.position = "top")+
  scale_colour_manual(values = brewer.pal(8,"Set1")[-6],
                      name = "")

# Plot discrimination
summary_hotdeck %>% 
  filter(method == "admit.chronic.age") %>% 
  ggplot(aes(x = split,
             y = discrimination))+
  geom_ribbon(aes(ymin = discrimination - discrimination_error,
                  ymax = discrimination + discrimination_error),
              fill = "Gray90")+
  geom_path()+
  geom_hline(linetype = "dashed",
             colour =  "red",
             yintercept = results$AUROC)+
  facet_wrap(~method)+
  theme_classic(20)

# Plot calibration
summary_hotdeck %>% 
    filter(method == "admit.chronic.age") %>% 
  ggplot(aes(x = split,
             y = calibration))+
  geom_ribbon(aes(ymin = calibration - calibration_error,
                  ymax = calibration + calibration_error),
              fill = "Gray90")+
  geom_path()+
  geom_hline(linetype = "dashed",
             colour =  "red",
             yintercept = results$chisq)+
  facet_wrap(~method)+
  theme_classic(20)
  
# Compare methods---------

# Combine
comparison_df <- rbind(
  results, summary_exclusion,
  summary_mean %>% filter(method != "median"),
  summary_hotdeck %>% filter(method == "admit.chronic.age") %>% 
    mutate(method = "hotdeck")
)

# Discrimination
discrim_plot <- comparison_df %>% 
  filter(method != "exclude") %>% 
  ggplot(aes(x = split,
             y = discrimination))+
  # geom_ribbon(aes(ymin = discrimination - discrimination_error,
  #                 ymax = discrimination + discrimination_error,
  #                 fill = method), alpha = 0.1)+
  geom_path(aes(colour = method), size = 1)+
  geom_hline(linetype = "dashed",
             colour =  "black",
             yintercept = results$discrimination)+
  theme_classic(20)+
  theme(legend.position = "top")+
  scale_colour_manual(values = brewer.pal(8,"Set1")[-6],
                      name = "")+
  scale_fill_manual(values = brewer.pal(8,"Set1")[-6],
                      name = "")

# Calibration
calib_plot <- comparison_df %>% 
  filter(method != "exclude") %>% 
  ggplot(aes(x = split,
             y = calibration))+
  geom_path(aes(colour = method), size = 1)+
  geom_hline(linetype = "dashed",
             colour =  "black",
             yintercept = results$calibration)+
  theme_classic(20)+
  theme(legend.position = "top")+
  scale_colour_manual(values = brewer.pal(8,"Set1")[-6],
                      name = "")+
  scale_fill_manual(values = brewer.pal(8,"Set1")[-6],
                    name = "")+
  scale_y_reverse()

# Assess optimality
opt_plot <- comparison_df %>% 
  filter(method != "exclude") %>% 
  ggplot(aes(x = discrimination,
             y = calibration,
             colour = method))+
  geom_point(size = 4)+
  geom_path()+
  theme_classic(20)+
  theme(legend.position = "top")+
  scale_colour_manual(values = brewer.pal(8,"Set1")[-6],
                      name = "")+
  scale_fill_manual(values = brewer.pal(8,"Set1")[-6],
                    name = "")+
  scale_y_reverse()

plot_grid(discrim_plot, calib_plot, opt_plot, nrow = 1)
