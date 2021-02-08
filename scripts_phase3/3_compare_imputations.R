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
pal <- brewer.pal(9, "Set1")[-6]

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
                      distance = 0,
                      mean_error = 0,
                      max = 0,
                      max_error = 0)

# Mean imputation-----
results_mean <- read_rds("data/impute/average.RDS")

# Extraction function
#x <- results_mean[[1]]
extraction_mean <- function(x)
{
  # Load packages
  require(magrittr)
  
  # Extract probabilities
  zero_probs <- x["zero_probs"][,1] %>% inverse_logit()
  int_mean_probs <- x["int_mean_probs"][,1] %>% inverse_logit()
  int_median_probs <- x["int_median_probs"][,1] %>% inverse_logit()
  
  # Measure distance
  zero_dist <- abs(zero_probs - probs)
  mean_dist <- abs(int_mean_probs - probs)
  median_dist <- abs(int_median_probs - probs)
  
  # Output
  output <- data.frame(method = c("assume_zero", "mean", "median"),
                       split = x["split"][1,1], n = x["n"][1,1],
                       mean_distance = c(mean(zero_dist),
                                         mean(mean_dist),
                                         mean(median_dist)),
                       max_distance = c(quantile(zero_dist, 0.95),
                                        quantile(mean_dist, 0.95),
                                        quantile(median_dist, 0.95)))
  return(output)
}

# Run extraction
psnice(value = 19)
cl <- makeCluster(ifelse(detectCores() <= 12,
                         detectCores() - 1,
                         12))
clusterExport(cl, varlist = c("probs",
                              "inverse_logit"))
metrics_mean <- parLapply(cl, 
                          results_mean,
                          extraction_mean)
stopCluster(cl)

# Summarise
summary_mean <- do.call(rbind.data.frame, metrics_mean) %>% 
  na.omit() %>% 
  group_by(method, split) %>% 
  summarise(distance = mean(mean_distance * 100),
            mean_error = sd(mean_distance * 100),
            max = mean(max_distance * 100),
            max_error = sd(max_distance * 100))

# Plot
summary_mean %>% 
  pivot_longer(c(3,5)) %>% 
  mutate(error = case_when(name == "distance" ~ mean_error,
                           name == "max" ~ max_error)) %>% 
  ggplot(aes(x = split,
             y = value))+
  geom_ribbon(aes(ymin = value - error,
                  ymax = value + error,
                  fill = method), alpha = 0.1)+
  geom_path(aes(colour = method))+
  geom_hline(linetype = "dashed",
             colour =  "black",
             yintercept = 0)+
  theme_classic(20)+
  theme(legend.position = "top")+
  facet_wrap(~name, scales = "free_y")+
  scale_colour_manual(values = pal,
                      name = "")+
  scale_fill_manual(values = pal,
                      name = "")+
  labs(x = "Proportion missing",
       y = "Prediction distance")


# Hot-deck imputation--------
results_hotdeck <- read_rds("data/impute/hotdeck.RDS")

# Extraction function
#x <- results_hotdeck[[1]]
extraction_hotdeck <- function(x)
{
  # Load packages
  require(magrittr)
  
  # List pooling methods
  pooling_methods <- c("admit.chronic.age",
                       "admit.age", "admit.chronic",
                       "age.chronic", "admit",
                       "age", "chronic")
  
  # Extract probabilities
  probs_admit.chronic.age <- x$probs[,1] %>% inverse_logit()
  probs_admit.age <- x$probs[,2] %>% inverse_logit()
  probs_admit.chronic <- x$probs[,3] %>% inverse_logit()
  probs_age.chronic <- x$probs[,4] %>% inverse_logit()
  probs_admit <- x$probs[,5] %>% inverse_logit()
  probs_age <- x$probs[,6] %>% inverse_logit()
  probs_chronic <- x$probs[,7] %>% inverse_logit()
  
  # Measure distance
  dist_admit.chronic.age <- abs(probs_admit.chronic.age - probs)
  dist_admit.age <- abs(probs_admit.age - probs)
  dist_admit.chronic <- abs(probs_admit.chronic - probs)
  dist_age.chronic <- abs(probs_age.chronic - probs)
  dist_admit <- abs(probs_admit - probs)
  dist_chronic <- abs(probs_age - probs)
  dist_age <- abs(probs_chronic - probs)

  # Output
  output <- data.frame(method = pooling_methods,
                       split = x$split[1], n = x$N[1],
                       mean_distance = c(mean(dist_admit.chronic.age),
                                         mean(dist_admit.age),
                                         mean(dist_admit.chronic),
                                         mean(dist_age.chronic),
                                         mean(dist_admit),
                                         mean(dist_chronic),
                                         mean(dist_age)),
                       max_distance = c(quantile(dist_admit.chronic.age, 0.95),
                                        quantile(dist_admit.age, 0.95),
                                        quantile(dist_admit.chronic, 0.95),
                                        quantile(dist_age.chronic, 0.95),
                                        quantile(dist_admit, 0.95),
                                        quantile(dist_chronic, 0.95),
                                        quantile(dist_age, 0.95)))
  return(output)
}

# Run extraction
psnice(value = 19)
cl <- makeCluster(ifelse(detectCores() <= 12,
                         detectCores() - 1,
                         12))
clusterExport(cl, varlist = c("probs", "inverse_logit"))
metrics_hotdeck <- parLapply(cl, 
                          results_hotdeck,
                          extraction_hotdeck)
stopCluster(cl)

# Summarise
summary_hotdeck <- do.call(rbind.data.frame, metrics_hotdeck) %>% 
  na.omit() %>% 
  group_by(method, split)  %>% 
  summarise(distance = mean(mean_distance * 100),
            mean_error = sd(mean_distance * 100),
            max = mean(max_distance * 100),
            max_error = sd(max_distance * 100))

# Plot -> All pooling methods
summary_hotdeck %>% 
  pivot_longer(c(3,5)) %>% 
  mutate(error = case_when(name == "distance" ~ mean_error,
                           name == "max" ~ max_error)) %>% 
  ggplot(aes(x = split,
             y = value))+
  geom_ribbon(aes(ymin = value - error,
                  ymax = value + error,
                  fill = method), alpha = 0.1)+
  geom_path(aes(colour = method))+
  geom_hline(linetype = "dashed",
             colour =  "black",
             yintercept = 0)+
  theme_classic(20)+
  theme(legend.position = "top")+
  facet_wrap(~name, scales = "free_y")+
  scale_colour_manual(values = pal,
                      name = "")+
  scale_fill_manual(values = pal,
                    name = "")+
  labs(x = "Proportion missing",
       y = "Prediction distance")

# Plot
summary_hotdeck %>% 
  pivot_longer(c(3,5)) %>% 
  mutate(error = case_when(name == "distance" ~ mean_error,
                           name == "max" ~ max_error)) %>% 
  filter(method == "admit.chronic.age") %>% 
  ggplot(aes(x = split,
             y = value))+
  geom_ribbon(aes(ymin = value - error,
                  ymax = value + error,
                  fill = method), alpha = 0.1)+
  geom_path(aes(colour = method))+
  geom_hline(linetype = "dashed",
             colour =  "black",
             yintercept = 0)+
  theme_classic(20)+
  theme(legend.position = "top")+
  facet_wrap(~name, scales = "free_y")+
  scale_colour_manual(values = pal,
                      name = "")+
  scale_fill_manual(values = pal,
                    name = "")+
  labs(x = "Proportion missing",
       y = "Prediction distance")

# Compare methods---------

# Combine
comparison_df <- rbind(
  results,
  summary_mean %>% filter(method != "median"),
  summary_hotdeck %>% filter(method == "admit.chronic.age") %>% 
    mutate(method = "hotdeck")
)

# Plot
comparison_df %>% 
  pivot_longer(c(3,5)) %>% 
  mutate(error = case_when(name == "distance" ~ mean_error,
                           name == "max" ~ max_error)) %>% 
  ggplot(aes(x = split,
             y = value))+
  geom_ribbon(aes(ymin = value - error,
                  ymax = value + error,
                  fill = method), alpha = 0.1)+
  geom_path(aes(colour = method))+
  geom_hline(linetype = "dashed",
             colour =  "black",
             yintercept = 0)+
  theme_classic(20)+
  theme(legend.position = "top")+
  facet_wrap(~name, scales = "free_y")+
  scale_colour_manual(values = pal,
                      name = "")+
  scale_fill_manual(values = pal,
                    name = "")+
  labs(x = "Proportion missing",
       y = "Prediction distance")


# Assess optimality
comparison_df %>% 
  ggplot(aes(x = distance,
             y = max,
             colour = method))+
  geom_errorbar(aes(ymin = max - max_error,
                    ymax = max + max_error))+
  geom_errorbarh(aes(xmin = distance - mean_error,
                    xmax = distance + mean_error))+
  geom_path()+
  geom_point(size = 4, aes(fill = method),
             colour = "black", shape = 21)+
  theme_classic(20)+
  theme(legend.position = "top")+
  scale_colour_manual(values = brewer.pal(8,"Set1")[-6],
                      name = "")+
  scale_fill_manual(values = brewer.pal(8,"Set1")[-6],
                    name = "")
