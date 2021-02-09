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
source("functions/calculate_apache_score.R")

# Load data
full_data <- read_csv("data/impute/complete_cases.csv") 
pal <- brewer.pal(9, "Set1")[-6]

# Complete cases------

# Score apache
full_data %<>% 
  mutate(apache_II_discharge = calculate_apache_scores(.))

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
                      calibration_error = 0,
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
                       AUC = c(auc_zero@y.values[[1]],
                               auc_mean@y.values[[1]],
                               auc_median@y.values[[1]]),
                       hoslem = c(cal_zero$statistic,
                                  cal_mean$statistic,
                                  cal_median$statistic),
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
clusterExport(cl, varlist = c("full_data", "inverse_logit",
                              "probs"))
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
            calibration_error = sd(hoslem),
            distance = mean(mean_distance * 100),
            mean_error = sd(mean_distance * 100),
            max = mean(max_distance * 100),
            max_error = sd(max_distance * 100))

# Plot
summary_mean %>% 
  pivot_longer(c(3,5,7,9)) %>% 
  mutate(error = case_when(name == "distance" ~ mean_error,
                           name == "max" ~ max_error,
                           name == "calibration" ~ calibration_error,
                           name == "discrimination" ~ discrimination_error)) %>% 
  ggplot(aes(x = split,
             y = value))+
  geom_ribbon(aes(ymin = value - error,
                  ymax = value + error,
                  fill = method), alpha = 0.1)+
  geom_path(aes(colour = method))+
  geom_hline(data = results %>% pivot_longer(c(3,5,7,9)),
             aes(yintercept = value), size = 1,
             linetype = "dashed")+
  theme_classic(20)+
  theme(legend.position = "top")+
  facet_wrap(~name, scales = "free_y")+
  scale_colour_manual(values = pal,
                      name = "")+
  scale_fill_manual(values = pal,
                    name = "")+
  labs(x = "Proportion missing")



# Compare methods---------

# Combine
comparison_df <- rbind(
  results,
  summary_mean %>% filter(method != "mean"),
  summary_hotdeck %>% filter(method == "admit.chronic.age") %>% 
    mutate(method = "hotdeck")
)

# Plot each dimension
comparison_df %>% 
  pivot_longer(c(3,5,7,9)) %>% 
  mutate(error = case_when(name == "distance" ~ mean_error,
                           name == "max" ~ max_error,
                           name == "calibration" ~ calibration_error,
                           name == "discrimination" ~ discrimination_error)) %>% 
  ggplot(aes(x = split,
             y = value))+
  geom_ribbon(aes(ymin = value - error,
                  ymax = value + error,
                  fill = method), alpha = 0.1)+
  geom_path(aes(colour = method))+
  geom_hline(data = results %>% pivot_longer(c(3,5,7,9)),
             aes(yintercept = value), size = 1,
             linetype = "dashed")+
  theme_classic(20)+
  theme(legend.position = "top")+
  facet_wrap(~name, scales = "free_y")+
  scale_colour_manual(values = pal,
                      name = "")+
  scale_fill_manual(values = pal,
                    name = "")+
  labs(x = "Proportion missing")

# Plot multidimensional