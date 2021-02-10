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
require(scales)

# Functions
source("functions/inverse_logit.R")
source("functions/calibration.R")
source("functions/calculate_apache_score.R")
source("functions/distance_4D.R")

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
                      max_error = 0,
                      dist_4D = 0,
                      error_4D = 0)

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
  
  # Measure 4-dimensional distance
  zero_4D <- distance_4D(results, auc_zero@y.values[[1]], cal_zero$statistic,
              mean(zero_dist), quantile(zero_dist, 0.95))
  mean_4D <- distance_4D(results, auc_mean@y.values[[1]], cal_mean$statistic,
              mean(mean_dist), quantile(mean_dist, 0.95))
  median_4D <- distance_4D(results, auc_median@y.values[[1]], cal_median$statistic,
              mean(median_dist), quantile(median_dist, 0.95))
  
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
                                        quantile(median_dist, 0.95)),
                       md_dist = c(zero_4D, mean_4D, median_4D))
  return(output)
}

# Run extraction
psnice(value = 19)
cl <- makeCluster(ifelse(detectCores() <= 12,
                         detectCores() - 1,
                         12))
clusterExport(cl, varlist = c("full_data", "inverse_logit",
                              "probs", "distance_4D",
                              "results"))
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
            max_error = sd(max_distance * 100),
            dist_4D = mean(md_dist),
            error_4D = sd(md_dist))

# Plot metrics
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

# Plot 4D
# summary_mean %>% 
#   ggplot(aes(x = split,
#              y = dist_4D))+
#   geom_ribbon(aes(ymin = dist_4D - error_4D,
#                   ymax = dist_4D + error_4D,
#                   fill = method), alpha = 0.1)+
#   geom_path(aes(colour = method))+
#   theme_classic(20)+
#   theme(legend.position = "top")+
#   scale_colour_manual(values = pal,
#                       name = "")+
#   scale_fill_manual(values = pal,
#                     name = "")+
#   labs(x = "Proportion missing",
#        y = "4D distance")

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
  
  # Measure 4-dimensional distance
  # admit.chronic.age
  admit.chronic.age_4D <- distance_4D(results, auc_admit.chronic.age@y.values[[1]], 
                                      cal_admit.chronic.age$statistic,
                                      mean(dist_admit.chronic.age), 
                                      quantile(dist_admit.chronic.age, 0.95))
  # admit.age
  admit.age_4D <- distance_4D(results, auc_admit.age@y.values[[1]], 
                              cal_admit.age$statistic,
                              mean(dist_admit.age), 
                              quantile(dist_admit.age, 0.95))
  # admit.chronic
  admit.chronic_4D <- distance_4D(results, auc_admit.chronic@y.values[[1]], 
                              cal_admit.chronic$statistic,
                              mean(dist_admit.chronic), 
                              quantile(dist_admit.chronic, 0.95))
  # age.chronic
  age.chronic_4D <- distance_4D(results, auc_age.chronic@y.values[[1]], 
                              cal_age.chronic$statistic,
                              mean(dist_age.chronic), 
                              quantile(dist_age.chronic, 0.95))
  # admit
  admit_4D <- distance_4D(results, auc_admit@y.values[[1]], 
                              cal_admit$statistic,
                              mean(dist_admit), 
                              quantile(dist_admit, 0.95))
  # chronic
  chronic_4D <- distance_4D(results, auc_chronic@y.values[[1]], 
                              cal_chronic$statistic,
                              mean(dist_chronic), 
                              quantile(dist_chronic, 0.95))
  # age
  age_4D <- distance_4D(results, auc_age@y.values[[1]], 
                              cal_age$statistic,
                              mean(dist_age), 
                              quantile(dist_age, 0.95))
  
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
                                  cal_age$statistic),
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
                                        quantile(dist_age, 0.95)),
                       md_dist = c(admit.chronic.age_4D,
                                   admit.age_4D,
                                   admit.chronic_4D,
                                   age.chronic_4D,
                                   admit_4D,
                                   chronic_4D,
                                   age_4D))
  return(output)
}

# Run extraction
psnice(value = 19)
cl <- makeCluster(ifelse(detectCores() <= 12,
                         detectCores() - 1,
                         12))
clusterExport(cl, varlist = c("full_data", "inverse_logit",
                              "probs", "distance_4D",
                              "results"))
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
            calibration_error = sd(hoslem),
            distance = mean(mean_distance * 100),
            mean_error = sd(mean_distance * 100),
            max = mean(max_distance * 100),
            max_error = sd(max_distance * 100),
            dist_4D = mean(md_dist),
            error_4D = sd(md_dist))

# Plot  -> All pooling methods
summary_hotdeck %>% 
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

# Plot 4D
# summary_hotdeck %>% 
#   ggplot(aes(x = split,
#              y = dist_4D))+
#   geom_ribbon(aes(ymin = dist_4D - error_4D,
#                   ymax = dist_4D + error_4D,
#                   fill = method), alpha = 0.1)+
#   geom_path(aes(colour = method))+
#   theme_classic(20)+
#   theme(legend.position = "top")+
#   scale_colour_manual(values = pal,
#                       name = "")+
#   scale_fill_manual(values = pal,
#                     name = "")+
#   labs(x = "Proportion missing",
#        y = "4D distance")

# Compare methods---------

# Combine
comparison_df <- rbind(
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
  geom_path(aes(colour = method), size = 1)+
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
# comparison_df %>% 
#   ggplot(aes(x = split,
#              y = dist_4D))+
#   geom_ribbon(aes(ymin = dist_4D - error_4D,
#                   ymax = dist_4D + error_4D,
#                   fill = method), alpha = 0.1)+
#   geom_path(aes(colour = method))+
#   theme_classic(20)+
#   theme(legend.position = "top")+
#   scale_colour_manual(values = pal,
#                       name = "")+
#   scale_fill_manual(values = pal,
#                     name = "")+
#   labs(x = "Proportion missing",
#        y = "4D distance")

# Multidimensional distance---------

# Rescale
rescale_df <- comparison_df %>%
  ungroup() %>%
  mutate(
    discrimination = rescale(discrimination,
                             from = c(
                               results$discrimination,
                               min(discrimination)
                             )),
    calibration = rescale(calibration,
                          from = c(results$calibration,
                                   max(calibration))),
    distance = rescale(distance,
                       from = c(results$distance,
                                max(distance))),
    max = rescale(max,
                  from = c(results$max,
                           max(max)))
  ) 

# Measure distance
md_dist <- distance_4D(data.frame(discrimination = 0,
                       calibration = 0,
                       distance = 0,
                       max = 0), rescale_df$discrimination,
            rescale_df$calibration,
            rescale_df$distance,
            rescale_df$max)

rescale_df %<>% 
  mutate(md_dist)


# Plot 
rescale_df %>%
  ggplot(aes(x = split,
             y = md_dist))+
  geom_path(aes(colour = method),
            size = 1)+
  geom_point(size = 4, aes(fill = method),
             colour = "black", shape = 21)+
  theme_classic(20)+
  theme(legend.position = "top")+
  scale_colour_manual(values = pal,
                      name = "")+
  scale_fill_manual(values = pal,
                    name = "")+
  labs(x = "Proportion missing",
       y = "Multidimensional distance")
