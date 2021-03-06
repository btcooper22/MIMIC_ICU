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
require(Metrics)

# Functions
source("functions/inverse_logit.R")
source("functions/calibration.R")
source("functions/calculate_apache_score.R")
source("functions/distance_4D.R")
source("functions/parallel_setup.R")

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
  require(Metrics)
  
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
  
  # Measure RMSE
  rmse(zero_probs, probs)
  
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
  filter(split < 0.6) %>% 
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
  filter(split < 0.6) %>% 
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

# PCA imputation----
results_PCA <- read_rds("data/impute/PCA.RDS")

# Extraction function
#x <- results_PCA[[1]]
extraction_PCA <- function(x)
{
  # Load packages
  require(magrittr)
  require(ROCR)
  require(ResourceSelection)
  require(foreach)
  
  # List number of components
  ncomps <- x[[3]]
  
  # Calculate discrimination
  auc_vector <- foreach(i = 1:length(ncomps),
                        .combine = "c") %do%
    {
      # Generate AUC
      auc_obj <- prediction(x[[4]][,i] %>% inverse_logit(),
                 full_data$readmission) %>% 
        performance(measure = "auc")
      
      # Return
      auc_obj@y.values[[1]]
    }

  # Calculate calibration
  cal_vector <- foreach(i = 1:length(ncomps),
                        .combine = "c") %do%
    {
      # Perform hosmer-lemeshow test
      cal_obj <- hoslem.test(full_data$readmission,
                             x[[4]][,i] %>% inverse_logit(), 10)
      
      # Return
      cal_obj$statistic
    }
  
  # Extract distances
  distance_frame <- foreach(i = 1:length(ncomps),
                        .combine = "rbind") %do%
    {
      # Extract probability
      probs_PCA <- x[[4]][,i] %>% inverse_logit()
      
      # Measure distance
      dist_PCA <- abs(probs_PCA - probs)
      
      # Return
      data.frame(mean = mean(dist_PCA),
                 max = quantile(dist_PCA, 0.95))
    }
  
  
  # Output
  output <- data.frame(method = ncomps,
                       split = x$split, n = x$n,
                       AUC = auc_vector,
                       hoslem = cal_vector,
                       mean_distance = distance_frame$mean,
                       max_distance = distance_frame$max)
  return(output)
}

# Run extraction
psnice(value = 19)
cl <- makeCluster(ifelse(detectCores() <= 12,
                         detectCores() - 1,
                         14))
clusterExport(cl, varlist = c("full_data", "inverse_logit",
                              "probs", "results"))
metrics_PCA <- parLapply(cl, 
                             results_PCA,
                             extraction_PCA)
stopCluster(cl)

# Summarise
summary_PCA <- do.call(rbind.data.frame, metrics_PCA) %>% 
  na.omit() %>% 
  group_by(method, split) %>% 
  summarise(discrimination = mean(AUC),
            discrimination_error = sd(AUC),
            calibration = mean(hoslem),
            calibration_error = sd(hoslem),
            distance = mean(mean_distance * 100),
            mean_error = sd(mean_distance * 100),
            max = mean(max_distance * 100),
            max_error = sd(max_distance * 100)) %>% 
  mutate(method = as.factor(method))

# Plot  -> All pooling methods
summary_PCA %>% 
  pivot_longer(c(3,5,7,9)) %>% 
  mutate(error = case_when(name == "distance" ~ mean_error,
                           name == "max" ~ max_error,
                           name == "calibration" ~ calibration_error,
                           name == "discrimination" ~ discrimination_error)) %>% 
  filter(split < 0.6) %>% 
  ggplot(aes(x = split,
             y = value))+
  # geom_ribbon(aes(ymin = value - error,
  #                 ymax = value + error,
  #                 fill = method), alpha = 0.1)+
  geom_path(aes(colour = method), size = 1)+
  geom_hline(data = results %>% pivot_longer(c(3,5,7,9)),
             aes(yintercept = value), size = 1,
             linetype = "dashed")+
  theme_classic(20)+
  theme(legend.position = "top")+
  facet_wrap(~name, scales = "free_y")+
  scale_colour_brewer(palette = "Set3")+
  labs(x = "Proportion missing")

# MICE imputation----
results_MICE <- read_rds("data/impute/MICE.RDS")

# Extraction function
#x <- results_MICE[[1]]
extraction_MICE <- function(x)
{
  # Load packages
  require(magrittr)
  require(ROCR)
  require(ResourceSelection)
  require(foreach)
  
  # List methods
  methods <- x[[3]]
  
  # Calculate discrimination
  auc_vector <- foreach(i = 1:length(methods),
                        .combine = "c") %do%
    {
      # Generate AUC
      auc_obj <- prediction(x[[4]][,i] %>% inverse_logit(),
                            full_data$readmission) %>% 
        performance(measure = "auc")
      
      # Return
      auc_obj@y.values[[1]]
    }
  
  # Calculate calibration
  cal_vector <- foreach(i = 1:length(methods),
                        .combine = "c") %do%
    {
      # Perform hosmer-lemeshow test
      cal_obj <- hoslem.test(full_data$readmission,
                             x[[4]][,i] %>% inverse_logit(), 10)
      
      # Return
      cal_obj$statistic
    }
  
  # Extract distances
  distance_frame <- foreach(i = 1:length(methods),
                            .combine = "rbind") %do%
    {
      # Extract probability
      probs_MICE <- x[[4]][,i] %>% inverse_logit()
      
      # Measure distance
      dist_MICE <- abs(probs_MICE - probs)
      
      # Return
      data.frame(mean = mean(dist_MICE),
                 max = quantile(dist_MICE, 0.95))
    }
  
  
  # Output
  output <- data.frame(method = methods,
                       split = x$split, n = x$n,
                       AUC = auc_vector,
                       hoslem = cal_vector,
                       mean_distance = distance_frame$mean,
                       max_distance = distance_frame$max)
  return(output)
}

# Run extraction
psnice(value = 19)
cl <- makeCluster(ifelse(detectCores() <= 12,
                         detectCores() - 1,
                         14))
clusterExport(cl, varlist = c("full_data", "inverse_logit",
                              "probs", "results"))
metrics_MICE <- parLapply(cl, 
                         results_MICE,
                         extraction_MICE)
stopCluster(cl)

# Summarise
summary_MICE <- do.call(rbind.data.frame, metrics_MICE) %>% 
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

# Plot  -> All pooling methods
summary_MICE %>% 
  pivot_longer(c(3,5,7,9)) %>% 
  mutate(error = case_when(name == "distance" ~ mean_error,
                           name == "max" ~ max_error,
                           name == "calibration" ~ calibration_error,
                           name == "discrimination" ~ discrimination_error)) %>% 
  filter(split < 0.6) %>% 
  ggplot(aes(x = split,
             y = value))+
  # geom_ribbon(aes(ymin = value - error,
  #                 ymax = value + error,
  #                 fill = method), alpha = 0.1)+
  geom_path(aes(colour = method), size = 1)+
  geom_hline(data = results %>% pivot_longer(c(3,5,7,9)),
             aes(yintercept = value), size = 1,
             linetype = "dashed")+
  theme_classic(20)+
  theme(legend.position = "top")+
  facet_wrap(~name, scales = "free_y")+
  scale_colour_brewer(palette = "Set1")+
  labs(x = "Proportion missing")

# amelia imputation----
results_amelia <- read_rds("data/impute/amelia.RDS")

# Extraction function
#x <- results_amelia[[1]]
extraction_amelia <- function(x)
{
  # Load packages
  require(magrittr)
  require(ROCR)
  require(ResourceSelection)
  
  # Calculate discrimination
  auc_obj <- prediction(x[[3]] %>% inverse_logit(),
                        full_data$readmission) %>%
    performance(measure = "auc")
  
  # Calculate calibration
  cal_obj <- hoslem.test(full_data$readmission,
                         x[[3]] %>% inverse_logit(), 10)
  
  
  # Extract distances
  probs_amelia <- x[[3]] %>% inverse_logit()
  
  # Measure distance
  dist_amelia <- abs(probs_amelia - probs)
  
  # Return
  data.frame(mean = mean(dist_amelia),
             max = quantile(dist_amelia, 0.95))
  
  
  # Output
  output <- data.frame(
    split = x$split,
    n = x$n,
    AUC = auc_obj@y.values[[1]],
    hoslem = cal_obj$statistic,
    mean_distance = mean(dist_amelia),
    max_distance = quantile(dist_amelia, 0.95)
  )
  return(output)
}

# Run extraction
psnice(value = 19)
cl <- makeCluster(ifelse(detectCores() <= 12,
                         detectCores() - 1,
                         14))
clusterExport(cl, varlist = c("full_data", "inverse_logit",
                              "probs", "results"))
metrics_amelia <- parLapply(cl, 
                          results_amelia,
                          extraction_amelia)
stopCluster(cl)

# Summarise
summary_amelia <- do.call(rbind.data.frame, metrics_amelia) %>% 
  na.omit() %>% 
  group_by(split) %>% 
  summarise(discrimination = mean(AUC),
            discrimination_error = sd(AUC),
            calibration = mean(hoslem),
            calibration_error = sd(hoslem),
            distance = mean(mean_distance * 100),
            mean_error = sd(mean_distance * 100),
            max = mean(max_distance * 100),
            max_error = sd(max_distance * 100))

# forest imputation----
results_forest <- read_rds("data/impute/forest.RDS")

# Extraction function
#x <- results_forest[[1]]
extraction_forest <- function(x)
{
  # Load packages
  require(magrittr)
  require(ROCR)
  require(ResourceSelection)
  
  # Calculate discrimination
  auc_obj <- prediction(x[[4]] %>% inverse_logit(),
                        full_data$readmission) %>%
    performance(measure = "auc")
  
  # Calculate calibration
  cal_obj <- hoslem.test(full_data$readmission,
                         x[[4]] %>% inverse_logit(), 10)
  
  
  # Extract distances
  probs_forest <- x[[4]] %>% inverse_logit()
  
  # Measure distance
  dist_forest <- abs(probs_forest - probs)
  
  # Return
  data.frame(mean = mean(dist_forest),
             max = quantile(dist_forest, 0.95))
  
  
  # Output
  output <- data.frame(
    split = x$split,
    n = x$n,
    AUC = auc_obj@y.values[[1]],
    hoslem = cal_obj$statistic,
    mean_distance = mean(dist_forest),
    max_distance = quantile(dist_forest, 0.95)
  )
  return(output)
}

# Run extraction
psnice(value = 19)
cl <- makeCluster(ifelse(detectCores() <= 12,
                         detectCores() - 1,
                         14))
clusterExport(cl, varlist = c("full_data", "inverse_logit",
                              "probs", "results"))
metrics_forest <- parLapply(cl, 
                            results_forest,
                            extraction_forest)
stopCluster(cl)

# Summarise
summary_forest <- do.call(rbind.data.frame, metrics_forest) %>% 
  na.omit() %>% 
  group_by(split) %>% 
  summarise(discrimination = mean(AUC),
            discrimination_error = sd(AUC),
            calibration = mean(hoslem),
            calibration_error = sd(hoslem),
            distance = mean(mean_distance * 100),
            mean_error = sd(mean_distance * 100),
            max = mean(max_distance * 100),
            max_error = sd(max_distance * 100))

# KNN imputation----
results_KNN <- read_rds("data/impute/KNN.RDS")

# Extraction function
#x <- results_KNN[[1]]
extraction_KNN <- function(x)
{
  # List number of K
  n_k <- x[[3]]
  
  if(!is.na(n_k[1]))
  {
    # Extract logit matrix
    logit_matrix <- data.frame(x[[4]])
    
    # Calculate discrimination
    auc_vector <- foreach(i = 1:length(n_k),
                          .combine = "c") %do%
      {
        # Generate AUC
        auc_obj <- prediction(logit_matrix[,i] %>% inverse_logit(),
                              full_data$readmission) %>% 
          performance(measure = "auc")
        
        # Return
        auc_obj@y.values[[1]]
      }
    
    # Calculate calibration
    cal_vector <- foreach(i = 1:length(n_k),
                          .combine = "c") %do%
      {
        # Perform hosmer-lemeshow test
        cal_obj <- hoslem.test(full_data$readmission,
                               logit_matrix[,i] %>% inverse_logit(), 10)
        
        # Return
        cal_obj$statistic
      }
    
    # Extract distances
    distance_frame <- foreach(i = 1:length(n_k),
                              .combine = "rbind") %do%
      {
        # Extract probability
        probs_KNN <- logit_matrix[,i] %>% inverse_logit()
        
        # Measure distance
        dist_KNN <- abs(probs_KNN - probs)
        
        # Return
        data.frame(mean = mean(dist_KNN),
                   max = quantile(dist_KNN, 0.95))
      }
    
    
    # Output
    output <- data.frame(K = n_k,
                         split = x$split, n = x$n,
                         AUC = auc_vector,
                         hoslem = cal_vector,
                         mean_distance = distance_frame$mean,
                         max_distance = distance_frame$max)
  } else
  {
    # Output
    output <- data.frame(K = NA,
                         split = x$split, n = x$n,
                         AUC = NA,
                         hoslem = NA,
                         mean_distance = NA,
                         max_distance = NA)
  }
  
  return(output)
}

# Run extraction
parallel_setup(15)
metrics_KNN <- foreach(i = 1:length(results_KNN),
        .combine = "rbind",
        .packages = c("magrittr", "ROCR",
                      "ResourceSelection",
                      "foreach")) %dopar%
  {
    extraction_KNN(results_KNN[[i]])
  }
stopImplicitCluster()

# Summarise
summary_KNN <- metrics_KNN %>% 
  na.omit() %>% 
  group_by(K, split) %>% 
  summarise(discrimination = mean(AUC),
            discrimination_error = sd(AUC),
            calibration = mean(hoslem),
            calibration_error = sd(hoslem),
            distance = mean(mean_distance * 100),
            mean_error = sd(mean_distance * 100),
            max = mean(max_distance * 100),
            max_error = sd(max_distance * 100)) %>% 
  mutate(K = as.factor(K))

# Plot  -> All K
summary_KNN %>% 
  pivot_longer(c(3,5,7,9)) %>% 
  mutate(error = case_when(name == "distance" ~ mean_error,
                           name == "max" ~ max_error,
                           name == "calibration" ~ calibration_error,
                           name == "discrimination" ~ discrimination_error)) %>% 
  filter(split < 0.6) %>% 
  ggplot(aes(x = split,
             y = value))+
  # geom_ribbon(aes(ymin = value - error,
  #                 ymax = value + error,
  #                 fill = method), alpha = 0.1)+
  geom_path(aes(colour = K), size = 1)+
  geom_hline(data = results %>% pivot_longer(c(3,5,7,9)),
             aes(yintercept = value), size = 1,
             linetype = "dashed")+
  theme_classic(20)+
  theme(legend.position = "top")+
  facet_wrap(~name, scales = "free_y")+
  scale_colour_manual(values = pal)+
  labs(x = "Proportion missing")

# Combine and write
comparison_df <- rbind(
  summary_mean %>% filter(method != "mean"),
  summary_hotdeck %>% filter(method == "admit.chronic.age") %>% 
    mutate(method = "hotdeck"),
  summary_PCA %>% filter(method == 5) %>% 
    mutate(method = "PCA"),
  summary_MICE %>% filter(method == "pmm") %>% 
    mutate(method = "MICE"),
  summary_amelia %>% 
    mutate(method = "amelia"),
  summary_forest %>% 
    mutate(method = "RF"),
  summary_KNN %>% filter(K == 3) %>% 
    na.omit() %>% 
    mutate(method = "KNN"),
) %>%   filter(split < 0.6)

write_csv(comparison_df, "data/impute/final_comparison.csv")



