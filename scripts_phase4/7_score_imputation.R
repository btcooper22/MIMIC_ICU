# Setup----

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
require(ResourceSelection)
require(ggplot2)
require(scales)
require(car)

# Load functions
source("functions/apache_score_mortality.R")
source("functions/inverse_logit.R")

# Load data
results <- read_rds("data/impute_mortality/average.RDS") %>% 
  c(read_rds("data/impute_mortality/MICE.RDS"),
    read_rds("data/impute_mortality/KNN.RDS"),
    read_rds("data/impute_mortality/forest.RDS"))
source("functions/data_loader_splitter.R")

# Complete cases----

# Build complete cases model
model_mortality <- glm(mort_inunit ~ apache_II,
                       data = apache_additional %>% filter(!is.na(apache_II)),
                       family = "binomial")

# Assess complete cases model
probs_complete <- 
data.frame(probs = predict(model_mortality, 
                           newdata = model_mortality$data),
           outcome = model_mortality$data$mort_inunit) %>% 
  mutate(probs = inverse_logit(probs))

# Discrimination
pred_complete <- prediction(probs_complete$probs,
                            probs_complete$outcome)
performance(pred_complete, measure = "auc")@y.values[[1]]

# Calibration
hoslem.test(probs_complete$outcome,
            probs_complete$probs, g = 10)

# Imputed datasets----

# Generate scores
results_scored <- foreach(i = 1:length(results)) %do%
  {
    print(names(results)[i])
    apache_score(cbind(results[[i]], apache_additional))
  }
names(results_scored) <- names(results)

# Set up parallel
n_cores <- 15
cl <- makeCluster(ifelse(detectCores() <= n_cores,
                          detectCores() - 1,
                          n_cores))
registerDoParallel(cl)
n_boot <- 10000

# Generate predictions and assess
results_final <- foreach(i = 1:length(results_scored),
        .combine = "rbind") %do%
  {
    print(names(results_scored)[i])
    
    # Extract scores
    scores_df <- results_scored[[i]]
    names(scores_df)[1] <- "apache_II"
    
    # Trim to only imputed
    scores_df %<>%
      mutate(old_score = apache_additional$apache_II) %>%
      filter(is.na(old_score))
    
    # Inner loop for bootstrap assessment
    results_boot <- foreach(j = 1:n_boot, .combine = "rbind",
            .packages = c("dplyr", "ROCR",
                          "ResourceSelection")) %dopar%
      {
        # Set splits
        scores_df$ID <- 1:nrow(scores_df)
        train_df <- scores_df %>% 
          group_by(mortality) %>% 
          slice_sample(prop = 0.75) %>% 
          ungroup()
        
        valid_df <- scores_df %>% 
          filter(ID %in% train_df$ID == FALSE)
        
        # Build model
        model_imputed <- glm(mortality ~ apache_II,
                             data = train_df,
                             family = "binomial")
        
        # Predict
        probs_df <- 
          data.frame(probs = predict(model_mortality, 
                                     newdata = valid_df),
                     outcome = valid_df$mortality) %>% 
          mutate(probs = inverse_logit(probs))
        
        # Discrimination
        pred <- prediction(probs_df$probs,
                           probs_df$outcome)
        auc <- performance(pred, measure = "auc")@y.values[[1]]
        
        # Calibration
        cal <- hoslem.test(probs_df$outcome,
                           probs_df$probs, g = 10)$statistic
        
        # Output
        data.frame(sample = j,
                   discrim = auc,
                   calib = cal %>% unname())
      }
    
    # Output
    results_boot %>%
      na.omit() %>% 
      summarise(discrimination = mean(discrim),
                discrim_error = sd(discrim),
                calibration = mean(calib),
                cal_error = sd(calib)) %>% 
      mutate(method = names(results_scored)[i])
  }
stopCluster(cl)
results_final

# Process names
names_split <- str_split(results_final$method, "_", simplify = TRUE)

# Quick plot
results_final %>% 
  mutate(class = names_split[,1],
         method = names_split[,2]) %>% 
  filter(method != "zero") %>% 
  ggplot(aes(x = discrimination, y = calibration,
             fill = class))+
  geom_point(size = 4, shape = 21)+
  theme_classic(20)+
  # geom_hline(yintercept = 27.339,
  #            linetype = "dashed")+
  # geom_vline(xintercept = 0.7191208,
  #            linetype = "dashed")+
  labs(y = "Calibration", x = "Discrimination")+
  scale_fill_brewer(palette = "Set1")+
  scale_y_reverse()

# Best from each method----

# Rescale
rescale_df <- results_final %>%
  mutate(
    discrimination = rescale(discrimination,
                             from = c(
                               max(discrimination),
                               min(discrimination)
                             )),
    calibration = rescale(calibration,
                          from = c(min(calibration),
                                   max(calibration))),
    class = names_split[,1],
    method = names_split[,2]) 

# Measure distance
dist_d <- (0 - rescale_df$discrimination) ^2 
dist_c <- (0 - rescale_df$calibration) ^2 
rescale_df$distance <- sqrt(dist_c + dist_d)

# Find best sub-method
best_methods <- rescale_df %>% 
  group_by(class) %>% 
  slice_min(distance) %>% 
  ungroup() %>% 
  mutate(names = paste(class, method,
                       sep = "_")) %>% 
  select(names) %>% deframe()

# Set up parallel
cl <- makeCluster(ifelse(detectCores() <= n_cores,
                         detectCores() - 1,
                         n_cores))
registerDoParallel(cl)

# Store bootstrap samples for best method
bootstrap_samples <- foreach(i = 1:length(best_methods),
                         .combine = "rbind") %do%
  {
    print(best_methods[i])
    
    # Extract scores
    scores_df <- results_scored[[best_methods[i]]]
    names(scores_df)[1] <- "apache_II"
    
    # Trim to only imputed
    scores_df %<>%
      mutate(old_score = apache_additional$apache_II) %>%
      filter(is.na(old_score))
    
    # Inner loop for bootstrap assessment
    results_boot <- foreach(j = 1:n_boot, .combine = "rbind",
                            .packages = c("dplyr", "ROCR",
                                          "ResourceSelection")) %dopar%
      {
        # Set splits
        scores_df$ID <- 1:nrow(scores_df)
        train_df <- scores_df %>% 
          group_by(mortality) %>% 
          slice_sample(prop = 0.75) %>% 
          ungroup()
        
        valid_df <- scores_df %>% 
          filter(ID %in% train_df$ID == FALSE)
        
        # Build model
        model_imputed <- glm(mortality ~ apache_II,
                             data = train_df,
                             family = "binomial")
        
        # Predict
        probs_df <- 
          data.frame(probs = predict(model_mortality, 
                                     newdata = valid_df),
                     outcome = valid_df$mortality) %>% 
          mutate(probs = inverse_logit(probs))
        
        # Discrimination
        pred <- prediction(probs_df$probs,
                           probs_df$outcome)
        auc <- performance(pred, measure = "auc")@y.values[[1]]
        
        # Calibration
        cal <- hoslem.test(probs_df$outcome,
                           probs_df$probs, g = 10)$statistic
        
        # Output
        data.frame(sample = j,
                   discrim = auc,
                   calib = cal %>% unname())
      }
    
    # Output
    results_boot %>%
      na.omit() %>% 
      mutate(method = best_methods[i])
  }
stopCluster(cl)

# Extract means and error
means_df <- results_final %>% 
  filter(method %in% best_methods)

# Plot
bootstrap_samples %>% 
  ggplot()+
  theme_classic(20)+
  labs(y = "Calibration", x = "Discrimination")+
  scale_colour_brewer(palette = "Set1",
                      name = "")+
  scale_fill_brewer(palette = "Set1",
                    name = "")+
  scale_y_reverse()+
  theme(legend.position = "top")+
  stat_ellipse(aes(x = discrim, y = calib,
                   colour = method), 
               size = 1, type = "norm")+
  geom_errorbar(data = means_df,
                aes(ymin = calibration - cal_error,
                    y = calibration,
                    ymax = calibration + cal_error,
                    x = discrimination,
                    colour = method))+
  geom_errorbarh(data = means_df,
                aes(xmin = discrimination - discrim_error,
                    y = calibration,
                    xmax = discrimination + discrim_error,
                    colour = method))+
  geom_point(data = means_df,
             aes(x = discrimination,
                 y = calibration,
                 fill = method),
             size = 4, shape = 21)

# Distributions
bootstrap_samples %>% 
  group_by(method) %>% 
  pivot_longer(2:3) %>% 
  ggplot(aes(x = value,
             fill = method))+
  theme_classic(20)+
  geom_density()+
  scale_fill_brewer(palette = "Set1")+
  theme(legend.position = "none")+
  facet_wrap(~method * name,
             scales = "free")
  