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
require(RColorBrewer)
pal <- c(brewer.pal(9, "Set1")[-6], "black")
require(bayestestR)
require(mice)
require(Amelia)

# Load functions
source("functions/apache_score_mortality.R")
source("functions/inverse_logit.R")

# Load results
results_long <- read_rds("data/impute_discharge/recent.RDS")

# Load data
apache_df <- read_csv("data/apache_discharge_missing.csv")
apache_scores_30d <- apache_df %>% 
  select(-row_id, -adm_id,
         -mort_30, -age,
         -chronic, -fractioninspiredoxygen,
         -apache_II, -missing_abg, -acute_renal_failure)
apache_additional_30d <- apache_df %>% 
  select(row_id, adm_id,
         mort_30, age,
         chronic, fractioninspiredoxygen,
         apache_II, missing_abg, acute_renal_failure)
rm(apache_df)

# Imputed datasets----

# Longitudinal
results_long_scored <- foreach(i = 1:length(results_long)) %do%
  {
    print(names(results_long)[i])
    
    # Build new missing_abg
    new_missing_abg <- is.na(results_long[[i]][["arterialpH"]]) &
             (is.na(results_long[[i]][["arterialoxygen"]]) |
                is.na(results_long[[i]][["arterialcarbon"]]))
    

    apache_additional_long <- apache_additional_30d %>% 
      mutate(missing_abg = new_missing_abg)
    
    apache_score(cbind(results_long[[i]], apache_additional_long))
  }
names(results_long_scored) <- names(results_long)

# Generate complete cases models
full_model_30d <- glm(mort_30 ~ apache_II,
                         data = apache_additional_30d,
                         family = "binomial")

# Find those still missing recent data
still_missing <- is.na(results_long_scored$recent_recent$apache_scores)
sum(still_missing)
round(mean(still_missing) * 100, 2)

# Set up parallel
n_cores <- 15
cl <- makeCluster(ifelse(detectCores() <= n_cores,
                         detectCores() - 1,
                         n_cores))
registerDoParallel(cl)

# Compare recent and mice/amelia in subset----
n_boot <- 100

# Load recent data
scores_df <- read_rds("data/impute_discharge/subset_recent.RDS")
names(scores_df)[1] <- "apache_II"

# Load mice subset
mice_subset <- read_rds("data/impute_discharge/subset_mice.RDS")

# Load amelia subset
amelia_subset <- read_rds("data/impute_discharge/subset_amelia.RDS")

# Bootstrap loop
results_boot_subset <- foreach(k = 1:n_boot, .combine = "rbind",
                               .packages = c("dplyr", "ROCR",
                                             "ResourceSelection",
                                             "tibble", "foreach",
                                             "mice", "Amelia")) %dopar%
  {
    # Determine validation split
    set.seed(k)
    validation_id <- scores_df %>% 
      rownames_to_column("row_id") %>% 
      slice_sample(prop = 0.1) %>% 
      ungroup() %>% 
      select(row_id) %>% 
      deframe() %>% 
      as.numeric()

    # Loop through M for mice
    mice_m <- foreach(m = 1:10, .combine = "rbind") %do%
      {
        # Select validation sample
        mice_complete <- complete(mice_subset$results, m)
        valid_df <- mice_complete[is.na(mice_subset$additional$apache_II),][validation_id,]
        
        # Score
        valid_df <- apache_score(cbind(valid_df, mice_subset$additional[validation_id,]))
        names(valid_df)[1] <- "apache_II"
        
        # Predict
        probs_df <- 
          data.frame(probs = predict(full_model_30d, 
                                     newdata = valid_df),
                     outcome = valid_df[,2]) %>% 
          mutate(probs = inverse_logit(probs)) %>% 
          na.omit()
        
        # Discrimination
        pred <- prediction(probs_df$probs,
                           probs_df[,2])
        auc <- performance(pred, measure = "auc")@y.values[[1]]
        
        # Calibration
        cal <- hoslem.test(probs_df[,2],
                           probs_df$probs, g = 10)$statistic
        
        # Output
        data.frame(imputation = m,
                   discrim = auc,
                   calib = cal %>% unname())
      }
    
    # Loop through M for amelia
    amelia_m <- foreach(m = 1:20, .combine = "rbind") %do%
      {
        # Select validation sample
        amelia_complete <- amelia_subset$results$imputations[[m]]
        valid_df <- amelia_complete[is.na(amelia_subset$additional$apache_II),][validation_id,]
        
        # Score
        valid_df <- apache_score(cbind(valid_df, amelia_subset$additional[validation_id,]))
        names(valid_df)[1] <- "apache_II"
        
        # Predict
        probs_df <- 
          data.frame(probs = predict(full_model_30d, 
                                     newdata = valid_df),
                     outcome = valid_df[,2]) %>% 
          mutate(probs = inverse_logit(probs)) %>% 
          na.omit()
        
        # Discrimination
        pred <- prediction(probs_df$probs,
                           probs_df[,2])
        auc <- performance(pred, measure = "auc")@y.values[[1]]
        
        # Calibration
        cal <- hoslem.test(probs_df[,2],
                           probs_df$probs, g = 10)$statistic
        
        # Output
        data.frame(imputation = m,
                   discrim = auc,
                   calib = cal %>% unname())
      }
    
    # Isolate longitudinal scores
    valid_df <- scores_df[validation_id,]
    
    # Predict
    probs_df <- 
      data.frame(probs = predict(full_model_30d, 
                                 newdata = valid_df),
                 outcome = valid_df$mort_30) %>% 
      mutate(probs = inverse_logit(probs)) %>% 
      na.omit()
    
    # Discrimination
    pred <- prediction(probs_df$probs,
                       probs_df$outcome)
    auc <- performance(pred, measure = "auc")@y.values[[1]]
    
    # Calibration
    cal <- hoslem.test(probs_df$outcome,
                       probs_df$probs, g = 10)$statistic
    
    # Take median AUC & cal over all M
    data.frame(sample = k,
               discrim_mice = median(mice_m$discrim, na.rm = T),
               calib_mice = median(mice_m$calib, na.rm = T),
               discrim_amelia = median(amelia_m$discrim, na.rm = T),
               calib_amelia = median(amelia_m$calib, na.rm = T),
               discrim_long = auc,
               calib_long = cal)
  }

# Output
results_boot_subset %>% 
  write_rds("models/boot_results_subset.RDS",
            compress = "gz")


# Hybrid of best mice/amelia and best recent----

# Load best longitudinal
recent_df <- results_long_scored[["recent_worsttotal"]]
names(recent_df)[1] <- "apache_II"

recent_df %<>%
  # Flag those still missing
  mutate(still_missing,
         old_score = apache_additional_30d$apache_II) %>% 
  rownames_to_column("row_id")

# Load mice hybrid
mice_hybrid <- read_rds("data/impute_discharge/hybrid_mice.RDS")

# Load amelia hybrid
amelia_hybrid <- read_rds("data/impute_discharge/hybrid_amelia.RDS")

# Bootstrap loop
results_boot_hybrid <- foreach(k = 1:n_boot, .combine = "rbind",
                               .packages = c("dplyr", "ROCR",
                                             "ResourceSelection",
                                             "tibble", "foreach",
                                             "mice", "Amelia")) %dopar%
  {
    # Determine validation split
    set.seed(k)
    validation_id <- recent_df %>% 
      filter(is.na(old_score)) %>% 
      slice_sample(prop = 0.1) %>% 
      ungroup() %>% 
      select(row_id) %>% 
      deframe() %>% 
      as.numeric()
    
    # Loop through M for mice
    mice_m <- foreach(m = 1:10, .combine = "rbind") %do%
      {
        # Select validation sample
        mice_complete <- complete(mice_hybrid, m)
        valid_df <- mice_complete[validation_id,]
        
        # Score
        valid_df <- apache_score(cbind(valid_df, apache_additional_30d[validation_id,]))
        names(valid_df)[1] <- "apache_II"
        
        # Predict
        probs_df <- 
          data.frame(probs = predict(full_model_30d, 
                                     newdata = valid_df),
                     outcome = valid_df[,2]) %>% 
          mutate(probs = inverse_logit(probs)) %>% 
          na.omit()
        
        # Discrimination
        pred <- prediction(probs_df$probs,
                           probs_df[,2])
        auc <- performance(pred, measure = "auc")@y.values[[1]]
        
        # Calibration
        cal <- hoslem.test(probs_df[,2],
                           probs_df$probs, g = 10)$statistic
        
        # Output
        data.frame(imputation = m,
                   discrim = auc,
                   calib = cal %>% unname())
      }
    
    # Loop through M for amelia
    amelia_m <- foreach(m = 1:20, .combine = "rbind") %do%
      {
        # Select validation sample
        amelia_complete <- amelia_hybrid$imputations[[m]]
        valid_df <- amelia_complete[validation_id,]
        
        # Score
        valid_df <- apache_score(cbind(valid_df, apache_additional_30d[validation_id,]))
        names(valid_df)[1] <- "apache_II"
        
        # Predict
        probs_df <- 
          data.frame(probs = predict(full_model_30d, 
                                     newdata = valid_df),
                     outcome = valid_df[,2]) %>% 
          mutate(probs = inverse_logit(probs)) %>% 
          na.omit()
        
        # Discrimination
        pred <- prediction(probs_df$probs,
                           probs_df[,2])
        auc <- performance(pred, measure = "auc")@y.values[[1]]
        
        # Calibration
        cal <- hoslem.test(probs_df[,2],
                           probs_df$probs, g = 10)$statistic
        
        # Output
        data.frame(imputation = m,
                   discrim = auc,
                   calib = cal %>% unname())
      }
    
    # Take median AUC & cal over all M
    data.frame(sample = k,
               discrim_mice = median(mice_m$discrim, na.rm = T),
               calib_mice = median(mice_m$calib, na.rm = T),
               discrim_amelia = median(amelia_m$discrim, na.rm = T),
               calib_amelia = median(amelia_m$calib, na.rm = T))
  }
  
# Write results
results_boot_hybrid %>% 
  write_rds("models/boot_results_hybrid.RDS",
            compress = "gz")