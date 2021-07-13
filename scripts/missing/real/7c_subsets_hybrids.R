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

# Load functions
source("functions/apache_score_mortality.R")
source("functions/inverse_logit.R")

# Load results
results_long <- read_rds("data/impute_discharge/recent.RDS")

results_multiple <- list(amelia_30d = read_rds("data/impute_discharge/amelia.RDS"),
                         mice_30d = read_rds("data/impute_discharge/MICE.RDS"))

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

# Multiple
cl <- makeCluster(12)
registerDoParallel(cl)
results_multiple_scored <- foreach(i = 1:length(results_multiple)) %do%
  {
    print(names(results_multiple[i]))
    
    # Middle loop for hyperparameters
    scored_imputations <- foreach(j = 1:length(results_multiple[[i]]),
                                  .packages = c("foreach", "mice",
                                                "Amelia", "dplyr")) %dopar%
    {
      df <- results_multiple[[i]][[j]]
      if(class(df) == "amelia")
      {
        print(df$m)
        # Determine additional
        if(grepl("inunit", names(results_multiple)[i]))
        {
          amelia_additional <- apache_additional_inunit
        }else
        {
          amelia_additional <- apache_additional_30d
        }
        
        # Inner loop through m
        foreach(k = 1:length(df$imputations)) %do%
          {
            apache_score(cbind(df$imputations[[k]],
                               amelia_additional))
          }
      }else
      {
        print(paste(unname(df$method[1]),
              unname(df$m), sep = "-"))
        # Determine additional
        if(grepl("inunit", names(results_multiple)[i]))
        {
          mice_additional <- apache_additional_inunit
        }else
        {
          mice_additional <- apache_additional_30d
        }
        
        # Inner loop through m
        foreach(k = 1:df$m) %do%
          {
            apache_score(cbind(mice::complete(df, k),
                               mice_additional))
          }
      }
    }
    
    # Sort names and output
    names(scored_imputations) <- names(results_multiple[[i]])
    scored_imputations
  }
names(results_multiple_scored) <- names(results_multiple)
stopCluster(cl)

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
n_boot <- 1000

# Generate predictions and assess
results_long_subset <- foreach(i = 1:length(results_long_scored),
        .combine = "rbind") %do%
  {
    print(names(results_long_scored)[i])

    # Extract scores
    scores_df <- results_long_scored[[i]]
    names(scores_df)[1] <- "apache_II"
    
    scores_df %<>%
      # Remove those still missing
      mutate(still_missing,
             old_score = apache_additional_30d$apache_II) %>% 
      filter(!still_missing) %>% 
      # Trim to only imputed
      filter(is.na(old_score))
    
    # Inner loop for bootstrap assessment
    results_boot_long <- foreach(j = 1:n_boot, .combine = "rbind",
                                   .packages = c("dplyr", "ROCR",
                                                 "ResourceSelection")) %dopar%
      {
        # Set splits
        set.seed(j)
        scores_df$ID <- 1:nrow(scores_df)
        valid_df <- scores_df %>%
          filter(!is.na(apache_II)) %>% 
          group_by(mort_30) %>% 
          slice_sample(prop = 0.1) %>% 
          ungroup()
        
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
        
        # Output
        data.frame(sample = j,
                   discrim = auc,
                   calib = cal %>% unname())
      }
    
    # Output
    results_boot_long %>%
          na.omit() %>% 
          summarise(discrimination = mean(discrim),
                    discrim_error = sd(discrim),
                    calibration = mean(calib),
                    cal_error = sd(calib)) %>% 
        mutate(method = names(results_long_scored)[i])
  }

# Second section for multiple imputations
results_mult_subset <- foreach(i = 1:length(results_multiple_scored),
                         .combine = "rbind") %do%
  {
    print(names(results_multiple_scored)[i])
    
    # Middle loop for hyperparameters
    scored_imputations <- foreach(j = 1:length(results_multiple_scored[[i]]),
                                  .combine = "rbind") %do%
      {
        df <- results_multiple_scored[[i]][[j]]

        # Amelia section
        if(grepl("amelia", names(results_multiple_scored)[i]))
        {
          print(length(df))
          # Determine additional
          if(grepl("inunit", names(results_multiple)[i]))
          {
            amelia_additional <- apache_additional_inunit 
            amelia_model <- full_model_inunit
          }else
          {
            amelia_additional <- apache_additional_30d %>% 
              filter(!still_missing)
            amelia_model <- full_model_30d
          }
          mort_name <- names(amelia_additional)[3]
          amelia_additional$row_id <- 1:nrow(amelia_additional)
          
          # Innermost loop for bootstrap assessment
          results_boot_amelia <- foreach(k = 1:n_boot, .combine = "rbind",
                                         .packages = c("dplyr", "ROCR",
                                                       "ResourceSelection",
                                                       "tibble", "foreach")) %dopar%
            {
              # Determine validation split
              set.seed(k)
              validation_id <- amelia_additional %>% 
                filter(is.na(apache_II)) %>% 
                select(all_of(c("row_id", mort_name))) %>% 
                group_by_at(mort_name) %>% 
                slice_sample(prop = 0.1) %>% 
                ungroup() %>% 
                select(row_id) %>% 
                deframe()
              
              # Innermost loop through m
              stats_m <- foreach(m = 1:length(df), .combine = "rbind") %do%
                {
                  # Select validation sample
                  valid_df <- df[[m]][validation_id,]
                  names(valid_df)[1] <- "apache_II"
                  
                  # Predict
                  probs_df <- 
                    data.frame(probs = predict(amelia_model, 
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
                         discrim = median(stats_m$discrim, na.rm = T),
                         calib = median(stats_m$calib, na.rm = T))
            }
  
          results_boot_amelia %>%
            na.omit() %>% 
            summarise(discrimination = mean(discrim),
                      discrim_error = sd(discrim),
                      calibration = mean(calib),
                      cal_error = sd(calib)) %>% 
            mutate(method = paste("amelia", length(df), sep = "_"))
        }else
        {
          print(names(results_multiple_scored[[i]])[j])
          # Determine additional
          if(grepl("inunit", names(results_multiple)[i]))
          {
            mice_additional <- apache_additional_inunit
            mice_model <- full_model_inunit
          }else
          {
            mice_additional <- apache_additional_30d %>% 
              filter(!still_missing)
            mice_model <- full_model_30d
          }
          mort_name <- names(mice_additional)[3]
          mice_additional$row_id <- 1:nrow(mice_additional)
          
          # Innermost loop for bootstrap assessment
          results_boot_mice <- foreach(k = 1:n_boot, .combine = "rbind",
                                         .packages = c("dplyr", "ROCR",
                                                       "ResourceSelection",
                                                       "tibble", "foreach")) %dopar%
            {
              # Determine validation split
              set.seed(k)
              validation_id <- mice_additional %>% 
                filter(is.na(apache_II)) %>% 
                select(all_of(c("row_id", mort_name))) %>% 
                group_by_at(mort_name) %>% 
                slice_sample(prop = 0.1) %>% 
                ungroup() %>% 
                select(row_id) %>% 
                deframe()
              
              # Innermost loop through m
              stats_m <- foreach(m = 1:length(df), .combine = "rbind") %do%
                {
                  # Select validation sample
                  valid_df <- df[[m]][validation_id,]
                  names(valid_df)[1] <- "apache_II"
                  
                  # Predict
                  probs_df <- 
                    data.frame(probs = predict(mice_model, 
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
                         discrim = median(stats_m$discrim, na.rm = T),
                         calib = median(stats_m$calib, na.rm = T))
            }
          
          results_boot_mice %>%
            na.omit() %>% 
            summarise(discrimination = mean(discrim),
                      discrim_error = sd(discrim),
                      calibration = mean(calib),
                      cal_error = sd(calib)) %>% 
            mutate(method = names(results_multiple_scored[[i]])[j])
          
        }
      }
    scored_imputations
  }

# Combine
results_subset <- results_long_subset %>% 
  mutate(type = "longitudinal") %>% 
  rbind(results_mult_subset %>% 
          mutate(type = "multiple"))

# Rescale
rescale_df_subset <- results_subset %>%
  mutate(
    discrimination = rescale(discrimination,
                             from = c(
                               max(discrimination),
                               min(discrimination))),
    calibration = rescale(calibration,
                          from = c(min(calibration),
                                   max(calibration)))) %>% 
  mutate(class = str_split(method, "_", simplify = TRUE)[,1])

# Measure distance
dist_d <- (0 - rescale_df_subset$discrimination) ^2 
dist_c <- (0 - rescale_df_subset$calibration) ^2 
rescale_df_subset$distance <- sqrt(dist_c + dist_d)

# Find best sub-method
best_methods <- rescale_df_subset %>% 
  group_by(class) %>% 
  slice_min(distance) %>% 
  ungroup()

# Store samples for best longitudinal
scores_df <- results_long_scored[["recent_worsttotal"]]
names(scores_df)[1] <- "apache_II"

scores_df %<>%
  # Remove those still missing
  mutate(still_missing,
         old_score = apache_additional_30d$apache_II) %>% 
  filter(!still_missing) %>% 
  # Trim to only imputed
  filter(is.na(old_score))

n_boot <- 10000

# Inner loop for bootstrap assessment
results_long_subset_full <- foreach(j = 1:n_boot, .combine = "rbind",
                             .packages = c("dplyr", "ROCR",
                                           "ResourceSelection")) %dopar%
  {
    # Set splits
    set.seed(j)
    scores_df$ID <- 1:nrow(scores_df)
    valid_df <- scores_df %>%
      filter(!is.na(apache_II)) %>% 
      group_by(mort_30) %>% 
      slice_sample(prop = 0.1) %>% 
      ungroup()
    
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
    
    # Output
    data.frame(sample = j,
               discrim = auc,
               calib = cal %>% unname())
  } %>% 
  mutate(type = "longitudinal")

# Make new multiple object
final_mult_list_subset <- list(amelia = list(a_20 = results_multiple_scored[["amelia_30d"]][["amelia_20"]]),
                               mice = list(m_q_10 = results_multiple_scored[["mice_30d"]][["MICE_quadratic-10"]]))

results_mult_subset_full <- foreach(i = 1:length(final_mult_list_subset),
                               .combine = "rbind") %do%
  {
    print(names(final_mult_list_subset)[i])
    
    # Middle loop for hyperparameters
    scored_imputations <- foreach(j = 1:length(final_mult_list_subset[[i]]),
                                  .combine = "rbind") %do%
      {
        df <- final_mult_list_subset[[i]][[j]]
        
        # Amelia section
        if(grepl("amelia", names(final_mult_list_subset)[i]))
        {
          print(length(df))
          # Determine additional
          if(grepl("inunit", names(results_multiple)[i]))
          {
            amelia_additional <- apache_additional_inunit 
            amelia_model <- full_model_inunit
          }else
          {
            amelia_additional <- apache_additional_30d %>% 
              filter(!still_missing)
            amelia_model <- full_model_30d
          }
          mort_name <- names(amelia_additional)[3]
          amelia_additional$row_id <- 1:nrow(amelia_additional)
          
          # Innermost loop for bootstrap assessment
          results_boot_amelia <- foreach(k = 1:n_boot, .combine = "rbind",
                                         .packages = c("dplyr", "ROCR",
                                                       "ResourceSelection",
                                                       "tibble", "foreach")) %dopar%
            {
              # Determine validation split
              set.seed(k)
              validation_id <- amelia_additional %>% 
                filter(is.na(apache_II)) %>% 
                select(all_of(c("row_id", mort_name))) %>% 
                group_by_at(mort_name) %>% 
                slice_sample(prop = 0.1) %>% 
                ungroup() %>% 
                select(row_id) %>% 
                deframe()
              
              # Innermost loop through m
              stats_m <- foreach(m = 1:length(df), .combine = "rbind") %do%
                {
                  # Select validation sample
                  valid_df <- df[[m]][validation_id,]
                  names(valid_df)[1] <- "apache_II"
                  
                  # Predict
                  probs_df <- 
                    data.frame(probs = predict(amelia_model, 
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
                         discrim = median(stats_m$discrim, na.rm = T),
                         calib = median(stats_m$calib, na.rm = T))
            }
          
          results_boot_amelia %>% 
            mutate(type = "amelia")
        }else
        {
          print(names(final_mult_list_subset[[i]])[j])
          # Determine additional
          if(grepl("inunit", names(results_multiple)[i]))
          {
            mice_additional <- apache_additional_inunit
            mice_model <- full_model_inunit
          }else
          {
            mice_additional <- apache_additional_30d %>% 
              filter(!still_missing)
            mice_model <- full_model_30d
          }
          mort_name <- names(mice_additional)[3]
          mice_additional$row_id <- 1:nrow(mice_additional)
          
          # Innermost loop for bootstrap assessment
          results_boot_mice <- foreach(k = 1:n_boot, .combine = "rbind",
                                       .packages = c("dplyr", "ROCR",
                                                     "ResourceSelection",
                                                     "tibble", "foreach")) %dopar%
            {
              # Determine validation split
              set.seed(k)
              validation_id <- mice_additional %>% 
                filter(is.na(apache_II)) %>% 
                select(all_of(c("row_id", mort_name))) %>% 
                group_by_at(mort_name) %>% 
                slice_sample(prop = 0.1) %>% 
                ungroup() %>% 
                select(row_id) %>% 
                deframe()
              
              # Innermost loop through m
              stats_m <- foreach(m = 1:length(df), .combine = "rbind") %do%
                {
                  # Select validation sample
                  valid_df <- df[[m]][validation_id,]
                  names(valid_df)[1] <- "apache_II"
                  
                  # Predict
                  probs_df <- 
                    data.frame(probs = predict(mice_model, 
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
                         discrim = median(stats_m$discrim, na.rm = T),
                         calib = median(stats_m$calib, na.rm = T))
            }
          
          results_boot_mice%>% 
            mutate(type = "mice")
        }
      }
    scored_imputations
  }

# Output
results_long_subset_full %>% 
  rbind(results_mult_subset_full) %>% 
  write_rds("models/boot_results_subset.RDS",
            compress = "gz")


# Hybrid of best mice/amelia and best recent----

# Extract best recent
recent_df <- results_long_scored[["recent_worsttotal"]]
names(recent_df)[1] <- "apache_II"

recent_df %<>%
  # Flag those still missing
  mutate(still_missing,
         old_score = apache_additional_30d$apache_II) %>% 
  rownames_to_column("row_id")

# Create mice list
mice_list <- foreach(i = 1:10, .combine = "cbind") %do%
  {
      results_multiple_scored[["mice_30d"]][["MICE_quadratic-10"]][[i]][,1]
  }

# Create mice hybrid
mice_hybrid <- cbind(recent_df, mice_list)

# Create amelia list
amelia_list <- foreach(i = 1:20, .combine = "cbind") %do%
  {
    results_multiple_scored[["amelia_30d"]][["amelia_20"]][[i]][,1]
  }

# Create amelia hybrid
amelia_hybrid <- cbind(recent_df, amelia_list)

# Bootstrap loop
n_boot <- 10000
results_boot_hybrid <- foreach(k = 1:n_boot, .combine = "rbind",
                               .packages = c("dplyr", "ROCR",
                                             "ResourceSelection",
                                             "tibble", "foreach")) %dopar%
  {
    # Determine validation split
    set.seed(k)
    validation_id <- recent_df %>% 
      filter(is.na(old_score)) %>% 
      slice_sample(prop = 0.1) %>% 
      ungroup() %>% 
      select(row_id) %>% 
      deframe()
    
    # Loop through M for mice
    mice_m <- foreach(m = 1:10, .combine = "rbind") %do%
      {
        # Select validation sample
        valid_df <- mice_hybrid[validation_id,]
        
        # Replace still missing values with imputed
        valid_df[is.na(valid_df$apache_II),2] <- valid_df[is.na(valid_df$apache_II), m + 5]
        
        # Predict
        probs_df <- 
          data.frame(probs = predict(full_model_30d, 
                                     newdata = valid_df),
                     outcome = valid_df[,3]) %>% 
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
        valid_df <- amelia_hybrid[validation_id,]
        
        # Replace still missing values with imputed
        valid_df[is.na(valid_df$apache_II),2] <- valid_df[is.na(valid_df$apache_II), m + 5]
        
        # Predict
        probs_df <- 
          data.frame(probs = predict(full_model_30d, 
                                     newdata = valid_df),
                     outcome = valid_df[,3]) %>% 
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