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
pal <- brewer.pal(9, "Set1")[-6]
require(bayestestR)

# Load functions
source("functions/apache_score_mortality.R")
source("functions/inverse_logit.R")

# Load results
results_inunit <- read_rds("data/impute_mortality/average.RDS") %>% 
  c(read_rds("data/impute_mortality/MICE.RDS"),
    read_rds("data/impute_mortality/KNN.RDS"),
    read_rds("data/impute_mortality/forest.RDS"),
    read_rds("data/impute_mortality/amelia.RDS"),
    read_rds("data/impute_mortality/PCA.RDS"))

results_30d <- read_rds("data/impute_discharge/average.RDS") %>% 
  c(read_rds("data/impute_discharge/MICE.RDS"),
    read_rds("data/impute_discharge/KNN.RDS"),
    read_rds("data/impute_discharge/forest.RDS"),
    read_rds("data/impute_discharge/amelia.RDS"),
    read_rds("data/impute_discharge/PCA.RDS"))

results_long <- read_rds("data/impute_discharge/recent.RDS")

# Load data
apache_df <- read_csv("data/apache_real_missing.csv")
apache_scores_inunit <- apache_df %>% 
  select(-row_id, -adm_id,
         -mort_inunit, -age,
         -chronic, -fractioninspiredoxygen,
         -apache_II, -missing_abg, -acute_renal_failure)
apache_additional_inunit <- apache_df %>% 
  select(row_id, adm_id,
         mort_inunit, age,
         chronic, fractioninspiredoxygen,
         apache_II, missing_abg, acute_renal_failure)
rm(apache_df)

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

# Generate scores
results_inunit_scored <- foreach(i = 1:length(results_inunit)) %do%
  {
    print(names(results_inunit)[i])
    apache_score(cbind(results_inunit[[i]], apache_additional_inunit))
  }
names(results_inunit_scored) <- names(results_inunit)

results_30d_scored <- foreach(i = 1:length(results_30d)) %do%
  {
    print(names(results_30d)[i])
    apache_score(cbind(results_30d[[i]], apache_additional_30d))
  }
names(results_30d_scored) <- names(results_30d)

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

results_30d_scored <- c(results_30d_scored, results_long_scored)

# Generate complete cases models
full_model_inunit <- glm(mort_inunit ~ apache_II,
                         data = apache_additional_inunit,
                         family = "binomial")

full_model_30d <- glm(mort_30 ~ apache_II,
                         data = apache_additional_30d,
                         family = "binomial")

# Set up parallel
n_cores <- 15
cl <- makeCluster(ifelse(detectCores() <= n_cores,
                          detectCores() - 1,
                          n_cores))
registerDoParallel(cl)
n_boot <- 10000

# Generate predictions and assess
results_final <- foreach(i = 1:length(results_30d_scored),
        .combine = "rbind") %do%
  {
    print(names(results_30d_scored)[i])
    
    # Smaller loop to allow longitudinal
    if(i <= length(results_inunit_scored))
    {
      # Extract scores
      scores_df <- results_inunit_scored[[i]]
      names(scores_df)[1] <- "apache_II"
      
      # Trim to only imputed
      scores_df %<>%
        mutate(old_score = apache_additional_inunit$apache_II) %>%
        filter(is.na(old_score))
      
      # Inner loop for bootstrap assessment
      results_boot_inunit <- foreach(j = 1:n_boot, .combine = "rbind",
                                     .packages = c("dplyr", "ROCR",
                                                   "ResourceSelection")) %dopar%
        {
          # Set splits
          set.seed(j)
          scores_df$ID <- 1:nrow(scores_df)
          valid_df <- scores_df %>% 
            filter(!is.na(apache_II)) %>% 
            group_by(mort_inunit) %>% 
            slice_sample(prop = 0.1) %>% 
            ungroup()
          
          # Predict
          probs_df <- 
            data.frame(probs = predict(full_model_inunit, 
                                       newdata = valid_df),
                       outcome = valid_df$mort_inunit) %>% 
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
    }

    # Repeat for 30-day mortality
    scores_df <- results_30d_scored[[i]]
    names(scores_df)[1] <- "apache_II"
    
    # Trim to only imputed
    scores_df %<>%
      mutate(old_score = apache_additional_30d$apache_II) %>%
      filter(is.na(old_score))
    
    # Inner loop for bootstrap assessment
    results_boot_30d <- foreach(j = 1:n_boot, .combine = "rbind",
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
    if(i <= length(results_inunit_scored))
    {
      results_boot_inunit %>%
        na.omit() %>% 
        summarise(discrimination = mean(discrim),
                  discrim_error = sd(discrim),
                  calibration = mean(calib),
                  cal_error = sd(calib)) %>% 
        mutate(mortality = "inunit",
               still_missing = sum(is.na(results_inunit_scored[[i]][[1]]))) %>% 
        rbind(
          results_boot_30d %>%
            na.omit() %>% 
            summarise(discrimination = mean(discrim),
                      discrim_error = sd(discrim),
                      calibration = mean(calib),
                      cal_error = sd(calib)) %>% 
            mutate(mortality = "30-day",
                   still_missing = sum(is.na(results_30d_scored[[i]][[1]])))
        ) %>% 
        mutate(method = names(results_inunit_scored)[i])
    }else
    {
        results_boot_30d %>%
          na.omit() %>% 
          summarise(discrimination = mean(discrim),
                    discrim_error = sd(discrim),
                    calibration = mean(calib),
                    cal_error = sd(calib)) %>% 
          mutate(mortality = "30-day",
                 still_missing = sum(is.na(results_30d_scored[[i]][[1]]))) %>% 
        mutate(method = names(results_30d_scored)[i])
    }

  }

# Construct hybrids----

# Extract best longitudinal
scores_df <- results_30d_scored[[61]]
names(scores_df)[1] <- "apache_II"

# Trim to only imputed
scores_df %<>%
  mutate(old_score = apache_additional_30d$apache_II) %>%
  filter(is.na(old_score))

# Extract best mice
scores_df_mice <- results_30d_scored[[11]]
names(scores_df_mice)[1] <- "apache_II"

# Trim to only imputed
scores_df_mice %<>%
  mutate(old_score = apache_additional_30d$apache_II) %>%
  filter(is.na(old_score))

# Find matching subsets
scores_df_mice <- scores_df_mice[which(is.na(scores_df$apache_II)),]
scores_df <- scores_df[which(!is.na(scores_df$apache_II)),]

scores_hybrid <- rbind(scores_df, scores_df_mice)

# Loop for bootstrap assessment
results_hybrid <- foreach(j = 1:n_boot, .combine = "rbind",
                            .packages = c("dplyr", "ROCR",
                                          "ResourceSelection")) %dopar%
  {
    # Set splits
    set.seed(j)
    scores_hybrid$ID <- 1:nrow(scores_hybrid)
    valid_df <- scores_hybrid %>%
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

# Process results
results_hybrid %<>% 
  na.omit() %>% 
  summarise(discrimination = mean(discrim),
            discrim_error = sd(discrim),
            calibration = mean(calib),
            cal_error = sd(calib)) %>% 
  mutate(mortality = "30-day",
         still_missing = 0) %>% 
  mutate(method = "hybrid_hybrid")

# Write and process----
results_final %>% rbind(results_hybrid)

results_final %>%
  write_rds("models/boot_results.RDS")

# Process names
names_split <- str_split(results_final$method, "_", simplify = TRUE)

# Quick plot
results_final %>% 
  mutate(class = names_split[,1],
         method = names_split[,2]) %>% 
  #filter(method != "zero") %>% 
  ggplot(aes(x = discrimination, y = calibration,
             fill = class))+
  geom_point(size = 4, shape = 21,
             alpha = 0.7)+
  theme_classic(20)+
  # geom_hline(yintercept = 27.339,
  #            linetype = "dashed")+
  # geom_vline(xintercept = 0.7191208,
  #            linetype = "dashed")+
  labs(y = "Calibration", x = "Discrimination")+
  scale_fill_brewer(palette = "Set1")+
  scale_y_reverse()+
  facet_wrap(~mortality,
             scales = "free")

# Best from each method----

# Rescale
rescale_df <- results_final %>%
  group_by(mortality) %>% 
  mutate(
    discrimination = rescale(discrimination,
                             from = c(
                               max(discrimination),
                               min(discrimination)
                             )),
    calibration = rescale(calibration,
                          from = c(min(calibration),
                                   max(calibration)))) %>% 
  ungroup() %>% 
  mutate(class = names_split[,1],
         method = names_split[,2])

# Measure distance
dist_d <- (0 - rescale_df$discrimination) ^2 
dist_c <- (0 - rescale_df$calibration) ^2 
rescale_df$distance <- sqrt(dist_c + dist_d)

# Find best sub-method
best_methods <- rescale_df %>% 
  group_by(class, mortality) %>% 
  slice_min(distance) %>% 
  ungroup() %>% 
  mutate(names = paste(class, method,
                       sep = "_")) %>% 
  select(mortality, names)

best_methods %<>% 
  rbind(data.frame(mortality = c("30-day", "inunit"),
                   names = "average_zero"))

# Store bootstrap samples for best method
bootstrap_samples <- foreach(i = 1:nrow(best_methods),
                         .combine = "rbind") %do%
  {
    print(best_methods[i,])
    mtype <- best_methods[i,1] %>% deframe
    
    # Extract scores
    if(mtype == "30-day")
    {
      scores_df <- results_30d_scored[[best_methods[i,2] %>% deframe()]] 
    }else
    {
      scores_df <- results_inunit_scored[[best_methods[i,2] %>% deframe()]] 
    }
    names(scores_df)[1:2] <- c("apache_II", "mortality")
    
    
    # Trim to only imputed
    if(mtype == "30-day")
    {
      scores_df %<>%
      mutate(old_score = apache_additional_30d$apache_II) %>%
      filter(is.na(old_score))
    }else
    {
      scores_df %<>%
      mutate(old_score = apache_additional_inunit$apache_II) %>%
      filter(is.na(old_score))
    }
    
    # Inner loop for bootstrap assessment
    results_boot <- foreach(j = 1:n_boot, .combine = "rbind",
                                   .packages = c("dplyr", "ROCR",
                                                 "ResourceSelection")) %dopar%
      {
        # Set splits
        scores_df$ID <- 1:nrow(scores_df)
        valid_df <- scores_df %>% 
          filter(!is.na(apache_II)) %>% 
          group_by(mortality) %>% 
          slice_sample(prop = 0.1) %>% 
          ungroup()
        
        # Predict
        if(mtype == "30-day")
        {
          probs_df <- 
            data.frame(probs = predict(full_model_30d, 
                                       newdata = valid_df),
                       outcome = valid_df$mortality) %>% 
            mutate(probs = inverse_logit(probs)) %>% 
            na.omit()
        }else
        {
          probs_df <- 
            data.frame(probs = predict(full_model_inunit, 
                                       newdata = valid_df),
                       outcome = valid_df$mortality) %>% 
            mutate(probs = inverse_logit(probs)) %>% 
            na.omit()
        }

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
      mutate(method = best_methods[i,2] %>% deframe(),
             mortality = mtype)
  }

# Post-process----

bootstrap_samples$method[bootstrap_samples$method == "average_zero"] <- "zero_zero"
bootstrap_samples %<>% 
  mutate(method = str_split(method, "_", simplify = TRUE)[,1])


# Extract means and error
names(best_methods)[2] <- "method"
means_df <- best_methods %>% 
  left_join(results_final) %>% 
  mutate(method = str_split(method, "_", simplify = TRUE)[,1])
means_df[14:15,2] <- "zero"

# Save
bootstrap_samples  %>%
  write_rds("models/boot_samples.RDS",
            compress = "gz")

# Plot
bootstrap_samples %>% 
  ggplot()+
  theme_classic(20)+
  labs(y = "Calibration", x = "Discrimination")+
  scale_colour_manual(values = pal,
                      name = "")+
  scale_fill_manual(values = pal,
                    name = "")+
  scale_y_reverse()+
  theme(legend.position = "top")+
  stat_ellipse(aes(x = discrim, y = calib,
                   colour = method),
               size = 2, type = "norm",
               level = 0.682)+
  # geom_errorbar(data = means_df,
  #               aes(ymin = calibration - cal_error,
  #                   y = calibration,
  #                   ymax = calibration + cal_error,
  #                   x = discrimination,
  #                   colour = method))+
  # geom_errorbarh(data = means_df,
  #               aes(xmin = discrimination - discrim_error,
  #                   y = calibration,
  #                   xmax = discrimination + discrim_error,
  #                   colour = method))+
  geom_point(data = means_df,
             aes(x = discrimination,
                 y = calibration,
                 fill = method),
             size = 4, shape = 21)+
  facet_wrap(~mortality)

# HDI
bootstrap_samples %>% 
  group_by(mortality, method) %>% 
  summarise(discrimination = hdi(discrim)[,2:3],
            calibration = hdi(calib)[,2:3])
stopCluster(cl)

# Compare subsets----

# Extract best longitudinal
scores_df <- results_30d_scored[[61]]
names(scores_df)[1] <- "apache_II"

# Trim to only imputed
scores_df %<>%
  mutate(old_score = apache_additional_30d$apache_II) %>%
  filter(is.na(old_score))

# Extract best mice
scores_df_mice <- results_30d_scored[[11]]
names(scores_df_mice)[1] <- "apache_II"

# Trim to only imputed
scores_df_mice %<>%
  mutate(old_score = apache_additional_30d$apache_II) %>%
  filter(is.na(old_score))

# Find matching subsets
scores_df_mice <- scores_df_mice[which(!is.na(scores_df$apache_II)),]
scores_df <- scores_df[which(!is.na(scores_df$apache_II)),]

# Loop for bootstrap assessment
results_boot_subset <- foreach(j = 1:1000, .combine = "rbind",
                        .packages = c("dplyr", "ROCR",
                                      "ResourceSelection")) %do%
  {
    print(j)
    set.seed(j)
    
    # Set splits
    scores_df$ID <- 1:nrow(scores_df)
    scores_df_mice$ID <- 1:nrow(scores_df_mice)
    
    valid_df_long <- scores_df %>% 
      filter(!is.na(apache_II)) %>% 
      group_by(mort_30) %>% 
      slice_sample(prop = 0.1) %>% 
      ungroup()
    
    valid_df_mice <- scores_df_mice %>% 
      filter(ID %in% valid_df_long$ID)
    
    # Predict
      probs_df_long <- 
        data.frame(probs = predict(full_model_30d, 
                                   newdata = valid_df_long),
                   outcome = valid_df_long$mort_30) %>% 
        mutate(probs = inverse_logit(probs)) %>% 
        na.omit()
      
      probs_df_mice <- 
        data.frame(probs = predict(full_model_30d, 
                                   newdata = valid_df_mice),
                   outcome = valid_df_mice$mort_30) %>% 
        mutate(probs = inverse_logit(probs)) %>% 
        na.omit()
    
    # Discrimination
    pred_long <- prediction(probs_df_long$probs,
                            probs_df_long$outcome)
    auc_long <- performance(pred_long, measure = "auc")@y.values[[1]]
    
    pred_mice <- prediction(probs_df_mice$probs,
                            probs_df_mice$outcome)
    auc_mice <- performance(pred_mice, measure = "auc")@y.values[[1]]
    
    # Calibration
    cal_long <- hoslem.test(probs_df_long$outcome,
                            probs_df_long$probs, g = 10)$statistic
    
    cal_mice <- hoslem.test(probs_df_mice$outcome,
                            probs_df_mice$probs, g = 10)$statistic
    
    if(is.na(cal_mice) | is.na(cal_long))
    {
      cal_long <- hoslem.test(probs_df_long$outcome,
                              probs_df_long$probs, g = 9)$statistic
      
      cal_mice <- hoslem.test(probs_df_mice$outcome,
                              probs_df_mice$probs, g = 9)$statistic
    }
    
    # Output
    data.frame(sample = j,
               discrim_long = auc_long,
               discrim_mice = auc_mice,
               calib_long = cal_long %>% unname(),
               calib_mice = cal_mice %>% unname())
  }


subset_plot <- results_boot_subset %>% 
  pivot_longer(2:5) %>% 
  mutate(measure = ifelse(grepl("discrim", name), "discrim", "calib"),
         method = ifelse(grepl("long", name), "longitudinal", "MICE")) %>% 
  pivot_wider(names_from = "measure",
             values_from = "value") %>% 
  group_by(sample, method) %>% 
  summarise(discrim = mean(discrim, na.rm = T),
            calib = mean(calib, na.rm = T)) %>% 
  ungroup()

subset_means <- subset_plot %>% 
  group_by(method) %>% 
  summarise(discrimination = mean(discrim, na.rm = T),
            calibration = mean(calib, na.rm = T))

subset_plot %>% 
  ggplot()+
  theme_classic(20)+
  labs(y = "Calibration", x = "Discrimination")+
  scale_colour_manual(values = pal,
                      name = "")+
  scale_fill_manual(values = pal,
                    name = "")+
  scale_y_reverse()+
  theme(legend.position = "top")+
  stat_ellipse(aes(x = discrim, y = calib,
                   colour = method),
               size = 2, type = "norm",
               level = 0.682)+
geom_point(data = subset_means,
           aes(x = discrimination,
               y = calibration,
               fill = method),
           size = 4, shape = 21)


  