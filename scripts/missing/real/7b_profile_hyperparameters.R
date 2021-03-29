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
n_boot <- 100

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

results_profile <- results_final %>% 
  mutate(class = names_split[,1],
         method = names_split[,2]) %>% 
  filter(class %in% c("MICE", "KNN",
                      "PCA", "RF",
                      "recent"))

# Plotting----------

# KNN
results_profile %>% 
  filter(class == "KNN") %>% 
  pivot_longer(c(1,3)) %>% 
  mutate(error = case_when(name == "discrimination" ~ discrim_error,
                           name == "calibration" ~ cal_error)) %>% 
  ggplot(aes(x = as.numeric(method), y = value))+
  geom_bar(stat = "identity",
           aes(fill = as.numeric(method)),
           colour = "black")+
  geom_errorbar(aes(ymin = value - error,
                    ymax = value + error),
                width = 0.25)+
  facet_wrap(~ mortality * name,
             scales = "free")+
  scale_fill_distiller(palette = "Spectral",
                       name = "KNN")+
  theme_classic(20)+
  labs(x = "", y = "")+
  theme(legend.position = "top")

# PCA
results_profile %>% 
  filter(class == "PCA") %>% 
  pivot_longer(c(1,3)) %>% 
  mutate(error = case_when(name == "discrimination" ~ discrim_error,
                           name == "calibration" ~ cal_error)) %>% 
  ggplot(aes(x = as.numeric(method), y = value))+
  geom_bar(stat = "identity",
           aes(fill = as.numeric(method)),
           colour = "black")+
  geom_errorbar(aes(ymin = value - error,
                    ymax = value + error),
                width = 0.25)+
  facet_wrap(~ mortality * name,
             scales = "free")+
  scale_fill_distiller(palette = "Spectral",
                       name = "PCA")+
  theme_classic(20)+
  labs(x = "", y = "")+
  theme(legend.position = "top")

# RF
results_profile %>% 
  filter(class == "RF") %>% 
  pivot_longer(c(1,3)) %>% 
  mutate(error = case_when(name == "discrimination" ~ discrim_error,
                           name == "calibration" ~ cal_error)) %>% 
  ggplot(aes(x = as.numeric(method), y = value))+
  geom_bar(stat = "identity",
           aes(fill = as.numeric(method)),
           colour = "black")+
  geom_errorbar(aes(ymin = value - error,
                    ymax = value + error),
                width = 0.25)+
  facet_wrap(~ mortality * name,
             scales = "free")+
  scale_fill_distiller(palette = "Spectral",
                       name = "RF")+
  theme_classic(20)+
  labs(x = "", y = "")+
  theme(legend.position = "top")


# MICE
results_profile %>% 
  filter(class == "MICE") %>% 
  pivot_longer(c(1,3)) %>% 
  mutate(error = case_when(name == "discrimination" ~ discrim_error,
                           name == "calibration" ~ cal_error)) %>% 
  ggplot(aes(x = method, y = value))+
  geom_bar(stat = "identity",
           aes(fill = method),
           colour = "black")+
  geom_errorbar(aes(ymin = value - error,
                    ymax = value + error),
                width = 0.25)+
  facet_wrap(~ mortality * name,
             scales = "free")+
  scale_fill_brewer(palette = "Set1",
                       name = "Mice")+
  theme_classic(20)+
  labs(x = "", y = "")+
  theme(legend.position = "top")

# Recent
results_profile %>% 
  filter(class == "recent") %>% 
  pivot_longer(c(1,3)) %>% 
  mutate(error = case_when(name == "discrimination" ~ discrim_error,
                           name == "calibration" ~ cal_error)) %>% 
  ggplot(aes(x = method, y = value))+
  geom_bar(stat = "identity",
           aes(fill = method),
           colour = "black")+
  geom_errorbar(aes(ymin = value - error,
                    ymax = value + error),
                width = 0.25)+
  facet_wrap(~ mortality * name,
             scales = "free")+
  scale_fill_brewer(palette = "Set1",
                    name = "Recent")+
  theme_classic(20)+
  labs(x = "", y = "")+
  theme(legend.position = "top")

