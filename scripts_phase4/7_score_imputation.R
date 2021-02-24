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

# Load functions
source("functions/apache_score_mortality.R")
source("functions/inverse_logit.R")

# Load data
results <- read_rds("data/impute_mortality/average.RDS") %>% 
  c(read_rds("data/impute_mortality/MICE.RDS"))
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
    print(i)
    apache_score(cbind(results[[i]], apache_additional))
  }
names(results_scored) <- names(results)

# Generate predictions and assess
results_final <- foreach(i = 1:length(results_scored),
        .combine = "rbind") %do%
  {
    # Extract scores
    scores_df <- results_scored[[i]]
    names(scores_df)[1] <- "apache_II"
    
    # Trim to only imputed
    # scores_df %<>%
    #   mutate(old_score = apache_additional$apache_II) %>%
    #   filter(is.na(old_score))
    
    # Build model
    # model_imputed <- glm(mortality ~ apache_II,
    #                        data = scores_df,
    #                        family = "binomial")
    
    # Predict
    probs_df <- 
      data.frame(probs = predict(model_mortality, 
                                 newdata = scores_df),
                 outcome = scores_df$mortality) %>% 
      mutate(probs = inverse_logit(probs))
    
    # Discrimination
    pred <- prediction(probs_df$probs,
                       probs_df$outcome)
    auc <- performance(pred, measure = "auc")@y.values[[1]]
    
    # Calibration
    cal <- hoslem.test(probs_df$outcome,
                probs_df$probs, g = 10)$statistic
    
    # Output
    data.frame(method = names(results_scored)[i],
               discrim = auc,
               calib = cal %>% unname())
  }

# Quick plot