# Packages
require(tidyverse)
require(readr)
require(magrittr)
require(lubridate)
require(doParallel)
require(tools)
require(ggplot2)
require(tidyr)
require(lmerTest)
require(Metrics)
require(equatiomatic)

# Read data
df <- read_csv("data/apache_data_full.csv") %>% 
  filter(type == "discharge",
         !is.na(sodium) &
           !is.na(potassium))

# Plot
df %>% 
  ggplot(aes(x = potassium,
             y = sodium))+
  geom_point()+
  geom_smooth(method = "lm")

# Model
glm(sodium ~ potassium,
    data = df) %>% summary()

# Set up parallel
psnice(value = 19)
registerDoParallel(ifelse(detectCores() <= 12,
                          detectCores() - 1,
                          12)
)


# Loop through charts and outputs
results <- foreach(i = 1:nrow(df), .combine = "rbind",
                   .packages = c("foreach", "readr",
                                 "dplyr", "lubridate",
                                 "tibble", "tidyr")) %dopar%
  {
    # Identify patient
    adm <- df$adm_id[i]
    
    # Load chart and output data
    labs_name <- paste("data/events/labevents_",
                       adm, ".csv", sep = "")
    
    if(file.exists(labs_name))
    {
      # Load data
      labs_df <- read_csv(labs_name)
      
      # Gather all serum sodium measurements
      sodium_measurements <- labs_df %>% 
        filter(itemid %in% c(50824, 50983)) %>% 
        select(charttime, valuenum) %>% 
        mutate(type = "sodium")
      
      # Gather all serum potassium measurements
      potassium_measurements <- labs_df %>% 
        filter(itemid %in% c(50822, 50971)) %>% 
        select(charttime, valuenum) %>% 
        mutate(type = "potassium")
      
      # Gather all serum chloride measurements
      chloride_measurements <- labs_df %>% 
        filter(itemid %in% c(50806, 50902)) %>% 
        select(charttime, valuenum) %>% 
        mutate(type = "chloride")
      
      # Combine, filter and pivot
      combined_df <- rbind(sodium_measurements,
                           potassium_measurements,
                           chloride_measurements) %>% 
        group_by(charttime) %>% 
        mutate(freq = n()) %>% 
        ungroup() %>% 
        filter(freq == 3) %>% 
        select(-freq) %>% 
        pivot_wider(names_from = "type",
                    values_from = "valuenum") %>% 
        select(-charttime)
      
      # Prepare full output
      combined_df$id <- adm
      combined_df
    }else
    {
      combined_df <- data.frame(sodium  = NA,
                                potassium  = NA,
                                chloride = NA,
                                id <- adm)
      combined_df
    }
  }
stopImplicitCluster()

# Process results
results$id %<>% as.factor()

# Create model
chloride_model <- lm(chloride ~ sodium * potassium,
                     data = results)
summary(chloride_model)

# Create mixed-effects model
chloride_model_lmer <- lmer(chloride ~ sodium * potassium + (1|id),
                            data = results)
summary(chloride_model_lmer)

coef(chloride_model)
chloride_coef <- coef(chloride_model_lmer)$id %>% colMeans()
chloride_coef

# Test and write
chloride_coef[1] + (chloride_coef[2] * 143) + 
  (chloride_coef[3] * 3.2) + (chloride_coef[4] * (143 * 3.2))

predict(chloride_model_lmer, newdata = data.frame(sodium = 143,
                                                  potassium = 3.2,
                                                  id = 165660))
write_rds(chloride_coef, "models/serum_chloride.RDS")

extract_eq(chloride_model)

# Cross-validate-----
# Create splits
splits_df <- data.frame(patient_id = unique(results$id),
                        split = floor(seq(0,10,length.out = length(unique(results$id)))) + 1)
splits_df$split[splits_df$split == 11] <- 10

# Set up parallel
cl <- makeCluster(10)
registerDoParallel(cl)

# Loop
xval <- foreach(i = 1:10, .combine = "rbind",
        .packages = c("dplyr", "magrittr", "lme4", "Metrics")) %dopar%
  {
    # Isolate validation patients
    validation_patients <- splits_df %>% 
      filter(split == i)
    
    # Isolate training patients
    training_patients <- splits_df %>% 
      filter(split != i)
    
    # Extract training values
    training_df <- results %>% 
      filter(id %in% training_patients$patient_id)
    
    # Extract validation values
    validation_df <- results %>% 
      filter(id %in% validation_patients$patient_id)
    
    # Train model
    chloride_model_train <- lmer(chloride ~ sodium * potassium + (1|id),
                                data = training_df)

    # Extract coefficients
    chloride_coef <- coef(chloride_model_lmer)$id %>% colMeans()
    
    # Predict from validation set
    validation_df %<>% 
      mutate(chloride_estimate = chloride_coef[1] + (chloride_coef[2] * sodium) + 
               (chloride_coef[3] * potassium) + (chloride_coef[4] * (sodium * potassium))) %>% 
      na.omit()
    
    # Calculate RMSE
    RMSE <- rmse(validation_df$chloride,validation_df$chloride_estimate)
    
    # Output
    data.frame(i, RMSE)
  }
stopCluster(cl)

mean(xval$RMSE)
sd(xval$RMSE)
