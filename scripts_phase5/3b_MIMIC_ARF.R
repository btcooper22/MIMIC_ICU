# Packages
require(tidyverse)
require(readr)
require(magrittr)
require(lubridate)
require(doParallel)
require(tools)
require(ggplot2)
require(caret)
require(e1071)

# Read data
df <- read_csv("data/apache_data_full.csv") %>% 
  filter(type == "discharge",
         !is.na(creatinine))

# Plot
df %>% 
  ggplot(aes(x = acute_renal_failure,
             y = creatinine))+
  geom_boxplot()

# Model
glm(acute_renal_failure ~ creatinine,
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
                                 "tibble")) %dopar%
  {
    # Identify patient
    adm <- df$adm_id[i]
    
    # Load chart and output data
    labs_name <- paste("data/events/labevents_",
                       adm, ".csv", sep = "")
    output_name <- paste("data/events/outputevents_",
                         adm, ".csv", sep = "") 
    
    if(file.exists(labs_name) &
       file.exists(output_name))
    {
    # Load data
    labs_df <- read_csv(labs_name)
    output_df <- read_csv(output_name)
    
    # Identify creatinine measurements
    creatinine <- labs_df %>% 
      filter(itemid %in% c(50912))
    
    # Identify urine output
    urine <- output_df %>% 
      # Carevue urine items
      filter(itemid %in% c(40055, 43175, 43175, 40069,
                           40094, 40715, 40473, 40085,
                           40057, 40056, 40405, 40428,
                           40086, 40096, 40651,
                           # Metavision urine items 
                           226559, 226560, 226561,
                           226584, 226563, 226564,
                           226565, 226567, 226557,
                           226558, 227488, 227489))
    
    # Loop through creatinine measurements
    output <- foreach(j = 1:nrow(creatinine),
            .combine = "rbind") %do%
      {
        # Extract measurement and time
        cre_val <- creatinine$valuenum[j]
        cre_time <- creatinine$charttime[j]
        
        # Create window for urine output
        ur_val <- urine %>% 
          filter(charttime < cre_time &
                   charttime > (cre_time - 86400)) %>% 
          # Sum urine output
          summarise(urine = sum(value, na.rm = T)) %>% 
          deframe()
        
        # Prepare output
        data.frame(cre_time,cre_val, ur_val)
      } %>% arrange(cre_time)
    
    # Prepare full output
    output$adm <- adm
    output$acute_renal_failure <- df$acute_renal_failure[i]
    output
    }else
    {
      output <- data.frame(cre_time = NA,
                           cre_val = NA,
                           ur_val = NA,
                           adm,
                           acute_renal_failure = df$acute_renal_failure[i])
      output
    }
  }
stopImplicitCluster()

results$adm %<>% as.factor()

results %<>% 
  filter(ur_val < 8000 &
           ur_val > 0) %>% 
  select(cre_val,
         ur_val,
         acute_renal_failure)

# Plot
results %>% 
  ggplot(aes(x = ur_val,
             y = cre_val))+
  geom_point(aes(colour = acute_renal_failure),
             alpha = 0.5)


# Create model
renal_model <- glm(acute_renal_failure ~ cre_val * ur_val,
                     data = results, family = "binomial",)
summary(renal_model)

# Test predictions to find threshold
results$pred <- predict(renal_model, newdata = results)

tuning_df <- foreach(i = seq(-2,2,0.001),
        .combine = "rbind") %do%
  {
    results$arf_pred <- results$pred > i
    
    conMat <- confusionMatrix(results$arf_pred %>% as.factor(),
                              results$acute_renal_failure %>% as.factor())
    
    data.frame(threshold = i,
               accuracy = conMat[["overall"]][["Accuracy"]],
               balanced_accuracy = conMat[["byClass"]][["Balanced Accuracy"]])
  }
plot(tuning_df$accuracy,
     tuning_df$balanced_accuracy)

# Measure distance from 1
tuning_df$dist <- sqrt(((tuning_df$accuracy - 1) ^2) + ((tuning_df$balanced_accuracy - 1)^2))

# Test and write
r_coef <- c(coefficients(renal_model), "threshold" = tuning_df$threshold[which.min(tuning_df$dist)])
r_coef[1] + (r_coef[2] * 2.8) + (r_coef[3] * 35) + (r_coef[4] * (2.8 * 35))
predict(renal_model, newdata = data.frame(cre_val = 2.8,
                                          ur_val = 35))
write_rds(r_coef, "models/acute_renal_failure.RDS")

