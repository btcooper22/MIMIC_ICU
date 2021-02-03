# Packages
require(dplyr)
require(magrittr)
require(readr)
require(tibble)
require(ggplot2)
require(ROCR)
require(foreach)
require(stringr)

source("functions/inverse_logit.R")
source("functions/calibration.R")

# Load data----- 
# Trim to APACHE-II
full_data <- read_csv("data/final_patients.csv") %>% 
  select(matches("_discharge$")) %>% 
  select(matches("^apache_")) %>% 
  cbind(read_csv("data/final_patients.csv") %>% 
          select(row_id, subject_id,
                 adm_id, readmission,
                 age),.)

# Binarise readmission
full_data %<>% 
  mutate(readmission = readmission == "Readmitted to ICU") %>% 
  # Remove cases where score comes from bicarbonate
  filter(!is.na(apache_respiratory_discharge)) %>% 
  select(-apache_bicarbonate_discharge)
  
# Confirm no missing values
table(full_data$readmission, useNA = "ifany")
summary(full_data)

# Create and assess APACHE model------
# Create model and predict
apache_model <- glm(readmission ~ apache_II_discharge,
                    data = full_data, family = "binomial")
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

# Simulate MCAR (missing completely at random)---------
# Set n
N <- 1

# Set splits
splits <- c(0.01, 0.05, seq(0.1, 0.8, 0.1))

# Pull APACHE matrix
apache_matrix <- full_data[,7:20]

# Outer loop over N
foreach(n = 1:N) %do%
  {
    # Inner loop for splits
    foreach(s = 1:length(splits)) %do%
      {
        # Define shape
        shape <- nrow(apache_matrix) * ncol(apache_matrix)
        
        # Draw randomisation array
        rand_array <- rbinom(shape, 1, splits[s]) %>% 
          as.logical()
        
        # Reshape array
        NA_array <- matrix(rand_array, 
               nrow(apache_matrix),
               ncol(apache_matrix))
        
        # Remove values
        data_out <- apache_matrix
        data_out[NA_array] <- NA
        
        # Write to file
        filename <- paste("data/MCAR/S", str_pad(splits[s], 4, "right", 0),
                          "_N", str_pad(n,3,"left",0), ".csv", sep = "")
        write.csv(data_out, filename,
                  row.names = FALSE)
      }
  }

# Basic methods -> Assume 0, Mean, hot-deck

# Intermediate methods -> MICE, KNN, MissForest

# Complex methods from competiton