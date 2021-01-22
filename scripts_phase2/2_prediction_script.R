# Packages
require(dplyr)
require(magrittr)
require(readr)
require(tibble)
require(brms)

source("functions/inverse_logit.R")

# Load data and models
results <- read_csv("data/patients_recalibrated.csv")
#hammer_frequentist <- read_rds("models/hammer_frequentist.rds")
hammer_bayes <- read_rds("models/hammer_bayes.rds")

# Extract bayesian posteriors
hammer_pst <- posterior_samples(hammer_bayes)
write_rds(hammer_pst, "models/hammer_posteriors.RDS")

# Set up prediction variables
male <- TRUE
general_surgery <- TRUE
cardiac_surgery <- FALSE
hyperglycemia <- FALSE
high_apache <- TRUE
fluid_balance_5L <- FALSE
ambulation <- TRUE
los_5 <- TRUE

# Function to select relevant posteriors
predict_posterior <- function(pst_df, sx = FALSE, gsurg = FALSE, csurg = FALSE, 
                              hyper = FALSE, apache = FALSE, fluid = FALSE, amb = FALSE,
                              los = FALSE)
{
  # Debug
  # sx <- male
  # gsurg <- general_surgery
  # csurg <- cardiac_surgery
  # hyper <- hyperglycemia
  # apache <- high_apache
  # fluid <- fluid_balance_5L
  # amb <- ambulation
  # los <- los_5
  
  pst_df <- hammer_pst
  
  # Combine variables
  varlist <- c(sx, gsurg, csurg, hyper, apache, fluid, amb, los)
  
  # Create coefficient list
  coef_names <- names(pst_df)[2:9]
  
  # Select relevant variables
  relevant_variables <- coef_names[varlist]
  
  # Extract from posterior data frame
  pred_df <- pst_df %>% 
    select(all_of(c("b_Intercept",relevant_variables)))
  
  # Create prediction posterior
  predictions <- pred_df %>% rowSums()
  
  # Transform to probabilities and output
  return(inverse_logit(predictions))
}

