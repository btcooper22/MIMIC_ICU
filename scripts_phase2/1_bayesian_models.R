# Packages
require(dplyr)
require(magrittr)
require(readr)
require(tibble)
require(brms)
require(glmnet)

# Fit frequentist models----------

# Load data
results <- read_csv("data/patients_recalibrated.csv")
train_df <- results %>% filter(type == "train")
validate_df <- results %>% filter(type == "validate")

# Fit hammer model
hammer_frequentist <- glm(readmission ~ sex + general_surgery + cardiac_surgery + hyperglycemia +
                            high_apache + fluid_balance_5L + ambulation + los_5,
                          data = train_df, family = "binomial")

# Vector of acceptable features
feature_id <- c(9:10, 12, 14, 16:25, 27:35, 51, 74)

# Binarise readmission column
lasso_df <- train_df %>%
  mutate(readmission = ifelse(readmission == TRUE, 1, 0))

# Select predictor and outcome vectors
x <- model.matrix(readmission~., lasso_df[,c(7,feature_id)])[,-1]
y <- lasso_df$readmission

# Find optimal lambda
cv <- cv.glmnet(x, y, alpha = 1)

# Fit LASSO cooper model
cooper_model <- glmnet(x, y, alpha = 1, lambda = cv$lambda.min)

# Fit Bayesian models----------
hammer_bayes <- brm(readmission ~ sex + general_surgery + cardiac_surgery + hyperglycemia +
                      high_apache + fluid_balance_5L + ambulation + los_5,
                    data = train_df, 
                    family = bernoulli(link = "logit"),
                    warmup = 500,
                    iter = 6250,
                    chains = 8,
                    cores = 8,
                    file = "models/hammer_bayes")

formu <- paste(names(train_df)[feature_id], collapse = " + ")

cooper_bayes <- brm(formula(paste("readmission ~", formu,
                                  sep = " ")),
                    data = train_df, 
                    prior = set_prior(lasso(df = cv$lambda.min)),
                    family = bernoulli(link = "logit"),
                    warmup = 500,
                    iter = 6250,
                    chains = 8,
                    cores = 8,
                    file = "models/cooper_bayes")

# Compare: Hammer----------

# Add into one data frame
fixef(hammer_bayes) %>% 
  as_tibble() %>% 
  rownames_to_column() %>% 
  select(Estimate, Q2.5, Q97.5)

