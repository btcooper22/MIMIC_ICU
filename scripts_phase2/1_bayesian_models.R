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
cooper_frequentist <- glmnet(x, y, alpha = 1, lambda = cv$lambda.min,
                             family = "binomial")

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

# Create Cooper formula
formu <- paste(names(train_df)[feature_id], collapse = " + ")

# Create scaled dataset
scale_df <- train_df %>% 
  select(all_of(c(7, feature_id))) %>% 
  mutate(apache_II = scale(apache_II),
         serum_glucose = scale(serum_glucose),
         respiratory_rate = scale(respiratory_rate),
         final_temp = scale(final_temp),
         final_SpO2 = scale(final_SpO2),
         final_lactate = scale(final_lactate),
         final_bp = scale(final_bp),
         apache_II_discharge = scale(apache_II_discharge))

cooper_bayes <- brm(formula(paste("readmission ~", formu,
                                  sep = " ")),
                    data = train_df, 
                    prior = set_prior(horseshoe(df = 1, par_ratio = 1)),
                    family = bernoulli(link = "logit"),
                    warmup = 500,
                    iter = 6250,
                    chains = 8,
                    cores = 8,
                    control = list(adapt_delta = 0.95))

# Compare: Hammer----------

# Add into one data frame
fixef(hammer_bayes, robust = TRUE) %>% 
  as.data.frame() %>%
  rownames_to_column("Variable") %>% 
  select(Variable, Estimate, Q2.5, Q97.5) %>% 
  as_tibble() %>% 
  rename(Bayesian = "Estimate") %>% 
  mutate(Frequentist = coefficients(hammer_frequentist)) %>% 
  select(Variable, Frequentist, Bayesian, Q2.5, Q97.5)

# Compare: Cooper--------

# Extract posteriors
posterior_df <- posterior_samples(cooper_bayes)

# Rescale
test <- posterior_df %>% 
  mutate(b_apache_II = (b_apache_II * attributes(scale_df$apache_II)[[3]]))

# Summarise
fixef(cooper_bayes)
cooper_frequentist$beta

