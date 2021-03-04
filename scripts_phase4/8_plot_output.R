# Packages
require(readr)
require(magrittr)
require(dplyr)
require(ROCR)
require(ggplot2)
require(ResourceSelection)
require(patchwork)

# Functions
source("functions/inverse_logit.R")
calibration <- function(.df)
{
  .df %>% 
    group_by(decile) %>%
    summarise(N = length(mortality),
              observed = (sum(mortality) / length(mortality)) * 100,
              predicted = mean(probs * 100),
              error = sd(probs * 100))
}

# Load data
data_inunit <- read_csv("data/apache_real_missing.csv")
data_30d <- read_csv("data/apache_discharge_missing.csv")

# Calculate mortality rates
sum(data_inunit$mort_inunit)
round(mean(data_inunit$mort_inunit) * 100, 2)

sum(data_30d$mort_30)
round(mean(data_30d$mort_30) * 100, 2)

# Measure missingess
data_inunit %>% 
  summarise(sum = sum(is.na(apache_II)),
            mean = mean(is.na(apache_II)) * 100)

data_30d %>% 
  summarise(sum = sum(is.na(apache_II)),
            mean = mean(is.na(apache_II)) * 100)

# Calibration and discrimination----

# Build models
model_inunit <- glm(mort_inunit ~ apache_II,
                    data = data_inunit, family = "binomial")

model_30d <- glm(mort_30 ~ apache_II,
                    data = data_30d, family = "binomial")

# Generate predictions
probs_inunit <- data.frame(probs = predict(model_inunit, newdata = data_inunit),
                           outcome = data_inunit$mort_inunit) %>% na.omit() %>% 
  mutate(probs = inverse_logit(probs))

probs_30d <- data.frame(probs = predict(model_30d, newdata = data_30d),
                           outcome = data_30d$mort_30) %>% na.omit() %>% 
  mutate(probs = inverse_logit(probs))

# Assess predictions
pred_inunit <- prediction(probs_inunit$probs,
                          probs_inunit$outcome)

pred_30d <- prediction(probs_30d$probs,
                          probs_30d$outcome)

# Generate AUC
performance(pred_inunit, measure = "auc")@y.values[[1]]
performance(pred_30d, measure = "auc")@y.values[[1]]

# Generate performance objects
performance_inunit <- performance(pred_inunit, "tpr", "fpr")
performance_30d <- performance(pred_30d, "tpr", "fpr")

# Plot AUC
plot_AUC <- data.frame(x = performance_inunit@x.values[[1]],
           y = performance_inunit@y.values[[1]],
           model = "In-unit mortality") %>% 
  rbind(
    data.frame(x = performance_30d@x.values[[1]],
               y = performance_30d@y.values[[1]],
               model = "30-day mortality")) %>% 
  ggplot(aes(x, y, colour = model))+
  geom_abline(slope = 1, intercept = 0,
              linetype = "dotted",
              size = 1)+
  geom_path(size = 1)+
  labs(x = "1 - Specificity",
       y = "Sensitivity")+
  theme_classic(20)+
  scale_color_manual(values = c("#984ea3", "#ff7f00"),
                     name = "")+
  theme(legend.position = "top")

# Split data into deciles
deciles_inunit <- tibble(
  mortality = probs_inunit$outcome,
  probs = probs_inunit$probs,
  decile = ntile(probs_inunit$probs, 10)
)

deciles_30d <- tibble(
  mortality = probs_30d$outcome,
  probs = probs_30d$probs,
  decile = ntile(probs_30d$probs, 10)
)

# Calculate calibration
calib_inunit <- calibration(deciles_inunit)
calib_30d <- calibration(deciles_30d)

# Quantify calibration
hoslem.test(probs_inunit$outcome,
            probs_inunit$probs, g = 10)
hoslem.test(probs_30d$outcome,
            probs_30d$probs, g = 10)

# Plot calibration
plot_cal <- rbind(calib_inunit %>% mutate(model = "In-unit mortality"), 
      calib_30d %>% mutate(model = "30-day mortality")) %>% 
  ggplot(aes(predicted, observed ,
             colour = model))+
  geom_abline(slope = 1, intercept = 0,
              size = 1, linetype = "dotted")+
  geom_path(size = 1)+
  geom_pointrange(aes(ymin = observed - error ,
                      ymax = observed + error,
                      y = observed,
                      x = predicted,
                      colour = model)) +
  theme_classic(20)+
  theme(legend.position = "top")+
  labs(x = "Predicted readmission",
       y = "Observed readmission")+
  scale_color_manual(values = c("#984ea3", "#ff7f00"),
                     name = "")+
  coord_cartesian(xlim = c(0,25),
                  ylim = c(0,25))

# Save plots
plot_AUC + plot_cal
# ggsave("writeup/presentation_figs/auc_cal.png",
#        width = 33.8, height = 17, units = "cm")
