# Packages
require(readr)
require(dplyr)
require(ggplot2)
require(scales)
require(tidyr)
require(RColorBrewer)
require(ResourceSelection)
require(ROCR)
require(magrittr)

source("functions/inverse_logit.R")
source("functions/calculate_apache_score.R")
source("functions/distance_4D.R")

# Load data----
comparison_df <- read_csv("data/impute/final_comparison.csv")
full_data <- read_csv("data/impute/complete_cases.csv") 
pal <- brewer.pal(9, "Set1")[-6]

# Score apache
full_data %<>% 
  mutate(apache_II_discharge = calculate_apache_scores(.))

# Create model and predict
apache_model <- glm(readmission ~ apache_II_discharge,
                    data = full_data, family = "binomial")
write_rds(apache_model, "models/apache_model.RDS")
probs <- predict(apache_model, newdata = full_data) %>% inverse_logit()

# Assess ROC
pred <- prediction(probs, full_data$readmission)
perf <- performance(pred, "tpr", "fpr")
auc <- performance(pred, measure = "auc")


# Assess calibration
hoslem_full <- hoslem.test(full_data$readmission,
                           probs, 10)

# Tabulate
results <- data.frame(method = "full_dataset",
                      split = 0,
                      discrimination = auc@y.values[[1]],
                      discrimination_error = 0,
                      calibration = hoslem_full$statistic,
                      calibration_error = 0,
                      distance = 0,
                      mean_error = 0,
                      max = 0,
                      max_error = 0,
                      dist_4D = 0,
                      error_4D = 0)

#Plot-----

# Plot each dimension
comparison_df %>% 
  pivot_longer(c(3,5,7,9)) %>% 
  mutate(error = case_when(name == "distance" ~ mean_error,
                           name == "max" ~ max_error,
                           name == "calibration" ~ calibration_error,
                           name == "discrimination" ~ discrimination_error)) %>% 
  ggplot(aes(x = split,
             y = value))+
  geom_ribbon(aes(ymin = value - error,
                  ymax = value + error,
                  fill = method), alpha = 0.1)+
  geom_path(aes(colour = method), size = 1)+
  geom_hline(data = results %>% pivot_longer(c(3,5,7,9)),
             aes(yintercept = value), size = 1,
             linetype = "dashed")+
  theme_classic(20)+
  theme(legend.position = "top")+
  facet_wrap(~name, scales = "free_y")+
  scale_colour_manual(values = pal,
                      name = "")+
  scale_fill_manual(values = pal,
                    name = "")+
  labs(x = "Proportion missing")

# Multidimensional distance---------

# Rescale
rescale_df <- comparison_df %>%
  ungroup() %>%
  mutate(
    discrimination = rescale(discrimination,
                             from = c(
                               results$discrimination,
                               min(discrimination)
                             )),
    calibration = rescale(calibration,
                          from = c(results$calibration,
                                   max(calibration))),
    distance = rescale(distance,
                       from = c(results$distance,
                                max(distance))),
    max = rescale(max,
                  from = c(results$max,
                           max(max)))
  ) 

# Measure distance
md_dist <- distance_4D(data.frame(discrimination = 0,
                                  calibration = 0,
                                  distance = 0,
                                  max = 0), rescale_df$discrimination,
                       rescale_df$calibration,
                       rescale_df$distance,
                       rescale_df$max)

rescale_df %<>% 
  mutate(md_dist)


# Plot 
rescale_df %>%
  ggplot(aes(x = split,
             y = md_dist))+
  geom_path(aes(colour = method),
            size = 1)+
  geom_point(size = 4, aes(fill = method),
             colour = "black", shape = 21)+
  theme_classic(20)+
  theme(legend.position = "top")+
  scale_colour_manual(values = pal,
                      name = "")+
  scale_fill_manual(values = pal,
                    name = "")+
  labs(x = "Proportion missing",
       y = "Multidimensional distance")

# Plots-----

# MICE, RF, PCA,

comparison_df %>% 
  pivot_longer(c(3,5,7,9)) %>% 
  filter(method == "median" |
           method == "MICE" |
           method == "RF" |
           method == "assume_zero") %>% 
  mutate(error = case_when(name == "distance" ~ mean_error,
                           name == "max" ~ max_error,
                           name == "calibration" ~ calibration_error,
                           name == "discrimination" ~ discrimination_error)) %>% 
  ggplot(aes(x = split,
             y = value))+
  geom_ribbon(aes(ymin = value - error,
                  ymax = value + error,
                  fill = method), alpha = 0.1)+
  geom_path(aes(colour = method), size = 1)+
  geom_hline(data = results %>% pivot_longer(c(3,5,7,9)),
             aes(yintercept = value), size = 1,
             linetype = "dashed")+
  theme_classic(20)+
  theme(legend.position = "none")+
  facet_wrap(~name, scales = "free_y")+
  scale_colour_manual(values = pal,
                      name = "")+
  scale_fill_manual(values = pal,
                    name = "")+
  labs(x = "Proportion missing",
       y = "")# + scale_y_reverse()
ggsave("writeup/presentation_figs/mice-RF.png",
       width = 33.8, height = 17, units = "cm")


rescale_df %>%
  filter(method == "median" |
           method == "MICE" |
           method == "RF" |
           method == "assume_zero") %>% 
  ggplot(aes(x = split,
             y = md_dist))+
  geom_path(aes(colour = method),
            size = 1)+
  geom_point(size = 4, aes(fill = method),
             colour = "black", shape = 21)+
  theme_classic(20)+
  theme(legend.position = "top")+
  scale_colour_manual(values = pal,
                      name = "")+
  scale_fill_manual(values = pal,
                    name = "")+
  labs(x = "Proportion missing",
       y = "Multidimensional distance")
