# Packages
require(readr)
require(magrittr)
require(dplyr)
require(ROCR)
require(ggplot2)
require(ResourceSelection)
require(patchwork)
require(forcats)
require(RColorBrewer)

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

# Calculate median APACHE-II
data_inunit %>% 
  filter(!is.na(apache_II)) %>% 
  summarise(median = median(apache_II),
            sd = sd(apache_II))

data_30d %>% 
  filter(!is.na(apache_II)) %>% 
  summarise(median = median(apache_II),
            sd = sd(apache_II))

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


# Measure missingness----

# Correct fio2
data_inunit$fractioninspiredoxygen[is.na(data_inunit$fractioninspiredoxygen)] <- 21
data_30d$fractioninspiredoxygen[is.na(data_30d$fractioninspiredoxygen)] <- 21

# Create missingness data frames
missing_inunit <-
  data.frame(
    temperature = is.na(data_inunit$temperature),
    map = is.na(data_inunit$systolicbp) |
      is.na(data_inunit$diastolicbp),
    pulse = is.na(data_inunit$pulse),
    respiratory = is.na(data_inunit$respiratory),
    sodium = is.na(data_inunit$sodium),
    potassium = is.na(data_inunit$potassium),
    creatinine = is.na(data_inunit$creatinine),
    haematocrit = is.na(data_inunit$haematocrit),
    whitebloodcount = is.na(data_inunit$whitebloodcount),
    gcs = is.na(data_inunit$glasgowcomaeye) |
      is.na(data_inunit$glasgowcomamotor) |
      is.na(data_inunit$glasgowcomaverbal),
    arterialpH = is.na(data_inunit$arterialpH) &
      is.na(data_inunit$bicarbonate),
    oxygenation = ifelse(
      data_inunit$fractioninspiredoxygen > 50,
      is.na(data_inunit$arterialoxygen) |
        is.na(data_inunit$arterialcarbon),
      is.na(data_inunit$arterialoxygen)
    ) &
      is.na(data_inunit$bicarbonate)
  )

missing_30d <-
  data.frame(
    temperature = is.na(data_30d$temperature),
    map = is.na(data_30d$systolicbp) |
      is.na(data_30d$diastolicbp),
    pulse = is.na(data_30d$pulse),
    respiratory = is.na(data_30d$respiratory),
    sodium = is.na(data_30d$sodium),
    potassium = is.na(data_30d$potassium),
    creatinine = is.na(data_30d$creatinine),
    haematocrit = is.na(data_30d$haematocrit),
    whitebloodcount = is.na(data_30d$whitebloodcount),
    gcs = is.na(data_30d$glasgowcomaeye) |
      is.na(data_30d$glasgowcomamotor) |
      is.na(data_30d$glasgowcomaverbal),
    arterialpH = is.na(data_30d$arterialpH) &
      is.na(data_30d$bicarbonate),
    oxygenation = ifelse(
      data_30d$fractioninspiredoxygen > 50,
      is.na(data_30d$arterialoxygen) |
        is.na(data_30d$arterialcarbon),
      is.na(data_30d$arterialoxygen)
    ) &
      is.na(data_30d$bicarbonate)
  )

# Frequency of missingess by variable
miss_table <- colSums(missing_inunit) %>%
  as.data.frame() %>% 
  mutate(model = "In-unit mortality") %>% 
  rbind(
    colSums(missing_30d) %>% 
      as.data.frame() %>% 
      mutate(model = "30-day mortality")) 
names(miss_table)[1] <- "Freq"
miss_table$varname <- rep(names(missing_30d), 2)
  
plot_miss1 <- miss_table %>% 
ggplot(aes(x = varname, y = Freq))+
  geom_bar(aes(fill = model),
           stat = "identity",
           position = "dodge",
           colour = "black") +
  theme_classic(20)+
  theme(legend.position = "top")+
  labs(x = "",
       y = "Frequency")+
  scale_fill_manual(values = c("#984ea3", "#ff7f00"),
                    name = "")

# Number of variables missing
plot_miss2 <- table(rowSums(missing_inunit)) %>% 
  as.data.frame() %>% 
  mutate(model = "In-unit mortality") %>% 
  rbind(table(rowSums(missing_30d)) %>% 
          as.data.frame() %>% 
          mutate(model = "30-day mortality")) %>%
  filter(Var1 != 0) %>% 
  ggplot(aes(x = Var1, y = Freq))+
  geom_bar(aes(fill = model),
           stat = "identity",
           position = "dodge",
           colour = "black") +
  theme_classic(20)+
  theme(legend.position = "top")+
  labs(x = "Number of variables missing",
       y = "Frequency")+
  scale_fill_manual(values = c("#984ea3", "#ff7f00"),
                     name = "")

ggsave("writeup/presentation_figs/misstogram_A.png",
       plot_miss2, width = 33.8, height = 17, units = "cm")

ggsave("writeup/presentation_figs/misstogram_B.png",
       plot_miss1, width = 33.8, height = 17, 
       units = "cm", scale = 1.5)

# Score imputations--------

pal <- brewer.pal(9, "Set1")[-6]

# Load data
results <- read_rds("models/boot_samples.RDS")

# Adjust factor levels
results %<>% 
  mutate(mortality = fct_recode(mortality, `30-day mortality` = "30-day",
                                `In-unit mortality` = "inunit"),
         method = fct_recode(method, Amelia = "amelia",
                             Median = "average",
                             Longitudinal = "recent",
                             `Random Forest` = "RF",
                             `Assume zero` = "zero"))

# Extract means
means_df <- results %>% 
  group_by(method, mortality) %>% 
  summarise(discrimination = mean(discrim),
            calibration = mean(calib)) %>% 
  ungroup()


# Zero only
p1 <- results %>% 
  filter(mortality == "30-day mortality",
         method == "Assume zero") %>% 
  ggplot()+
  theme_classic(20)+
  labs(y = "Calibration", x = "Discrimination")+
  scale_colour_manual(values = c("#999999"),
                      name = "")+
  scale_fill_manual(values = c("#999999"),
                    name = "")+
  scale_y_reverse()+
  theme(legend.position = "top")+
  stat_ellipse(aes(x = discrim, y = calib,
                   colour = method),
               size = 2, type = "norm",
               level = 0.682)+
  geom_point(data = means_df %>% 
               filter(mortality == "30-day mortality",
                      method == "Assume zero"),
             aes(x = discrimination,
                 y = calibration,
                 fill = method),
             size = 4, shape = 21)+
  coord_cartesian(xlim = c(0.65,1),
                  ylim = c(22, 0))

p2 <- results %>% 
  filter(mortality == "In-unit mortality",
         method == "Assume zero") %>% 
  ggplot()+
  theme_classic(20)+
  labs(y = "Calibration", x = "Discrimination")+
  scale_colour_manual(values = c("#999999"),
                      name = "")+
  scale_fill_manual(values = c("#999999"),
                    name = "")+
  scale_y_reverse()+
  theme(legend.position = "top")+
  stat_ellipse(aes(x = discrim, y = calib,
                   colour = method),
               size = 2, type = "norm",
               level = 0.682)+
  geom_point(data = means_df %>% 
               filter(mortality == "In-unit mortality",
                      method == "Assume zero"),
             aes(x = discrimination,
                 y = calibration,
                 fill = method),
             size = 4, shape = 21)+
  coord_cartesian(xlim = c(0.25,0.83),
                  ylim = c(30, 0))

p1 + p2
ggsave("writeup/presentation_figs/ellipse_1.png",
       width = 33.8, height = 14, units = "cm")

# Zero + average
p1 <- results %>% 
  filter(mortality == "30-day mortality",
         method == "Assume zero" |
           method == "Median") %>% 
  ggplot()+
  theme_classic(20)+
  labs(y = "Calibration", x = "Discrimination")+
  scale_colour_manual(values = c("black","#999999"),
                      name = "")+
  scale_fill_manual(values = c("black","#999999"),
                    name = "")+
  scale_y_reverse()+
  theme(legend.position = "top")+
  stat_ellipse(aes(x = discrim, y = calib,
                   colour = method),
               size = 2, type = "norm",
               level = 0.682)+
  geom_point(data = means_df %>% 
               filter(mortality == "30-day mortality",
                      method == "Assume zero" |
                        method == "Median"),
             aes(x = discrimination,
                 y = calibration,
                 fill = method),
             size = 4, shape = 21)+
  coord_cartesian(xlim = c(0.65,1),
                  ylim = c(22, 0))

p2 <- results %>% 
  filter(mortality == "In-unit mortality",
         method == "Assume zero" |
           method == "Median") %>% 
  ggplot()+
  theme_classic(20)+
  labs(y = "Calibration", x = "Discrimination")+
  scale_colour_manual(values = c("black","#999999"),
                      name = "")+
  scale_fill_manual(values = c("black","#999999"),
                    name = "")+
  scale_y_reverse()+
  theme(legend.position = "top")+
  stat_ellipse(aes(x = discrim, y = calib,
                   colour = method),
               size = 2, type = "norm",
               level = 0.682)+
  geom_point(data = means_df %>% 
               filter(mortality == "In-unit mortality",
                      method == "Assume zero"|
                        method == "Median"),
             aes(x = discrimination,
                 y = calibration,
                 fill = method),
             size = 4, shape = 21)+
  coord_cartesian(xlim = c(0.25,0.83),
                  ylim = c(30, 0))

p1 + p2
ggsave("writeup/presentation_figs/ellipse_2.png",
       width = 33.8, height = 14, units = "cm")


# All except recent & mice
p1 <- results %>% 
  filter(mortality == "30-day mortality",
         method %in% c("MICE", "Recent") == FALSE) %>% 
  ggplot()+
  theme_classic(20)+
  labs(y = "Calibration", x = "Discrimination")+
  scale_colour_manual(values = c("#ff7f00", "black","#984ea3",
                                 "#a65628","#4daf4a","#999999"),
                      name = "")+
  scale_fill_manual(values = c("#ff7f00", "black","#984ea3",
                               "#a65628","#4daf4a","#999999"),
                    name = "")+
  scale_y_reverse()+
  theme(legend.position = "top")+
  stat_ellipse(aes(x = discrim, y = calib,
                   colour = method),
               size = 2, type = "norm",
               level = 0.682)+
  geom_point(data = means_df %>% 
               filter(mortality == "30-day mortality",
                      method %in% c("MICE", "Recent") == FALSE),
             aes(x = discrimination,
                 y = calibration,
                 fill = method),
             size = 4, shape = 21)+
  coord_cartesian(xlim = c(0.65,1),
                  ylim = c(22, 0))

p2 <- results %>% 
  filter(mortality == "In-unit mortality",
         method %in% c("MICE", "Recent") == FALSE) %>% 
  ggplot()+
  theme_classic(20)+
  labs(y = "Calibration", x = "Discrimination")+
  scale_colour_manual(values = c("#ff7f00", "black","#984ea3",
                                 "#a65628","#4daf4a","#999999"),
                      name = "")+
  scale_fill_manual(values = c("#ff7f00", "black","#984ea3",
                               "#a65628","#4daf4a","#999999"),
                    name = "")+
  scale_y_reverse()+
  theme(legend.position = "top")+
  stat_ellipse(aes(x = discrim, y = calib,
                   colour = method),
               size = 2, type = "norm",
               level = 0.682)+
  geom_point(data = means_df %>% 
               filter(mortality == "In-unit mortality",
                      method %in% c("MICE", "Recent") == FALSE),
             aes(x = discrimination,
                 y = calibration,
                 fill = method),
             size = 4, shape = 21)+
  coord_cartesian(xlim = c(0.25,0.83),
                  ylim = c(30, 0))

p1 + p2
ggsave("writeup/presentation_figs/ellipse_3.png",
       width = 33.8, height = 14, units = "cm")


# Add mice
p1 <- results %>% 
  filter(mortality == "30-day mortality",
         method != "Recent") %>% 
  ggplot()+
  theme_classic(20)+
  labs(y = "Calibration", x = "Discrimination")+
  scale_colour_manual(values = c("#ff7f00", "black","#984ea3", "#377eb8",
                                 "#a65628","#4daf4a","#999999"),
                      name = "")+
  scale_fill_manual(values = c("#ff7f00", "black","#984ea3","#377eb8",
                               "#a65628","#4daf4a","#999999"),
                    name = "")+
  scale_y_reverse()+
  theme(legend.position = "top")+
  stat_ellipse(aes(x = discrim, y = calib,
                   colour = method),
               size = 2, type = "norm",
               level = 0.682)+
  geom_point(data = means_df %>% 
               filter(mortality == "30-day mortality",
                      method != "Recent"),
             aes(x = discrimination,
                 y = calibration,
                 fill = method),
             size = 4, shape = 21)+
  coord_cartesian(xlim = c(0.65,1),
                  ylim = c(22, 0))

p2 <- results %>% 
  filter(mortality == "In-unit mortality",
         method != "Recent") %>% 
  ggplot()+
  theme_classic(20)+
  labs(y = "Calibration", x = "Discrimination")+
  scale_colour_manual(values = c("#ff7f00", "black","#984ea3", "#377eb8",
                                 "#a65628","#4daf4a","#999999"),
                      name = "")+
  scale_fill_manual(values = c("#ff7f00", "black","#984ea3","#377eb8",
                               "#a65628","#4daf4a","#999999"),
                    name = "")+
  scale_y_reverse()+
  theme(legend.position = "top")+
  stat_ellipse(aes(x = discrim, y = calib,
                   colour = method),
               size = 2, type = "norm",
               level = 0.682)+
  geom_point(data = means_df %>% 
               filter(mortality == "In-unit mortality",
                      method != "Recent"),
             aes(x = discrimination,
                 y = calibration,
                 fill = method),
             size = 4, shape = 21)+
  coord_cartesian(xlim = c(0.25,0.83),
                  ylim = c(30, 0))

p1 + p2
ggsave("writeup/presentation_figs/ellipse_4.png",
       width = 33.8, height = 14, units = "cm")

# All
p1 <- results %>% 
  filter(mortality == "30-day mortality") %>% 
  ggplot()+
  theme_classic(20)+
  labs(y = "Calibration", x = "Discrimination")+
  scale_colour_manual(values = c("#ff7f00", "black","#984ea3","#377eb8",
                                 "#a65628","#e41a1c","#4daf4a","#999999"),
                      name = "")+
  scale_fill_manual(values = c("#ff7f00", "black","#984ea3","#377eb8",
                               "#a65628","#e41a1c","#4daf4a","#999999"),
                    name = "")+
  scale_y_reverse()+
  theme(legend.position = "top")+
  stat_ellipse(aes(x = discrim, y = calib,
                   colour = method),
               size = 2, type = "norm",
               level = 0.682)+
  geom_point(data = means_df %>% 
               filter(mortality == "30-day mortality"),
             aes(x = discrimination,
                 y = calibration,
                 fill = method),
             size = 4, shape = 21)+
  coord_cartesian(xlim = c(0.65,1),
                  ylim = c(22, 0))

p2 <- results %>% 
  filter(mortality == "In-unit mortality") %>% 
  ggplot()+
  theme_classic(20)+
  labs(y = "Calibration", x = "Discrimination")+
  scale_colour_manual(values = c("#ff7f00", "black","#984ea3", "#377eb8",
                                 "#a65628","#4daf4a","#999999"),
                      name = "")+
  scale_fill_manual(values = c("#ff7f00", "black","#984ea3","#377eb8",
                               "#a65628","#4daf4a","#999999"),
                    name = "")+
  scale_y_reverse()+
  theme(legend.position = "top")+
  stat_ellipse(aes(x = discrim, y = calib,
                   colour = method),
               size = 2, type = "norm",
               level = 0.682)+
  geom_point(data = means_df %>% 
               filter(mortality == "In-unit mortality"),
             aes(x = discrimination,
                 y = calibration,
                 fill = method),
             size = 4, shape = 21)+
  coord_cartesian(xlim = c(0.25,0.83),
                  ylim = c(30, 0))

p1 + p2
ggsave("writeup/presentation_figs/ellipse_5.png",
       width = 33.8, height = 14, units = "cm")


# Subsets and combinations----

# Load data
results_subset <- read_rds("models/boot_results_subset.RDS") 

# Adjust factor levels
results_subset %<>% 
  mutate(type = fct_recode(type, Amelia = "amelia",
                             Longitudinal = "longitudinal",
                             MICE = "mice"))

# Extract means
means_df <- results_subset %>% 
  group_by(type) %>% 
  summarise(discrimination = mean(discrim),
            calibration = mean(calib)) %>% 
  ungroup()


# Plot
results_subset %>% 
  # filter(mortality == "30-day mortality",
  #        method == "Assume zero") %>% 
  ggplot()+
  theme_classic(20)+
  labs(y = "Calibration", x = "Discrimination")+
  scale_colour_manual(values = c("#e41a1c", "#a65628", 
                                 "#984ea3"), name = "")+
  scale_fill_manual(values = c("#e41a1c", "#a65628", 
                               "#984ea3"), name = "")+
  scale_y_reverse()+
  theme(legend.position = "top")+
  stat_ellipse(aes(x = discrim, y = calib,
                   colour = type),
               size = 2, type = "norm",
               level = 0.682,
               alpha = 0.8)+
  geom_point(data = means_df ,
             aes(x = discrimination,
                 y = calibration,
                 fill = type),
             size = 4, shape = 21)
