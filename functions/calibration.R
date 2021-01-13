calibration <- function(.df, probs_name)
{
  .df %>% 
    summarise(N = length(readmission),
              observed = (sum(readmission) / length(readmission)) * 100,
              predicted = mean(get(probs_name) * 100),
              error = sd(get(probs_name) * 100))
}