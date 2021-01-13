calibration <- function(.df, probs_name, quantiles_name)
{
  .df %>% 
  select(any_of(c("readmission", probs_name, quantiles_name))) %>% 
    group_by_at(quantiles_name) %>%
    summarise(N = length(readmission),
              observed = (sum(readmission) / length(readmission)) * 100,
              predicted = mean(get(probs_name) * 100),
              error = sd(get(probs_name) * 100))
}