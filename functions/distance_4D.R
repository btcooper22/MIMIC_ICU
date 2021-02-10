distance_4D <- function(comparison_point, new_discrimination,
                        new_calibration, new_distance, new_max)
{
  # Debug
  # comparison_point <- results
  # new_discrimination <-  auc_zero@y.values[[1]]
  # new_calibration <-  cal_zero$statistic
  # new_distance <-  mean(zero_dist)
  # new_max  <-  quantile(zero_dist, 0.95)
  
  # Measure distances
  discrimination_dist <- (comparison_point$discrimination - new_discrimination)^2 %>% unname()
  calibration_dist <- (comparison_point$calibration - new_calibration)^2 %>% unname()
  distance_dist <- (comparison_point$distance - new_distance)^2 %>% unname()
  max_dist <- (comparison_point$max - new_max)^2 %>% unname()
  
  # Sum and output
  return(
  sqrt(discrimination_dist + 
         calibration_dist + 
         distance_dist + 
         max_dist)
  )
}