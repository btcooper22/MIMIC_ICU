brier_extraction <- function(obs, pred)
{
  # Debug
  # obs <- patients$readmission == "Readmitted to ICU"
  # pred <- probs_hammer
  
  # Calculate brier score
  brier_score <- BrierScore(obs, pred)
  
  # Return details
  return(
    tibble(
      score = brier_score$bs,
      reliability = brier_score$rel,
      resolution = brier_score$res
    )
  )
}