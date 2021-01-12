inverse_logit <- function(logit){
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
}