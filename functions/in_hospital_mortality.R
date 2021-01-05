in_hospital_mortality <- function(df) #add_inhospital_mortality_to_icustays
{
  mortality <- !is.na(df$deathtime) & ((df$admittime <= df$deathtime) & (df$dischtime >= df$deathtime))
  mortality <-  mortality | !is.na(df$deathtime) & 
    !is.na(df$dod) & 
    ((df$intime <= df$dod) & (df$outtime >= df$dod))
  return(mortality)
}