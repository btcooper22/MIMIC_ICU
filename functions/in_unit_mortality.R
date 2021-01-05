in_unit_mortality <- function(df) #add_inunit_mortality_to_icustays
{
  mortality <- !is.na(df$deathtime) & ((df$intime <= df$deathtime) & (df$outtime >= df$deathtime))
  mortality = mortality | !is.na(df$deathtime) & 
    !is.na(df$dod) & 
    ((df$intime <= df$dod) & (df$outtime >= df$dod))
  return(mortality)
}