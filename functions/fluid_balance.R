# Adapted from concepts in https://github.com/MIT-LCP/mimic-code/tree/master/concepts

# https://github.com/MIT-LCP/mimic-code/blob/master/concepts/fluid_balance/colloid_bolus.sql
# https://github.com/MIT-LCP/mimic-code/blob/master/concepts/fluid_balance/crystalloid_bolus.sql
# https://github.com/MIT-LCP/mimic-code/blob/master/concepts/fluid_balance/urine_output.sql

fluid_balance <- function(icustay_df)
{
  # Debug
  icustay_df <- mimic_preproc$stays %>% 
    filter(hadm_id == adm) 
  
  # Extract database source
  dbsource <- icustay_df$dbsource
  
  # Extract inputs for metavision
  if(dbsource == "metavision")
  {
    # Read input events
    inputfile <- paste("data/events/inputeventsmv_", icustay_df$hadm_id, ".csv", sep = "")
    input_events <- read_csv(inputfile)
    
    # Select colloid and crystalloid bolus items
    fluids_in <- input_events %>% 
      filter(itemid %in% c(220864, 220862, 225174,
                           225795, 225796, 220861,
                           220863, 225158, 225828,
                           225797, 225159, 225161,
                           225823, 225825, 225827,
                           225941, 226089, 225944)) %>% 
      filter(statusdescription != "Rewritten")
    

    # Find plasma transfusions
    ffp_in <- input_events %>% 
      filter(itemid == 220970)
    
    # Find RBC transfusions
    rbc_in <- input_events %>% 
      filter(itemid == 225168)
    
    # Add to total fluids
    fluids_in %<>% rbind(ffp_in, rbc_in)
    
  }else
    # Extract inputs for carevue
  {
    
  }
  
  # Read output events
  outputfile <- paste("data/events/outputevents_", icustay_df$hadm_id, ".csv", sep = "")
  output_events <- read_csv(outputfile)
  
  # Select urine events
  urine_out <- output_events %>% 
    # Carevue urine items
    filter(itemid %in% c(40055, 43175, 43175, 40069,
                         40094, 40715, 40473, 40085,
                         40057, 40056, 40405, 40428,
                         40086, 40096, 40651,
                         # Metavision urine items 
                         226559, 226560, 226561,
                         226584, 226563, 226564,
                         226565, 226567, 226557,
                         226558, 227488, 227489))
  
  # Correct for input of GU irrigants
  urine_out %<>% 
    mutate(value = ifelse(itemid == 227488, -value, value))
  
  # Standardise units to ml and round
  urine_out %<>% 
    mutate(amount = ifelse(valueuom == "L", 
                           value * 1000, value) %>% 
             round())
  
  # Combine
  fluid_balance_df <- rbind(
    # Input events
    fluids_in %>% 
      select(starttime, amount) %>% 
      mutate(eventtype = "input") %>% 
      rename("value" = amount,
             "time" = starttime),
    # Output events
    urine_out %>% 
      select(charttime, value) %>% 
      mutate(eventtype = "output",
             value = -value) %>% 
      rename("time" = charttime)
  )
  
  # Sort and calculate balance
  fluid_balance_df %<>% 
    arrange(time) %>% 
    mutate(balance = cumsum(value))
  
  # Output max
  return(max(fluid_balance_df$balance))
}

# EXAMPLE EVENTS
# 100012 - metavision, fresh frozen plasma
# 100119 - carevue

events_carevue <- dir("data/events", pattern = "inputeventscv")
for(E in 1:length(events_carevue))
{
  if(any(read.csv(paste("data/events/", events_carevue[E], sep = ""))$itemid %in% 220970))
  {
    print(events_carevue[E])
    break
  }
}