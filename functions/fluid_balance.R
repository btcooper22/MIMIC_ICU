# Adapted from concepts in https://github.com/MIT-LCP/mimic-code/tree/master/concepts

# https://github.com/MIT-LCP/mimic-code/blob/master/concepts/fluid_balance/colloid_bolus.sql
# https://github.com/MIT-LCP/mimic-code/blob/master/concepts/fluid_balance/crystalloid_bolus.sql
# https://github.com/MIT-LCP/mimic-code/blob/master/concepts/fluid_balance/urine_output.sql

fluid_balance <- function(icustay_df, chart_df, readmission)
{
  # Debug
  # icustay_df <- mimic_preproc$stays %>%
  #  filter(hadm_id == adm)
  # chart_df <- charts
  #readmission <- outcomes$readmission[i]
  
  # Find icu stay
  if(readmission)
  {
    icustay_df <- icustay_df[1,]
  }
  
  # Extract database source
  dbsource <- icustay_df$dbsource
  
  # Extract inputs for metavision
  if(dbsource == "metavision")
  {
    # Read input events
    inputfile <- paste("data/events/inputeventsmv_", icustay_df$hadm_id, ".csv", sep = "")
    input_events <- read_csv(inputfile, col_types = c("ddddTTddcdcTdddccccddcdddcccTdd"))
    
    # Select colloid and crystalloid bolus items and round
    fluids_in <- input_events %>% 
      filter(itemid %in% c(220864, 220862, 225174,
                           225795, 225796, 220861,
                           220863, 225158, 225828,
                           225797, 225159, 225161,
                           225823, 225825, 225827,
                           225941, 226089, 225944)) %>% 
      filter(statusdescription != "Rewritten") %>% 
      mutate(amount = ifelse(amountuom == "L", 
                             amount * 1000, amount) %>% 
               round())

    # Find plasma transfusions
    ffp_in <- input_events %>% 
      filter(itemid == 220970)
    
    # Find RBC transfusions
    rbc_in <- input_events %>% 
      filter(itemid == 225168)
    
    # Add to total fluids
    fluids_in %<>% rbind(ffp_in, rbc_in)
    
  }else # Extract inputs for carevue
  {
    # Read input events
    inputfile <- paste("data/events/inputeventscv_", icustay_df$hadm_id, ".csv", sep = "")
    input_events <- read_csv(inputfile)
    
    # Select bolus items and round
    fluids_in <- input_events %>% 
      # Colloid bolus
      filter(itemid %in% c(30008, 30009, 42832, 40548,
                           45403, 44203, 30181, 46564,
                           43237, 43353, 30012, 46313,
                           30011, 30016, 42975, 42944,
                           46336, 46729, 40033, 45410,
                           42731)) %>% 
      filter(amount > 100,
             amount < 2000)
    
    fluids_in %<>% 
      rbind(input_events %>% 
        # Crystalloid bolus
        filter(itemid %in% c(
          30015, 30018, 30020, 30021,
          30058, 30060, 30061, 30063,
          30065, 30143, 30159, 30160,
          30169, 30190, 40850, 41491,
          42639, 42187, 43819, 41430,
          40712, 44160, 42383, 42297,
          42453, 40872, 41915, 41490,
          46501, 45045, 41984, 41371,
          41582, 41322, 40778, 41896,
          41428, 43936, 44200, 41619,
          40424, 41457, 41581, 42844,
          42429, 41356, 40532, 42548,
          44184, 44521, 44741, 44126,
          44110, 44633, 44983, 44815,
          43986, 45079, 46781, 45155,
          43909, 41467, 44367, 41743,
          40423, 44263, 42749, 45480,
          44491, 41695, 46169, 41580,
          41392, 45989, 45137, 45154,
          44053, 41416, 44761, 41237,
          44426, 43975, 44894, 41380,
          42671))) %>% 
      filter(!is.na(amount),
             amount > 248,
             amount <= 2000,
             amountuom == "ml") %>%
      mutate(amount = round(amount))
    
    # Extract colloids from charted events
    charted_fluids <- chart_df %>% 
      filter(ITEMID %in% c(2510, 3087, 6937,
                           3087, 3088)) %>% 
      filter(is.na(VALUENUM)) %>%
      mutate(VALUENUM = round(VALUENUM))
    
    # Find plasma transfusions
    ffp_in <- input_events %>% 
      filter(itemid %in% c(30005, 30180)) %>% 
      filter(amount > 0)
    
    # Find RBC transfusions
    rbc_in <- input_events %>% 
      filter(itemid %in% c(30179, 30001, 30004)) %>% 
      filter(amount > 0)
    
    # Add to total fluids
    fluids_in %<>% 
      rbind(ffp_in, rbc_in) %>% 
      select(charttime, amount) %>% 
      rename("starttime" = charttime) %>% 
      rbind(
        charted_fluids %>% 
          select(CHARTTIME, VALUENUM) %>% 
          rename("starttime" = CHARTTIME,
                 "amount" = VALUENUM)
      )
    
  }
  
  # Read output events
  outputfile <- paste("data/events/outputevents_", icustay_df$hadm_id, ".csv", sep = "")
  if(file.exists(outputfile))
  {
    output_events <- read_csv(outputfile, col_types = c("ddddTddcTdlll"))
    
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
      filter(!is.na(value)) %>% 
      mutate(balance = cumsum(value))
    
    # Output max
    return(max(fluid_balance_df$balance))
  }else
  {
    return(NA)
  }

}

# EXAMPLE EVENTS
# 100012 - metavision, fresh frozen plasma
# 101019 - carevue, ffp
# 100295 - carvue, fluid bolus
