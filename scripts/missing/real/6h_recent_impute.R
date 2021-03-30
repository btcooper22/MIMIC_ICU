# Packages
require(dplyr)
require(magrittr)
require(readr)
require(foreach)
require(stringr)
require(tidyr)
require(tibble)
require(lubridate)
require(doParallel)
require(tools)

# Functions
implausability_filter <- function(val, good_range)
{
  val[val < good_range[1]] <- NA
  val[val > good_range[2]] <- NA
  return(val)
}

worst_value <- function(val, bad_direction)
{
  if(all(is.na(val)))
  {
    return(NA)
  }else 
  {
    if(bad_direction == "max")
    {
      return(max(val, na.rm = TRUE))
    }else if(bad_direction == "min")
    {
      return(min(val, na.rm = TRUE))
    } else
    {
      return(val[which.max(abs(val - as.numeric(bad_direction)))])
    }
  }
}

# Load data
apache_df <- read_csv("data/apache_discharge_missing.csv") %>% 
  filter(is.na(apache_II))
apache_scores <- apache_df %>% 
  select(-row_id, -adm_id,
         -mort_30, -age,
         -chronic, -fractioninspiredoxygen,
         -apache_II, -missing_abg, -acute_renal_failure)
apache_additional <- apache_df %>% 
  select(row_id, adm_id,
         mort_30, age,
         chronic, fractioninspiredoxygen,
         apache_II, missing_abg, acute_renal_failure)
rm(apache_df)

mimic_icustays <- read_rds("data/mimic_preprocessed_missing.RDS")$stays

# Preparation----

# Create database list
db_list <- data.frame(
  variable = names(apache_scores),
  source = c(
    rep("chart", 6),
    rep("lab", 5),
    rep("chart", 3),
    "lab",
    rep("chart", 2)
  )
)

# Create and join dictionary
itemid <- data.frame(
  variable = "temperature",
  itemid = c(676, 677, 678, 679, 223761, 223762),
  min = 30, max = 40, worst = "37.2"
) %>% 
  rbind(
    data.frame(variable = "systolicbp",
               itemid = c(6, 51, 455, 6701, 220050, 220179),
               min = 70, max = 220, worst = "max"),
    data.frame(variable = "diastolicbp",
               itemid = c(8364, 8368, 8441, 
                          8555, 220051, 220180),
               min = 40, max = 140, worst = "max"),
    data.frame(variable = "pulse",
               itemid = c(211, 220045),
               min = 20, max = 200, worst = "90"),
    data.frame(variable = "respiratory",
               itemid = c(614, 615, 618, 219, 619, 653, 1884,
                          8113, 1635, 3603, 224688, 224689,
                          224690, 220210),
               min = 0, max = 60, worst = "max"),
    data.frame(variable = "arterialpH",
               itemid = c(1126, 780,4753, 223830),
               min = 7, max = 8, worst = "7.55"),
    data.frame(variable = "sodium",
               itemid = c(50824, 50983),
               min = 100, max = 160, worst = "140"),
    data.frame(variable = "potassium",
               itemid = c(50822, 50971),
               min = 0, max = 8.5, worst = "4.5"),
    data.frame(variable = "creatinine",
               itemid = c(50912),
               min = 0, max = 6, worst = "max"),
    data.frame(variable = "haematocrit",
               itemid = c(51221, 50810),
               min = 15, max = 60, worst = "44"),
    data.frame(variable = "whitebloodcount",
               itemid = c(51300, 51301),
               min = 0, max = 50, worst = "max"),
    data.frame(variable = "glasgowcomaeye",
               itemid = c(184, 220739),
               min = 0, max = 4, worst = "min"),
    data.frame(variable = "glasgowcomamotor",
               itemid = c(454, 223901),
               min = 0, max = 6, worst = "min"),
    data.frame(variable = "glasgowcomaverbal",
               itemid = c(723, 223900),
               min = 0, max = 6, worst = "min"),
    data.frame(variable = "bicarbonate",
               itemid = c(50803, 50882),
               min = 0, max = 60, worst = "27"),
    data.frame(variable = "arterialoxygen",
               itemid = c(4203, 3785, 3837, 3828, 220224, 779),
               min = 20, max = 250, worst = "min"),
    data.frame(variable = "arterialcarbon",
               itemid = c(220235, 777, 778, 3784, 3835, 3836),
               min = 8, max = 200, worst = "40")
  )

db_list %<>% 
  left_join(itemid)

# Main loop----

# Prepare parallel
n_cores <- 15
cl <- makeCluster(ifelse(detectCores() <= n_cores,
                         detectCores() - 1,
                         n_cores))
registerDoParallel(cl)

imputed_results <- foreach(i = 1:nrow(apache_scores), .combine = "rbind",
                           .packages = c("dplyr", "tibble",
                                         "lubridate", "foreach")) %dopar%
  {
    # Flag which missing
    missing_vars <- names(apache_scores)[is.na(apache_scores[i,])] %>% 
      enframe() %>% 
      rename(variable = "value") %>% 
      left_join(db_list)
    
    # Find ICU stay ID
    in_out_times <- mimic_icustays %>% 
      filter(hadm_id == apache_additional$adm_id[i]) %>% 
      select(intime, outtime,
             icustay_id) %>% 
      slice_min(outtime) 
    
    # Extract discharge time
    discharge_time <- in_out_times %>% 
      select(outtime) %>% deframe() %>% force_tz("UTC")
    
    # Loop through chart variables
    if("chart" %in% missing_vars$source)
    {
      # Load chart data
      chart_df <- read.csv(paste("data/events/chartevents_",
                                 apache_additional$adm_id[i],
                                 ".csv", sep = ""))
      # Loop
      varnames <- unique(missing_vars %>% 
                           filter(source == "chart") %>% 
                           select(variable) %>% deframe())
      imputed_vars <- foreach(v = 1:length(varnames), .combine = "rbind") %do%
        {
          # Extract variable information
          var_df <- missing_vars %>% 
            filter(variable == varnames[v])
          
          # Filter chart
          chart_filtered <- chart_df %>% 
            # Select item of interest
            filter(ITEMID %in% var_df$itemid) %>% 
            # Filter to before discharge
            filter(CHARTTIME < discharge_time) %>% 
            # Correct temp if in F
            mutate(VALUENUM = ifelse(ITEMID %in% c(678, 679, 223761),
                                  (VALUENUM - 32) * 5/9,
                                  VALUENUM)) %>% 
            # Apply implausability filter
            mutate(VALUENUM = implausability_filter(VALUENUM, c(var_df$min[1],
                                                                var_df$max[1]))) %>% 
            # Sort by date 
            arrange(CHARTTIME)
          
          # Extract measurements
          chart_measurements <- chart_filtered %>% 
            select(VALUENUM) %>% 
              deframe() %>% rev()
          
          # Count measurements
          n_measure <- length(chart_measurements)
          
          # Extract values
          data.frame(variable = varnames[v],
                     recent = chart_measurements[!is.na(chart_measurements)][1],
                     mean5 = mean(chart_measurements[1:min(c(n_measure, 5))], na.rm = T),
                     median5 = median(chart_measurements[1:min(c(n_measure, 5))], na.rm = T),
                     worst5 = worst_value(chart_measurements[1:min(c(n_measure, 5))],
                                          var_df$worst[1]),
                     mean10 = mean(chart_measurements[1:min(c(n_measure, 10))], na.rm = T),
                     median10 = median(chart_measurements[1:min(c(n_measure, 10))], na.rm = T),
                     worst10 = worst_value(chart_measurements[1:min(c(n_measure, 10))],
                                          var_df$worst[1]),
                     meantotal = mean(chart_measurements, na.rm = T),
                     mediantotal = median(chart_measurements, na.rm = T),
                     worsttotal = worst_value(chart_measurements,var_df$worst[1]))
        }
      
      # Also extract FiO2 for each case
      fio2_filtered <- chart_df %>% 
        filter(ITEMID %in% c(223835, 2981, 3420)) %>% 
        # Filter to before discharge
        filter(CHARTTIME < discharge_time) %>% 
        # Apply implausability filter
        mutate(VALUENUM = implausability_filter(VALUENUM, c(21, 100))) %>% 
        # Sort by date 
        arrange(CHARTTIME) %>% 
        mutate(VALUENUM = ifelse(is.na(VALUENUM), 21, VALUENUM))
      
      # Extract measurements
      fio2_measurements <- fio2_filtered %>% 
        select(VALUENUM) %>% 
        deframe() %>% rev()
      
      # Count measurements
      n_measure <- length(fio2_measurements)
      
      # Extract values
      output_chart <- rbind(imputed_vars,
            data.frame(variable = "fractioninspiredoxygen",
                       recent = fio2_measurements[!is.na(fio2_measurements)][1],
                       mean5 = mean(fio2_measurements[1:min(c(n_measure, 5))], na.rm = T),
                       median5 = median(fio2_measurements[1:min(c(n_measure, 5))], na.rm = T),
                       worst5 = worst_value(fio2_measurements[1:min(c(n_measure, 5))],
                                            "max"),
                       mean10 = mean(fio2_measurements[1:min(c(n_measure, 10))], na.rm = T),
                       median10 = median(fio2_measurements[1:min(c(n_measure, 10))], na.rm = T),
                       worst10 = worst_value(fio2_measurements[1:min(c(n_measure, 10))],
                                             "max"),
                       meantotal = mean(fio2_measurements, na.rm = T),
                       mediantotal = median(fio2_measurements, na.rm = T),
                       worsttotal = worst_value(fio2_measurements, "max")))
    }
    
    # Loop through lab variables
    if("lab" %in% missing_vars$source)
    {
      # Load lab data
      lab_df <- read.csv(paste("data/events/labevents_",
                                 apache_additional$adm_id[i],
                                 ".csv", sep = ""))
      # Loop
      varnames <- unique(missing_vars %>% 
                           filter(source == "lab") %>% 
                           select(variable) %>% deframe())
      imputed_vars <- foreach(v = 1:length(varnames), .combine = "rbind") %do%
        {
          # Extract variable information
          var_df <- missing_vars %>% 
            filter(variable == varnames[v])
          
          # Filter lab
          lab_filtered <- lab_df %>% 
            # Select item of interest
            filter(itemid %in% var_df$itemid) %>% 
            # Filter to before discharge
            filter(charttime < discharge_time) %>% 
            # Apply implausability filter
            mutate(valuenum = implausability_filter(valuenum, c(var_df$min[1],
                                                                var_df$max[1]))) %>% 
            # Sort by date 
            arrange(charttime)
          
          # Flip NA
          lab_filtered$valuenum[is.na(lab_filtered$valuenum)] <- 
            as.numeric(lab_filtered$value)[is.na(lab_filtered$valuenum)]

          # Extract measurements
          lab_measurements <- lab_filtered %>% 
            select(valuenum) %>% 
            deframe() %>% rev()
          
          # Count measurements
          n_measure <- length(lab_measurements)
          
          # Extract values
          data.frame(variable = varnames[v],
                     recent = lab_measurements[!is.na(lab_measurements)][1],
                     mean5 = mean(lab_measurements[1:min(c(n_measure, 5))], na.rm = T),
                     median5 = median(lab_measurements[1:min(c(n_measure, 5))], na.rm = T),
                     worst5 = worst_value(lab_measurements[1:min(c(n_measure, 5))],
                                          var_df$worst[1]),
                     mean10 = mean(lab_measurements[1:min(c(n_measure, 10))], na.rm = T),
                     median10 = median(lab_measurements[1:min(c(n_measure, 10))], na.rm = T),
                     worst10 = worst_value(lab_measurements[1:min(c(n_measure, 10))],
                                           var_df$worst[1]),
                     meantotal = mean(lab_measurements, na.rm = T),
                     mediantotal = median(lab_measurements, na.rm = T),
                     worsttotal = worst_value(lab_measurements,var_df$worst[1]))
        }
      
      output_lab <- imputed_vars
    }
    
    # Combine
    if("chart" %in% missing_vars$source &
       "lab" %in% missing_vars$source)
    {
      output <- rbind(output_chart, output_lab)
    }else if("chart" %in% missing_vars$source)
    {
      output <- output_chart
    }else
    {
      output <- output_lab
    }
    
    # Add additional metadata
    output %>% 
      mutate(adm = apache_additional$adm_id[i],
             row_id = i) %>% 
      relocate(row_id, adm)
  }
stopCluster(cl)

# Reformulate----

# Load normal data
rm(apache_additional, apache_scores)
source("functions/data_loader_splitter.R")
apache_df <- cbind(apache_additional, apache_scores)

# Pivot imputed data
imputed_results %<>% 
  pivot_longer(4:13, "method")

# Loop through methods
methods_list <- unique(imputed_results$method)

results <- foreach(i = 1:length(methods_list)) %do%
  {
    print(i)
    # Create subset
    imputed_subset <- imputed_results %>% 
      filter(method == methods_list[i])
    
    # Identify individuals
    individuals <- unique(imputed_subset$adm)
    
    # Set aside dataset
    new_df <- apache_df
    
    # Loop through individuals
    for(j in 1:nrow(new_df))
      {
      #print(j)
        # Check if imputation done
        if(new_df$adm_id[j] %in% individuals)
        {
          # Extract individual imputations
          ind_imp <- imputed_subset %>% 
            filter(adm == new_df$adm_id[j])
          
          # Fill values
          for(k in 1:nrow(ind_imp))
          {
            new_df[j,which(names(new_df) == ind_imp$variable[k])] <- ind_imp$value[k]
          }
        }
    }
    new_df %>%
      select(-row_id, -adm_id,
             -mort_30, -age,
             -chronic, -fractioninspiredoxygen,
             -apache_II, -missing_abg, -acute_renal_failure)
  }

names(results) <- paste("recent", methods_list, sep = "_")

write_rds(results, "data/impute_discharge/recent.RDS",
          compress = "gz")

