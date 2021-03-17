# Load packages
require(tidyverse)
require(magrittr)
require(lubridate)
require(doParallel)
require(tools)

# Load mimic dataset
source("functions/mimic_load.R")

# Load outcomes
outcomes <- read_csv("data/outcomes_mortality.csv")

# chartevents----------

header <- names(read_csv("C:/Users/benco/charts/chartevents_aa", n_max = 1))
files <- dir("C:/Users/benco/charts/", pattern = "_")

# Make database of hadm_ids in which split
if(!file.exists("data/chart_database_mortality.csv"))
{
  # Prepare parallel options
  psnice(value = 19)
  registerDoParallel(ifelse(detectCores() <= 12,
                            detectCores() - 1,
                            12)
  )

  chart_database <- foreach(i = 1:length(files), .combine = "rbind",
                            .packages = "readr") %dopar%
    {
      # Load in chart
      print(files[i])
      chart_search <- read_csv(paste("C:/Users/benco/charts/", files[i], sep = ""),
                               col_types = "ccccccccccccccc", col_names = header)
      
      # List admissions
      adm_list <- unique(chart_search$HADM_ID)
      
      # Remove files
      rm(chart_search)
      gc()
      
      # Add to output
      data.frame(chart_split = files[i], hadm_id = adm_list)
    }
  write_csv(chart_database, "data/chart_database_mortality.csv")
  stopImplicitCluster()
  gc()
}else
{
  chart_database <- read_csv("data/chart_database_mortality.csv")
}

# Filter outcomes by missing chart data
outcomes %<>% 
  filter(outcomes$hadm_id %in% chart_database$hadm_id == TRUE)
write_csv(outcomes, "data/outcomes_mortality.csv")

# Create filenames
chartevents_filenames <- data.frame(filename = paste("data/events/chartevents_", outcomes$hadm_id,
                                                   ".csv", sep = ""), ID = 1:nrow(outcomes)) %>% 
  mutate(done = file.exists(filename)) %>% 
  filter(done == FALSE)


# Prepare parallel options
psnice(value = 19)
registerDoParallel(ifelse(detectCores() <= 10,
                          detectCores() - 1,
                          10)
)

# Extract chart events
print(noquote("Extracting chart events"))
foreach(i = chartevents_filenames$ID, .packages = c("readr", "dplyr",
                               "magrittr", "foreach")) %dopar%
{
  # Find subject info
  print(i)
  subj <- outcomes$subject_id[i]
  adm <- outcomes$hadm_id[i]
  
  # Determine filename
  filename <- paste("data/events/chartevents_", outcomes$hadm_id[i],
                    ".csv", sep = "")
  
  if(!file.exists(filename))
  {
    # Find correct chart split
    chart_n <- chart_database$chart_split[chart_database$hadm_id==adm]
    
    # Load chart split and filter
    chart_df <- foreach(cn = 1:length(chart_n), .combine = "rbind") %do% 
      {
        print(cn)
        chart_in <- read_csv(paste("C:/Users/benco/charts/", chart_n[cn], sep = ""),
                 col_types = "ccccccccccccccc", col_names = header)
        
        # Subset
        chart_subset <- chart_in %>% 
          filter(HADM_ID == adm)
        
        # Clear
        rm(chart_in)
        gc()
        
        # Output
        chart_subset
      }

    # Write to file
    chart_df %>% 
      write_csv(filename)
    
    rm(chart_df)
    gc()
  }
}
stopImplicitCluster()

# labevents------------

# Load all lab events
labevents <- mimic$labevents %>% 
  collect()

# Prepare parallel options
psnice(value = 19)
registerDoParallel(ifelse(detectCores() <= 6,
                          detectCores() - 1,
                          6)
)

# Create filenames
labevents_filenames <- data.frame(filename = paste("data/events/labevents_", outcomes$hadm_id,
                  ".csv", sep = ""), ID = 1:nrow(outcomes)) %>% 
  mutate(done = file.exists(filename)) %>% 
  filter(done == FALSE)

print(noquote("Extracting lab events"))
foreach(i = labevents_filenames$ID, .packages = c("dplyr","magrittr",
                                 "readr")) %dopar%
{
 # Find labevents
  labs <- labevents %>% 
    # Select admission
    filter(hadm_id == outcomes$hadm_id[i]) 
  
  # Write to file
  filename <- paste("data/events/labevents_", outcomes$hadm_id[i],
                    ".csv", sep = "")
  if(!file.exists(filename))
  {
    labs %>% 
      write_csv(filename)
  }
}
stopImplicitCluster()
rm(labevents)
gc()


# inputevents_mv-------------
inputevents_mv <- mimic$inputevents_mv %>% 
  collect()

# Prepare parallel options
psnice(value = 19)
registerDoParallel(ifelse(detectCores() <= 6,
                          detectCores() - 1,
                          6)
)

# Create filenames
inputevents_mv_filenames <- data.frame(filename = paste("data/events/inputeventsmv_", outcomes$hadm_id,
                                                   ".csv", sep = ""), ID = 1:nrow(outcomes)) %>% 
  mutate(done = file.exists(filename)) %>% 
  filter(done == FALSE)


print(noquote("Extracting metavision input events"))
foreach(i = inputevents_mv_filenames$ID, .packages = c("dplyr","magrittr",
                                            "readr")) %dopar%
  {
    # Find labevents
    inputs <- inputevents_mv %>% 
      # Select admission
      filter(hadm_id == outcomes$hadm_id[i]) 
    
    if(nrow(inputs) > 0)
    {
      # Write to file
      filename <- paste("data/events/inputeventsmv_", outcomes$hadm_id[i],
                        ".csv", sep = "")
      if(!file.exists(filename))
      {
        inputs %>% 
          write_csv(filename)
      }
    }
  }
stopImplicitCluster()
rm(inputevents_mv)
gc()

# inputevents_cv-------------
inputevents_cv <- mimic$inputevents_cv %>% 
  collect()

# Prepare parallel options
psnice(value = 19)
registerDoParallel(ifelse(detectCores() <= 6,
                          detectCores() - 1,
                          6)
)

# Create filenames
inputevents_cv_filenames <- data.frame(filename = paste("data/events/inputeventscv_", outcomes$hadm_id,
                                                        ".csv", sep = ""), ID = 1:nrow(outcomes)) %>% 
  mutate(done = file.exists(filename)) %>% 
  filter(done == FALSE)


print(noquote("Extracting carevue input events"))
foreach(i = inputevents_cv_filenames$ID, .packages = c("dplyr","magrittr",
                                            "readr")) %dopar%
  {
    # Find labevents
    inputs <- inputevents_cv %>% 
      # Select admission
      filter(hadm_id == outcomes$hadm_id[i]) 
    
    if(nrow(inputs) > 0)
    {
      # Write to file
      filename <- paste("data/events/inputeventscv_", outcomes$hadm_id[i],
                        ".csv", sep = "")
      if(!file.exists(filename))
      {
        inputs %>% 
          write_csv(filename)
      }
    }
  }
stopImplicitCluster()
rm(inputevents_cv)
gc()

# outputevents -------------
outputevents <- mimic$outputevents %>% 
  collect()

# Prepare parallel options
psnice(value = 19)
registerDoParallel(ifelse(detectCores() <= 12,
                          detectCores() - 1,
                          12)
)

# Create filenames
outputevents_filenames <- data.frame(filename = paste("data/events/outputevents_", outcomes$hadm_id,
                                                        ".csv", sep = ""), ID = 1:nrow(outcomes)) %>% 
  mutate(done = file.exists(filename)) %>% 
  filter(done == FALSE)


print(noquote("Extracting output events"))
foreach(i = outputevents_filenames$ID, .packages = c("dplyr","magrittr",
                                            "readr")) %dopar%
  {
    # Find labevents
    outputs <- outputevents %>% 
      # Select admission
      filter(hadm_id == outcomes$hadm_id[i]) 
    
    if(nrow(outputs) > 0)
    {
      # Write to file
      filename <- paste("data/events/outputevents_", outcomes$hadm_id[i],
                        ".csv", sep = "")
      if(!file.exists(filename))
      {
        outputs %>% 
          write_csv(filename)
      }
    }
  }
stopImplicitCluster()
rm(outputevents)
gc()

