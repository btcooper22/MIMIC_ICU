# Packages
require(dplyr)
require(magrittr)
require(readr)
require(tibble)
require(foreach)
require(stringr)
require(doParallel)
require(tools)
require(tidyr)
require(forcats)

# Functions
source("functions/filename_metadata.R")
source("functions/make_pool.R")
source("functions/parallel_setup.R")
source("functions/calculate_apache_score.R")

# Load data
full_data <- read_csv("data/impute/complete_cases.csv") 

# Load MIMIC database
source("functions/mimic_load.R")

# Find admissions data
admissions <- mimic$admissions %>% 
  select(hadm_id, diagnosis,
         admission_type) %>% 
  collect() %>% 
  rename(adm_id = "hadm_id")

# Join to database
full_data %<>% 
  left_join(admissions)

# Create combinations list
pooling_vars <- list(
  c("admission_type", "apache_chronic_discharge",
    "apache_age_discharge"),
  c("admission_type", "apache_age_discharge"),
  c("admission_type", "apache_chronic_discharge"),
  c("apache_age_discharge", "apache_chronic_discharge"),
  "admission_type",
  "apache_age_discharge",
  "apache_chronic_discharge"
)

# Main loop through files-------

# Load model for prediction
apache_model <- read_rds("models/apache_model.RDS")

# Locate files
files <- dir("data/MCAR/", full.names = TRUE)

# Set up parallel
ptm <- proc.time()
parallel_setup(14)

output <- foreach(f = 1:length(files), .packages = c("dplyr", "tidyr","stringr",
                                         "tibble","magrittr", "forcats",
                                         "foreach", "readr")) %dopar%
  {
    # Extract metadata
    metadata <- filename_metadata(files[f])
    
    # Create filename
    filename <- paste("data/impute/hotdeck/hotdeck_S",
                      str_pad(metadata[2], 4, "right", 0),"_N",
                      str_pad(metadata[1],3,"left",0), ".RDS", sep = "")
    
    if(!file.exists(filename))
    {
    # Load file
    mcar_df <- read.csv(files[f]) %>% 
    # Join pooling data
    cbind(full_data %>% 
            select(admission_type))
    
    # Find NA columns
    mcar_NA <- mcar_df %>% 
      select_if(function(x) any(is.na(x))) %>% 
      # Add identifying data again
      cbind(mcar_df[,13:16])
    
    # Loop through NA columns
    imputed_NA <- foreach(i = 1:(ncol(mcar_NA)-4)) %do%
      {
        # Find column name
        col_name <- names(mcar_NA)[i]
        
        # Identify data
        NA_variable <- subset(mcar_NA, is.na(mcar_NA[,i]))
        
        # Identify complete cases
        present_variable <- subset(mcar_NA, !is.na(mcar_NA[,i]))
      
        # Loop through pooling options
        replace_values <- foreach(j = 1:length(pooling_vars)) %do%
        {
          # Debug
          print(paste(i,j))
          
          # Generate unique pooling groups
          if(length(pooling_vars[[j]]) > 1)
          {
            NA_variable$pooling_group_var <- apply(NA_variable[,pooling_vars[[j]]],
                                                   1, paste, collapse = "_" )
          }else
          {
            NA_variable$pooling_group_var <- NA_variable[,pooling_vars[[j]]]
          }
          
          # Loop through pool groups
          foreach(k = 1:length(unique(NA_variable$pooling_group_var)),
                  .combine = "c") %do%
            {
              # Extract identifier
              identifier <- unique(NA_variable$pooling_group_var)[k]
              
              # Generate candidate pool
              candidate_pool <- make_pool(pooling_vars[[j]],
                                          NA_variable %>% 
                                            filter(pooling_group_var == identifier),
                                          present_variable)
              
              # Sample from candidate pool
              sample(candidate_pool[,col_name],
                     sum(NA_variable$pooling_group_var == identifier),
                     replace = TRUE)
            }
        }
        
        # Replace NAs for each method
        imputed_col <- foreach(j = 1:length(replace_values)) %do%
          {
            fullcol <- mcar_NA[,i]
            fullcol[is.na(fullcol)] <- replace_values[[j]]
            fullcol
          }
        
        # Output
        imputed_col
      }
    
    # Reattach
    full_imputed <- foreach(i = 1:length(pooling_vars),
                            .combine = "cbind") %do%
      {
        full_df <- foreach(j = 1:(ncol(mcar_df)-4),
                .combine = "cbind") %do%
          {
            # i is methods, j is columns
            if(any(is.na(mcar_df[,j])))
            {
              # Replace with imputed
              col_values <- imputed_NA[[j]][[i]]
            }else
            {
              # Do not replace
              col_values <- mcar_df[,j]
            }
            
            # Output
            col_values
          }
        
        # Attach remaining values
        full_df %<>% as.data.frame() %>% 
          cbind(mcar_df[,13:15])
        names(full_df) <- names(mcar_df)[1:15]
        
        # Sum and predict
        probs <-  full_df %>% 
          mutate(apache_II_discharge = calculate_apache_scores(.)) %>% 
          predict(apache_model, newdata = .)
        
        # Output
        probs
      }
    
    # Create output
    output <- list(pooling = pooling_vars,
         probs = full_imputed)
    
    # Write to file
    write_rds(output, filename,
              compress = "gz")
    }else
    {
      output <- read_rds(filename)
    }
    
    # Output
    output$N <- metadata[1]
    output$split <- metadata[2]
    output
  }
stopImplicitCluster() 
proc.time() - ptm # 54 min

# Compile all
write_rds(output, "data/impute/hotdeck.RDS",
          compress = "gz")
