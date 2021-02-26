# Packages
require(dplyr)
require(magrittr)
require(readr)
require(foreach)
require(caret)
require(stringr)
require(RANN)

# Load data-------
source("functions/data_loader_splitter.R")

# Set number of K to profile
k_list <- 2:15

# Loop through K
imputed_dfs <- foreach(k = 1:length(k_list),
                       .errorhandling = "remove") %do%
  {
    print(k)
    
    # Prepare KNN imputation
    KNN_preproc <-
      preProcess(apache_scores %>% as.data.frame(),
                 method = "knnImpute",
                 k = k_list[k],
                 knnSummary = mean)
    
    # Perform imputation
    KNN_results <- KNN_preproc %>%
      predict(apache_scores,
              na.action = na.pass)
    
    # De-normalise data
    procNames <- data.frame(
      col = names(KNN_preproc$mean),
      mean = KNN_preproc$mean,
      sd = KNN_preproc$std
    )
    
    for (i in procNames$col) {
      KNN_results[i] <-
        KNN_results[i] * KNN_preproc$std[i] + KNN_preproc$mean[i]
    }
    
    # Return imputed data
    list(k = k_list[k],
         KNN_results)
  }

# Disentangle and write
all_KNN <- foreach(i = 1:length(imputed_dfs)) %do%
  {
    imputed_dfs[[i]][[2]]
  }

names(all_KNN) <- paste("KNN", k_list, sep = "_")

write_rds(all_KNN, "data/impute/KNN.RDS",
          compress = "gz")