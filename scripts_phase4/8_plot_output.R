# Packages
require(readr)
require(magrittr)
require(dplyr)

# Load data
data_inunit <- read_csv("data/apache_real_missing.csv")
data_30d <- read_csv("data/apache_discharge_missing.csv")

# Calculate mortality rates
sum(data_inunit$mort_inunit)
round(mean(data_inunit$mort_inunit) * 100, 2)

sum(data_30d$mort_30)
round(mean(data_30d$mort_30) * 100, 2)

# Measure missingess
data_inunit %>% 
  summarise(sum = sum(is.na(apache_II)),
            mean = mean(is.na(apache_II)) * 100)

data_30d %>% 
  summarise(sum = sum(is.na(apache_II)),
            mean = mean(is.na(apache_II)) * 100)

# Calibration and discrimination----
