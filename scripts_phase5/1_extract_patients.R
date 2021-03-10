# Packages
require(bigrquery)
require(DBI)
require(dplyr)

# Find project name
bq_projects()

# Find datasets
list_datasets("prd-proj-decovid-358375")

# Find tables
list_tables("prd-proj-decovid-358375","icnarc_analytics")

# Set ICNARC table as connection
con <- dbConnect(
  bigrquery::bigquery(),
  project = "prd-proj-decovid-358375",
  dataset = "icnarc_analytics"
)
dbListTables(con)

# Read in ICNARC table
icnarc <- tbl(con, "ICNARC") %>% collect()
names(icnarc)
