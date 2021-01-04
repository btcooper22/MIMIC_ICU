require(DBI)
require(dplyr)
print(noquote("Loading MIMIC-III dataset..."))

con <- DBI::dbConnect(RPostgres::Postgres(),
                      user = "postgres",
                      password = "postgres",
                      dbname = "mimic")
dbSendQuery(con, "set search_path to mimiciii;")

# Load mimic tables
table_list <- c("admissions", "callout", "caregivers", "chartevents", "cptevents",
                "d_cpt", "d_icd_diagnoses", "d_icd_procedures", "d_items", "d_labitems",
                "datetimeevents", "diagnoses_icd", "drgcodes", "icustays", "inputevents_cv",
                "inputevents_mv", "labevents", "microbiologyevents", "noteevents", "outputevents",
                "patients", "prescriptions", "procedureevents_mv", "procedures_icd", "services",
                "transfers")
# Create blank list
mimic <- list()

# Loading loop
for(i in 1:length(table_list))
{
  # Generate name
  obj_name <- paste("mimic", table_list[i], sep = "_")
  
  # Pull table
  loaded_table <- tbl(con, table_list[i])
  
  # Insert into list
  mimic[[i]] <- loaded_table
  
  # Rename
  names(mimic)[i] <- table_list[i]
}
rm(i, obj_name, table_list, 
   loaded_table, con)