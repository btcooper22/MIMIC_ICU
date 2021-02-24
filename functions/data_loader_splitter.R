apache_df <- read_csv("data/apache_real_missing.csv")
apache_scores <- apache_df %>% 
  select(-row_id, -adm_id,
         -mort_inunit, -age,
         -chronic, -fractioninspiredoxygen,
         -apache_II, -missing_abg)
apache_additional <- apache_df %>% 
  select(row_id, adm_id,
         mort_inunit, age,
         chronic, fractioninspiredoxygen,
         apache_II, missing_abg)
rm(apache_df)