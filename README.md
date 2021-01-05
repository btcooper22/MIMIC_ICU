# MIMIC_ICU
Code to define and compare existing ICU readmission risk scores within the MIMIC-III database. This consists of four core scripts.

## extract_patients

This script uses the MIMIC-III database and filters by eligibility criteria to produce a dataset of `usable' ICU patients, which can be redefined at any time. For data security reasons, this dataset file will not be uploaded to github. Currently extracts patients who were admitted under or moved to a surgical service.

## preprocess_data

For all patients defined in *extract_patients* this script pre-processes all MIMIC data concerning their demographic information and ICU stay. This process largely follows the workflow laid out by Lin et al. 2018 [^1] and available in Python code at https://github.com/Jeffreylin0925/MIMIC-III_ICU_Readmission_Analysis and https://github.com/YerevaNN/mimic3-benchmarks. It outputs all demographic, stay/transfer and input (diagnosis/procedure/prescription) data for selected patients. Events data (chart, lab or output) are not processed, as they are both very large, and not needed in their entirity for subsequent analyses. This will be added in if needed.

## define_variables

Independent of which patients are included, this script defines and calculates the input and outcome variables required by each of the three risk scores, and outputs this to another file, also not uploaded.

## compare_scores

This script takes the variables calculated above, and uses them to generate the three risk scores. Scores are then compared in the standard method.

[^1]: Lin Y-W, Zhou Y, Faghri F, Shaw MJ, Campbell RH (2019) Analysis and prediction of unplanned intensive care unit readmission using recurrent neural networks with long short-term memory. PLoS ONE 14(7): e0218942. https://doi.org/10.1371/journal.pone.0218942