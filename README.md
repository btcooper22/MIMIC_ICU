# MIMIC_ICU
Code to define and compare existing ICU readmission risk scores within the MIMIC-III database. This consists of three core scripts.

## extract_patients

This script uses the MIMIC-III database and filters by eligibility criteria to produce a dataset of `usable' ICU patients, which can be redefined at any time. For data security reasons, this dataset file will not be uploaded to github.

## define_variables

Independant of which patients are included, this script defines and calculates the variables required by each of the three risk scores, and outputs this to another file, also not uploaded.

## compare_scores

This script takes the variables calculated above, and uses them to generate the three risk scores. Scores are then compared in the standard method.