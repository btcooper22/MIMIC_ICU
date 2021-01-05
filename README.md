# MIMIC_ICU
Code to define and compare existing ICU readmission risk scores within the MIMIC-III database. This consists of five core modules.

## extract_patients

This module uses the MIMIC-III database and filters by eligibility criteria to produce a dataset of `usable' ICU patients, which can be redefined at any time. For data security reasons, this dataset file will not be uploaded to github. Currently extracts patients who were admitted under or moved to a surgical service.

## preprocess_data

For all patients defined in *extract_patients* this module pre-processes all MIMIC data concerning their demographic information and ICU stay. This process largely follows the workflow laid out by Lin et al. 2018 [^1] and available in Python code at https://github.com/Jeffreylin0925/MIMIC-III_ICU_Readmission_Analysis and https://github.com/YerevaNN/mimic3-benchmarks. It outputs all demographic, stay/transfer and input (diagnosis/procedure/prescription) data for selected patients. Events data (chart, lab or output) are not processed, as they are both very large, and not needed in their entirity for subsequent analyses. This will be added in if needed.

## define_outcomes

From the data, this module define the measured outcomes needed for the three risk scores. The scores of Hammer et al. 2020 [^2], Martin et al. 2019 [^3] and Frost et al. 2010 [^4] all use readmission within the same hospital admission, which I restrict to 'within the same hospitalisation event as the surgery', excluding prior or subsequent hospitalisations.

## extract_predictors

For all patients specified in earlier modules, this module extracts the predictors needed for each risk score. Patients with missing data for some scores but not others are still included - this is dealt with in the *compare_scores* module.

## compare_scores

This module takes the variables calculated above, and uses them to generate the three risk scores. Scores are then compared in the standard method.

[^1]: Lin Y-W, Zhou Y, Faghri F, Shaw MJ, Campbell RH (2019) Analysis and prediction of unplanned intensive care unit readmission using recurrent neural networks with long short-term memory. PLoS ONE 14(7): e0218942. https://doi.org/10.1371/journal.pone.0218942

[^2]: Hammer M, Grabitz SD, Teja B, (2020) A Tool to Predict Readmission to the Intensive Care Unit in Surgical Critical Care Patients—The RISC Score. Journal of Intensive Care Medicine  https://doi.org/10.1177/0885066620949164

[^3]: Martin LA, Kilpatrick JA, Al-Dulaimi R, Mone MC, Tonna JE, Barton RG, Brooke BS (2019) Predicting ICU readmission among surgical ICU patients: Development and validation of a clinical nomogram. Surgery 165 (2): 373-380. https://doi.org/10.1016/j.surg.2018.06.053

[^4]: Frost SA, Tam V, Alexandrou E, Hunt L, Salamonson Y, Davidson PM, Parr MJ, Hillman KM (2010) Readmission to intensive care: development of a nomogram for individualising risk. Crit Care Resusc. 12(2): 83-89. PMID: 20513215