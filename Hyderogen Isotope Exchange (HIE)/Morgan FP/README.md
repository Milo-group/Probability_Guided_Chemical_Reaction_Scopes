This folder contains Morgan ECFP features for the Hydrogen Isotope Exchange (HIE) case. Near-constant fingerprint bits are removed before modeling. SMILES and training labels are taken from **RDKit Fragments - Published**.

Run `python morgan_fp_filter.py` in this folder before the R scripts.

**Files list**

`morgan_fp_filter.py`  
Builds Morgan fingerprints and writes filtered training and prediction feature tables.

`morgan_feature_config.R`  
Shared settings for loading the Morgan feature CSVs in the R scripts.

`FP Model Construction_HIE.R`  
Model search, similarity under-sampling, confusion matrices, probability heatmaps, and predictions.

`FP Cross validation_HIE.R`  
Cross-validation, leave-one-out analysis, and PDF model report.

`Training_morgan_filtered.csv`  
Training data with `Name`, `class`, and retained `MFP_*` columns (generated).

`Prediction_morgan_filtered.csv`  
Prediction substrates with the same feature columns (generated).

`CrossValidation.ModelReport.HIE.MorganFP.pdf`  
Cross-validation model report (generated).
