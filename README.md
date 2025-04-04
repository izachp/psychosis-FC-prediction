Data and code for preprint '[Functional Coupling and Longitudinal Outcome Prediction in First-Episode Psychosis](https://www.medrxiv.org/content/10.1101/2025.04.01.25325005v1)'

This work uses multiple cross-validated prediction algorithms to evaluate whether resting-state functional coupling (FC), acquired via fMRI, could predict the longitudinal changes in symptoms and functioning of first-episode psychosis patients in the [STAGES clinical trial](https://academic.oup.com/schizbullopen/article/1/1/sgaa015/5810294). *Short answer: it could not.*

**This repository includes:**
- Clinical outcome data
- Preprocessed and parcellated fMRI timeseries
- Code for creating FC matrices
- Code for running the three prediction algorithms
- Code for assessing prediction performance & significance and generating results figures

**NOTE:** Since the STAGES dataset used here is not open-access, we are unable to provide raw MRI data for patients. The clinical outcome data has also been z-scored, which may minimally affect any results generated.

Details on how to use this repository are provided in `scripts`.

## File Descriptions
1. `bprs_data.csv`: Contains z-scores corresponding to each patient's total on the Brief Psychiatric Rating Scale, assessed at baseline, 6 months, or 12 months (fkTimepoint = 1, 4, and 5, respectively).
2. `sofas_data.csv`: Contains z-scores corresponding to each patient's total on the Social and Occupational Functioning Assessment Scale, assessed at baseline, 6 months, or 12 months (fkTimepoint = 1, 4, and 5, respectively).
3. `low_signal_regions.txt`: Contains names and parcel IDs of regions excluded from the CPM and KRR models due to low signal.
4. `subjects_ses-1.txt`: List of subjects with usable FC data at baseline.
5. `subjects_ses-2.txt`: List of subjects with usable FC data at 3 months.

## Citation
[PREPRINT] Pope IZ, Chopra S, Holmes A, Francey SM, O’Donoghue B, Cropley VL, et al. Functional Coupling and Longitudinal Outcome Prediction in First-Episode Psychosis. medRxiv. 2025 Apr 4. doi:10.1101/2025.04.01.25325005v1

If you use any of the prediction algorithms implemented here, please cite their work:

**Connectome-based predictive modelling**

Shen X, Finn ES, Scheinost D, Rosenberg MD, Chun MM, Papademetris X, et al. Using connectome-based predictive modeling to predict individual behavior from brain connectivity. Nat Protoc. 2017 Mar;12(3):506–18.

**Kernel Ridge Regression**

Li J, Kong R, Liégeois R, Orban C, Tan Y, Sun N, et al. Global signal regression strengthens association between resting-state functional connectivity and behavior. NeuroImage. 2019 Aug 1;196:126–41.

**Multilayer Meta-matching**

Chen P, An L, Wulan N, Zhang C, Zhang S, Ooi LQR, et al. Multilayer meta-matching: Translating phenotypic prediction models from multiple datasets to small data. Imaging Neuroscience. 2024 Jul 17;2:1–22.
