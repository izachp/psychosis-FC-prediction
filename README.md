Data and code for submitted manuscript 'Functional coupling and longitudinal outcome prediction in first-episode psychosis'

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

