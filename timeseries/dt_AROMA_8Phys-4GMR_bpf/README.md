These folders contain individual fMRI time series, preprocessed using the 'dt_AROMA_8Phys-4GMR_bpf' detailed in the Supplement. Briefly, each voxel underwent:

1. Linear detrending
2. Removal of motion-related signals using the automated ICA-AROMA algorithm
3. Removal of mean signals from grey matter, white matter, and cerebrospinal fluid tissues, as well as their squares, temporal derivatives, and squares of temporal derivatives
4. Band-pass filtering at 0.008-0.08 Hz

Time series were then calculated by taking the mean signal for each parcel, weighting each voxel's contribution by its probability of being grey matter, as estimated by fMRIPrep.

## Folder descriptions
1. `ses-1` contains baseline fMRI data for the 328-region parcellation (CPM & KRR)
2. `ses-2` contains 3-month fMRI data for the 328-region parcellation (CPM & KRR)
3. `ses-1_419-reg` contains baseline fMRI data for the 419-region parcellation (multilayer meta-matching)
