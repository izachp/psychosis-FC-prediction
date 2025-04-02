These folders contain individual fMRI time series, preprocessed using the 'dt_AROMA_8Phys-4GMR_bpf' detailed in the Supplement. Briefly, each voxel underwent:

1. Linear detrending
2. Removal of motion-related signals using the automated ICA-AROMA algorithm
3. Removal of mean signals from grey matter, white matter, and cerebrospinal fluid tissues, as well as their squares, temporal derivatives, and squares of temporal derivatives
4. Band-pass filtering at 0.008-0.08 Hz

Time series were then calculated by taking the mean signal for each parcel, weighting each voxel's contribution by its probability of being grey matter, as estimated by fMRIPrep.

`ses-1 ` contains baseline fMRI data, parcellated into 328 regions (Schaefer300 + Scale II Melbourne Subcortex - low signal regions), and is used for CPM & KRR
`ses-2` contains 3-month fMRI data, parcellated into 328 regions (Schaefer300 + Scale II Melbourne Subcortex - low signal regions), and is used for CPM & KRR
`ses-1_419-reg` contains baseline fMRI, parcellated into 419 regions (Schaefer400 + FreeSurfer's subcortical segmentation), and is used for multilayer meta-matching
