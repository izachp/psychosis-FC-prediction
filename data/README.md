## Folder descriptions
- `timeseries` contains patients' fMRI time series, preprocessed using either the `dt_AROMA_8Phys-4GMR_bpf` or `dt_AROMA_8Phys_bpf` pipeline and parcellated using either the `419_parc` or `328_parc` scheme (details provided in supplementary information). Our main analysis used the former pipeline and parcellation. `ses-1` and `ses-2` contain data from baseline and 3-month scans, respectively.

## File descriptions
- `low_signal_regions.txt`: The four regions of the globi pallidi with consistently low fMRI signal, which are automatically excluded when computing FC from `328_parc` data, as done in some exploratory models.
- `subjects_ses-<ses>.txt`: List of patients with usable data from baseline or 3-month scans.
- All other files contain *z*-scored clinical outcome data. Our main analyses used `<sofas/bprs>_change_scores.txt`.