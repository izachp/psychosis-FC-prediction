Data and code supporting '[Functional Coupling and Longitudinal Outcome Prediction in First-Episode Psychosis](https://doi.org/10.1016/j.bpsgos.2025.100589)'

We used resting-state functional magnetic resonance imaging (fMRI) to evaluate whether baseline functional coupling (FC), or 3-month changes (ΔFC), could predict 6- and 12-month changes in symptoms and functioning of first-episode psychosis patients in the [STAGES randomised control trial](https://doi.org/10.1093/schizbullopen/sgaa015). *Short answer: it could not.*

Three cross-validated algorithms were chosen:
- [Connectome-based predictive modelling (CPM)](https://doi.org/10.1038/nprot.2016.178)
- [Kernel ridge regression (KRR)](https://doi.org/10.1016/j.neuroimage.2019.04.016)
- [Multilayer meta-matching](https://doi.org/10.1162/imag_a_00233)

**This repository includes:**
- Clinical outcome data
- Preprocessed and parcellated fMRI timeseries
- Code for creating FC (and ΔFC) matrices
- Code for running prediction models
- Code for assessing prediction performance/significance and generating results figures

**NOTE:** Due to ethics constraints on the STAGES dataset, we are unable to provide fMRI images. The clinical outcome data has also been z-scored (within timepoints where applicable), which may create minor discrepancies when using this repository to reproduce our results.

## Setup
- Python dependencies are listed in `pyproject.toml`, and can be set up with [UV](https://docs.astral.sh/uv/) by running `uv sync` from the project root in the terminal.
- The R scripts for making FC matrices and generating results require the `tidyverse`, `ggplot2`, and `ggpubr` libraries. 
- The KRR analyses require [these MATLAB scripts from the CBIG repository](https://github.com/ThomasYeoLab/CBIG/tree/master/utilities/matlab/predictive_models/KernelRidgeRegression) to be downloaded.
- The meta-matching analyses also require [these files](https://github.com/ThomasYeoLab/Meta_matching_models/tree/main/rs-fMRI) to be downloaded, along with [these pre-trained model files](https://github.com/ThomasYeoLab/Meta_matching_models/releases/tag/v2.0-rsfMRI). Unzip the 656MB file from the latter link and place the `models` directory inside the `v2.0` directory from the former link.

Details on how to reproduce the paper's results are provided in `scripts`, and information is also provided in `data`.

## Citation
Pope IZ, Chopra S, Holmes A, Francey SM, O’Donoghue B, Cropley VL, et al. Functional Coupling and Longitudinal Outcome Prediction in First-Episode Psychosis. Biological Psychiatry Global Open Science. 2025 Aug 11. doi:10.1016/j.bpsgos.2025.100589

If you use any of the prediction algorithms implemented here, please cite their work:

**Connectome-based predictive modelling**

Shen X, Finn ES, Scheinost D, Rosenberg MD, Chun MM, Papademetris X, et al. Using connectome-based predictive modeling to predict individual behavior from brain connectivity. Nat Protoc. 2017 Mar;12(3):506–18.

**Kernel Ridge Regression**

Li J, Kong R, Liégeois R, Orban C, Tan Y, Sun N, et al. Global signal regression strengthens association between resting-state functional connectivity and behavior. NeuroImage. 2019 Aug 1;196:126–41.

**Multilayer Meta-matching**

Chen P, An L, Wulan N, Zhang C, Zhang S, Ooi LQR, et al. Multilayer meta-matching: Translating phenotypic prediction models from multiple datasets to small data. Imaging Neuroscience. 2024 Jul 17;2:1–22.
