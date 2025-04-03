These scripts were run using MATLAB R2019b, R-4.4.0, and Python 3.9.18. The prediction algorithm scripts also make use of parallel processing and batch scripting on a HPC system, which should be adjusted according to your needs.

## Dependencies
MATLAB
- Kernel Ridge Regression requires separately downloading [these scripts from the CBIG repository](https://github.com/ThomasYeoLab/CBIG/tree/master/utilities/matlab/predictive_models/KernelRidgeRegression) and [these hyperparameter optimisation values](https://github.com/ThomasYeoLab/CBIG/blob/master/stable_projects/preprocessing/Li2019_GSR/KernelRidgeRegression/lambda_set.mat). Credit goes to Li et al. (2019) for their work here, please cite them if using this implementation!

R
- Packages for making FC: `pheatmap`, `RcolorBrewer`, `dplyr`, `tidyverse`, `stringr`
- Packages for generating results: `ggplot2`, `readxl`, `dplyr`, `ggpubr`, `writexl`, `gtools`, `tidyverse`
 
Python
- Packages for our usage of KRR: `pandas`, `numpy`
- Multilayer meta-matching requires separately downloading [these scripts from the CBIG repository](https://github.com/ThomasYeoLab/Meta_matching_models/tree/main/rs-fMRI), with setup and package requirements explained in their v2.0 README. Again, all credit to Chen et al. (2024) for this tool, please cite them if using!
- Some extra packages are required for our usage: `scipy.stats`, `pandas`, `concurrent.futures`

## How to run
1. `make_FC.R`: This script creates a `conmats` directory with the same structure as that of `timeseries`, then populates it with individual FC and ΔFC matrices in `.txt` format, as well as `.png` visualisations.
   - Please specify the path to this repository in line 8.
2. `submit_CPM.sh`: This script submits eight HPC cluster jobs via `run_CPM.sh`, which itself calls upon `CPM.py`.
   - Please specify the path to this repository in line 25, then adjust lines 3-19 in `run_CPM.sh` according to your needs.
   - Upon running, a `cpm` directory will be created and gradually populated with results.
   - Also note that `CPM.py` begins by saving `cpm/<scale>/<timepoint>/<fc_type>/change_observed.xlsx`, so check that all eight of these outputs have been saved correctly before proceeding with the next steps, as the KRR & meta-matching scripts require these tables as inputs.
   - If you wish to reproduce our supplementary CPM models, simply adjust line 29 to use p < 0.05 or p < 0.001.
3. `submit_KRR.sh`: This script submits eight HPC cluster jobs via `run_KRR.sh`, which itself calls upon `make_KRR_inputs.m`, `run_KRR_repeats.m`, and `run_KRR_perms.m`.
   - Please specify the path to this repository in line 23 and the CBIG repository in line 20, then adjust lines 3-18 in `run_KRR.sh` according to your needs.
   - Upon running, `krr` directory will be created and gradually populated with results.
4. `submit_meta-matching.sh`: This script submits four HPC cluster jobs via `run_meta-matching.sh`, which itself calls upon `meta_matching.py`.
   - Please specify the path to this repository in line 19 and CBIG's meta-matching repository in line 16, then adjust lines 3-19 in `run_meta-matching.sh`.
   - Upon running, a `meta-matching` directory will be created and gradually populated with results.
5. `make_results_figures.R`: This script creates a `figures` directory, then populates it with the following figures:
   - Box-plots of performance (predicted-observed *r*) for all models (Figures 2-4 in the paper)
   - Distributions of observed clinical outcomes (Figure S2 in the supplement)
   - Prediction performances (*r_mean*) for all models, superimposed against family-wise error (FWE)-corrected null distributions (Figures S3-S5 in the supplement)
   - This script will also save models' stats (*r_mean*, *p*, *p_FWE*) in the `cpm`, `krr`, and `meta-matching` directories (Tables 1-3 in the paper)
