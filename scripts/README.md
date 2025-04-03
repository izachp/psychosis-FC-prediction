These scripts were run using MATLAB R2019b, R-4.4.0, and Python 3.9.18. The prediction algorithm scripts also make use of parallel processing and batch scripting on a HPC system, which may need to be adjusted according to your needs.

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
