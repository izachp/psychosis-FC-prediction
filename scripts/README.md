These scripts were run using MATLAB R2019b, R-4.4.0, and Python 3.9.18. The prediction algorithm scripts also make use of parallel processing and `sbatch` jobs on a high-performance computing cluster ([MASSIVE](doi.org/10.3389/fninf.2014.00030)), which should be adjusted according to your setup.

Some scripts need to be adjusted near the top:
- `make_FC.R`, `submit_<algorithm>.sh`, and `make_figures.R` require paths to this repository, the CBIG repository, and/or the meta-matching files
- `make_FC.R`, `submit_<algorithm>.sh`, and `make_figures.R` contain options for null modelling and exploratory models (alternate preprocessing, parcellation, outcome, group, or CPM feature selection threshold)
- `run_<algorithm>.sh` requires `sbatch` job settings

**NOTE:** Null modelling is disabled by default, as this can take a lot of time and compute.

## How to run
1. `make_FC.R`: Populates a `data/conmats` directory with FC and ΔFC matrices for each patient.
2. `submit_<algorithm>.sh`: Submits `sbatch` jobs of `run_<algorithm>.sh`, which itself calls upon algorithm-specific scripts.
   - The three algorithms can be run in any order
   - The results for each prediction model will appear in a `results/<algorithm>` directory.
   - Job logs will appear in a `./slurmout` directory
3. `make_figures.R`: Populates a `results/figures` directory with:
   - Box plots of performance (predicted-observed *r*) for all models
   - Distributions of observed outcomes (Figure S2)
   - Mean prediction performances (*r_mean*) for all models, superimposed against family-wise error (FWE)-corrected null distributions
   - Stats (*r_mean*, *p*, *p_FWE*) will be saved in a `results/stats` directory
   - Also included is a *t*-test comparing outcomes between patients with usable 3-month FC data and those without
