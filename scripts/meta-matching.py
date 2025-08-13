import os
import sys
import random
from scipy import stats
import torch
import pickle
import sklearn
from sklearn.model_selection import train_test_split
import numpy as np
import warnings
warnings.filterwarnings("ignore")
import pandas as pd
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

# ---------------------------------------------- SETUP -----------------------------------------
# Parse arguments from the bash script
repo = sys.argv[1]            # Path to this repo
meta = sys.argv[2]            # Path to `Meta_matching_models/rs-fMRI` from the CBIG repo
scale = sys.argv[3]           # 'bprs' or 'sofas'
timepoint = sys.argv[4]       # '4' or '5', (6/12 months post-baseline), supps used 'lm_slopes'
preproc = sys.argv[5]         # 'dt_AROMA_8Phys-4GMR_bpf', supps used 'dt_AROMA_8Phys_bpf'
outcome_type = sys.argv[6]    # 'change_scores', supps used 'lm_slopes'
group = sys.argv[7]           # 'all', supps used '_placebo' and '_medication'

# Parameters
splits = 100
permutations = 1000
# ----------------------------------------------------------------------------------------------

# Get meta-matching v2.0 model files
sys.path.append(meta)
sys.path.append(f'{meta}/utils')

from CBIG_model_pytorch import check_models_v20
from CBIG_model_pytorch import demean_norm
from CBIG_model_pytorch import stacking
from CBIG_model_pytorch import multilayer_metamatching_infer

# Check whether meta-matching v2.0 model files exist and are up-to-date
check_models_v20(f'{meta}/v2.0/models')

def run_splits(y_input, x_input):
    """
    Runs multilayer meta-matching across many splits of data into the folds.
    Returns the predicted-observed r-values for each subject, for each split.
    """
    
    # Initialise results dataframe
    r_vals = np.full(splits, np.nan, dtype='float64')
    
    # Run meta-matching model & store results for each iteration
    for n in range(splits):
        
        # Set seed
        random.seed(n)
        np.random.seed(n)
        torch.manual_seed(n)
        torch.cuda.manual_seed(n)
        torch.cuda.manual_seed_all(n)
    
        # Split data into training and testing
        x_train, x_test, y_train, y_test = train_test_split(x_input, y_input, test_size=0.25, random_state=n)
        n_subj_train, n_subj_test = x_train.shape[0], x_test.shape[0]
        
        # Inference using all DNN + RR models from all source datasets
        dataset_names = {'extra-large': 'UKBB', 'large': 'ABCD', 'medium': ['GSP', 'HBN', 'eNKI']}
        y_train_pred, y_names = multilayer_metamatching_infer(x_train, y_train, f'{meta}/v2.0/models', dataset_names)
        y_test_pred, _ = multilayer_metamatching_infer(x_test, y_test, f'{meta}/v2.0/models', dataset_names)
        
        # Stacking
        y_test_final=np.zeros((y_test_pred.shape[0], y_train.shape[1]))
        for i in range(y_train.shape[1]):
            # For each test phenotype, perform stacking by developing a KRR model
            y_test_temp, _ = stacking(y_train_pred, y_test_pred, y_train[:,i].view(), [0.00001, 0.0001, 0.001, 0.004, 0.007, 0.01, 0.04, 0.07, 0.1, 0.4, 0.7, 1, 1.5, 2, 2.5, 3, 3.5, 4, 5, 10, 15, 20])
            y_test_final[:,i] = y_test_temp.flatten()
        
        r_vals[n] = stats.pearsonr(y_test_final[:,0], y_test[:,0])[0]
        
    return r_vals
  
def run_permutation(i):
    """
    Runs multilayer meta-matching across many splits, using randomly permuted outcome data.
    Returns the mean r-value across all splits.
    This function is called in parallel for each permutation.
    """
    # Permute behavioural data
    np.random.seed(i)
    perm_outcomes = np.random.permutation(outcomes)
    
    try:
        # Run CPM on permuted outcomes, retain r_mean
        r_vals = run_splits(perm_outcomes, all_fc_data)
        r_mean = np.mean(r_vals)
        
        return i, r_mean
      
    except Exception as e:
        print(f"Error in permutation {i+1}: {e}")
        return i, None
    
# -------------------------------------------------- PREPARE DATA --------------------------------------------------
print(f"""
Preparing data for meta-matching
Scale: {scale}
Timepoint: {timepoint}
Preprocessing: {preproc}
Outcome type: {outcome_type}
Group: {group}
Splits: {splits}
Permutations: {permutations}
""")

# Get list of subjects with usable baseline fMRI data
subj_list = pd.read_csv(f'{repo}/data/subjects_ses-1.txt', header=None).squeeze().tolist()
        
# Load clinical data
if group == 'all':
    group = ''

outcomes = pd.read_csv(f'{repo}/data/{scale}_{outcome_type}{group}.txt', sep='\t')

if outcome_type == 'change_scores':
    outcomes = outcomes[outcomes['timepoint'] == int(timepoint)]

# Keep subjects with both fMRI & clinical data
subj_list = np.array([value for value in subj_list if value in outcomes['BIDS_ID'].values])
outcomes = outcomes[outcomes['BIDS_ID'].isin(subj_list)]

outcomes = outcomes['outcome'].values
outcomes = outcomes.reshape(-1, 1)

# Initialise group FC matrix
all_fc_data = np.zeros((len(subj_list), int(419*418/2)))

# Retain only lower triangle of FC matrix, and normalise
for i, subject in enumerate(subj_list):
    conmat = np.loadtxt(os.path.join(f'{repo}/data/conmats/{preproc}/419_parc/ses-1/conmat_{subject}_ses-1.txt'))
    all_fc_data[i, :] = conmat[np.tril_indices(419, k=-1)]
all_fc_data = demean_norm(all_fc_data)

# ---------------------------------------- RUN MODEL ----------------------------------------
# Time and a Word
start_time = time.time()
print('Starting!')

# Run multilayer meta-matching & save prediction performance for every split
r_vals = run_splits(outcomes, all_fc_data)
np.savetxt(f'{repo}/results/meta-matching/{scale}/{timepoint}/r_vals_{preproc}{group}.txt',
           r_vals, fmt='%.5f')

print(f'Elapsed time for model: {time.time() - start_time:.2f} seconds')
                
# ----------------------------------------- RUN NULL MODEL IF SPECIFIED -----------------------------------------
if sys.argv[8] == '1':
  
    print(f'Running {permutations} permutations of null model!')
    start_time = time.time()
  
    # Lower splits to reduce compute
    splits = 20

    # Generate permutations
    perm = np.zeros((permutations, len(subj_list)), dtype=np.int32)
    rng = np.random.default_rng(seed=1)

    for i in range(permutations):
        perm[i] = rng.permutation(len(subj_list))

    # Initialise array to store null r_mean values
    null_r = np.full((permutations, 2), np.nan, dtype='float64')
    
    # Parallel processing
    with ProcessPoolExecutor(max_workers=int(sys.argv[9])) as executor:
        futures = [
            executor.submit(
                run_permutation, i
            ) for i in range(permutations)
        ]
        
        # Save results
        for future in as_completed(futures):
            i, r_mean = future.result()
            if r_mean is not None:
                null_r[i, 0] = i + 1  # Store permutation index
                null_r[i, 1] = r_mean
                np.savetxt(f'{repo}/results/meta-matching/{scale}/{timepoint}/null_r_means_{preproc}{group}.txt',
                           null_r, fmt='%.5f')

    print(f"Execution time: {time.time() - start_time:.2f} seconds")
