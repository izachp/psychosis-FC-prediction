# This script has been adapted from a Jupyter notebook written by Emily Finn, available at github.com/esfinn/cpm_tutorial

import numpy as np
from scipy import stats
import pandas as pd
from pathlib import Path
import os
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
import sys

# ---------------------------------------------- SETUP ----------------------------------------------
# Parse arguments from the bash script
repo = sys.argv[1]
scale = sys.argv[2]           # 'bprs' or 'sofas'
timepoint = sys.argv[3]       # '4' or '5', (6/12 months post-baseline), supps used 'lm_slopes'
fc_type = sys.argv[4]         # 'fc_baseline' or 'fc_change_bl3m'
p_thresh = float(sys.argv[5]) # 0.01, supps used p<0.05 and p<0.001
group = sys.argv[6]           # 'all', supps used '_placebo' and '_medication'
preproc = sys.argv[7]         # 'dt_AROMA_8Phys-4GMR_bpf', supps used 'dt_AROMA_8Phys_bpf'
parc = sys.argv[8]            # '419_parc', supps used '328_parc'
outcome_type = sys.argv[9]    # 'change_scores', supps used 'lm_slopes'

# Parameters
k = 4
splits = 100
permutations = 1020   # Excess accounts for some permutations failing (SVD doesn't converge in np.polyfit)
# ---------------------------------------------------------------------------------------------------

def read_in_matrices(subj_list, file_suffix, data_dir):
    """
    Reads in a set of individual-subject connectivity matrices stored in data_dir,
    
    Returns a dataframe that is subjects x edges (by vectorizing the upper triangle of each FC matrix).
    
    Assumes:
    - each matrix is stored in a separate file beginning with the subject ID, and
    - matrices are symmetric (squareform); i.e., for a parcellation with 268 nodes, matrices should be 268 x 268
    """
    
    all_fc_data = {}
            
    for subj in subj_list:
      
        file = [f for f in os.listdir(data_dir) if subj in f and file_suffix in f]

        # read it in and make sure it's symmetric and has reasonable dimensions
        tmp = np.loadtxt(data_dir / file[0])
        assert tmp.shape[0]==tmp.shape[1]>1, "Matrix seems to have incorrect dimensions: {}".format(tmp.shape)
        
        # take just the upper triangle and store it in a dictionary
        all_fc_data[subj] = tmp[np.triu_indices_from(tmp, k=1)]
        
    # Convert dictionary into dataframe
    all_fc_data = pd.DataFrame.from_dict(all_fc_data, orient='index')
    
    return all_fc_data

def mk_kfold_indices(rng, n_subs):
    """
    Splits list of subjects into k folds for cross-validation.
    """
    
    n_subs_per_fold = n_subs // k
    indices = np.repeat(np.arange(k), n_subs_per_fold)
    remainder = n_subs % k
    if remainder > 0:
        indices = np.concatenate([indices, np.arange(remainder)])
    rng.shuffle(indices)
    return indices

def split_train_test(indices, test_fold, subj_list):
    """
    For a subj list, k-fold indices, and given fold, returns lists of train_subs and test_subs
    """

    train_mask = indices != test_fold
    test_mask = indices == test_fold
    return subj_list[train_mask], subj_list[test_mask]

def get_train_test_data(behav_data, train_subs, test_subs):
    """
    Extracts requested FC and behavioral data for a list of train_subs and test_subs
    """

    train_vcts = all_fc_data.loc[train_subs, :].values
    test_vcts = all_fc_data.loc[test_subs, :].values
    train_behav = behav_data.loc[train_subs, 'outcome'].values

    return train_vcts, train_behav, test_vcts
  
def select_features(train_vcts, train_behav):
    """
    Runs the CPM feature selection step: 
    - correlates each edge with behavior, and returns a mask of edges that are correlated above some threshold, one for each tail (positive and negative)
    """

    # assert train_vcts.index.equals(train_behav.index), "Row indices of FC vcts and behavior don't match!"

    # Vectorized correlation and p-value calculation
    n = train_behav.shape[0]
    cov = np.dot(train_behav - train_behav.mean(), train_vcts - train_vcts.mean(axis=0)) / (n - 1)
    corr = cov / (np.std(train_behav, ddof=1) * np.std(train_vcts, axis=0, ddof=1))
    t_stat = corr * np.sqrt((n - 2) / (1 - corr**2))
    p_val = 2 * (1 - stats.t.cdf(np.abs(t_stat), df=n - 2))
    
    # Define positive and negative masks
    mask_dict = {
        "pos": (p_val < p_thresh) & (corr > 0),
        "neg": (p_val < p_thresh) & (corr < 0)
    }
    
    return mask_dict
  
def build_model(train_vcts, mask_dict, train_behav):
    """
    Builds a CPM model:
    - takes a feature mask, sums all edges in the mask for each subject, and uses simple linear regression to relate summed network strength to behavior
    """

    # assert train_vcts.index.equals(train_behav.index), "Row indices of FC vcts and behavior don't match!"

    model_dict = {}

    # Loop through pos and neg tails
    for tail, mask in mask_dict.items():
        X = train_vcts[:, mask].sum(axis=1)
        # Use np.polyfit for vectorized regression
        slope, intercept = np.polyfit(X, train_behav, 1)
        model_dict[tail] = (slope, intercept)

    return model_dict
  
def apply_model(test_vcts, mask_dict, model_dict):
    """
    Applies a previously trained linear regression model to a test set to generate predictions of behavior.
    """

    behav_pred = {}
    for tail, mask in mask_dict.items():
        X = test_vcts[:, mask].sum(axis=1)
        slope, intercept = model_dict[tail]
        behav_pred[tail] = slope * X + intercept
    return behav_pred

def cpm_wrapper(behav_data, rng):
    """
    Runs CPM (pos & neg models) for a single split of the data into folds, using each for testing once.
    Returns the predicted (pos & neg models) and observed behavior for each subject in the testing set.
    """

    # assert all_fc_data.index.equals(behav_data.index), "Row (subject) indices of FC vcts and behavior don't match!"
    
    indices = mk_kfold_indices(rng, len(subj_list))
    behav_obs_pred = np.full((len(subj_list), 3), np.nan)
    
    for fold in range(k):
        train_subs, test_subs = split_train_test(indices, fold, subj_list)
        train_vcts, train_behav, test_vcts = get_train_test_data(behav_data, train_subs, test_subs)
        mask_dict = select_features(train_vcts, train_behav)
        model_dict = build_model(train_vcts, mask_dict, train_behav)
        behav_pred = apply_model(test_vcts, mask_dict, model_dict)
        test_idx = np.where(indices == fold)[0]
        behav_obs_pred[test_idx, 0] = behav_pred["pos"]
        behav_obs_pred[test_idx, 1] = behav_pred["neg"]
        
    behav_obs_pred[:, 2] = behav_data['outcome'].values
    return behav_obs_pred

def run_splits(behav_data):
    """
    Runs CPM across many random splits of subjects across the folds.
    Returns the predicted-observed r-values for each subject, for each split of the pos and neg models.
    """

    # Set random seed for reproducibility
    rng = np.random.default_rng(seed=1)
    
    # Initialise results dataframe
    r_vals = np.full((splits, 2), np.nan)
    
    # Run CPM & store r-values for each split
    for i in range(splits):
        
        behav_obs_pred = cpm_wrapper(behav_data, rng)
        
        r_vals[i, 0] = stats.pearsonr(behav_obs_pred[:,2], behav_obs_pred[:,0])[0]
        r_vals[i, 1] = stats.pearsonr(behav_obs_pred[:,2], behav_obs_pred[:,1])[0]

    return r_vals
  
def run_permutation(i):
    """
    Runs CPM across many splits, using randomly permuted outcome data.
    Returns the mean r-value for the positive and negative models across all splits.
    This function is called in parallel for each permutation.
    """
    
    # Shuffle clinical outcomes according to a permutation
    perm_outcomes = outcomes.take(perm[i])
    perm_outcomes.index = subj_list
    
    try:
        # Run CPM on permuted outcomes, retain only mean correlation
        r_vals = run_splits(perm_outcomes)
        r_mean = np.mean(r_vals, axis=0)
        
        return i, r_mean
      
    except np.linalg.LinAlgError as e:
        print(f"Error in permutation {i+1}: {e}")
        return i, None

# -------------------------- PREPARE DATA --------------------------
print(f"""
Preparing data for CPM
Scale: {scale}
Timepoint: {timepoint}
FC type: {fc_type}
Preprocessing: {preproc}
Parcellation: {parc}
Outcome type: {outcome_type}
Feature selection threshold: p < {p_thresh}
Group: {group}
Splits: {splits}
Folds: {k}
""")

# List subjects with fMRI data
if fc_type == 'fc_baseline':
    subj_list = pd.read_csv(f'{repo}/data/subjects_ses-1.txt', header=None).squeeze().tolist()

elif fc_type == 'fc_change_bl3m':
    subj_list = pd.read_csv(f'{repo}/data/subjects_ses-2.txt', header=None).squeeze().tolist()

else:
    print('Invalid FC type. Please enter "fc_baseline" or "fc_change_bl3m"')
    sys.exit()

# Load clinical data
if group == 'all':
    group = ''

outcomes = pd.read_csv(f'{repo}/data/{scale}_{outcome_type}{group}.txt', sep='\t')

if outcome_type == 'change_scores':
    outcomes = outcomes[outcomes['timepoint'] == int(timepoint)]
    outcomes = outcomes.drop(columns=['timepoint'])

# Keep subjects with both fMRI & clinical data
subj_list = np.array([value for value in subj_list if value in outcomes['BIDS_ID'].values])
outcomes = outcomes[outcomes['BIDS_ID'].isin(subj_list)]
outcomes = outcomes.set_index('BIDS_ID')

# Save observed clinical outcomes
outcomes.to_csv(f'{repo}/results/cpm/{scale}/{timepoint}/{fc_type}/outcomes_observed{group}.csv', index=True)

# Load FC
if fc_type == 'fc_baseline':
    all_fc_data = read_in_matrices(subj_list, '_ses-1', Path(f'{repo}/data/conmats/{preproc}/{parc}/ses-1'))
  
elif fc_type == 'fc_change_bl3m':
    all_fc_data = read_in_matrices(subj_list, '_change', Path(f'{repo}/data/conmats/{preproc}/{parc}/change'))
  
else:
    print('Invalid FC type. Please enter "fc_baseline" or "fc_change_bl3m"')
    sys.exit()

n_edges = all_fc_data.shape[1]

# --------------------------- RUN MODEL ---------------------------

# Time and a Word
start_time = time.time()
print('Starting!')

# Run all splits of CPM & save r-values
r_vals = run_splits(outcomes)
np.savetxt(f'{repo}/results/cpm/{scale}/{timepoint}/{fc_type}/r_vals_{preproc}_{parc}_{p_thresh}{group}.txt', r_vals, fmt='%.5f')

print(f"Execution time: {time.time() - start_time:.2f} seconds")

# ------------- RUN NULL MODEL IF SPECIFIED BY USER --------------
if sys.argv[10] == '1':
  
    print(f'Running null model with {permutations} permutations...')
    start_time = time.time()
    
    # Generate permutations
    perm = np.zeros((permutations, len(subj_list)), dtype=np.int32)
    rng = np.random.default_rng(seed=1)

    for i in range(permutations):
        perm[i] = rng.permutation(len(subj_list))

    # Initialise array to store null model r-values
    null_r = np.empty((permutations, 3), dtype='float64')
    
    # Run permutations in parallel
    with ProcessPoolExecutor(max_workers=int(sys.argv[11])) as executor:
        futures = [
            executor.submit(
                run_permutation, i
            ) for i in range(permutations)
        ]
        
        for future in as_completed(futures):
            i, r_mean = future.result()
            if r_mean is not None:
                null_r[i, 0] = i + 1  # Store permutation index
                null_r[i, 1] = r_mean[0]  # Store mean r for positive model
                null_r[i, 2] = r_mean[1]  # Store mean r for negative model
                np.savetxt(f'{repo}/results/cpm/{scale}/{timepoint}/{fc_type}/null_r_means_{preproc}_{parc}_{p_thresh}{group}.txt', null_r, fmt='%.5f')

    print(f"Execution time: {time.time() - start_time:.2f} seconds")
