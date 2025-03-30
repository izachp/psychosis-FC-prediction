# This script has been adapted from a script written by Emily Finn, available at github.com/esfinn/cpm_tutorial

import numpy as np
import scipy as sp
import pandas as pd
from pathlib import Path
import os
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import openpyxl
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
import sys

repo = sys.argv[1]

# Select model and set parameters
scale = sys.argv[2]           # 'bprs' or 'sofas'
timepoint = sys.argv[3]       # '4' or '5', corresponding to 6 and 12 months post-baseline
fc_type = sys.argv[4]         # 'fc_baseline' or 'fc_change_bl3m'

p_thresh = sys.argv[5]        # Supplementary results also used p<0.05 and p<0.001
k = 4
splits = 100
permutations = 1020   # Excess accounts for some permutations failing (SVD doesn't converge in np.polyfit)

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

def mk_kfold_indices(rng):
    """
    Splits list of subjects into k folds for cross-validation.
    """
    
    n_subs = len(subj_list)
    n_subs_per_fold = n_subs//k # floor integer for n_subs_per_fold

    indices = [[fold_no]*n_subs_per_fold for fold_no in range(k)] # generate repmat list of indices
    remainder = n_subs % k # figure out how many subs are left over
    remainder_inds = list(range(remainder))
    indices = [item for sublist in indices for item in sublist]    
    [indices.append(ind) for ind in remainder_inds] # add indices for remainder subs

    assert len(indices)==n_subs, "Length of indices list does not equal number of subjects, something went wrong"

    rng.shuffle(indices) # shuffles in place

    return np.array(indices)

def split_train_test(indices, test_fold):
    """
    For a subj list, k-fold indices, and given fold, returns lists of train_subs and test_subs
    """

    train_inds = np.where(indices!=test_fold)
    test_inds = np.where(indices==test_fold)

    train_subs = []
    for sub in subj_list[train_inds]:
        train_subs.append(sub)

    test_subs = []
    for sub in subj_list[test_inds]:
        test_subs.append(sub)

    return (train_subs, test_subs)
  
def get_train_test_data(behav_data, train_subs, test_subs):

    """
    Extracts requested FC and behavioral data for a list of train_subs and test_subs
    """

    train_vcts = all_fc_data.loc[train_subs, :]
    test_vcts = all_fc_data.loc[test_subs, :]
    
    train_behav = behav_data.loc[train_subs, 'outcome']

    return (train_vcts, train_behav, test_vcts)
  
def select_features(train_vcts, train_behav, train_covars):
    
    """
    Runs the CPM feature selection step: 
    - correlates each edge with behavior, and returns a mask of edges that are correlated above some threshold, one for each tail (positive and negative)
    """

    assert train_vcts.index.equals(train_behav.index), "Row indices of FC vcts and behavior don't match!"

    mask_dict = {}
    
    # Correlate all edges with behav vector
    cov = np.dot(train_behav.T - train_behav.mean(), train_vcts - train_vcts.mean(axis=0)) / (train_behav.shape[0]-1)
    corr = cov / np.sqrt(np.var(train_behav, ddof=1) * np.var(train_vcts, axis=0, ddof=1))

    # Get p-values of correlation at each edge (two-tailed)
    n = train_behav.shape[0]
    t_stat = corr * np.sqrt((n - 2) / (1 - corr**2))
    p_val = 2 * sp.stats.t.sf(np.abs(t_stat), df=n-2)
    
    # Define positive and negative masks
    mask_dict["pos"] = (p_val < p_thresh) & (corr > 0)
    mask_dict["neg"] = (p_val < p_thresh) & (corr < 0)
    
    return mask_dict
  
def build_model(train_vcts, mask_dict, train_behav):
    """
    Builds a CPM model:
    - takes a feature mask, sums all edges in the mask for each subject, and uses simple linear regression to relate summed network strength to behavior
    """

    assert train_vcts.index.equals(train_behav.index), "Row indices of FC vcts and behavior don't match!"

    model_dict = {}

    # Loop through pos and neg tails
    t = 0
    for tail, mask in mask_dict.items():
        X = train_vcts.values[:, mask].sum(axis=1)
        y = train_behav
        (slope, intercept) = np.polyfit(X, y, 1)
        model_dict[tail] = (slope, intercept)
        t+=1

    return model_dict
  
def apply_model(test_vcts, mask_dict, model_dict):
    """
    Applies a previously trained linear regression model to a test set to generate predictions of behavior.
    """

    behav_pred = {}

    # Loop through pos and neg tails
    t = 0
    for tail, mask in mask_dict.items():
        X = test_vcts.loc[:, mask].sum(axis=1)

        slope, intercept = model_dict[tail]
        behav_pred[tail] = slope*X + intercept
        t+=1

    return behav_pred
  
def cpm_wrapper(behav_data, rng):

    assert all_fc_data.index.equals(behav_data.index), "Row (subject) indices of FC vcts and behavior don't match!"
    
    indices = mk_kfold_indices(rng)
    
    # Initialize df for storing observed and predicted behavior
    col_list = []
    for tail in ["pos", "neg"]:
        col_list.append("Proportional change predicted (" + tail + ")")
    col_list.append("Proportional change observed")
    behav_obs_pred = pd.DataFrame(index=subj_list, columns = col_list)
    
    for fold in range(k):
        train_subs, test_subs = split_train_test(indices, test_fold=fold)
        train_vcts, train_behav, train_covars, test_vcts = get_train_test_data(behav_data, train_subs, test_subs)
        mask_dict = select_features(train_vcts, train_behav, train_covars)
        model_dict = build_model(train_vcts, mask_dict, train_behav)
        behav_pred = apply_model(test_vcts, mask_dict, model_dict)
        for tail, predictions in behav_pred.items():
            behav_obs_pred.loc[test_subs, "Proportional change predicted (" + tail + ")"] = predictions
            
    behav_obs_pred.loc[subj_list, "Proportional change observed"] = behav_data['outcome']
    
    return behav_obs_pred

def run_splits(behav_data):

    # Set random seed for reproducibility
    rng = np.random.default_rng(seed=1)
    
    # Initialise results dataframe
    r_vals = pd.DataFrame(index=range(splits), columns=['pos', 'neg'], dtype='float64')
    
    # Run CPM & store r-values for each split
    for i in range(splits):
        
        behav_obs_pred = cpm_wrapper(behav_data, rng)
        
        r_vals.iloc[i, 0] = sp.stats.pearsonr(behav_obs_pred.iloc[:,2], behav_obs_pred.iloc[:,0])[0]
        r_vals.iloc[i, 1] = sp.stats.pearsonr(behav_obs_pred.iloc[:,2], behav_obs_pred.iloc[:,1])[0]

    r_vals = r_vals.round(5)
    
    return(r_vals)

# ----------- PREPARE DATA --------------

# Load clinical data for scale and timepoint of interest
data = pd.read_csv(repo + '/' + scale + '_data.csv')
end_scores = data[data['fkTimePointID'] == int(timepoint)]

# Read in subjects list to account for missing data / motion exclusion
if fc_type == 'fc_baseline':
    with open(repo + '/subjects_ses-1.txt') as f:
        subj_list = f.read().splitlines()
        
elif fc_type == 'fc_change_bl3m':
    with open(repo + '/subjects_ses-2.txt') as f:
        subj_list = f.read().splitlines()
        
else:
    print('Invalid FC type. Please enter "fc_baseline" or "fc_change_bl3m"')
    sys.exit()

# Select only subjects with complete clinical data
subj_list = np.array([value for value in subj_list if value in end_scores['BIDS_ID'].values])

# Get clinical scores at baseline and follow-up
end_scores = end_scores[end_scores['BIDS_ID'].isin(subj_list)]
end_scores = end_scores[scale + '_score'].values
base_scores = data[data['fkTimePointID'] == 1 & data['BIDS_ID'].isin(subj_list)]
base_scores = base_scores[scale + '_score'].values

# Calculate clinical outcomes as proportional change
outcomes = (end_scores - base_scores) / base_scores
outcomes = outcomes[~np.isnan(outcomes)]

outcomes = pd.DataFrame({'BIDS_ID': subj_list, 'outcome': outcomes}, index=subj_list)
outcomes = outcomes.set_index('BIDS_ID')

# Save observed clinical outcomes
outcomes.to_excel(repo + '/cpm/' + scale + '/' + timepoint + '/' + fc_type + '/change_observed.xlsx')

# Load functional connectivity data
if fc_type == 'fc_baseline':
    all_fc_data = read_in_matrices(subj_list, '_ses-1', Path(repo + '/conmats/dt_AROMA_8Phys-4GMR_bpf/ses-1'))
  
elif fc_type == 'fc_change_bl3m':
    all_fc_data = read_in_matrices(subj_list, '_change', Path(repo + '/conmats/dt_AROMA_8Phys-4GMR_bpf/change'))
  
else:
    print('Invalid FC type. Please enter "fc_baseline" or "fc_change_bl3m"')
    sys.exit()

n_edges = all_fc_data.shape[1]

# ----------- RUN MODEL --------------

# Time and a Word
start_time = time.time()
print(f'Running {splits} splits of CPM! \nScale: {scale} \nTimepoint: {timepoint} \nFC type: {fc_type}')

# Run splits of CPM & save r-values
r_vals = run_splits(outcomes)
r_vals.to_excel(repo + '/cpm/' + scale + '/' + timepoint + '/' + fc_type + '/r_vals.xlsx', index=False)

end_time = time.time()
print(f"Execution time: {end_time - start_time:.2f} seconds")

# ----------- NULL MODEL --------------

# Generate permutations
perm = np.zeros((permutations, len(subj_list)), dtype=np.int8)
rng = np.random.default_rng(seed=1)

for i in range(permutations):
    perm[i] = rng.permutation(len(subj_list))

# Initialise dataframe to store null model r-values
null_r = pd.DataFrame(index=range(permutations), columns=['permutation', 'pos', 'neg'], dtype='float64')

def run_permutation(i):
    
    # Shuffle clinical outcomes according to a permutation
    perm_outcomes = outcomes.take(perm[i])
    perm_outcomes.index = subj_list
    
    try:
        # Run CPM on permuted outcomes, retain r_mean
        r_vals = run_splits(perm_outcomes)
        r_mean = r_vals.mean()
        
        return i, r_mean
      
    except np.linalg.LinAlgError as e:     # This protects null_r when a permutation fails (SVD doesn't converge in np.polyfit)
      print(f"Error in permutation {i+1}: {e}")
      return i, None

def main():
    
    with ProcessPoolExecutor(max_workers=sys.argv[6]) as executor:
        futures = [
            executor.submit(
                run_permutation, i
            ) for i in range(permutations)
        ]
        
        for future in as_completed(futures):
            i, r_mean = future.result()
            if r_mean is not None:
                null_r.loc[i] = [i+1, r_mean['pos'], r_mean['neg']]
                print(f'Finished permutation {i+1}')
                null_r.to_excel(repo + '/cpm/' + scale + '/' + timepoint + '/' + fc_type + '/null_r_means.xlsx', index=False)

# ------------ RUN NULL MODEL --------------
print('Running null model...')
start_time = time.time()

if __name__ == '__main__':
    main()

end_time = time.time()
print(f"Execution time: {end_time - start_time:.2f} seconds")
