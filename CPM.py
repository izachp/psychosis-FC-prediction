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

# NOTE: I'ven't got the covariate-corrected model working yet, so this script only runs the standard model.

# Select model and set parameters
scale = sys.argv[1]                                        # 'bprs' or 'sofas'
timepoint = sys.argv[2]                                    # '3', '4', or '5', corresponding to 3, 6, and 12 months post-baseline
fc_type = sys.argv[3]                                      # 'fc_baseline' or 'fc_change_bl3m'
p_thresh = 0.01
corr_type = 'standard'   # 'standard' or 'corrected'
covariates = ['Age_bl', 'sex', 'Mean_FD_Jenk']
k = 4
iterations = 100
permutations = 1000
permutations = round(permutations*1.05)  # Add excess to account for some permutations failing (SVD doesn't converge in np.polyfit)

# Get functional connectivity matrices at baseline
def read_in_matrices(subj_list, file_suffix, data_dir, zscore=False):
    """
    Reads in a set of individual-subject connectivity matrices stored in data_dir,
    
    Returns a dataframe that is subjects x edges (by vectorizing the upper triangle of each FC matrix).
    
    Assumes:
    - each matrix is stored in a separate file beginning with the subject ID, and
    - matrices are symmetric (squareform); i.e., for a parcellation with 268 nodes, matrices should be 268 x 268
    """
    
    all_fc_data = {}
            
    for subj in subj_list:
        # try to find this subject's matrix
        if file_suffix:
            file = [f for f in os.listdir(data_dir) if subj in f and file_suffix in f]
        else:
            file = [f for f in os.listdir(data_dir) if subj in f]
            
        # make sure there is one and only one file    
        if len(file) ==0:
            raise ValueError("No data found for subject {}".format(subj))
        if len(file) >1:
            raise ValueError("More than one matrix found for subject {}! Specify a suffix?".format(subj))
        
        # read it in and make sure it's symmetric and has reasonable dimensions
        tmp = np.loadtxt(data_dir / file[0])
        assert tmp.shape[0]==tmp.shape[1]>1, "Matrix seems to have incorrect dimensions: {}".format(tmp.shape)
        
        # take just the upper triangle and store it in a dictionary
        if ~zscore:
            all_fc_data[subj] = tmp[np.triu_indices_from(tmp, k=1)]
        if zscore:
            all_fc_data[subj] = sp.stats.zscore(tmp[np.triu_indices_from(tmp, k=1)])
        
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
    
    train_behav = behav_data.loc[train_subs, 'change_score']
    train_covars = behav_data.loc[train_subs, covariates]

    return (train_vcts, train_behav, train_covars, test_vcts)
  
def select_features(train_vcts, train_behav, train_covars):
    
    """
    Runs the CPM feature selection step: 
    - correlates each edge with behavior, and returns a mask of edges that are correlated above some threshold, one for each tail (positive and negative)
    """

    assert train_vcts.index.equals(train_behav.index), "Row indices of FC vcts and behavior don't match!"

    mask_dict = {}
    
    # Correlate all edges with behav vector
    if corr_type == 'standard':
        cov = np.dot(train_behav.T - train_behav.mean(), train_vcts - train_vcts.mean(axis=0)) / (train_behav.shape[0]-1)
        corr = cov / np.sqrt(np.var(train_behav, ddof=1) * np.var(train_vcts, axis=0, ddof=1))

        # Get p-values of correlation at each edge (two-tailed)
        n = train_behav.shape[0]
        t_stat = corr * np.sqrt((n - 2) / (1 - corr**2))
        p_val = 2 * sp.stats.t.sf(np.abs(t_stat), df=n-2)
        
    elif corr_type == 'corrected':
        '''
        # Regress out covariates
        # Make list of regression weights to use in test set
    
    
        # OLD STUFF
  
        # Partial correlation
        # Combine ranked data
        ranked_train_all = pd.concat([ranked_train_vcts, ranked_train_behav, ranked_train_covars], axis=1)
        
        print(ranked_train_all.info())
        print('Done combining')
        # Compute the covariance matrix of the combined data
        cov_matrix = np.cov(ranked_train_all.iloc[:100], rowvar=False)  # This crashes RStudio! TESTING FIRST 100 COLUMNS ONLY
        
        # COMMENT OUT FOR TESTING
            
        print('Done cov')
        # Compute the inverse of the covariance matrix
        inv_cov_matrix = np.linalg.inv(cov_matrix)
        print('Done inv')
        # Initialize arrays to store correlations and p-values
        corr = np.zeros(n_edges)
        p_val = np.ones(n_edges)
        
        
        # Extract indices for target variables
        target_index = data.columns.get_loc(ranked_train_behav.name)
        covar_indices = [data.columns.get_loc(cov) for cov in covariates]
        
        for i, edge in enumerate(ranked_train_vcts.columns):
            edge_index = data.columns.get_loc(edge)
            
           
            # Partial correlation computation
            r_xy = -inv_cov_matrix[edge_index, target_index] / np.sqrt(
                inv_cov_matrix[edge_index, edge_index] * inv_cov_matrix[target_index, target_index]
            )
            
            # Store the correlation
            corr[i] = r_xy
            
            # Compute the p-value (using Fisher's z-transformation)
            n = data.shape[0]
            z = np.arctanh(r_xy)
            se = 1 / np.sqrt(n - len(covar_indices) - 3)
            z_score = z / se
            p_val[i] = 2 * (1 - sp.stats.norm.cdf(np.abs(z_score)))
      
        ## EVEN OLDER STUFF
        
        brainbehav = pd.concat([train_vcts, train_behav, train_covars], axis=1)
    
        corr = np.zeros(train_vcts.shape[1])
        p_val = np.ones(train_vcts.shape[1])
        
        for i, edge in enumerate(train_vcts.columns):
  
            stats = partial_corr(brainbehav, x=edge, y='change_score', covar=covariates, method='spearman')
            corr[i] = stats['r']
            p_val[i] = stats['p-val']
        '''
    else:
        raise ValueError("Invalid correlation type. Please enter 'standard' or 'corrected'")
    
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
    
    # Initialize array for storing feature masks
    all_masks = {}
    all_masks["pos"] = np.zeros((k, n_edges))
    all_masks["neg"] = np.zeros((k, n_edges))
    
    for fold in range(k):
        train_subs, test_subs = split_train_test(indices, test_fold=fold)
        train_vcts, train_behav, train_covars, test_vcts = get_train_test_data(behav_data, train_subs, test_subs)
        mask_dict = select_features(train_vcts, train_behav, train_covars)
        all_masks["pos"][fold,:] = mask_dict["pos"]
        all_masks["neg"][fold,:] = mask_dict["neg"]
        model_dict = build_model(train_vcts, mask_dict, train_behav)
        behav_pred = apply_model(test_vcts, mask_dict, model_dict)
        for tail, predictions in behav_pred.items():
            behav_obs_pred.loc[test_subs, "Proportional change predicted (" + tail + ")"] = predictions
            
    behav_obs_pred.loc[subj_list, "Proportional change observed"] = behav_data['change_score']
    
    return behav_obs_pred, all_masks

def run_cpm(behav_data):

    # Set random seed for reproducibility
    rng = np.random.default_rng(seed=1)
    
    # Initialise results dataframes
    r_vals = pd.DataFrame(index=range(iterations), columns=['pos', 'neg'], dtype='float64')                 # Stores r values for pos & neg models for each iteration
    mask_counts = pd.DataFrame(index=range(iterations), columns=['pos', 'neg'], dtype='int16')              # Stores number of pos/neg edges for each iteration
    
    pos_edges = np.zeros((n_edges, 2), dtype=np.int32)                                              # Stores edge number and count of inclusion in masks /folds*iterations
    pos_edges[:,0] = np.arange(n_edges)+1
    neg_edges = np.zeros((n_edges, 2), dtype=np.int32)
    neg_edges[:,0] = np.arange(n_edges)+1
    
    pos_preds = pd.DataFrame(index=subj_list, columns=range(iterations), dtype='float64')                          # Stores predicted change scores for each iteration (scatter plot)
    neg_preds = pd.DataFrame(index=subj_list, columns=range(iterations), dtype='float64')
    
    # Run CPM & store results for each iteration
    for i in range(iterations):
        
        behav_obs_pred, all_masks = cpm_wrapper(behav_data, rng)
        
        mask_counts.iloc[i, 0] = all_masks['pos'].sum().astype(np.int16)
        mask_counts.iloc[i, 1] = all_masks['neg'].sum().astype(np.int16)
        
        pos_edges[:, 1] += all_masks['pos'].sum(axis=0).astype(np.int32)
        neg_edges[:, 1] += all_masks['neg'].sum(axis=0).astype(np.int32)
        
        r_vals.iloc[i, 0] = sp.stats.pearsonr(behav_obs_pred.iloc[:,2], behav_obs_pred.iloc[:,0])[0]
        r_vals.iloc[i, 1] = sp.stats.pearsonr(behav_obs_pred.iloc[:,2], behav_obs_pred.iloc[:,1])[0]
        
        pos_preds.iloc[:,i] = behav_obs_pred.iloc[:,0].astype(np.float64)
        neg_preds.iloc[:,i] = behav_obs_pred.iloc[:,1].astype(np.float64)

    # Remove unused edges from mask
    pos_edges = pos_edges[pos_edges[:,1] != 0]
    neg_edges = neg_edges[neg_edges[:,1] != 0]
    
    # Titles
    pos_edges = pd.DataFrame(pos_edges, columns=['edge', 'count'])
    neg_edges = pd.DataFrame(neg_edges, columns=['edge', 'count'])
    
    # Round floats to save space
    r_vals = r_vals.round(5)
    pos_preds = pos_preds.round(2)
    neg_preds = neg_preds.round(2)
    
    return(r_vals, mask_counts, pos_edges, neg_edges, pos_preds, neg_preds)

# ----------- PREPARE DATA --------------

# Load behavioural data
data = pd.read_csv('~/STAGES_clinical/' + scale + '_data.csv')
end_scores = data[data['fkTimePointID'] == int(timepoint)]

# Read in subjects list to account for missing data / motion exclusion
if fc_type == 'fc_baseline':
    with open('/home/ipop0003/kg98/isaac/scripts/subjects_ses-1.txt') as f:
        subj_list = f.read().splitlines()
        
elif fc_type == 'fc_change_bl3m':
    with open('/home/ipop0003/kg98/isaac/scripts/subjects_ses-2.txt') as f:
        subj_list = f.read().splitlines()
        
else:
    print('Invalid FC type. Please enter "fc_baseline" or "fc_change_bl3m"')
    sys.exit()
    
# Equalise subjects across models with same scale and FC type, allowing FWE correction by assigning the same null permutation indices
if fc_type == 'fc_change_bl3m' and timepoint == '4':
    subj_list.remove('sub-038')

# Select only subjects with complete clinical data
subj_list = np.array([value for value in subj_list if value in end_scores['BIDS_ID'].values])

# Get clinical scores at baseline and follow-up
end_scores = end_scores[end_scores['BIDS_ID'].isin(subj_list)]
end_scores = end_scores[scale + '_score'].values
base_scores = data[data['fkTimePointID'] == 1 & data['BIDS_ID'].isin(subj_list)]
base_scores = base_scores[scale + '_score'].values

# Calculate proportional change scores
change_scores = (end_scores - base_scores) / base_scores
change_scores = change_scores[~np.isnan(change_scores)]

# Get confounder data for partial correlation - age, sex, head motion
demog_data = pd.read_excel('~/STAGES_clinical/STAGES_demog_clin_massive.xlsx', sheet_name='Demogs')
demog_data = demog_data[['BIDS_ID', 'Age_bl', 'sex']]
demog_data = demog_data.drop_duplicates(subset='BIDS_ID')

motion_data = pd.read_excel('~/kg98/isaac/scripts/motion_exclusion.xlsx')
motion_data = motion_data[motion_data['session'] == 1]
motion_data = motion_data[['BIDS_ID', 'Mean_FD_Jenk']]

all_behav_data = pd.DataFrame({'BIDS_ID': subj_list, 'change_score': change_scores}, index=subj_list)
all_behav_data = all_behav_data.merge(demog_data, on='BIDS_ID', how='inner')
all_behav_data = all_behav_data.merge(motion_data, on='BIDS_ID', how='inner')
all_behav_data = all_behav_data.set_index('BIDS_ID')

# Save behavioural data
filedir = '~/kg98/isaac/cpm/' + scale + '/' + timepoint + '/' + fc_type + '/'
all_behav_data.to_excel(filedir + 'change_observed.xlsx')

# Load functional connectivity data
if fc_type == 'fc_baseline':
    all_fc_data = read_in_matrices(subj_list, '_ses-1', Path('/home/ipop0003/kg98/isaac/conmats/dt_AROMA_8Phys-4GMR_bpf/ses-1'))
  
elif fc_type == 'fc_change_bl3m':
    all_fc_data = read_in_matrices(subj_list, '_change', Path('/home/ipop0003/kg98/isaac/conmats/dt_AROMA_8Phys-4GMR_bpf/change'))
  
else:
    print('Invalid FC type. Please enter "fc_baseline" or "fc_change_bl3m"')
    sys.exit()

n_edges = all_fc_data.shape[1]

# ----------- RUN ITERATIONS OF CPM --------------

# Time and a Word
start_time = time.time()
print(f'Running 100 iterations of CPM! \nScale: {scale} \nTimepoint: {timepoint} \nFC type: {fc_type}')

# Run 100 iterations of CPM
r_vals, mask_counts, pos_edges, neg_edges, pos_preds, neg_preds = run_cpm(all_behav_data)

# Save results as .xlsx files
r_vals.to_excel(filedir + 'r_vals.xlsx', index=False)
mask_counts.to_excel(filedir + 'mask_counts.xlsx', index=False)
  
with pd.ExcelWriter(filedir + 'edge_incidence.xlsx') as writer:
    pos_edges.to_excel(writer, sheet_name='Positive', index=False)
    neg_edges.to_excel(writer, sheet_name='Negative', index=False)

with pd.ExcelWriter(filedir + 'change_predictions.xlsx') as writer:
    pos_preds.to_excel(writer, sheet_name='Positive', header=False)
    neg_preds.to_excel(writer, sheet_name='Negative', header=False)

end_time = time.time()
print('Elapsed time: ' + str(round((end_time - start_time),2)) + ' seconds')

# ----------- NULL MODEL --------------

# Generate permutations
perm = np.zeros((permutations, len(subj_list)), dtype=np.int8)
rng = np.random.default_rng(seed=1)

for i in range(permutations):
    perm[i] = rng.permutation(len(subj_list))

null_r = pd.DataFrame(index=range(permutations), columns=['permutation', 'pos', 'neg'], dtype='float64')

# Permute behavioural data and rerun CPM
def run_permutation(i):
    
    perm_all_behav = all_behav_data.take(perm[i])
    perm_all_behav.index = subj_list
    
    try: 
        r_vals, mask_counts, pos_edges, neg_edges, pos_preds, neg_preds = run_cpm(perm_all_behav)
        
        # # Save results as .xlsx files TEMPORARILY COMMENTED OUT (ISAAC)
        # filedir = '~/kg98/isaac/cpm/' + scale + '/' + timepoint + '/' + fc_type #+ '/' + cpm_kwargs['corr_type']
        # 
        # r_vals.to_excel(filedir + '/null_r_vals/permutation_' + str(i+1) + '.xlsx', index=False)
        # mask_counts.to_excel(filedir + '/null_mask_counts/permutation' + str(i+1) +'.xlsx', index=False)
        # 
        # with pd.ExcelWriter(filedir + '/null_edge_incidence/permutation_' + str(i+1) + '.xlsx') as writer:
        #     pos_edges.to_excel(writer, sheet_name='Positive', index=False)
        #     neg_edges.to_excel(writer, sheet_name='Negative', index=False)
            
        # Take mean of each column in r_vals
        r_vals = r_vals.mean()
        
        return i, r_vals
      
    # This protects null_r when a permutation fails (SVD doesn't converge in np.polyfit)
    except np.linalg.LinAlgError as e:
      print(f"Error in permutation {i+1}: {e}")
      return i, None

def main():
    
    with ProcessPoolExecutor(max_workers=24) as executor:
        futures = [
            executor.submit(
                run_permutation, i
            ) for i in range(permutations)
        ]
        
        for future in as_completed(futures):
            i, r_vals = future.result()
            if r_vals is not None:
                null_r.loc[i] = [i+1, r_vals['pos'], r_vals['neg']]
                print(f'Finished permutation {i+1}')
                null_r.to_excel('~/kg98/isaac/cpm/' + scale + '/' + timepoint + '/' + fc_type + '/null_r_means.xlsx', index=False)

# ------------ RUN NULL MODEL --------------
print('Running null model...')
start_time = time.time()

if __name__ == '__main__':
    main()

end_time = time.time()
print(f"Execution time: {end_time - start_time:.2f} seconds")
