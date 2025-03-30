import os
import sys
import random
import scipy
import torch
import pickle
import sklearn
import numpy as np
import warnings

warnings.filterwarnings("ignore")
sys.path.append(sys.argv[1]) 

from CBIG_model_pytorch import check_models_v20
from CBIG_model_pytorch import demean_norm
from sklearn.model_selection import train_test_split
from CBIG_model_pytorch import stacking
from CBIG_model_pytorch import multilayer_metamatching_infer
from scipy.stats.stats import pearsonr
import pandas as pd
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

# Set path to this repo & meta-matching repo
repo = sys.argv[1]
meta = sys.argv[2] + '/'

path_v20 = os.path.join(meta, 'v2.0')
model_v20_path = os.path.join(path_v20, 'models')
sys.path.append(os.path.join(meta, "utils"))

# Check whether meta-matching v2.0 model files exist and are up-to-date
check_models_v20(model_v20_path)

# Set model parameters
scale = sys.argv[3]
timepoint = sys.argv[4]
splits = 100
permutations = 1000

# get subject list
data = pd.read_excel(os.path.join(repo, 'cpm/' + scale + '/' + timepoint + '/fc_baseline/change_observed.xlsx'))
subjects = data['BIDS_ID'].values

# Get change scores
y_input = data['change_score'].values
y_input = y_input.reshape(-1, 1)

# Initialise group FC matrix
x_input = np.zeros((len(subjects), int(419*418/2)))

# Retain only lower triangle of FC matrix, and normalise
for i, subject in enumerate(subjects):
    conmat = np.loadtxt(os.path.join(repo, 'conmats/dt_AROMA_8Phys-4GMR_bpf/ses-1_419-reg/conmat_' + subject + '_ses-1.txt'))
    x_input[i, :] = conmat[np.tril_indices(419, k=-1)]
x_input = demean_norm(x_input)

def run_model(y_input=y_input):
    
    # Initialise results dataframe
    r_vals = pd.DataFrame(index=range(splits), columns=['r value'], dtype='float64')
    
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
        y_train_pred, y_names = multilayer_metamatching_infer(x_train, y_train, model_v20_path, dataset_names)
        y_test_pred, _ = multilayer_metamatching_infer(x_test, y_test, model_v20_path, dataset_names)
        
        # Stacking
        y_test_final=np.zeros((y_test_pred.shape[0], y_train.shape[1]))
        for i in range(y_train.shape[1]):
            # For each test phenotype, perform stacking by developing a KRR model
            y_test_temp, _ = stacking(y_train_pred, y_test_pred, y_train[:,i].view(), [0.00001, 0.0001, 0.001, 0.004, 0.007, 0.01, 0.04, 0.07, 0.1, 0.4, 0.7, 1, 1.5, 2, 2.5, 3, 3.5, 4, 5, 10, 15, 20])
            y_test_final[:,i] = y_test_temp.flatten()
        
        r_vals.iloc[n,:] = scipy.stats.stats.pearsonr(y_test_final[:,0], y_test[:,0])[0]
    r_vals = r_vals.round(5)
        
    return r_vals

# ----------- RUN MODEL --------------
        
# Time and a Word
start_time = time.time()
print(f'Running {splits} splits of multilayer meta-matching! \nScale: {scale} \nTimepoint: {timepoint}')

# Run splits of multilayer meta-matching model
r_vals = run_model()

# Save results
r_vals.to_excel(repo + '/meta-matching/' + scale + '/' + timepoint + '/r_vals.xlsx', index=False)

end_time = time.time()
print(f'Elapsed time for model: {end_time - start_time:.2f} seconds')

# ----------- NULL MODEL --------------
# Lower splits to reduce compute
splits = 20

# Initialise dataframe to store null model results
null_r = pd.DataFrame(index=range(permutations), columns=['permutation', 'r_val'], dtype='float64')
  
def run_permutation(i):
    try:
        # Permute behavioural data
        np.random.seed(i)
        perm_y_input = np.random.permutation(y_input)

        # Run model with permuted behavioural data
        r_vals = run_model(perm_y_input)
            
        # Calculate r_mean
        r_mean = r_vals['r value'].mean()
        
        print(f"Finished permutation {i}, r_mean: {r_mean}")
        return i, r_mean
      
    except Exception as e:
        print(f"Error in permutation {i}: {e}")
        return i, None

def main():
    with ProcessPoolExecutor(max_workers=sys.argv[5]) as executor:
        futures = [
            executor.submit(
                run_permutation, i
            ) for i in range(permutations)
        ]
        
        for future in as_completed(futures):
            try:
                i, r_mean = future.result()
                null_r.loc[i] = [i + 1, r_mean]
                null_r.to_excel(repo + '/meta-matching/' + scale + '/' + timepoint + '/null_r_means.xlsx', index=False)
                print(f"Permutation {i} complete")
            except Exception as e:
                print(f"Error: {e}")
                
# ------------ RUN NULL MODEL --------------
print(f'Running {permutations} permutations of null model!')
start_time = time.time()

if __name__ == '__main__':
    main()

end_time = time.time()
print(f"Execution time: {end_time - start_time:.2f} seconds")
