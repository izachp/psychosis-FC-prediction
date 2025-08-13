import os
import numpy as np
import sys

permutations = 1000

indir = sys.argv[1]
preproc = sys.argv[2]
parc = sys.argv[3]
group = sys.argv[4]

# Get all files in null directory
if group == 'all':
  group = ''

files = os.listdir(f'{indir}/{preproc}_{parc}{group}/nulls/')

# Ensure the null model is complete
if len(files) < permutations:
  print(f'Null model in {indir} ({preproc}, {parc}, {group}) has only completed {len(files)} permutations, expected {permutations}.')
  sys.exit(1)

# Sort files by number in filename
files = sorted(files, key=lambda x: int(x.split('perm')[1].split('.')[0]))

# Initialise an array to hold the mean r values
null_r = np.full((permutations,), np.nan, dtype='float64')

# Read each file and compute the mean of the r values
for i, file in enumerate(files):
  
  # Get null r values
  r_vals = np.loadtxt(f'{indir}/{preproc}_{parc}{group}/nulls/{file}', dtype='float64')
  
  # Put mean of null r values in dataframe
  null_r[i] = np.mean(r_vals)
  
np.savetxt(f'{indir}/null_r_means_{preproc}_{parc}{group}.txt', null_r, fmt='%.5f')
