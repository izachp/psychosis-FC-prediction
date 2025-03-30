import os
import pandas as pd
import numpy as np
import sys

permutations = 1000

# Get all files in null directory
outdir = sys.argv[1]
files = os.listdir(outdir + '/nulls/')

# Ensure the null model is complete
if len(files) < permutations:
  print('Null model in ' + outdir + ' has not completed ' + permutations + ' permutations')
  sys.exit(1)

# Sort files by number in filename (e.g: permutation_50.xlsx)
files = sorted(files, key=lambda x: int(x.split('perm')[1].split('.')[0]))

# Create a dataframe to store the means
null_r = pd.DataFrame(index=range(permutations), columns=['r_mean'], dtype='float64')

# Read in files
t = 0
for file in files:
  
  # Get null r values
  df = pd.read_excel(outdir + '/nulls/' + file, header=None)
  
  # Put mean of null r values in dataframe
  null_r.loc[t, 'r_mean'] = df.mean().values[0]
  
  t += 1
  
# Save to excel
null_r.to_excel(outdir + '/null_r_means.xlsx', index=False, header=False)
