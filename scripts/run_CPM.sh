#!/bin/env bash

#SBATCH --job-name=CPM
#SBATCH --account=
#SBATCH --time=03:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=4096
#SBATCH --mail-user=
#SBATCH --mail-type=FAIL
#SBATCH --output=

# NOTE: Above settings may need to be adjusted

# Activate python environment with required packages
source /path/to/environment

# Define number of parallel processes to use for null modelling
export max_workers=24

# Run model
python $repo/scripts/CPM.py \
       $repo \
       $scale $timepoint $fc_type \
       $p_thresh \
       $max_workers
