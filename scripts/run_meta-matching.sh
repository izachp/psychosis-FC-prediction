#!/bin/env bash

#SBATCH --job-name=Meta-matching
#SBATCH --account=
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=20000
#SBATCH --mail-user=
#SBATCH --mail-type=FAIL
#SBATCH --output=

# NOTE: Above settings may need to be adjusted

# Activate python environment with all required packages for meta-matching, as outlined in the README
source /path/to/meta-matching/environment

# Define number of parallel processes to use for null modelling
export max_workers=24

# Run model
python $repo/scripts/meta-matching.py \
       $repo $meta \
       $scale $timepoint \
       $max_workers
