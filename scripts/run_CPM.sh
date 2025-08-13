#!/bin/bash

#SBATCH --job-name=Pope_CPM
#SBATCH --account=
#SBATCH --time=12
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=4G
#SBATCH --mail-user=
#SBATCH --mail-type=FAIL
#SBATCH --output=slurmout/job-%j.out

# NOTE: Above settings may need to be adjusted, especially if running the null models
# On our system, null modelling took 18 CPUs, 4G per CPU, 6 hours
# If results stop outputting but the job keeps running, try increasing mem-per-cpu

# Activate python environment with required packages
source $repo/.venv/bin/activate

# Define number of parallel processes to use for null modelling
# (set this equal to --cpus-per-task above)
export max_workers=18

# Run model
python $repo/scripts/CPM.py \
       $repo \
       $scale $timepoint $fc_type \
       $feature_selection_threshold \
       $group $preprocessing $parcellation $outcome_type \
       $run_null_model $max_workers\
