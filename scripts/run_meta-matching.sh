#!/bin/env bash

#SBATCH --job-name=Pope_Meta
#SBATCH --account=
#SBATCH --time=12
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=6G
#SBATCH --mail-user=
#SBATCH --mail-type=FAIL
#SBATCH --output=slurmout/job-%j.out

# NOTE: Above settings may need to be adjusted, especially if running the null models
# On our system, null modelling used 6 CPUs, 14G per CPU, 10 hours
# If results stop outputting but the job keeps running, try increasing mem-per-cpu

source $repo/.venv/bin/activate

# Define number of parallel processes to use for null modelling
# (set this equal to --cpus-per-task above)
export max_workers=6

# Run model
python $repo/scripts/meta-matching.py \
       $repo $meta \
       $scale $timepoint \
       $preprocessing $outcome_type $group \
       $run_null_model $max_workers
