#!/bin/env bash

# This script runs Multilayer Meta-Matching (Chen et al. 2024) using functional coupling and clinical outcome data from STAGES. A total of 4 models are run, which differ in terms of the:
#
#     1. Clinical scale predicted:
#        - Brief Psychiatric Rating Scale (BPRS)
#        - Social and Occupational Functioning Assessment Scale (SOFAS)
#
#     2. Timeframe of clinical change predicted:
#        - Baseline to 6 months
#        - Baseline to 12 months
#
# For each model, 100 iterations of 4-fold cross-validated multilayer meta-matching are run, along with 1000 null permutations (each with 20 iterations) for significance testing.

# Specify path to meta-matching repo
export meta=/path/to/Meta_matching_models/rs-fMRI/utils

# Specify path to this repo
export repo=/path/to/repo
mkdir $repo/meta-matching

# Loop through models, making results directories and running models
for scale in sofas bprs; do

    export scale=$scale
    mkdir $repo/meta-matching/$scale

    for timepoint in 4 5; do

	export timepoint=$timepoint
	mkdir $repo/meta-matching/$scale/$timepoint

	# Check if model has already been run
	if [ -f $repo/meta-matching/$scale/$timepoint/null_r_means.xlsx ]; then

	    echo 'All results for '$scale $timepoint' already exist!'	

	else sbatch $repo/scripts/run_meta-matching.sh

	fi

    done

done
