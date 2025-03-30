#!/bin/env bash

# This script runs Connectome-based Predictive Modelling (Shen et al. 2017) using functional coupling and clinical outcome data from STAGES. A total of 16 models are run, which differ in terms of the:
#
#     1. Clinical scale predicted:
#        - Brief Psychiatric Rating Scale (BPRS)
#        - Social and Occupational Functioning Assessment Scale (SOFAS)
#
#     2. Timeframe of clinical change predicted:
#        - Baseline to 6 months
#        - Baseline to 12 months
#
#     3. FC type used to predict:
#        - Baseline FC
#        - Change in FC (baseline to 3 months)
#
#     4. Feature selection type used:
#        - Positive (r > 0 for FC vs change scores in training set)
#        - Negative (r < 0 for FC vs change scores in training set)
#        * Note: Both are run within each call of CPM.py
#
# For each model, 100 splits of 4-fold cross-validation CPM are run, along with 1000 null permutations (each with 100 splits) for significance testing.

# Specify path to this repo
export repo=/path/to/repo
mkdir $repo/cpm

# Set feature selection threshold (Supps also test p=0.001 & p=0.05)
export $p_thresh=0.01

# Loop through models, making results directories and running CPM
for scale in sofas bprs; do

    export scale=$scale
    mkdir $repo/cpm/$scale

    for timepoint in 4 5; do

	export timepoint=$timepoint
	mkdir $repo/cpm/$scale/$timepoint

	for fc_type in fc_baseline fc_change_bl3m; do

	    export fc_type=$fc_type
	    mkdir $repo/cpm/$scale/$timepoint/$fc_type

	    # Check if model has already been run
	    if [ -f $repo/cpm/$scale/$timepoint/$fc_type/null_r_means.xlsx ]; then
		echo 'CPM results for '$scale $timepoint $fc_type' already exist!'
		
	    else sbatch $repo/scripts/run_CPM.sh

	    fi
		
	done

    done

done
