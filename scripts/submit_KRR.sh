#!/bin/env bash

# This script runs Kernel Ridge Regression (Li et al. 2019) using functional coupling and clinical outcome data from STAGES. A total of 8 models are run, which differ in terms of the:
#
#     1. Clinical scale predicted:
#        - Brief Psychiatric Rating Scale (BPRS)
#        - Social and Occupational Functioning Assessment Scale (SOFAS)
#
#     2. Timeframe of clinical change predicted:
#        - Baseline to 6 months
#        - Baseline to 12 months
#
#     3. Functional connectivity type used:
#        - Baseline FC
#        - Change in FC (baseline to 3 months)
#
# For each model, 100 splits of 4-fold cross-validation KRR are run, along with 1000 null permutations (each with 50 splits) for significance testing.

# Specify path to CBIG repo with required KRR scripts
export cbig=/path/to/cbig

# Specify path to this repo
export repo=/path/to/repo
mkdir $repo/krr

# Loop through models, making directories/inputs and running models
for scale in bprs sofas; do

    export scale=$scale
    mkdir $repo/krr/$scale

    for timepoint in 4 5; do

	export timepoint=$timepoint
	mkdir $repo/krr/$scale/$timepoint

	for fc_type in fc_baseline fc_change_bl3m; do

	    export fc_type=$fc_type
	    mkdir $repo/krr/$scale/$timepoint/$fc_type

	    # Check if model has already been run
	    if [ -f $repo/krr/$scale/$timepoint/$fc_type/null_r_means.xlsx ]; then
		echo 'KRR results for '$scale $timepoint $fc_type' already exist!'
		
	    else sbatch $repo/scripts/run_KRR.sh

	    fi
		
	done

    done

done
