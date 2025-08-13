#!/bin/env bash

# This script runs Multilayer Meta-Matching (Chen et al. 2024) using functional coupling and clinical outcome data from STAGES.
# A total of 4 models are run in the main analysis, which differ in terms of the:
#
#     1. Clinical scale predicted:
#        - Brief Psychiatric Rating Scale (BPRS)
#        - Social and Occupational Functioning Assessment Scale (SOFAS)
#
#     2. Timeframe of clinical change predicted:
#        - Baseline to 6 months
#        - Baseline to 12 months
#
# For each model, 100 iterations of 4-fold cross-validated meta-matching are run, along with 1000 null permutations (each with 20 iterations) for significance testing.

# ------------------------------------------ SETUP ------------------------------------------
# Specify path to meta-matching repo
export meta=path/to/Meta_matching_models/rs-fMRI

# Specify path to this repo
export repo=path/to/psychosis-FC-prediction

# Set finer details of model
# (supplementary results options commented; sometimes only one can be swapped in at a time)
export preprocessing=dt_AROMA_8Phys-4GMR_bpf # or dt_AROMA_8Phys_bpf
export outcome_type=change_scores            # or lm_slopes
export group=all                             # or _placebo/_medication
scales=('sofas' 'bprs')                      # or bprs_pos

export run_null_model=0                      # 1 to run, but hours until done is no fun
# -------------------------------------------------------------------------------------------

# Loop through models to make results directories and run meta-matching
for scale in ${scales[@]}; do
    export scale=$scale

    if [ $outcome_type = change_scores ]; then

		if [ $group = all ]; then

			for timepoint in 4 5; do
				export timepoint=$timepoint

				# Check if model has already been run
				if [ -f $repo/results/meta-matching/$scale/$timepoint/null_r_means.xlsx ]; then
					echo 'Meta-matching results for '$scale $timepoint $preprocessing' already exist!'

				else mkdir -p $repo/results/meta-matching/$scale/$timepoint
					# Submit job
					sbatch $repo/scripts/run_meta-matching.sh

				fi

			done

			# Handle group option
		elif [[ $group = _placebo || $group = _medication ]]; then
			export timepoint=4

			if [ -f $repo/results/meta-matching/$scale/$timepoint/null_r_means_${preprocessing}$group.txt ]; then
				echo 'Meta-matching results for '$group $scale $timepoint $preprocessing' already exist!'
			
			else mkdir -p $repo/results/meta-matching/$scale/$timepoint
				# Submit job
				sbatch $repo/scripts/run_meta-matching.sh
			fi
			
		else echo 'Error: group must be all, _placebo, or _medication'
		fi

	# Handle lm_slopes option
    elif [ $outcome_type = lm_slopes ]; then
		export timepoint=lm_slopes

		if [ -f $repo/results/meta-matching/$scale/$timepoint/null_r_means_${preprocessing}.txt ]; then
			echo 'Meta-matching results for '$scale $timepoint $preprocessing' already exist!'
			
		else mkdir -p $repo/results/meta-matching/$scale/$timepoint
			# Submit job
			sbatch $repo/scripts/run_meta-matching.sh
		fi

    else echo 'Error: outcome_type must be change_scores or lm_slopes'
    fi

done
