#!/bin/env bash

# This script runs Connectome-based Predictive Modelling (Shen et al. 2017) using functional coupling and clinical outcome data from STAGES.
# A total of 16 models are run in the main analysis, which differ in terms of the:
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
# For each model, 100 splits of 4-fold cross-validated CPM are run, along with 1000 null permutations (each with 100 splits) for significance testing.

# ------------------------------------------ SETUP ------------------------------------------
# Specify path to this repo
export repo=path/to/psychosis-FC-prediction

# Set finer details of model
# (supplementary results options commented; sometimes only one can be swapped in at a time)
export feature_selection_threshold=0.01      # or 0.05/0.001
export preprocessing=dt_AROMA_8Phys-4GMR_bpf # or dt_AROMA_8Phys_bpf
export parcellation=419_parc                 # or 328_parc
export outcome_type=change_scores            # or lm_slopes
export group=all                             # or _placebo/_medication
scales=('sofas' 'bprs')                      # or bprs_pos

export run_null_model=0 					 # 1 to run, but hours until done is no fun
# -------------------------------------------------------------------------------------------

# Loop through models to make results directories and run CPM
for scale in ${scales[@]}; do
    export scale=$scale

    if [ $outcome_type = change_scores ]; then
	
		if [ $group = all ]; then
			
			for timepoint in 4 5; do
				export timepoint=$timepoint 

				for fc_type in fc_baseline fc_change_bl3m; do
					export fc_type=$fc_type

					# Check if model has already been run
					if [ -f $repo/results/cpm/$scale/$timepoint/$fc_type/null_r_means_${preprocessing}_${parcellation}_${feature_selection_threshold}.txt ]; then
					
						echo 'CPM results for '$scale $timepoint $fc_type $preprocessing $parcellation $feature_selection_threshold' already exist!'
					
					else mkdir -p $repo/results/cpm/$scale/$timepoint/$fc_type
						# Submit job
						sbatch $repo/scripts/run_CPM.sh
					fi
					
				done

			done

		# Handle group option
		elif [[ $group = _placebo || $group = _medication ]]; then
			export timepoint=4
			
			for fc_type in fc_baseline fc_change_bl3m; do
				export fc_type=$fc_type

				if [ -f $repo/results/cpm/$scale/$timepoint/$fc_type/null_r_means_${preprocessing}_${parcellation}_${feature_selection_threshold}$group.txt ]; then
					
					echo 'CPM results for '$group $scale $timepoint $fc_type $preprocessing $parcellation $feature_selection_threshold' already exist!'
					
				else mkdir -p $repo/results/cpm/$scale/$timepoint/$fc_type
					# Submit job
					sbatch $repo/scripts/run_CPM.sh
				fi
			
			done
		else echo 'Error: group must be all, _placebo, or _medication'
		fi
	
	# Handle lm_slopes option
    elif [ $outcome_type = lm_slopes ]; then
		export timepoint=lm_slopes
		export fc_type=fc_baseline

		if [ -f $repo/results/cpm/$scale/$timepoint/$fc_type/null_r_means_${preprocessing}_${parcellation}_${feature_selection_threshold}.txt ]; then
			
			echo 'CPM results for '$scale $timepoint $fc_type $preprocessing $parcellation $feature_selection_threshold' already exist!'
			
		else mkdir -p $repo/results/cpm/$scale/$timepoint/$fc_type
			# Submit job
			sbatch $repo/scripts/run_CPM.sh
		fi

    else echo 'Error: outcome_type must be change_scores or lm_slopes'
    fi

done
