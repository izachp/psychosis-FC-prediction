#!/bin/env bash

# This script runs Kernel Ridge Regression (Li et al. 2019) using functional coupling and clinical outcome data from STAGES.
# A total of 8 models are run in the main analysis, which differ in terms of the:
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
# For each model, 100 splits of 4-fold cross-validation KRR are run, along with 1000 null permutations (each with 50 splits) for significance testing.

# ------------------------------------------ SETUP ------------------------------------------
# Specify path to CBIG repo
export cbig=path/to/CBIG

# Specify path to this repo
export repo=path/to/psychosis-FC-prediction

# Set finer details of model
# (supplementary results options commented; sometimes only one can be swapped in at a time)
export preprocessing=dt_AROMA_8Phys-4GMR_bpf # or dt_AROMA_8Phys_bpf
export parcellation=419_parc                 # or 328_parc
export outcome_type=change_scores            # or lm_slopes
export group=all                             # or _placebo/_medication
scales=('sofas' 'bprs')                      # or bprs_pos

export run_null_model=0                      # 1 to run, but hours until done is no fun
# -------------------------------------------------------------------------------------------

# Loop through models to make results directories and run KRR
for scale in ${scales[@]}; do
    export scale=$scale

    if [ $outcome_type = change_scores ]; then

		if [ $group = all ]; then
			
			for timepoint in 4 5; do	
				export timepoint=$timepoint

				for fc_type in fc_baseline fc_change_bl3m; do
					export fc_type=$fc_type
					export indir=$repo/results/krr/$scale/$timepoint/$fc_type

					# Check if model has already been run
					if [ -f $indir/null_r_means_${preprocessing}_${parcellation}.txt ]; then
						echo 'KRR results for '$scale $timepoint $fc_type $preprocessing $parcellation' already exist!'
					
					else mkdir -p $indir
						# Submit job
						sbatch $repo/scripts/run_KRR.sh
					fi
					
				done

			done

			# Handle group option
		elif [[ $group = _placebo || $group = _medication ]]; then
			export timepoint=4
			
			for fc_type in fc_baseline fc_change_bl3m; do
				export fc_type=$fc_type
				export indir=$repo/results/krr/$scale/$timepoint/$fc_type

				if [ -f $indir/null_r_means_${preprocessing}_${parcellation}$group.txt ]; then
					echo 'KRR results for '$group $scale $timepoint $fc_type $preprocessing $parcellation' already exist!'
					
				else mkdir -p $indir
					# Submit job
					sbatch $repo/scripts/run_KRR.sh
				fi
			
			done
			
		else echo 'Error: group must be all, _placebo, or _medication'
		fi
	
	# Handle lm_slopes option
    elif [ $outcome_type = lm_slopes ]; then
		export timepoint=lm_slopes
		export fc_type=fc_baseline
		export indir=$repo/results/krr/$scale/$timepoint/$fc_type

		if [ -f $indir/null_r_means_${preprocessing}_${parcellation}.txt ]; then
			
			echo 'KRR results for '$scale $timepoint $fc_type $preprocessing $parcellation' already exist!'
			
		else mkdir -p $indir
			# Submit job
			sbatch $repo/scripts/run_KRR.sh
		fi

    else echo 'Error: outcome_type must be change_scores or lm_slopes'
    fi

done