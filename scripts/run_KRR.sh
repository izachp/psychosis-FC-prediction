#!/bin/bash

#SBATCH --job-name=Pope_KRR
#SBATCH --account=
#SBATCH --time=30
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=6G
#SBATCH --cpus-per-task=6
#SBATCH --mail-user=
#SBATCH --mail-type=FAIL
#SBATCH --output=slurmout/job-%j.out

# NOTE: Above settings may need to be adjusted, especially if running the null models
# On our system, null modelling used 16 CPUs, 6G per CPU, 18 hours
# If results stop outputting but the job keeps running, try increasing mem-per-cpu

# Set max CPUs engaged at once, for parallel processing during null modelling
# Must be less than --cpus-per-task, probably half or a third
export max_cpus=8

module purge
module load matlab/r2019b

cd $repo/scripts

# First wrangle the FC & outcome data into a format that KRR can use
echo 'Making inputs for '$scale $timepoint $fc_type $preprocessing $parcellation $outcome_type $group
#matlab -nodisplay -nojvm -nosplash -nodesktop -r "make_KRR_inputs('$repo', '$scale', '$timepoint', '$fc_type', '$preprocessing', '$parcellation', '$outcome_type', '$group'); exit;"

echo 'Running KRR!'
#matlab -nodisplay -nojvm -nosplash -nodesktop -r "KRR('$indir', '$cbig', '$preprocessing', '$parcellation', '$group'); exit;" > /dev/null 2>&1 # Suppress output

# Run null model if specified
if [ $run_null_model = '1' ]; then

    echo 'Running null model!'

    # Loop through permutations, running MATLAB in parallel
    for perm in {1..1000}; do

	mkdir -p $indir/${preprocessing}_${parcellation}/nulls/$perm
	
	matlab -nodisplay -nojvm -nosplash -nodesktop -r "KRR_null('$indir', '$cbig', '$perm', '$preprocessing', '$parcellation', '$group'); exit;" > /dev/null 2>&1 & # Run in background & suppress output to log

	# Halt looping if at max capacity
	echo 'Commenced permutation '$perm
	if [[ $(jobs -r -p | wc -l) -ge $max_cpus ]]; then
            wait -n
	fi
    done

    wait

    # Make a single file with the 1000 null r_mean values
    source $repo/.venv/bin/activate
    python combine_krr_nulls.py $indir $preprocessing $parcellation $group

fi

echo 'All done!'



 

