#!/bin/bash

#SBATCH --job-name=KRR
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --mem=20000
#SBATCH --cpus-per-task=24
#SBATCH --mail-user=
#SBATCH --mail-type=FAIL
#SBATCH --output=

# NOTE: Above settings may need to be adjusted

# Activate python environment for combining all null results at the end
source /path/to/environment

# Set max CPUs engaged at once, for parallel processing during null modelling
export max_cpus=6

module purge
module load matlab/r2019b

cd $repo/scripts

# First wrangle the FC & outcome data into a format that KRR can use
echo 'Making inputs!'
matlab -nodisplay -r "make_KRR_inputs('$repo', '$scale', '$timepoint', '$fc_type'); exit;"

# Set directory
outdir=$repo/krr/$scale/$timepoint/$fc_type

echo 'Running KRR!'
matlab -nodisplay -r "run_KRR_repeats('$outdir', '$cbig'); exit;"

# ---------- NULL MODEL ----------

echo 'Running KRR null model for '$scale' at timepoint '$timepoint' using '$fc_type' inputs'

mkdir $outdir/nulls

# Loop through permutations, running MATLAB in parallel
for perm in {1..1000}; do

    mkdir $outdir/$perm
    
    matlab -nodisplay -r "run_KRR_perms('$outdir', '$cbig', '$perm'); exit;" > /dev/null 2>&1 & # Suppress output to log

    # Halt looping if at max capacity
    echo 'Permutation '$perm' complete'
    if [[ $(jobs -r -p | wc -l) -ge $max_cpus ]]; then
        wait -n
    fi
done

wait

# Make a single file with the 1000 null r_mean values
python $repo/scripts/combine_krr_nulls.py $outdir

echo 'All done!'



 

