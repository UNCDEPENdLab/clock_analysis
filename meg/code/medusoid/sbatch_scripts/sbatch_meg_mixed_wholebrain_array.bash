#!/bin/bash

#SBATCH -p smp
#SBATCH -N 1
#SBATCH --mem 20g
#SBATCH -n 1
#SBATCH -t 23:00:00
#SBATCH --mail-user=ayd1@pitt.edu
#SBATCH --mail-type=ALL
#SBATCH -c 1
# [ -z "$incrementby" ] && echo "No incrementby env variable passed in" && exit 1
[ -z "$epoch" ] && echo "No epoch env variable passed in" && exit 1
[ -z "$regressor" ] && echo "No regressor env variable passed in" && exit 1
[ -z "$SLURM_ARRAY_TASK_ID" ] && echo "No file start env variable passed in" && exit 1
[ -z "$base_num" ] && echo "No base file number env variable passed in" && exit 1

declare -i SLURM_ARRAY_TASK_ID
declare -i base_num
export sourcefilestart=$((SLURM_ARRAY_TASK_ID+base_num))
#fid=$( basename "$sourcefile" )
#export filename=$fid

cd ../ #sbatch scripts live in subdirectory
R CMD BATCH --no-save --no-restore meg_medusa_mixed_by_wholebrain_crc.R time_freq_${epoch}_${regressor}_${sourcefilestart}_$(date +%s).Rout
