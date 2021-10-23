#!/bin/bash

#SBATCH -p smp
#SBATCH -N 1
#SBATCH --mem 20g
#SBATCH -n 1
#SBATCH -t 5:00:00
#SBATCH --mail-user=ayd1@pitt.edu
#SBATCH --mail-type=ALL
#SBATCH -c 32
[ -z "$sourcefile" ] && echo "No sourcefile env variable passed in" && exit 1
[ -z "$epoch" ] && echo "No epoch env variable passed in" && exit 1

# for debug only:
# sourcefile="/bgfs/adombrovski/tfr_rds1/RT/frontal_r_group_freq_split_f_04.204.rds"

fid=$( basename "$sourcefile" )
export filename=$fid
export encode=TRUE
export rt_predict=TRUE
export domain=tf
export group_sensors=FALSE

cd ../ #sbatch scripts live in subdirectory
R CMD BATCH --no-save --no-restore meg_medusa_mixed_by_combined_crc.R time_freq_${epoch}_${fid}_$(date +%s).Rout
