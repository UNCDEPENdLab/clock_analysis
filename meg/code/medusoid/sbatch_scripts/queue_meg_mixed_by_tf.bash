#!/bin/bash

epochs="clock"
basedir="/bgfs/adombrovski/tfr_rds1"



for e in $epochs; do
declare -a flist
flist=$( find ${basedir}/${e} -type f -iname "*freq_split*" )
#flist=$( find ${basedir}/${e} -type f -iname "*temporal_l_group_freq_split_f_03.536*" )
for ss in $flist; do
# alignment=$e
sbatch --export=epoch=${e},sourcefile="$ss",encode=TRUE,rt_predict=TRUE,domain=tf,group_sensors=FALSE sbatch_meg_mixed_by_tf_one_file.bash
done

done
