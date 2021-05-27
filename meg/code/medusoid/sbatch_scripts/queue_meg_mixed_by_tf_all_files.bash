#!/bin/bash

epochs="RT clock"
basedir="/bgfs/adombrovski/tfr_rds1/"
encode=FALSE
rt_predict=TRUE
alignment=RT
domain=tf
group_sensors=FALSE
for e in $epochs; do
slist=$( find ${basedir}/${e}/grouped_tf -type f -iname "*_group.RDS" )
arr=($slist)
for ss in ${arr(1)}; do
encode=TRUE
rt_predict=TRUE
alignment=${e}
sbatch --export=epoch=${e},sourcefile="$ss" sbatch_meg_mixed_by_tf_one_file.bash
done
done