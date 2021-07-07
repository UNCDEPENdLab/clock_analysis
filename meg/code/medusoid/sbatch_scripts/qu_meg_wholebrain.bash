#!/bin/bash
epochs="RT"
basedir="/bgfs/adombrovski/tfr_rds1"
# Can set this to whatever you want the count to be
countby=10
for e in $epochs; do
#################################
# OLD CODE FOR GRABBING ALL FILES
#declare -a flist
#flist=$( find ${basedir}/${e} -type f -iname "*freq_split*" )
#flist=$( find ${basedir}/${e} -type f -iname "*temporal_l_group_freq_split_f_03.536*" )
#for ss in $flist; do
#################################
# Initialized at 1, will be incremented in lo
counter=1
# Contact the above two strings to get the correct dir
epochdir=$basedir/${e}
echo $basedir/${e}
#echo $epochdir
# Get the filecount
filecount=`ls $epochdir/*freq_t*  | wc -l`
# filecount=12
# Iterate over the files in a chunked manner
until [ $counter -gt $filecount ]
do
# alignment=$e
#sbatch --export=epoch=${e},sourcefile="$ss",encode=TRUE,rt_predict=TRUE,domain=tf,group_sensors=FALSE sbatch_meg_mixed_by_tf_one_file.bash
############################

# Run the sbatch for this chunk
sbatch --export=epoch=${e},sourcefilestart="$counter",incrementby="$countby" sbatch_meg_mixed_wholebrain.bash
############################
# Print statement for testing
echo $counter
# increment the counter
((counter=counter+$countby))
done
done