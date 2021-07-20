#!/bin/bash

epochs="RT"
basedir="/bgfs/adombrovski/tfr_rds1"
outfname="meg_mixed_by_tf_ddf_wholebrain_entropy_change_rs_"
# Can set this to whatever you want the count to be
countby=10
step=5 # has to be a divisor of countby (10)
for e in $epochs; do
    #################################
    # OLD CODE FOR GRABBING ALL FILES
    #declare -a flist
    #flist=$( find ${basedir}/${e} -type f -iname "*freq_split*" )
    #flist=$( find ${basedir}/${e} -type f -iname "*temporal_l_group_freq_split_f_03.536*" )
    #for ss in $flist; do
    #################################
    # Initialized at 1, will be incremented in loop, should not need modified
    counter=1
    # Contact the above two strings to get the correct dir
    epochdir="$basedir/${e}"
    #echo $basedir/${e}
    #echo $epochdir
    # Get the filecount
    filecount=`ls $epochdir/*freq_t* | wc -l`
    # Iterate over the files in a chunked manner
    #c=1
    until [ $counter -gt $filecount ]
    do
      	# if the expected output file for this iteraction does not exist
        fpath="$epochdir/$outfname${e}$counter"
        #echo $fpath
        if ! test -f "$fpath"
        then
            ############################
            # Run the sbatch for this chunk
            sbatch --export=epoch=${e},sourcefilestart="$counter",incrementby="$step" sbatch_meg_mixed_wholebrain.bash
            sbatch --export=epoch=${e},sourcefilestart="$(($counter+$step))",incrementby="$step" sbatch_meg_mixed_wholebrain.bash

            # Print statement for testing
            #echo $fpath
            echo "$counter"
            echo "$(($counter+$step))"
            ############################
        # exit the if statement
        fi
    # increment the counter
    ((counter=counter+$countby))
        #((c=c+1))
    done
done
