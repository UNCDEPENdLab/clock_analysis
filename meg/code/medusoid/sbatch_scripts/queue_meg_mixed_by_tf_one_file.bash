
#!/bin/bash

#SBATCH -p general
#SBATCH -N 1
#SBATCH --mem 64g
#SBATCH -n 8
#SBATCH -t 2:00:00

[ -z "$sourcefile" ] && echo "No sourcefile env variable passed in" && exit 1
[ -z "$epoch" ] && echo "No epoch env variable passed in" && exit 1

sid=$( basename "$sourcefile" )
sid=${sid/.rds/}
cd ../ #sbatch scripts live in subdirectory

R CMD BATCH --no-save --no-restore time_freq.R tf_rout_files/time_freq_${epoch}_${sid}.Rout