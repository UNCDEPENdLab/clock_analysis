#!/bin/bash

#SBATCH -p general
#SBATCH -N 1
#SBATCH --mem 212g
#SBATCH -n 36
#SBATCH -t 24:00:00
#SBATCH --mail-type=all
#SBATCH --mail-user=mnhallq@email.unc.edu

module use /proj/mnhallqlab/sw/modules
module load r/4.0.3_depend

#export epoch=RT

R CMD BATCH --no-save --no-restore meg_medusa_mixed_by_clock.R meg_medusa_mixed_by_clock_${alignment}_rtp${rt_predict}_enc${encode}.Rout
