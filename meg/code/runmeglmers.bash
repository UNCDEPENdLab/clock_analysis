#!/usr/bin/env sh

#PBS -A mnh5174_a_g_hc_default
#PBS -l nodes=4:ppn=22:himem
#PBS -l walltime=54:00:00
#PBS -j oe
#PBS -m abe
#PBS -W group_list=mnh5174_collab

#env
cd $PBS_O_WORKDIR
module use /gpfs/group/mnh5174/default/sw/modules

module load gcc/5.3.1
module load r/3.5.0

export PATH

R CMD BATCH --no-save --no-restore clock_MEG_lmer.R
