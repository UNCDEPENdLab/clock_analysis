#!/usr/bin/env sh

#PBS -A mnh5174_a_g_hc_default
#PBS -l nodes=1:ppn=20
#PBS -l walltime=10:00:00
#PBS -j oe
#PBS -M michael.hallquist@psu.edu
#PBS -m abe
#PBS -W group_list=mnh5174_collab

cd $PBS_O_WORKDIR

export G=/gpfs/group/mnh5174/default

module use $G/sw/modules

module load r/3.5.0
module load fsl/5.0.11
module load afni/18.1.15

export PATH

#the fsl_lvl1_pipeline_file environment variable must be passed in through qsub, which is picked up by the R script

R CMD BATCH --no-save --no-restore execute_fsl_lvl1_pipeline.R
