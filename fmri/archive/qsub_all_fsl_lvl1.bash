#!/usr/bin/env sh

#PBS -l nodes=1:ppn=40
#PBS -l walltime=60:00:00
#PBS -A mnh5174_collab
#PBS -j oe
#PBS -M michael.hallquist@psu.edu
#PBS -m abe

env
cd $PBS_O_WORKDIR

module load R
module load fsl

R CMD BATCH run_fsl_lvl1.R
