#!/usr/bin/env sh

#PBS -A m5m_a_g_sc_default
#PBS -l nodes=1:ppn=4
#PBS -l walltime=120:00:00
#PBS -l pmem=8gb
#PBS -j oe
#PBS -M michael.hallquist@psu.edu
#PBS -m abe
#PBS -W group_list=mnh5174_collab 

cd $PBS_O_WORKDIR

module use /gpfs/group/mnh5174/default/sw/modules
module load r/3.6.0

export PATH

R CMD BATCH --no-save --no-restore $to_run
