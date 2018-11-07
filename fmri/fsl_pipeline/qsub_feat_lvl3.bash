#!/usr/bin/env sh

#PBS -A mnh5174_a_g_hc_default
#PBS -l nodes=1:ppn=30
#PBS -l walltime=72:00:00
#PBS -j oe
#PBS -M michael.hallquist@psu.edu
#PBS -m abe

export G=/gpfs/group/mnh5174/default

module use $G/sw/modules

#env
cd $PBS_O_WORKDIR

module load fsl/5.0.11

if [ ! -r "$torun" ]; then
    echo "Cannot find fsf file: $torun"
    exit 1
fi

#use local fork-based feat
${G}/lab_resources/bin/feat_parallel "$torun"
