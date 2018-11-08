#!/usr/bin/env sh

#PBS -A mnh5174_a_g_hc_default
#PBS -l nodes=1:ppn=30
#PBS -l walltime=72:00:00
#PBS -j oe
#PBS -M michael.hallquist@psu.edu
#PBS -m abe

export G=/gpfs/group/mnh5174/default

module use $G/sw/modules

ni_tools="$G/lab_resources"

#location of MRI template directory for preprocessing
MRI_STDDIR="${ni_tools}/standard"

#add preprocessing scripts that may be called in this pipeline
PATH="${ni_tools}/c3d-1.1.0-Linux-x86_64/bin:${ni_tools}/fmri_processing_scripts:${ni_tools}/fmri_processing_scripts/autopreproc:${ni_tools}/bin:${PATH}"

export PATH MRI_STDDIR

#env
cd $PBS_O_WORKDIR

module load fsl/5.0.11

if [ ! -r "$torun" ]; then
    echo "Cannot find fsf file: $torun"
    exit 1
fi

#use local fork-based feat to run all slices simultaneously
${G}/lab_resources/bin/feat_parallel "$torun"

#agglomerate into single AFNI file
${ni_tools}/fmri_processing_scripts/feat_lvl2_to_afni.R --gfeat-dir ${torun/.fsf/.gfeat} --no_subjstats --no_varcope --stat_outfile ${torun/.fsf/_gfeat_stats}
