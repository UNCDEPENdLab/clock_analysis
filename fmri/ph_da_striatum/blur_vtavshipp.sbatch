#!/bin/bash

#SBATCH -N 1
#SBATCH -n 16
#SBATCH --time=8:00:00
#SBATCH --mem-per-cpu=15G
cd $SLURM_SUBMIT_DIR

module use /proj/mnhallqlab/sw/modules
#module load afni/21.3.04

module load fsl/6.0.4
module load afni/23.0.7

bash smooth_vta_vs_hippo_in_mask.bash
