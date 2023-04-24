# !/bin/bash
# create mask of just VTA, VS, Hippo

set -ex
3dcalc -overwrite -a "masks/bilateral_striatum_tight_7Networks_2.3mm.nii.gz<3,8>" -expr "step(a)" -prefix "masks/bilateral_vs_2.3mm.nii.gz"


fslmaths masks/pauli_da_midbrain_2.3mm.nii.gz \
  -add masks/bilateral_vs_2.3mm.nii.gz \
  -add /proj/mnhallqlab/projects/clock_analysis/fmri/hippo_voxelwise/masks/long_axis_l_2.3mm.nii.gz \
  -add /proj/mnhallqlab/projects/clock_analysis/fmri/hippo_voxelwise/masks/long_axis_r_2.3mm.nii.gz \
  -bin \
  masks/vtavship_mask_2.3mm.nii.gz
