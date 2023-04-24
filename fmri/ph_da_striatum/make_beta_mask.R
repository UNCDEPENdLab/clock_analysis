# make integer mask of regions to extract for VS, VTA, Hipp project

library(RNifti)

# simple script to bin hippocampal long axis betas into 12 bins for storage and analysis
library(data.table)
library(dplyr)

hip_l <- readNifti("/proj/mnhallqlab/projects/clock_analysis/fmri/hippo_voxelwise/masks/long_axis_l_2.3mm.nii.gz")
hip_r <- readNifti("/proj/mnhallqlab/projects/clock_analysis/fmri/hippo_voxelwise/masks/long_axis_r_2.3mm.nii.gz")

hl_pos <- which(hip_l > 0)
hl_vals <- hip_l[hl_pos]

hr_pos <- which(hip_r > 0)
hr_vals <- hip_r[hip_r > 0]

hb_vals <- c(hl_vals, hr_vals)


nbins <- 12
bin_cuts <- seq(min(hb_vals) - 1e-5, max(hb_vals) + 1e-5, length.out = nbins + 1)

# create 12 bins
hl_bins <- as.integer(cut(hl_vals, bin_cuts))
hr_bins <- as.integer(cut(hr_vals, bin_cuts))

hip_l[hl_pos] <- hl_bins
hip_r[hr_pos] <- hr_bins

writeNifti(hip_l, "masks/long_axis_l_12bins_2.3mm.nii.gz")
writeNifti(hip_r, "masks/long_axis_r_12bins_2.3mm.nii.gz")

# striatum and VTA
# 1 = L VS
# 2 = R VS
# 3 = L VTA
# 4 = R VTA

striatum <- readNifti("masks/bilateral_striatum_tight_7Networks_2.3mm.nii.gz")
striatum[striatum %in% c(1, 2, 4:7, 9, 10)] <- 0 # not relevant
striatum[striatum == 3] <- 1
striatum[striatum == 8] <- 2

vta <- readNifti("masks/pauli_da_midbrain_2.3mm.nii.gz")
vta_pos <- which(vta > 0, arr.ind=TRUE)
vta_loc <- voxelToWorld(vta_pos, vta)
r_vox <- which(vta_loc[, 1] > 0)
l_vox <- which(vta_loc[, 1] < 0)

vta[vta_pos[l_vox, ]] <- 3
vta[vta_pos[r_vox, ]] <- 4

vta_str <- striatum + vta
writeNifti(vta_str, "masks/vta_vs_2.3mm.nii.gz")

