library(dependlab)

masks <- c("l_accumbens_thr20.nii", "l_ant_hippo_p75_p92.nii", "l_post_hippo_p08_p25.nii", "r_accumbens_thr20.nii", "r_ant_hippo_p75_p92.nii", "r_post_hippo_p08_p25.nii", "vmpfc_clust1_z5.7_2009c.nii")
masks <- file.path("/gpfs/group/mnh5174/default/clock_analysis/fmri/hippo_voxelwise/hippo_dcm/masks", masks)
sdirs <- list.dirs("/gpfs/group/mnh5174/default/MMClock/MR_Proc", recursive=FALSE)
spm_dirs <- sapply(sdirs, function(x) {
  expect <- file.path(x, "mni_5mm_aroma", "sceptic-clock-feedback-v_entropy-pe_max-preconvolve_fse_groupfixed", "SPM.mat")
  if (file.exists(expect)) {
    return(dirname(expect))
  } else {
    return(NULL)
  }
})

spm_dirs <- unname(unlist(spm_dirs))

#aa <- spm_extract_anatomical_rois(spm_dirs, masks[-7], session=1, extent=0, adjust_F_index=1, contrast_index=NULL, ncores=4, spm_path="/gpfs/group/mnh5174/default/lab_resources/spm12", matlab_path="/opt/aci/sw/matlab/R2017b/bin")

#vmpfc with thresholding
aa <- spm_extract_anatomical_rois(spm_dirs, masks[7], session=1, extent=0, adjust_F_index=1, contrast_index=4, ncores=4, threshold=0.05,
  spm_path="/gpfs/group/mnh5174/default/lab_resources/spm12", matlab_path="/opt/aci/sw/matlab/R2017b/bin")

#contrasts
#1 effects of interest F
#2 clock
#3 clock_nopmod (first and last trials)
#4 clock x entropy pmod
#5 feedback
#6 feedback_nopmod
#7 feedback x pe
