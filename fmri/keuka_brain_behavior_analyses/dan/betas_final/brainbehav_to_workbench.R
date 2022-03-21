library(reticulate)
library(RNifti)
library(glue)
#write.csv(to_plot %>% select(-rhs, -effect), file="fmri_brainbehavior_parcel_betas_ptfce.csv", row.names=FALSE)
use_python("/usr/local/bin/python3")

nib <- import("nibabel")
np <- import("numpy")
pd <- import("pandas")

# loc <- RNifti::readNifti("Schaefer2018_200Parcels_7Networks_order.dscalar.nii")
vol <- RNifti::readNifti("Schaefer2018_200Parcels_7Networks_order_fonov_1mm_ants.nii.gz")

template <- nib$load('Schaefer2018_200Parcels_7Networks_order.dscalar.nii')
template_data <- template$get_fdata() #do operations here

#to_plot <- data.table::fread("fmri_brainbehavior_parcel_betas_ptfce.csv")
#to_plot <- data.table::fread("fmri_brainbehavior_parcel_betas_200.csv")
to_plot <- data.table::fread("fmri_brainbehavior_parcel_entropy_betas_200.csv")
terms <- grep("fmri_beta", unique(to_plot$term), value=TRUE)

# make roi label image (isn't this what the original dscalar is anyway?)
# new_data <- np$zeros(dim(template_data))
# 
# for (idx in unique(roi_df$mask_value)) {
#   new_data[template_data==idx] <- roi_df$mask_value[roi_df$mask_value == idx]
# }
# 
# new_cii <- nib$cifti2$Cifti2Image(new_data, template$header)
# new_cii$to_filename(paste0('Schaefer_ptfce_rois.dscalar.nii'))

#tmp <- niftiHeader(vol)
#tmp$dim[1] <- 4

# create both surface and volume images
for (tt in terms) {
  new_data <- np$zeros(dim(template_data))
  new_stat <- vol
  new_stat <- new_stat*0
  new_p <- new_stat # p-value
  
  roi_df <- to_plot[term == tt, ]
  
  for (idx in unique(roi_df$mask_value)) {
    new_data[template_data==idx] <- roi_df$statistic[roi_df$mask_value == idx]
    new_stat[vol == idx] <- roi_df$statistic[roi_df$mask_value == idx]
    new_p[vol == idx] <- 1 - roi_df$p.value[roi_df$mask_value == idx] # store 1-p for threshold
  }
  
  new_cii <- nib$cifti2$Cifti2Image(new_data, template$header)
  new_cii$to_filename(file.path("workbench", paste0(make.names(tt), 'entropy_200.dscalar.nii')))
  
  #writeNifti(new_stat, file=file.path("volumetric", paste0(make.names(tt), '.nii.gz')))
  out_f <- file.path("volumetric", paste0(make.names(tt), 'entropy_200.nii.gz'))
  stat_f <- tempfile(fileext = ".nii.gz")
  p_f <- tempfile(fileext = ".nii.gz")
  writeNifti(new_stat, file=stat_f)
  writeNifti(new_p, file=p_f)
  system(glue("/Users/hallquist/abin/3dTcat -overwrite -prefix {out_f} {stat_f} {p_f}"))
  system(glue("/Users/hallquist/abin/3drefit -relabel_all_str 'tstat 1-p' {out_f}"))
  unlink(c(stat_f, p_f))
}

