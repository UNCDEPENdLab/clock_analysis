library(dplyr)
library(data.table)
library(glue)
ddf <- readRDS("/Volumes/GoogleDrive/My Drive/SCEPTIC_fMRI/dan_medusa/dan_parcels_clock_encode_medusa_fmri.rds")


analysis_dir <- "~/Data_Analysis/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/medusa_final"
#analysis_dir <- "/proj/mnhallqlab/projects/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/medusa_final"
setwd(analysis_dir)

source("../get_trial_data.R")
source("../betas_final/parcel_brain_behavior_functions.R")

mask_nifti <- "Schaefer2018_200Parcels_7Networks_order_fonov_1mm_ants.nii.gz"
mask_cifti <- "Schaefer2018_200Parcels_7Networks_order.dscalar.nii"

to_plot <- ddf$coef_df_reml %>%
  filter(effect=="fixed")

# fill_mask_with_stats(mask_nifti = mask_nifti, mask_cifti = mask_cifti, mask_col = "atlas_value",
#                      stat_dt = to_plot, subbrik_cols = c("estimate", "statistic"),
#                      subbrik_labels = c("beta", "t"),
#                      split_on = c("model_name", "term", "evt_time"), afni_dir = "~/abin", out_dir = getwd())


fill_mask_with_stats(mask_nifti = mask_nifti, mask_cifti = mask_cifti, mask_col = "atlas_value",
                     stat_dt = to_plot, stat_cols = c("estimate", "statistic"),
                     subbrik_labels = c("beta", "t"), stack_along = "evt_time",
                     split_on = c("model_name", "term"), afni_dir = "~/abin", 
                     img_prefix = NULL, out_dir = file.path(analysis_dir, "clock_encode"))

ddf <- readRDS("/Volumes/GoogleDrive/My Drive/SCEPTIC_fMRI/dan_medusa/dan_parcels_rt_encode_medusa_fmri.rds")

to_plot <- ddf$coef_df_reml %>%
  filter(effect=="fixed")

fill_mask_with_stats(mask_nifti = mask_nifti, mask_cifti = mask_cifti, mask_col = "atlas_value",
                     stat_dt = to_plot, stat_cols = c("estimate", "statistic"),
                     subbrik_labels = c("beta", "t"), stack_along = "evt_time",
                     split_on = c("model_name", "term"), afni_dir = "~/abin", 
                     img_prefix = NULL, out_dir = file.path(analysis_dir, "rt_encode"))


pdiff <- abs(to_plot$p.value - .005)
to_plot$statistic[order(pdiff)]
