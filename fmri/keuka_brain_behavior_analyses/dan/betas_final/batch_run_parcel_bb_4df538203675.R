if (!require(pacman)) { install.packages('pacman'); library(pacman) }
pacman::p_load(fmri.pipeline, tidyverse, data.table, sfsmisc)
if (file.exists('parcel_input_snapshot.RData')) {
  load('parcel_input_snapshot.RData')
} else {
  stop('Cannot load input environment object: parcel_input_snapshot.RData')
}
to_plot <- mixed_by_betas('/Volumes/GoogleDrive/My Drive/SCEPTIC_fMRI/wholebrain_betas/L1m-echange/Schaefer_444_final_2009c_2.3mm_cope_l2.csv.gz', labels_df, trial_df, mask_file = 'Schaefer_444_final_2.3mm.nii.gz',
rhs_form = fmri.pipeline:::named_list(int, slo), ncores = 4, afni_dir = '/Users/alexdombrovski/abin',
split_on = c('l1_cope_name', 'l2_cope_name', 'mask_value'))
