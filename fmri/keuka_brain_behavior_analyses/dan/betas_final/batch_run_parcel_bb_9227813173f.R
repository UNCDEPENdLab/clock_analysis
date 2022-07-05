if (!require(pacman)) { install.packages('pacman'); library(pacman) }
pacman::p_load(fmri.pipeline, tidyverse, data.table, sfsmisc)
if (file.exists('parcel_input_snapshot.RData')) {
  load('parcel_input_snapshot.RData')
} else {
  stop('Cannot load input environment object: parcel_input_snapshot.RData')
}
to_plot <- mixed_by_betas('/Volumes/GoogleDrive/My Drive/SCEPTIC_fMRI/wholebrain_betas/L1m-abspe_by_rew/Schaefer2018_200Parcels_7Networks_order_fonov_2.3mm_ants_cope_l2.csv.gz', labels_df, trial_df, mask_file = 'Schaefer_444_final_2.3mm.nii.gz',
rhs_form = fmri.pipeline:::named_list(int, slo), ncores = 2, afni_dir = '/Users/alexdombrovski/abin',
out_prefix = 'Schaefer_400_DAN_manual_labels_',
    
                   #                         rhs_model_formulae = fmri.pipeline:::named_list(int, slo), ncores = 16, afni_dir = '/proj/mnhallqlab/sw/afni',
                   #                         calculate = c('parameter_estimates_reml', 'fit_statistics'),
                   #                         beta_level = 2L, focal_contrast = 'overall', emtrends_spec = emtrends_spec,
                   split_on = c('l1_cope_name', 'l2_cope_name', 'mask_value'))
