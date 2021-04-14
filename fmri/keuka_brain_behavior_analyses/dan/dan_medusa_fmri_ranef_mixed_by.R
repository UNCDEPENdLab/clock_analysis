library(tidyverse)
library(lme4)


medusa_dir = "~/Box/SCEPTIC_fMRI/dan_medusa/"
cache_dir = "~/Box/SCEPTIC_fMRI/dan_medusa/cache"
repo_directory <- "~/code/clock_analysis"
gc()

# load MEDUSA deconvolved data
# load(file.path(cache_dir, 'rt_dan_tall_ts.Rdata'))
# saveRDS(rt_comb, file = file.path(cache_dir, 'rt_dan_tall_ts.RDS'))
message("Reading in decons")
rt_comb <- readRDS(file.path(cache_dir, 'rt_dan_tall_ts.RDS'))

message("Reading in behavioral variables")
df <- readRDS(file.path(repo_directory, "fmri/keuka_brain_behavior_analyses/dan/fmri_trial_df_medusa.RDS"))
message("Merging")
d <- merge(df, rt_comb, by = c("id", "run", "run_trial"))
gc()
# mixed_by call
source("~/code/fmri.pipeline/R/mixed_by.R")
splits = c("stream", "side", "evt_time")
encode_formula = formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + scale(rt_vmax_lag)  + scale(rt_vmax_change) + 
                           v_entropy_wi + v_entropy_wi_change + kld3_lag  + v_max_wi  + scale(abs_pe) + outcome + (outcome + scale(abs_pe)|id))
message("Running mixed_by")
ddf <- mixed_by(d, outcomes = "decon_interp", rhs_model_formulae = encode_formula , split_on = splits,
                ncores = 8, refit_on_nonconvergence = 5, padjust_by = NULL,
                tidy_args = "ran_vals")
setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/dan/plots/rt_decode')
encode_results_fname = "rt_encode_output_streams_mixed_by_abs_pe_ranef.RDS"
saveRDS(file = encode_results_fname, ddf)
