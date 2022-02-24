library(tidyverse)
library(lme4)
library(data.table)

#medusa_dir = "~/Box/SCEPTIC_fMRI/dan_medusa/"
cache_dir <- "/Users/hallquist/OneDrive - University of North Carolina at Chapel Hill/medusa"
#cache_dir = "~/Box/SCEPTIC_fMRI/dan_medusa/cache"
#repo_directory <- "~/code/clock_analysis"
repo_directory <- "~/Data_Analysis/clock_analysis"
gc()

# load MEDUSA deconvolved data
# load(file.path(cache_dir, 'rt_dan_tall_ts.Rdata'))
# saveRDS(rt_comb, file = file.path(cache_dir, 'rt_dan_tall_ts.RDS'))
message("Reading in decons")
rt_comb <- readRDS(file.path(cache_dir, 'rt_dan_tall_ts.RDS')) %>%
  dplyr::select(-score_csv, -rt_vmax_change_next) %>% setDT()
#intersect(names(trial_df), names(rt_comb))

message("Reading in behavioral variables")
trial_df <- readRDS(file.path(repo_directory, "fmri/keuka_brain_behavior_analyses/dan/fmri_trial_df_medusa.RDS")) %>%
  dplyr::select(id, run, run_trial, trial_neg_inv_sc, rt_csv_sc, rt_lag_sc,
                v_entropy_wi, v_entropy_wi_change, kld3_lag, v_max_wi, abs_pe, outcome)

#for simplicity, change all matrix columns (wi-centering and scaling) to numeric
trial_df <- as.data.frame(lapply(trial_df, function(x) {
  if (inherits(x, "matrix")) { x <- as.vector(x) }
  return(x)
}))
setDT(trial_df)

message("Merging")
d <- merge(trial_df, rt_comb, by = c("id", "run", "run_trial"))
rm(rt_comb)
gc()

# cor(trial_df$v_entropy_wi, trial_df$v_entropy_wi_change, use="pairwise")
# vv <- trial_df %>% select(v_entropy_wi, v_entropy_wi_change) %>%
#   mutate(el1=lag(v_entropy_wi, 1), el2=lag(v_entropy_wi, 2)) %>% na.omit()
# ccf(vv$v_entropy_wi, vv$v_entropy_wi_change)
#
# cor(vv$v_entropy_wi_change, vv$el1)
# cor(vv$v_entropy_wi_change, vv$el2)
# acf(vv$v_entropy_wi)

# mixed_by call
source("~/Data_Analysis/r_packages/fmri.pipeline/R/mixed_by.R")

splits = c("stream", "side", "evt_time")
encode_formula <- formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + scale(rt_vmax_lag)  + scale(rt_vmax_change) +
                           v_entropy_wi + v_entropy_wi_change + kld3_lag  + v_max_wi  + scale(abs_pe) + outcome +
                           (outcome + scale(abs_pe) + v_entropy_wi | id) )

encode_formula_change <- formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + scale(rt_vmax_lag)  + scale(rt_vmax_change) +
                                   v_entropy_wi + v_entropy_wi_change + kld3_lag  + v_max_wi  + scale(abs_pe) + outcome +
                                   (outcome + scale(abs_pe) + v_entropy_wi + v_entropy_wi_change | id) )

encode_formula_kld <- formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + scale(rt_vmax_lag)  + scale(rt_vmax_change) +
                                   v_entropy_wi + v_entropy_wi_change + kld3_lag  + v_max_wi  + scale(abs_pe) + outcome +
                                   (outcome + scale(abs_pe) + v_entropy_wi + kld3_lag | id) )


# message("Running mixed_by")
ddf <- mixed_by(d, outcomes = "decon_interp", rhs_model_formulae = list(
  entropy_pe=encode_formula,
  entropy_change=encode_formula_change,
  entropy_kld=encode_formula_kld), split_on = splits,
  ncores = 16, refit_on_nonconvergence = 5, padjust_by = NULL,
  tidy_args = "ran_vals")

saveRDS(ddf, file="encode_ranefs_fmri.rds")
# setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/dan/plots/rt_decode')
# encode_results_fname = "rt_encode_output_streams_mixed_by_abs_pe_ranef.RDS"
# saveRDS(file = encode_results_fname, ddf)

