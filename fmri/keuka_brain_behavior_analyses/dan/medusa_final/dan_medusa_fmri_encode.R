library(tidyverse)
library(lme4)
library(data.table)

out_dir <- "/Volumes/GoogleDrive/My Drive/SCEPTIC_fMRI/dan_medusa"
#repo_directory <- "~/code/clock_analysis"
repo_directory <- "~/Data_Analysis/clock_analysis"

visuomotor_long <- TRUE # what we want to load in load_medusa_data_dan.R
source(file.path(repo_directory, "fmri/keuka_brain_behavior_analyses/dan/load_medusa_data_dan.R"))

# mixed_by call
source("~/Data_Analysis/r_packages/fmri.pipeline/R/mixed_by.R")

# helper function to compile list of formulae
named_list <- function(...) {
  vnames <- as.character(match.call())[-1]
  return(setNames(list(...), vnames))
}

# no need to analyze time points with tons of missingness -- subset to times with plenty of data
# also drop out any missing decon data since that is the DV
clock_visuomotor_long %>% group_by(evt_time) %>%
  summarise(isna=sum(is.na(decon_interp)))

clock_visuomotor_long <- clock_visuomotor_long %>% filter(evt_time >= -4 & evt_time <= 6 & !is.na(decon_interp))

clock_visuomotor_long_online %>% group_by(evt_time) %>%
  summarise(isna=sum(is.na(decon_interp)))

clock_visuomotor_long_online <- clock_visuomotor_long_online %>% filter(!is.na(decon_interp) & evt_time <= 4)

rt_visuomotor_long %>% group_by(evt_time) %>%
  summarise(isna=sum(is.na(decon_interp)))

rt_visuomotor_long <- rt_visuomotor_long %>% filter(evt_time >= -4 & evt_time <= 6 & !is.na(decon_interp))

alignment <- "clock"
#alignment <- "rt"
#alignment <- "clock_online"

setDT(trial_df)

message("Merging")
if (alignment=="clock") {
  # subset to columns of interest
  trial_df <- trial_df %>%
    dplyr::select(id, run, run_trial, trial_neg_inv_sc, rt_csv_sc, rt_lag_sc, pe_max_lag,
                  v_entropy_wi, v_entropy_wi_change_lag, kld3_lag, v_max_wi, abs_pe_lag, outcome_lag) %>%
    mutate(log_kld3_lag = log(kld3_lag + .00001))
  
  d <- merge(trial_df, clock_visuomotor_long, by = c("id", "run", "run_trial"))
  d <- d %>% tidyr::separate(visuomotor_side, into=c("vm_gradient", "side"), sep="_")  
} else if (alignment == "clock_online") {
  # subset to columns of interest
  trial_df <- trial_df %>%
    dplyr::select(id, run, run_trial, trial_neg_inv_sc, rt_csv_sc, rt_lag_sc, pe_max_lag,
                  v_entropy_wi, v_entropy_wi_change_lag, kld3_lag, v_max_wi, abs_pe_lag, outcome_lag) %>%
    mutate(log_kld3_lag = log(kld3_lag + .00001))
  
  d <- merge(trial_df, clock_visuomotor_long_online, by = c("id", "run", "run_trial"))
  d <- d %>% tidyr::separate(visuomotor_side, into=c("vm_gradient", "side"), sep="_") 
} else if (alignment == "rt") {
  # subset to columns of interest
  trial_df <- trial_df %>%
    dplyr::select(id, run, run_trial, trial_neg_inv_sc, rt_csv_sc, rt_lag_sc, pe_max,
                  v_entropy_wi, v_entropy_wi_change, kld3, v_max_wi, abs_pe, outcome) %>%
    mutate(log_kld3 = log(kld3 + .00001))
  
  d <- merge(trial_df, rt_visuomotor_long, by = c("id", "run", "run_trial"))
  d <- d %>% tidyr::separate(visuomotor_side, into=c("vm_gradient", "side"), sep="_")
}

rm(rt_visuomotor_long)
rm(clock_visuomotor_long)
gc()

# cor(trial_df$v_entropy_wi, trial_df$v_entropy_wi_change, use="pairwise")
# vv <- trial_df %>% select(v_entropy_wi, v_entropy_wi_change) %>%
#   mutate(el1=lag(v_entropy_wi, 1), el2=lag(v_entropy_wi, 2)) %>% na.omit()
# ccf(vv$v_entropy_wi, vv$v_entropy_wi_change)
#
# cor(vv$v_entropy_wi_change, vv$el1)
# cor(vv$v_entropy_wi_change, vv$el2)
# acf(vv$v_entropy_wi)


if (alignment == "clock" || alignment == "clock_online") {
  # basal analysis
  # vmax_wi: targeted analysis to demonstrate MT+ (vmax-positive) versus DAN (vmax-negative)
  enc_clock_base <- formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + v_max_wi +
                           v_entropy_wi + v_entropy_wi_change_lag + abs_pe_lag + outcome_lag +
                           (1 | id) )
  
  enc_clock_rslope <- formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + v_max_wi +
                             v_entropy_wi + v_entropy_wi_change_lag + abs_pe_lag + outcome_lag +
                             (abs_pe_lag + v_entropy_wi + v_entropy_wi_change_lag + v_max_wi | id) )
  
  
  # clock-aligned using kld3_lag
  enc_clock_kld <- formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + v_max_wi + 
                          v_entropy_wi + v_entropy_wi_change_lag + abs_pe_lag + outcome_lag + kld3_lag +
                          (1 | id) )
  
  # log version of kld given the strong skew
  enc_clock_logkld <- formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + v_max_wi + 
                             v_entropy_wi + v_entropy_wi_change_lag + abs_pe_lag + outcome_lag + log_kld3_lag +
                             (1 | id) )

  # signed PE
  enc_clock_pe <- formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + v_max_wi + 
                         v_entropy_wi + v_entropy_wi_change_lag + pe_max_lag +
                         (1 | id) )
  
  # abs_pe x outcome interaction
  enc_clock_int <- formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + v_max_wi +
                          v_entropy_wi + v_entropy_wi_change_lag + abs_pe_lag * outcome_lag +
                          (1 | id) )
  
  # trial interactions
  enc_clock_trial <- formula(~ v_max_wi*run_trial + v_entropy_wi*run_trial + v_entropy_wi_change_lag*run_trial + abs_pe_lag*run_trial +
                               rt_csv_sc*run_trial + rt_lag_sc + outcome_lag +
                               (1 | id) )
  
  flist <- named_list(enc_clock_base, enc_clock_rslope, enc_clock_kld, enc_clock_logkld, enc_clock_pe, enc_clock_int, enc_clock_trial)
} else if (alignment == "rt") {
  # basal analysis
  enc_rt_base <- formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + v_max_wi +
                           v_entropy_wi + v_entropy_wi_change + abs_pe + outcome +
                           (1 | id) )
  # killed: rt_vmax
  
  # vmax_wi: targeted analysis to demonstrate MT+ (vmax-positive) versus DAN (vmax-negative)
  enc_rt_rslope <- formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + v_max_wi +
                             v_entropy_wi + v_entropy_wi_change + abs_pe + outcome +
                             (abs_pe + v_entropy_wi + v_entropy_wi_change + v_max_wi | id) )
  
  enc_rt_kld <- formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + v_max_wi + 
                          v_entropy_wi + v_entropy_wi_change + abs_pe + outcome + kld3 +
                          (1 | id) )
  
  enc_rt_logkld <- formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + v_max_wi + 
                          v_entropy_wi + v_entropy_wi_change + abs_pe + outcome + log_kld3 +
                          (1 | id) )

  # signed PE
  enc_rt_pe <- formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + v_max_wi + 
                         v_entropy_wi + v_entropy_wi_change + pe_max +
                         (1 | id) )
  
  # abs_pe x outcome interaction
  enc_rt_int <- formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + v_max_wi +
                           v_entropy_wi + v_entropy_wi_change + abs_pe * outcome +
                           (1 | id) )
  
  # trial interactions
  enc_rt_trial <- formula(~ rt_csv_sc*run_trial + rt_lag_sc + v_max_wi*run_trial +
                           v_entropy_wi*run_trial + v_entropy_wi_change*run_trial + abs_pe*run_trial + outcome +
                           (1 | id) )
  
  flist <- named_list(enc_rt_base, enc_rt_rslope, enc_rt_kld, enc_rt_logkld, enc_rt_pe, enc_rt_int, enc_rt_trial)
}


splits <- c("vm_gradient", "side", "evt_time")

message("Running mixed_by")
ddf <- mixed_by(d, outcomes = "decon_interp", rhs_model_formulae = flist,
                split_on = splits, scale_predictors = c("abs_pe", "abs_pe_lag", "pe_max", "pe_max_lag", "run_trial"),
                tidy_args = list(effects = c("fixed", "ran_vals", "ran_pars", "ran_coefs"), conf.int = TRUE), 
                calculate = c("parameter_estimates_reml"), ncores = 16, refit_on_nonconvergence = 5, padjust_by = "term")

saveRDS(ddf, file=file.path(out_dir, paste0(alignment, "_encode_medusa_fmri.rds")))
