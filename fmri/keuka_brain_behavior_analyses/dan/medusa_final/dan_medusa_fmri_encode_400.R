library(tidyverse)
library(lme4)
library(data.table)
library(readxl)
library(fmri.pipeline) # mixed_by call: use fmri.pipeline installation

out_dir <- "~/OneDrive - University of Pittsburgh/Documents/SCEPTIC_fMRI/dan_medusa"
decon_dir <- "~/OneDrive - University of Pittsburgh/Documents/SCEPTIC_fMRI/dan_medusa"
schaefer_dir <- "~/code/schaefer_wb_parcellation"
# schaefer_dir <- "~/Data_Analysis/schaefer_wb_parcellation"
repo_directory <- "~/code/clock_analysis"
# repo_directory <- "~/Data_Analysis/clock_analysis"

# move to new 400 labels
labels <- setDT(read_excel("~/OneDrive - University of Pittsburgh/Documents/SCEPTIC_fMRI/dan_medusa/schaefer_400_remap/MNH DAN Labels 400 Good Only 47 parcels.xlsx")) %>%
  rename(roi_num7=roi7_400) %>%
  select(roi_num7, mnh_label_400, network7_400, network17_400, parcel_group, hemi)

# sanity checks -- network labels from different inputs match
#all.equal(labels$network17_400, labels$network17)
#all.equal(labels$network7_400, labels$network7)

visuomotor_long <- TRUE # what we want to load in load_medusa_data_dan.R
# visuomotor_long <- FALSE # no 4 visuomotor
# tall_only <- TRUE # get parcelwise
schaefer_400 <- TRUE # indicate that we want the 400-parcel data

cleanup <- T # whether to remove big dataframes from memory or keep them for model tweaks
wm <- TRUE # compare with working memory model

source(file.path(repo_directory, "fmri/keuka_brain_behavior_analyses/dan/load_medusa_data_dan.R"))

if (isTRUE(tall_only)) {
  # rename for internal consistency
  rt_comb <- rt_comb %>%
    select(id, run, run_trial, evt_time, visuomotor_side, decon_interp, label)
    
  rt_visuomotor_long <- rt_comb
  clock_visuomotor_long <- clock_comb %>%
    select(id, run, run_trial, evt_time, visuomotor_side, decon_interp, label)
  rm(rt_comb, clock_comb)
}

# uh oh, some mismatches in original -- corrected by 5Jul2022 reprocessing of event alignment to match final N=47 400 DAN
# setdiff(unique(labels$roi_num7), unique(rt_visuomotor_long$roi_num7))
# setdiff(unique(rt_visuomotor_long$roi_num7), unique(dan_labels$roi_num7))

emm = T # extract EMMEANS estimates, e.g. for hi/lo abs(PE)

#alignment <- "clock"
#alignment <- "clock_online"
alignment <- "rt"

message("Merging")
if (alignment=="clock") {
  # subset to columns of interest
  trial_df <- trial_df %>%
    dplyr::select(id, run, run_trial, trial_neg_inv_sc, rt_csv_sc, rt_lag_sc, pe_max_lag,
                  v_entropy_wi, v_entropy_wi_change_lag, kld3_lag, v_max_wi, abs_pe_lag, outcome_lag) %>%
    mutate(log_kld3_lag = log(kld3_lag + .00001))
  
  d <- merge(trial_df, clock_visuomotor_long, by = c("id", "run", "run_trial"))
  d <- d %>% tidyr::separate(visuomotor_side, into=c("vm_gradient17", "side"), sep="_")
} else if (alignment == "clock_online") {
  # subset to columns of interest
  trial_df <- trial_df %>%
    dplyr::select(id, run, run_trial, trial_neg_inv_sc, rt_csv_sc, rt_lag_sc, pe_max_lag,
                  v_entropy_wi, v_entropy_wi_change_lag, kld3_lag, v_max_wi, abs_pe_lag, outcome_lag) %>%
    mutate(log_kld3_lag = log(kld3_lag + .00001))
  
  d <- merge(trial_df, clock_visuomotor_long_online, by = c("id", "run", "run_trial"))
  d <- d %>% tidyr::separate(visuomotor_side, into=c("vm_gradient17", "side"), sep="_")
} else if (alignment == "rt") {
  # subset to columns of interest
  trial_df <- trial_df %>%
    dplyr::select(id, run, run_trial, trial, trial_neg_inv_sc, rt_csv_sc, rt_lag_sc, pe_max, rew_om_c, abs_pe_c, abspexrew,
                  v_entropy_wi, v_entropy_wi_change, kld3, v_max_wi, abs_pe, outcome, v_entropy_wi_full, v_entropy_wi_change_full) %>%
    mutate(log_kld3 = log(kld3 + .00001)) #%>% select(-trial)  # keep trial for merging
  if (wm) {wm_df <- readRDS(file.path(repo_directory, "fmri/keuka_brain_behavior_analyses/dan/wm_entropy/mmclock_wm_entropy.rds"))
  trial_df <- trial_df %>% inner_join(wm_df, by = c("id", "run", "trial")) %>% group_by(id, run) %>% arrange(id, run, trial) %>% mutate(
    wm_entropy_wi = scale(choice_entropy + rewom_entropy),
    wm_choice_entropy_wi  = scale(choice_entropy),
    wm_reward_entropy_wi  = scale(reward_entropy),
    wm_entropy_wi_change = lead(wm_entropy_wi) - wm_entropy_wi,
    wm_choice_entropy_wi_change = lead(wm_choice_entropy_wi) - wm_choice_entropy_wi,
    wm_reward_entropy_wi_change = lead(wm_reward_entropy_wi) - wm_reward_entropy_wi
  ) %>% ungroup()
  # inspect correlations: nothing noteworthy
  cormat <- psych::corr.test(trial_df %>% select(contains("entropy_wi")))
  pdf("sceptic_wm_entropy_correlations.pdf", height = 10, width = 10)
  corrplot::corrplot(corr = cormat$r, p.mat = cormat$p, order = "hclust", 
                     number.digits = 2, addCoef.col="white")
  dev.off()
  d}
  d <- merge(trial_df, rt_visuomotor_long, by = c("id", "run", "run_trial"))
  d <- d %>% tidyr::separate(visuomotor_side, into=c("vm_gradient17", "side"), sep="_")
  if (wm) {d <- d %>% filter(!is.na(wm_entropy_wi_change))}
}

if (cleanup) {
# clean up a few variables to reduce memory burden
rm(rt_visuomotor_long)
rm(clock_visuomotor_long)}
#rm(clock_visuomotor_long_online)
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
  
  # log version of kld given the strong skew
  enc_clock_logkld_rslope <- formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + v_max_wi + 
                                       v_entropy_wi + v_entropy_wi_change_lag + abs_pe_lag + outcome_lag + log_kld3_lag +
                                       (abs_pe_lag + v_entropy_wi + v_entropy_wi_change_lag + v_max_wi | id) )

  # signed PE
  enc_clock_pe <- formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + v_max_wi + 
                         v_entropy_wi + v_entropy_wi_change_lag + pe_max_lag +
                         (1 | id) )
  
  # abs_pe x outcome interaction
  enc_clock_int <- formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + v_max_wi +
                          v_entropy_wi + v_entropy_wi_change_lag + abs_pe_lag * outcome_lag +
                          (1 | id) )
  
  # centered reward and abs pe to compute interaction
  # enc_clock_int_cent <- formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + v_max_wi +
  #                            v_entropy_wi + v_entropy_wi_change_lag + abs_pe_c + rew_om_c + abspexrew +
  #                            (1 | id) )
  
  # trial interactions
  enc_clock_trial <- formula(~ v_max_wi*run_trial + v_entropy_wi*run_trial + v_entropy_wi_change_lag*run_trial + abs_pe_lag*run_trial +
                               rt_csv_sc*run_trial + rt_lag_sc + outcome_lag +
                               (1 | id) )
  
  flist <- fmri.pipeline:::named_list(enc_clock_base, enc_clock_rslope, enc_clock_kld, enc_clock_logkld,
                                      enc_clock_logkld_rslope, enc_clock_pe, enc_clock_int, enc_clock_trial)
} else if (alignment == "rt") {
  # basal analysis
  enc_rt_base <- formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + v_max_wi +
                           v_entropy_wi + v_entropy_wi_change + abs_pe + outcome +
                           (1 | id) )
  # basal analysis
  enc_rt_base_full <- formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + v_max_wi +
                           v_entropy_wi_full + v_entropy_wi_change_full + abs_pe + outcome +
                           (1 | id) )
  # basal model without entropy, only entropy change
  enc_rt_change <- formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + v_max_wi +
                           v_entropy_wi_change + abs_pe + outcome +
                           (1 | id) )
  
  # working memory model with a single, summed entropy variable
  enc_rt_wm <- formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + v_max_wi + 
                                wm_entropy_wi + wm_entropy_wi_change + abs_pe + outcome + 
                                (1 | id) )
  

  # separate choice and reward entropies
  enc_rt_wm_exp <- formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + v_max_wi + 
                         wm_choice_entropy_wi + wm_choice_entropy_wi_change + wm_reward_entropy_wi + wm_reward_entropy_wi_change +
                           abs_pe + outcome + 
                         (1 | id) )
  
  enc_rt_wm_exp_rslope <- formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + v_max_wi + 
                             wm_choice_entropy_wi + wm_choice_entropy_wi_change + wm_reward_entropy_wi + wm_reward_entropy_wi_change +
                             abs_pe + outcome + 
                             (wm_choice_entropy_wi + wm_choice_entropy_wi_change + wm_reward_entropy_wi + wm_reward_entropy_wi_change| id) )
  

  enc_rt_wm_exp_sceptic <- formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + v_max_wi + 
                                     v_entropy_wi + v_entropy_wi_change + 
                             wm_choice_entropy_wi + wm_choice_entropy_wi_change + wm_reward_entropy_wi + wm_reward_entropy_wi_change +
                             abs_pe + outcome + 
                             (1 | id) )
  # to be comprehensive, all random slopes
  enc_rt_wm_exp_sceptic_rslope <- formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + v_max_wi + 
                                     v_entropy_wi + v_entropy_wi_change + 
                                     wm_choice_entropy_wi + wm_choice_entropy_wi_change + wm_reward_entropy_wi + wm_reward_entropy_wi_change +
                                     abs_pe + outcome + 
                                     (v_entropy_wi + v_entropy_wi_change + 
                                        wm_choice_entropy_wi + wm_choice_entropy_wi_change + wm_reward_entropy_wi + wm_reward_entropy_wi_change| id) )

  enc_rt_wm_exp_sceptic_full <- formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + v_max_wi + 
                                     v_entropy_wi_full + v_entropy_wi_change_full + 
                                     wm_choice_entropy_wi + wm_choice_entropy_wi_change + wm_reward_entropy_wi + wm_reward_entropy_wi_change +
                                     abs_pe + outcome + 
                                     (1 | id) )
  
  # what the hell, pile WM entropies on top of SCEPTIC entropy
  
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
  
  enc_rt_wm_logkld <- formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + v_max_wi + 
                             wm_entropy_wi + wm_entropy_wi_change + abs_pe + outcome + log_kld3 +
                             (1 | id) ) # perhaps an unfair comparison since logkld may absorb WM choice entropy variance
  
  enc_rt_logkld_rslope <- formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + v_max_wi + 
                                    v_entropy_wi + v_entropy_wi_change + abs_pe + outcome + log_kld3 +
                                    (abs_pe + v_entropy_wi + v_entropy_wi_change + v_max_wi | id) )

  # signed PE
  enc_rt_pe <- formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + v_max_wi + 
                         v_entropy_wi + v_entropy_wi_change + pe_max +
                         (1 | id) )
  
  # abs_pe x outcome interaction
  enc_rt_int <- formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + v_max_wi +
                           v_entropy_wi + v_entropy_wi_change + abs_pe * outcome +
                           (1 | id) )
  
  enc_rt_int_cent <- formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + v_max_wi +
                                  v_entropy_wi + v_entropy_wi_change + abs_pe_c + rew_om_c + abspexrew +
                                  (1 | id) )
  
  # centered outcome approach to get absPE at average of reward/omissions
  enc_rt_int <- formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + v_max_wi +
                          v_entropy_wi + v_entropy_wi_change + abs_pe * outcome +
                          (1 | id) )
  
  
  # trial interactions
  enc_rt_trial <- formula(~ rt_csv_sc*run_trial + rt_lag_sc + v_max_wi*run_trial +
                           v_entropy_wi*run_trial + v_entropy_wi_change*run_trial + abs_pe*run_trial + outcome +
                           (1 | id) )
  
  # flist <- fmri.pipeline:::named_list(enc_rt_base, enc_rt_rslope, enc_rt_kld, enc_rt_logkld, 
  #                                     enc_rt_logkld_rslope, enc_rt_pe, enc_rt_int, enc_rt_int_cent, enc_rt_trial)
  flist <- fmri.pipeline:::named_list(enc_rt_wm_exp_rslope, enc_rt_wm_exp_sceptic_rslope)
  
}


# for parcelwise only, split on both vm_gradient and label so both are preserved in outputs
# splits <- c("vm_gradient17", "label", "side", "evt_time")

splits <- c("vm_gradient17", "side", "evt_time")

message("Running mixed_by")
ddf <- mixed_by(d, outcomes = "decon_interp", rhs_model_formulae = flist,
                split_on = splits, scale_predictors = c("abs_pe", "abs_pe_lag", "pe_max", "pe_max_lag", "run_trial"),
                tidy_args = list(effects = c("fixed", "ran_vals", "ran_pars", "ran_coefs"), conf.int = TRUE), 
                calculate = c("parameter_estimates_reml", "fit_statistics"), ncores = 9, refit_on_nonconvergence = 5, padjust_by = "term",
                emtrends_spec = list(
                  abspe = list(outcome = "decon_interp", model_name = "enc_rt_int", var = "abs_pe", specs = c("outcome"))
                ))

saveRDS(ddf, file=file.path(out_dir, paste0(alignment, "_400_final47_encode_medusa_fmri.rds")))
# saveRDS(ddf, file=file.path(out_dir, paste0(alignment, "_400_final47_encode_medusa_fmri_parcelwise.rds")))
# saveRDS(ddf, file=file.path(out_dir, paste0(alignment, "_400_final47_encode_medusa_fmri_wm_rslope.rds")))

table(ddf$coef_df_reml$vm_gradient17)
table(ddf$coef_df_reml$label)

# dd <- ddf$emtrends_list$abspe
# dd <- dd[,c(-4, -5)]
# ggplot(dd, aes(x=evt_time, y=abs_pe.trend, ymin=abs_pe.trend-std.error, ymax=abs_pe.trend+std.error, color=vm_gradient17)) + 
#   geom_pointrange() + geom_line() + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
#   facet_grid(outcome...1 ~ side)
