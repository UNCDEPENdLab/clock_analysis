library(tidyverse)
library(lme4)
library(data.table)
library(readxl)
library(fmri.pipeline)
out_dir <- "/Volumes/GoogleDrive/My Drive/SCEPTIC_fMRI/explore_medusa"
repo_directory <- "~/code/clock_analysis"
# repo_directory <- "~/Data_Analysis/clock_analysis"

# visuomotor_long <- TRUE # what we want to load in load_medusa_data_dan.R
# source(file.path(repo_directory, "fmri/keuka_brain_behavior_analyses/dan/load_medusa_data_dan.R"))



emm = T # extract EMMEANS estimates, e.g. for hi/lo abs(PE)

# mixed_by call
source("~/code/fmri.pipeline/R/mixed_by.R")

# helper function to compile list of formulae
named_list <- function(...) {
  vnames <- as.character(match.call())[-1]
  return(setNames(list(...), vnames))
}

# # no need to analyze time points with tons of missingness -- subset to times with plenty of data
# # also drop out any missing decon data since that is the DV
# clock_visuomotor_long %>% group_by(evt_time) %>%
#   summarise(isna=sum(is.na(decon_interp)))
# 
# clock_visuomotor_long <- clock_visuomotor_long %>% filter(evt_time >= -4 & evt_time <= 6 & !is.na(decon_interp))
# 
# clock_visuomotor_long_online %>% group_by(evt_time) %>%
#   summarise(isna=sum(is.na(decon_interp)))
# 
# clock_visuomotor_long_online <- clock_visuomotor_long_online %>% filter(!is.na(decon_interp) & evt_time <= 4)
# 
# rt_visuomotor_long %>% group_by(evt_time) %>%
#   summarise(isna=sum(is.na(decon_interp)))
# 
# rt_visuomotor_long <- rt_visuomotor_long %>% filter(evt_time >= -4 & evt_time <= 6 & !is.na(decon_interp))

# move to new 400 labels
labels <- setDT(read_excel("/Volumes/GoogleDrive/My Drive/SCEPTIC_fMRI/dan_medusa/schaefer_400_remap/MNH DAN Labels 400 Good Only 47 parcels.xlsx")) %>%
  rename(roi_num7=roi7_400) %>%
  select(roi_num7, mnh_label_400, network7_400, network17_400, parcel_group, hemi) %>% mutate(atlas_value = as.integer(roi_num7))

# rt_visuomotor_long <- fread("~/data/explore_medusa_400/exp_decon_rt_aligned.csv.gz") 
# forty_seven <- unique(labels$atlas_value)
# rt_visuomotor_long_400_47 <- rt_visuomotor_long %>% filter(atlas_value %in% forty_seven)
# fwrite(rt_visuomotor_long_400_47, "explore_decon_rt_aligned_47_400.csv")  
setwd("~/data/explore_medusa_400/")

files <-  gsub("//", "/", list.files(pattern = "interpolated", full.names = T))
message(paste0("Found ", length(files), " files."))
csl <- lapply(files, function(x) {
  print(x)
  df <- fread(x) 
  df$id <- as.integer(stringi::stri_extract_first(x, regex = "\\d+"))
  # if (class(df)=="list") {
  df$run <- stringi::stri_extract_last(x, regex = "\\d+")
  # } else if (ncol(df)==3) {
  # df <- df$fit_df
  # }
  return(df)
})
rt_visuomotor_long <- data.table::rbindlist(csl)

setwd(out_dir)
forty_seven <- unique(labels$atlas_value)
rt_visuomotor_long_400_47 <- rt_visuomotor_long %>% filter(atlas_value %in% forty_seven) %>% mutate(run = as.integer(run))
rt_visuomotor_long_400_47 <- rt_visuomotor_long_400_47 %>% inner_join(labels, by = "atlas_value")
# save
setwd(file.path(paste0(out_dir, "/data")))
saveRDS(rt_visuomotor_long, file = "explore_rt_decon_all_444_parcels.rds")
saveRDS(rt_visuomotor_long_400_47, file = "explore_rt_decon_dan_400_47.rds")
# rt_visuomotor_long <- rt_visuomotor_long %>% inner_join(dan_labels, by = "atlas_value") %>% filter(!is.na(Stream))

setwd("~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/dan")

source("get_trial_data.R")
source("medusa_final/plot_medusa.R") # careful -- plot_medusa sets out_dir, need to reset
out_dir <- "/Volumes/GoogleDrive/My Drive/SCEPTIC_fMRI/explore_medusa"
source("~/code/fmri.pipeline/R/mixed_by.R")

trial_df <- setDT(get_trial_data(dataset = "explore", repo_directory = "~/OneDrive - University of Pittsburgh/Documents/SCEPTIC_fMRI/EXPLORE_Medusa/"))  

trial_df <- trial_df %>% dplyr::select(id, run, run_trial, trial, trial_neg_inv_sc, rt_csv_sc, rt_lag_sc, pe_max, rew_om_c, abs_pe_c, abspexrew,
                v_entropy_wi, v_entropy_wi_change, kld3, v_max_wi, abs_pe, outcome) %>%
  mutate(log_kld3 = log(kld3 + .00001)) #,
         # id = as.character(id))
# rt_visuomotor_long <- rt_visuomotor_long %>% inner_join(trial_df)
# rt_visuomotor_long <- rt_visuomotor_long %>% filter(evt_time > -3 & evt_time < 6)

#alignment <- "clock"
#alignment <- "clock_online"
alignment <- "rt"

# setDT(trial_df)

# save objects and add subject_df
setwd(file.path(paste0(out_dir, "/data")))
sub_df <- readRDS("./explore_n146.rds") %>% select(-ipipds_total, -neoffi_total)

message("Merging")
if (alignment=="clock") {
  # subset to columns of interest
  trial_df <- trial_df %>%
    dplyr::select(id, run, run_trial, trial_neg_inv_sc, rt_csv_sc, rt_lag_sc, pe_max_lag,
                  v_entropy_wi, v_entropy_wi_change_lag, kld3_lag, v_max_wi, abs_pe_lag, outcome_lag) %>%
    mutate(log_kld3_lag = log(kld3_lag + .00001))
  trial_df <- inner_join(trial_df, sub_df, by = "id")
  d <- merge(trial_df, clock_visuomotor_long, by = c("id", "run", "trial"))
  d <- d %>% tidyr::separate(visuomotor_side, into=c("vm_gradient", "side"), sep="_")  
} else if (alignment == "clock_online") {
  # subset to columns of interest
  trial_df <- trial_df %>%
    dplyr::select(id, run, run_trial, trial_neg_inv_sc, rt_csv_sc, rt_lag_sc, pe_max_lag,
                  v_entropy_wi, v_entropy_wi_change_lag, kld3_lag, v_max_wi, abs_pe_lag, outcome_lag) %>%
    mutate(log_kld3_lag = log(kld3_lag + .00001))
  trial_df <- inner_join(trial_df, sub_df, by = "id")
  d <- merge(trial_df, clock_visuomotor_long_online, by = c("id", "run", "trial"))
  d <- d %>% tidyr::separate(visuomotor_side, into=c("vm_gradient", "side"), sep="_") 
} else if (alignment == "rt") {
  # subset to columns of interest
  trial_df <- trial_df %>%
    dplyr::select(id, run, run_trial, trial, trial_neg_inv_sc, rt_csv_sc, rt_lag_sc, pe_max, rew_om_c, abs_pe_c, abspexrew,
                  v_entropy_wi, v_entropy_wi_change, kld3, v_max_wi, abs_pe, outcome) %>%
    mutate(log_kld3 = log(kld3 + .00001))
  trial_df <- inner_join(trial_df, sub_df, by = "id")
  d <- merge(trial_df, rt_visuomotor_long_400_47, by = c("id", "run", "trial"))
  d <- d %>% rename(side = "hemi")
  # d <- d %>% tidyr::separate(visuomotor_side, into=c("vm_gradient", "side"), sep="_")
}

# rm(rt_visuomotor_long)
# rm(clock_visuomotor_long)
# rm(clock_visuomotor_long_online)
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
  
  # centered reward and abs pe to compute interaction
  # enc_clock_int_cent <- formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + v_max_wi +
  #                            v_entropy_wi + v_entropy_wi_change_lag + abs_pe_c + rew_om_c + abspexrew +
  #                            (1 | id) )
  
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
  # explore interactions of abs_pe and v_entropy_wi_change with cognitive variables, psychopathology and group
  enc_rt_base_age_edu_grp <- formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + v_max_wi + v_entropy_wi + 
                                            v_entropy_wi_change*education_yrs + 
                                            v_entropy_wi_change*age + 
                                            v_entropy_wi_change*Group + 
                                       abs_pe*education_yrs + 
                                       abs_pe*age + 
                                       abs_pe*Group + 
                                            outcome +
                           (1 | id) )
  enc_rt_base_age_edu_grp_mmse <- formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + v_max_wi + v_entropy_wi + 
                                       v_entropy_wi_change*education_yrs + 
                                       v_entropy_wi_change*age + 
                                       v_entropy_wi_change*Group + 
                                         v_entropy_wi_change*mmse + 
                                       abs_pe*education_yrs + 
                                       abs_pe*age + 
                                       abs_pe*Group +
                                         abs_pe*mmse + 
                                       outcome +
                                       (1 | id) )
  
  enc_rt_base_age_edu_upps_anx <- formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + v_max_wi + v_entropy_wi + 
                                            v_entropy_wi_change*education_yrs + 
                                            v_entropy_wi_change*age + 
                                            v_entropy_wi_change*uppsp_negative_urgency + 
                                              v_entropy_wi_change*uppsp_lack_of_premeditation + 
                                            v_entropy_wi_change*anxiety_dx + 
                                            abs_pe*education_yrs + 
                                            abs_pe*age + 
                                            abs_pe*uppsp_negative_urgency +
                                            abs_pe*uppsp_lack_of_premeditation + 
                                              abs_pe*anxiety_dx + 
                                            outcome +
                                            (1 | id) )
  
  enc_rt_mmse_grp_upps_anx_rslope <- formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + v_max_wi + v_entropy_wi + 
                                            v_entropy_wi_change*mmse + 
                                            v_entropy_wi_change*uppsp_negative_urgency + 
                                            v_entropy_wi_change*uppsp_lack_of_premeditation + 
                                            v_entropy_wi_change*anxiety_dx + 
                                            v_entropy_wi_change*Group + 
                                            abs_pe*mmse + 
                                            abs_pe*uppsp_negative_urgency +
                                            abs_pe*uppsp_lack_of_premeditation + 
                                            abs_pe*anxiety_dx + 
                                            abs_pe*Group + 
                                            outcome +
                                            (v_entropy_wi_change + abs_pe | id) )
  
    
  
  # vmax_wi: targeted analysis to demonstrate MT+ (vmax-positive) versus DAN (vmax-negative)
  enc_rt_rslope <- formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + v_max_wi +
                             v_entropy_wi + v_entropy_wi_change + abs_pe + outcome +
                             (abs_pe + v_entropy_wi + v_entropy_wi_change + v_max_wi | id) )
  enc_rt_rslope_kld <- formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + v_max_wi +
                             v_entropy_wi + v_entropy_wi_change + abs_pe + outcome + kld3 +
                             (abs_pe + v_entropy_wi + v_entropy_wi_change + v_max_wi | id) )
  
  enc_rt_rslope_logkld <- formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + v_max_wi +
                                 v_entropy_wi + v_entropy_wi_change + abs_pe + outcome + log_kld3 +
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
  
  # signed PE
  enc_rt_pe_kld <- formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + v_max_wi + 
                         v_entropy_wi + v_entropy_wi_change + pe_max + kld3 +
                         (1 | id) )
  # signed PE w/o covariates
  enc_rt_pe_only <- formula(~ pe_max +
                         (1 | id) )

  enc_rt_pe_rslope <- formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + v_max_wi + 
                         v_entropy_wi + v_entropy_wi_change + pe_max +
                         (pe_max | id) )
  
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
  flist <- named_list(enc_rt_mmse_grp_upps_anx_rslope)
  # flist <- named_list(enc_rt_rslope, enc_rt_kld, enc_rt_pe)
  # flist <- named_list(enc_rt_base, enc_rt_rslope, enc_rt_kld, enc_rt_logkld, enc_rt_pe, enc_rt_int, enc_rt_int_cent, enc_rt_trial)
}


splits <- c("parcel_group", "side", "evt_time")

message("Running mixed_by")
ddf <- mixed_by(d, outcomes = "decon_mean", rhs_model_formulae = flist,
                split_on = splits, scale_predictors = c("abs_pe", "abs_pe_lag", "pe_max", "pe_max_lag", "run_trial",
                                                        names(sub_df)[8:25]),
                tidy_args = list(effects = c("fixed", "ran_vals", "ran_pars", "ran_coefs"), conf.int = TRUE), 
                calculate = c("parameter_estimates_reml"), ncores = 9, refit_on_nonconvergence = 5, padjust_by = "term"#,
                # emtrends_spec = list(
                #   abspe = list(outcome = "decon_mean", model_name = "enc_rt_int", var = "abs_pe", specs = c("outcome"))
                #)
                )

saveRDS(ddf, file=file.path(out_dir, paste0(alignment, "_", splits[1], "_grp_mmse_imp_anx_rslope_explore_400_47_encode_July20_2022.rds")))

# aggregate all results into one RDS
setwd(out_dir)
# files <-  gsub("//", "/", list.files(pattern = ".*400.*47.*rds", full.names = T))
# message(paste0("Found ", length(files), " files."))
# csl <- lapply(files, function(x) {
#   print(x)
#   df <- readRDS(x) %>% Filter(Negate(is.null),.)
#   df <- df$coef_df_reml
#   # df$id <- as.integer(stringi::stri_extract_first(x, regex = "\\d+"))
#   # # if (class(df)=="list") {
#   # df$run <- stringi::stri_extract_last(x, regex = "\\d+")
#   # } else if (ncol(df)==3) {
#   # df <- df$fit_df
#   # }
#   return(df)
# })
# ddf_all <- data.table::rbindlist(csl)
# saveRDS(ddf_all, file=file.path(out_dir, paste0(alignment, "_", splits[1], "_all_models_explore_400_47_encode_July19_2022.rds")))
# 
ddf <- ddf$coef_df_reml %>% dplyr::filter(evt_time <= 5) %>%
  filter(effect=="fixed") %>%
  group_by(term, model_name) %>%
  mutate(p_FDR=p.adjust(p.value, method="fdr")) %>%
  ungroup() %>% setDT()

plot_medusa(ddf, x="evt_time", y="estimate", ymin="estimate - std.error", ymax="estimate + std.error", color="parcel_group", facet_by="side",
            out_dir=file.path(out_dir, "explore_400_47"),  p.value="padj_BY_term")


dd <- ddf$emtrends_list$abspe
dd <- dd[,c(-4, -5)]
ggplot(dd, aes(x=evt_time, y=abs_pe.trend, ymin=abs_pe.trend-std.error, ymax=abs_pe.trend+std.error, color=vm_gradient)) + 
  geom_pointrange() + geom_line() + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
  facet_grid(outcome...1 ~ side)
