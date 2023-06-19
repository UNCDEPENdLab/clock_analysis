library(tidyverse)
library(lme4)
library(data.table)
library(readxl)
library(fmri.pipeline)
out_dir <- "~/OneDrive - University of Pittsburgh/Documents/SCEPTIC_fMRI/explore_medusa"
repo_directory <- "~/code/clock_analysis"
# repo_directory <- "~/Data_Analysis/clock_analysis"

reprocess = F
emm = T # extract EMMEANS estimates, e.g. for hi/lo abs(PE)
rerun_old_models = F
alignment = "rt"
bilateral = T # combine R/L, include as crossed effect with random slopes
# mixed_by call
# source("~/code/fmri.pipeline/R/mixed_by.R")

# helper function to compile list of formulae
named_list <- function(...) {
  vnames <- as.character(match.call())[-1]
  return(setNames(list(...), vnames))
}

# setwd(file.path(paste0(out_dir,"/dan_beta_1")))
setwd(file.path(paste0(out_dir,"/data")))
if (reprocess) {
  # move to new 400 labels
  labels <- setDT(read_excel(file.path(paste0(repo_directory, "/fmri/keuka_brain_behavior_analyses/dan/MNH DAN Labels 400 Good Only 47 parcels.xlsx")))) %>%
    rename(roi_num7=roi7_400) %>%
    select(roi_num7, parcel_group, hemi) %>% mutate(atlas_value = as.integer(roi_num7))
  # 
  # 
  # files <-  gsub("//", "/", list.files(pattern = "interpolated", full.names = T))
  # message(paste0("Found ", length(files), " files."))
  # csl <- lapply(files, function(x) {
  #   print(x)
  #   df <- fread(x) 
  #   df$id <- as.integer(stringi::stri_extract_first(x, regex = "\\d+"))
  #   # if (class(df)=="list") {
  #   df$run <- stringi::stri_extract_last(x, regex = "\\d+")
  #   # } else if (ncol(df)==3) {
  #   # df <- df$fit_df
  #   # }
  #   return(df)
  # })
  # rt_visuomotor_long <- data.table::rbindlist(csl)
  # setwd(out_dir)
  rt_visuomotor_long <- read_csv("rt_aligned_dan_beta_1.csv.gz", 
                                 col_types = cols(id = col_integer()))
  rt_visuomotor_long <- rt_visuomotor_long %>% mutate(run = as.integer(parse_number(run)))
  forty_seven <- unique(labels$atlas_value)
  rt_visuomotor_long_400_47 <- rt_visuomotor_long %>% filter(atlas_value %in% forty_seven) %>% inner_join(labels, by = "atlas_value")
  
  # clock_visuomotor_long <- read_csv("clock_aligned_444_dan.csv.gz", 
  #                                col_types = cols(id = col_integer()))
  # 
  # ids <- unique(clock_visuomotor_long$id)
  # 
  # # inspect clock-aligned more closely
  # test <- clock_visuomotor_long %>% filter(id == ids[2] & run == "run2") 
  pdf("check_rt_decon_alignment.pdf", height = 30, width = 60)
  ggplot(rt_visuomotor_long_400_47, aes(evt_time, decon_mean)) + geom_smooth() + facet_wrap(id~run)
  dev.off()
  # 
  # sdf <- clock_visuomotor_long %>% filter(evt_time < 6.6) %>% group_by(id, run, trial) %>% summarise(max_time_mean = evt_time[which.max(decon_mean)],
  #                                                                         max_time_median = evt_time[which.max(decon_mean)]) %>%
  #   group_by(id, run) %>% summarise(mean_peak = mean(max_time_mean))
  #   hist(sdf$mean_peak)
  # pdf("hist_max_clock_decon.pdf", height = 5, width = 8)
  # ggplot(sdf, aes(mean_peak, color = run)) + geom_histogram() 
  # dev.off()
  
  rsdf <- rt_visuomotor_long_400_47 %>% filter(evt_time < 5 & evt_time > -5) %>% group_by(id, run, trial) %>% summarise(max_time_mean = evt_time[which.max(decon_mean)],
                                                                                                                        max_time_median = evt_time[which.max(decon_mean)]) %>%
    group_by(id, run) %>% summarise(mean_peak = mean(max_time_mean))
  hist(rsdf$mean_peak, 30)
  pdf("hist_max_rt_decon.pdf", height = 5, width = 8)
  ggplot(rsdf, aes(mean_peak, color = as.character(run))) + geom_histogram() 
  dev.off()
  
  # save
  setwd(file.path(paste0(out_dir, "/data")))
  saveRDS(rt_visuomotor_long, file = "explore_rt_decon_all_444_parcels_beta1.rds")
  # saveRDS(rt_visuomotor_long_400_47, file = "explore_rt_decon_dan_400_47_beta_1.rds")
} else {
  # rt_visuomotor_long_400_47 <-  readRDS("explore_rt_decon_dan_400_47_beta1.rds")
  rt_visuomotor_long_400_47 <-  setDT(readRDS("explore_rt_decon_all_444_parcels_beta1.rds"))
  labels <- setDT(read_excel(file.path(paste0(repo_directory, "/fmri/keuka_brain_behavior_analyses/dan/MNH DAN Labels 400 Good Only 47 parcels.xlsx")))) %>%
    rename(roi_num7=roi7_400) %>%
    select(roi_num7, parcel_group, hemi) %>% mutate(atlas_value = as.integer(roi_num7))
  }

# get lags

rt_visuomotor_long_400_47 <- rt_visuomotor_long_400_47 %>% select(!any_of(c("decon_sd", "decon_median"))) %>% group_by(id, run, trial, atlas_value) %>% arrange(id, atlas_value, run, trial, evt_time) %>%
  mutate(decon_lag = lag(decon_mean)) %>% ungroup()


rt_visuomotor_long_400_47 <- rt_visuomotor_long_400_47 %>%  inner_join(labels, by = "atlas_value") #%>% filter(!is.na(Stream))

setwd(file.path(paste0(repo_directory, "/fmri/keuka_brain_behavior_analyses/dan/")))

source("get_trial_data.R")
source("medusa_final/plot_medusa.R") # careful -- plot_medusa sets out_dir, need to reset
# source("~/code/fmri.pipeline/R/mixed_by.R")

trial_df <- setDT(get_trial_data(dataset = "explore", repo_directory))  

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
sub_df <- readRDS("./explore_n146.rds") 
to_scale <- names(sub_df[8:57] %>% select_if(is.numeric))
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
    mutate(log_kld3 = log(kld3 + .00001),
           id = as.integer(id))
  trial_df <- inner_join(trial_df, sub_df, by = "id")
  d <- merge(trial_df, rt_visuomotor_long_400_47, by = c("id", "run", "trial"))
  d <- d %>% rename(side = "hemi")
  # d <- d %>% tidyr::separate(visuomotor_side, into=c("vm_gradient", "side"), sep="_")
}

## check data alignment:
# ggplot(d, aes(evt_time, decon_mean)) + geom_smooth()

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
  clock_base <- formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + v_max_wi +
                          v_entropy_wi + v_entropy_wi_change_lag + abs_pe_lag + outcome_lag +
                          (1 | id) )
  
  clock_rslope <- formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + v_max_wi +
                            v_entropy_wi + v_entropy_wi_change_lag + abs_pe_lag + outcome_lag +
                            (abs_pe_lag + v_entropy_wi + v_entropy_wi_change_lag + v_max_wi | id) )
  
  
  # clock-aligned using kld3_lag
  clock_kld <- formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + v_max_wi + 
                         v_entropy_wi + v_entropy_wi_change_lag + abs_pe_lag + outcome_lag + kld3_lag +
                         (1 | id) )
  
  # log version of kld given the strong skew
  clock_logkld <- formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + v_max_wi + 
                            v_entropy_wi + v_entropy_wi_change_lag + abs_pe_lag + outcome_lag + log_kld3_lag +
                            (1 | id) )
  
  # signed PE
  clock_pe <- formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + v_max_wi + 
                        v_entropy_wi + v_entropy_wi_change_lag + pe_max_lag +
                        (1 | id) )
  
  # abs_pe x outcome interaction
  clock_int <- formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + v_max_wi +
                         v_entropy_wi + v_entropy_wi_change_lag + abs_pe_lag * outcome_lag +
                         (1 | id) )
  
  # centered reward and abs pe to compute interaction
  # clock_int_cent <- formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + v_max_wi +
  #                            v_entropy_wi + v_entropy_wi_change_lag + abs_pe_c + rew_om_c + abspexrew +
  #                            (1 | id) )
  
  # trial interactions
  clock_trial <- formula(~ v_max_wi*run_trial + v_entropy_wi*run_trial + v_entropy_wi_change_lag*run_trial + abs_pe_lag*run_trial +
                           rt_csv_sc*run_trial + rt_lag_sc + outcome_lag +
                           (1 | id) )
  
  flist <- named_list(clock_base, clock_rslope, clock_kld, clock_logkld, clock_pe, clock_int, clock_trial)
} else if (alignment == "rt") {
  # main model
  rt_rslope_logkld3_bl <- formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + v_max_wi + v_entropy_wi + 
                                    v_entropy_wi_change +
                                    abs_pe + 
                                    outcome + log_kld3 +
                                    (v_entropy_wi_change + abs_pe | id) + (1 | side) )
  rt_N_C_age_ed_mmse <- formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + v_max_wi + v_entropy_wi + 
                                  v_entropy_wi_change*neo_neuroticism +
                                  v_entropy_wi_change*neo_conscientiousness +
                                  abs_pe*neo_neuroticism +
                                  abs_pe*neo_conscientiousness +
                                  outcome + log_kld3 +
                                  v_entropy_wi_change*education_yrs + 
                                  v_entropy_wi_change*age + 
                                  v_entropy_wi_change*mmse + 
                                  abs_pe*education_yrs + 
                                  abs_pe*age + 
                                  abs_pe*mmse + 
                                  (1| id) )
  
  rt_group_att_age_ed_mmse <- formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + v_max_wi + v_entropy_wi + 
                                  v_entropy_wi_change*Group_a +
                                  abs_pe*Group_a +
                                  outcome + log_kld3 +
                                  v_entropy_wi_change*education_yrs + 
                                  v_entropy_wi_change*age + 
                                  v_entropy_wi_change*mmse + 
                                  abs_pe*education_yrs + 
                                  abs_pe*age + 
                                  abs_pe*mmse + 
                                  (1| id) )
  rt_group_leth_age_ed_mmse <- formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + v_max_wi + v_entropy_wi + 
                                        v_entropy_wi_change*group_leth +
                                        abs_pe*group_leth +
                                        outcome + log_kld3 +
                                        v_entropy_wi_change*education_yrs + 
                                        v_entropy_wi_change*age + 
                                        v_entropy_wi_change*mmse + 
                                        abs_pe*education_yrs + 
                                        abs_pe*age + 
                                        abs_pe*mmse + 
                                        (1| id) )
  # all previously run models parked here for compactness
  if (rerun_old_models) {
    rt_upps_all_subsc_rslope <- formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + v_max_wi + v_entropy_wi + 
                                          v_entropy_wi_change*uppsp_lack_of_perseveration +
                                          v_entropy_wi_change*uppsp_lack_of_premeditation +
                                          v_entropy_wi_change*uppsp_negative_urgency +
                                          v_entropy_wi_change*uppsp_positive_urgency +
                                          abs_pe*uppsp_lack_of_perseveration +
                                          abs_pe*uppsp_lack_of_premeditation +
                                          abs_pe*uppsp_positive_urgency +
                                          abs_pe*uppsp_negative_urgency + 
                                          outcome +
                                          (abs_pe + v_entropy_wi_change| id) )
    
    
    # basal analysis
    rt_int_only <- formula(~   (1 | id) )
    
    rt_base <- formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + v_max_wi +
                         v_entropy_wi + v_entropy_wi_change + abs_pe + outcome +
                         (1 | id) )
    
    rt_upps_meds <- formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + v_max_wi + v_entropy_wi + 
                              v_entropy_wi_change*athf +
                              v_entropy_wi_change*opioid +
                              v_entropy_wi_change*sedhyp +
                              v_entropy_wi_change*antipsychotic +
                              v_entropy_wi_change*uppsp_negative_urgency + 
                              abs_pe*athf +
                              abs_pe*opioid +
                              abs_pe*sedhyp +
                              abs_pe*antipsychotic +
                              abs_pe*uppsp_negative_urgency + 
                              outcome +
                              (1 | id) )
    # explore interactions of abs_pe and v_entropy_wi_change with cognitive variables, psychopathology and group
    rt_base_age_edu_grp <- formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + v_max_wi + v_entropy_wi + 
                                     v_entropy_wi_change*education_yrs + 
                                     v_entropy_wi_change*age + 
                                     v_entropy_wi_change*Group + 
                                     abs_pe*education_yrs + 
                                     abs_pe*age + 
                                     abs_pe*Group + 
                                     outcome +
                                     (1 | id) )
    
    rt_base_age_edu_grp_mmse <- formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + v_max_wi + v_entropy_wi + 
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
    
    rt_base_age_edu_upps_anx <- formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + v_max_wi + v_entropy_wi + 
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
    rt_mmse_grp_upps_anx_rslope <- formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + v_max_wi + v_entropy_wi + 
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
    
    rt_mmse_grp_neg_urge_demo_rslope <- formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + v_max_wi + v_entropy_wi + 
                                                  v_entropy_wi_change*mmse + 
                                                  v_entropy_wi_change*uppsp_negative_urgency + 
                                                  abs_pe*mmse + 
                                                  abs_pe*uppsp_negative_urgency +
                                                  outcome +
                                                  (v_entropy_wi_change + abs_pe | id) )
    
    rt_mmse_grp_neg_urge_only_rslope <- formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + v_max_wi + v_entropy_wi + 
                                                  v_entropy_wi_change*uppsp_negative_urgency + 
                                                  abs_pe*uppsp_negative_urgency +
                                                  outcome +
                                                  (v_entropy_wi_change + abs_pe | id) )
    rt_mmse_grp_neg_urge_only_rint <- formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + v_max_wi + v_entropy_wi + 
                                                v_entropy_wi_change*uppsp_negative_urgency + 
                                                abs_pe*uppsp_negative_urgency +
                                                outcome +
                                                (1| id) )
    
    
    rt_upps_all_subsc <- formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + v_max_wi + v_entropy_wi + 
                                   v_entropy_wi_change*uppsp_lack_of_perseveration +
                                   v_entropy_wi_change*uppsp_lack_of_premeditation +
                                   v_entropy_wi_change*uppsp_negative_urgency +
                                   v_entropy_wi_change*uppsp_positive_urgency +
                                   abs_pe*uppsp_lack_of_perseveration +
                                   abs_pe*uppsp_lack_of_premeditation +
                                   abs_pe*uppsp_positive_urgency +
                                   abs_pe*uppsp_negative_urgency + 
                                   outcome +
                                   (1 | id) )
    # signed PE
    rt_pe <- formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + v_max_wi + 
                       v_entropy_wi + v_entropy_wi_change + pe_max +
                       (1 | id) )
    
    # signed PE
    rt_pe_kld <- formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + v_max_wi + 
                           v_entropy_wi + v_entropy_wi_change + pe_max + kld3 +
                           (1 | id) )
    # signed PE w/o covariates
    rt_pe_only <- formula(~ pe_max +
                            (1 | id) )
    
    rt_pe_rslope <- formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + v_max_wi + 
                              v_entropy_wi + v_entropy_wi_change + pe_max +
                              (pe_max | id) )
    
    # abs_pe x outcome interaction
    rt_int <- formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + v_max_wi +
                        v_entropy_wi + v_entropy_wi_change + abs_pe * outcome +
                        (1 | id) )
    
    rt_int_cent <- formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + v_max_wi +
                             v_entropy_wi + v_entropy_wi_change + abs_pe_c + rew_om_c + abspexrew +
                             (1 | id) )
    
    # centered outcome approach to get absPE at average of reward/omissions
    rt_int <- formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + v_max_wi +
                        v_entropy_wi + v_entropy_wi_change + abs_pe * outcome +
                        (1 | id) )
    
  }
  
  # flist <- named_list(rt_N_C_age_ed_mmse, rt_group_att_age_ed_mmse, rt_group_leth_age_ed_mmse)
  # flist <- named_list(rt_int_only)
  # flist <- named_list(rt_int_only)

  # flist <- named_list(rt_upps_all_subsc_rslope)
  # flist <- named_list(rt_rslope, rt_kld, rt_pe)
  # flist <- named_list(rt_base, rt_rslope, rt_kld, rt_logkld, rt_pe, rt_int, rt_int_cent, rt_trial)
  rt_base_lag <- formula(~ decon_lag + trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + v_max_wi +
                       v_entropy_wi + v_entropy_wi_change + abs_pe + outcome +
                       (1 | id) )
  
  flist <- named_list(rt_base_lag)
  d <- d %>% filter(!is.na(decon_lag))

}
splits <- c("parcel_group", "evt_time") # bilateral
# splits <- c("parcel_group", "side", "evt_time")
message(paste0("Running mixed_by for ", print(flist)))
ddf <- mixed_by(d, outcomes = "decon_mean", rhs_model_formulae = flist,
                split_on = splits, scale_predictors = c("abs_pe", "abs_pe_lag", "pe_max", "pe_max_lag", "run_trial", "decon_lag",
                                                        to_scale),
                tidy_args = list(effects = c("fixed", "ran_vals", "ran_pars", "ran_coefs"), conf.int = TRUE), 
                calculate = c("parameter_estimates_reml"), ncores = detectCores() - 1, refit_on_nonconvergence = 5, padjust_by = "term"#,
                # emtrends_spec = list(
                #   abspe = list(outcome = "decon_mean", model_name = "rt_int", var = "abs_pe", specs = c("outcome"))
                #)
)

saveRDS(ddf, file=file.path(out_dir, paste0(alignment, "_", flist[1], "May_30_2023.rds")))

# saveRDS(ddf, file=file.path(out_dir, paste0(alignment, "_", splits[1], "_att_leth_N_C_Nov_29_2022.rds")))

# saveRDS(ddf, file=file.path(out_dir, paste0(alignment, "_", splits[1], "_rt_int_Sept_2022.rds")))


# saveRDS(ddf, file=file.path(out_dir, paste0(alignment, "_", splits[1], "_upps_subscales_rslope_explore_400_47_encode_Aug4_2022.rds")))
# ddf <- readRDS()
# 
df <- ddf$coef_df_reml %>% dplyr::filter(evt_time > -3 & evt_time < 5) %>%
  filter(effect=="fixed") %>%
  group_by(term, model_name) %>%
  mutate(p_FDR=p.adjust(p.value, method="fdr"),
         parcel_group = factor(parcel_group, order = T, levels = c("MT+", "PPCcaudal", "PPCrostral", "Premotor"), labels = c("MT+", "Caudal PPC", "Rostral PPC", "Premotor")))  %>% 
  ungroup() %>% setDT()

plot_medusa(df, x="evt_time", y="estimate", ymin="estimate - std.error", ymax="estimate + std.error", color="parcel_group", #facet_by="side",
            out_dir=file.path(out_dir, "dan_beta_1"),  p.value="padj_BY_term")

# plot manually to make sure it's not in the plotting code
colors <- RColorBrewer::brewer.pal(4, "Dark2") %>% setNames(c("1" = "MT+","2" = "Premotor","3" = "Rostral PPC","4" = "Caudal PPC"))
ggplot(df %>% filter(term == "(Intercept)"), aes(evt_time, estimate, color=parcel_group)) +
  geom_line(size=1, position=position_dodge(width=0.4)) + 
  geom_point(aes(size=as.numeric(1/p_FDR)), position=position_dodge(width=0.4)) +
  scale_color_manual(values = colors) +
  geom_hline(yintercept = 0, size=1.5, alpha=0.6) +
  geom_vline(xintercept = 0, size=1.5, alpha=0.6) #+
scale_size_manual(values=c(0.5, 0.8, 1.1, 1.4)) + theme_bw(base_size=15)

# aggregate all results into one RDS
setwd(out_dir)
files <-  gsub("//", "/", list.files(pattern = ".*400.*47.*rds", full.names = T))
message(paste0("Found ", length(files), " files."))
csl <- lapply(files, function(x) {
  print(x)
  df <- readRDS(x) %>% Filter(Negate(is.null),.)
  df <- df$coef_df_reml
  # df$id <- as.integer(stringi::stri_extract_first(x, regex = "\\d+"))
  # # if (class(df)=="list") {
  # df$run <- stringi::stri_extract_last(x, regex = "\\d+")
  # } else if (ncol(df)==3) {
  # df <- df$fit_df
  # }
  return(df)
})
ddf_all <- data.table::rbindlist(csl)
saveRDS(ddf_all, file=file.path(out_dir, paste0(alignment, "_", splits[1], "_all_models_explore_400_47_encode_Aug4_2022.rds")))

# 
# dd <- ddf$emtrends_list$abspe
# dd <- dd[,c(-4, -5)]
# ggplot(dd, aes(x=evt_time, y=abs_pe.trend, ymin=abs_pe.trend-std.error, ymax=abs_pe.trend+std.error, color=vm_gradient)) + 
#   geom_pointrange() + geom_line() + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
#   facet_grid(outcome...1 ~ side)
