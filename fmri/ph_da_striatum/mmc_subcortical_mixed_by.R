# 2022-03-31 AndyP
# mixed_by model building for AIC comparison
# Last modified by Angela I Sep 2023 for Explore group comparisons
# Mixed by runs an MLM at each time point
# ddf is the output structure; reml and ml are two different convergence methods (use reml - contains AIC and BIC and log likelihood)
# ddf$coeff_df_reml - has all fixed and random effects; value is the estimate and estimated error is std.error

rm(list=ls())

library(stringr)
library(pracma)
library(wesanderson)
library(tidyverse)
library(data.table) 
library(fmri.pipeline)

wm = FALSE # run working memory analyses a la DAN paper
align_to = "RT" #options are clock or response
# AI - these flags allow you to collapse across hemispheres and also to be able to split into networks; default is to do by ROI number (atlas_value)
do_subregion = TRUE #this will collapse across subregions and hemispheres
do_structure = FALSE  #collapses across coarser groups and hemispheres
do_network = FALSE #this will collapse across networks and hemispheres
do_hippocampusAP = FALSE #run hippocampus anterior vs. posterior analysis
#models = c(5, 5.1, 5.2, 5.5, 5.6, 5.7, 7, 7.5, 7.6) #models to run, clock
#models = c(5, 5.1, 5.2, 5.3, 5.4, 5.5, 5.6, 5.7, 5.8, 5.9, 7, 7.5, 7.6) #models to run, feedback
models = c(4) # AD: running only the basic models
region_to_run = "subcortical" #options are cortical (vmPFC) and subcortical
do_vPFC_HC = FALSE
split_hipp_amyg = TRUE #if want to separate out hippocampus and amygdala

subcortical_cache_dir = '~/OneDrive - University of Pittsburgh/Documents/SCEPTIC_fMRI/mmclock_subcortical_medusa/data'  # where the decons are
subcortical_output_dir = '~/OneDrive - University of Pittsburgh/Documents/SCEPTIC_fMRI/mmclock_subcortical_medusa/'  # where the decons are
repo_directory = "/Users/alexdombrovski/code/clock_analysis"
map_dir = "~/code/schaefer_wb_parcellation/labels/"
source("/Users/alexdombrovski/code/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/get_trial_data.R")
ncores = 22 #for MacPro, total = 24
# source("~/code/fmri.pipeline/R/mixed_by.R")

setwd(subcortical_output_dir)

# add trial-level behavioral variables ----
trial_df <- get_trial_data(repo_directory) %>%
  select(-asc_trial) %>%
  mutate(across(c(id, run, trial, run_trial), as.integer)) %>%
  dplyr::select(id, run, run_trial, trial, trial_neg_inv_sc, rt_csv_sc, rt_lag_sc, pe_max, rew_om_c, abs_pe_c, abspexrew,
                v_entropy_wi, v_entropy_wi_change, kld3, v_max_wi, abs_pe, reward, reward_lag, v_entropy_wi_full, v_entropy_wi_change_full, iti_ideal) %>%
  mutate(log_kld3 = log(kld3 + .00001)) %>% 
  group_by(id, run) %>% arrange(run_trial) %>%
  mutate(v_max_lead_sc = lead(v_max_wi),
         iti_sc = scale(iti_ideal),
         rt_csv_lag_sc = lag(rt_csv_sc),
         reward_rec = if_else(reward=="omission", -0.5, 0.5),
         reward_lag_rec = if_else(reward_lag=="omission", -0.5, 0.5),
         abs_pe_max_sc = scale(abs(pe_max)),
         pe_max_lag = lag(pe_max),
         pe_max_lag_sc = scale(pe_max_lag),
         abs_pe_max_lag_sc = scale(abs(pe_max_lag)),
         rt_csv_lag_sc = lag(rt_csv_sc),
         iti_prev_sc = lag(iti_sc)) %>% ungroup


# read in decons
message("Loading medusa data from cache: ", subcortical_cache_dir)
if(align_to == "RT"){
  print("aligning to response")
  subcortical_df = fread(file.path(subcortical_cache_dir,  'Schaefer_444_Subcortical_2.3mm_rt_long_decon_aligned.csv.gz')) #response-aligned
  vta_df <- readRDS("data/VTA_medusa_rt_aligned.rds") %>% rename(subregion = atlas_value) %>% mutate(atlas_value = 445)
  shared_cols <- intersect(names(vta_df), names(subcortical_df))
  subcortical_df <- rbind(subcortical_df %>% select(all_of(shared_cols)), vta_df %>% select(all_of(shared_cols)))
  # cortical_df = fread(file.path(cortical_cache_dir,  'rt_aligned_200_vmpfc.csv.gz')) #response aligned
} else if(align_to == "clock"){
  print("aligning to clock")
  subcortical_df = fread(file.path(subcortical_cache_dir,  'Schaefer_444_Subcortical_2.3mm_clock_long_decon_aligned.csv.gz')) # clock-aligned
  # cortical_df = fread(file.path(cortical_cache_dir,  'clock_aligned_200_vmpfc.csv.gz')) #
}

if (wm) {wm_df <- readRDS(file.path(repo_directory, "fmri/keuka_brain_behavior_analyses/dan/wm_entropy/mmclock_wm_entropy.rds"))
trial_df <- trial_df %>% inner_join(wm_df, by = c("id", "run", "trial")) %>% group_by(id, run) %>% arrange(id, run, trial) %>% mutate(
  wm_entropy_wi = scale(choice_entropy + rewom_entropy),
  wm_choice_entropy_wi  = scale(choice_entropy),
  wm_reward_entropy_wi  = scale(reward_entropy),
  wm_entropy_wi_change = lead(wm_entropy_wi) - wm_entropy_wi,
  wm_choice_entropy_wi_change = lead(wm_choice_entropy_wi) - wm_choice_entropy_wi,
  wm_reward_entropy_wi_change = lead(wm_reward_entropy_wi) - wm_reward_entropy_wi
) %>% ungroup()
}
d <- inner_join(trial_df, subcortical_df, by = c("id", "run", "trial")) 


#Add vmpfc to the subcortical df
# vmpfc_roi_list <- c(55,56,65,66,67,84,86,88,89,159,160,161,170,171,191,192,194)
#vmpfc_roi_list <- append(vmpfc_roi_list,seq(201,230,1))
# cortical_df <- cortical_df %>% filter(atlas_value %in% vmpfc_roi_list)


#fix how run is coded in subcortical_df
#Get labels file; change to the adjusted one (region 191 changed from Default to Limbic for symmetry)??
#labels <- read_csv("~/Documents/RESEARCH/ANALYSES/fMRI/Dec2022/extracted_values/12Dec2022/region_labels_244_adj.csv") %>% mutate(roi_num = as.numeric(roi_num)) %>% inner_join(read_csv("~/Documents/RESEARCH/ANALYSES/fMRI/Dec2022/extracted_values/12Dec2022/region_lookup_244.csv"), by = "roi_num")
setwd(map_dir)
labels <- read_csv("Schaefer_444_region_labels.csv") %>% mutate(roi_num = as.numeric(roi_num)) %>% 
  inner_join(read_csv("Schaefer_444_region_lookup.csv"), by = "roi_num") %>% mutate(atlas_value = as.numeric(roi_num)) 
# label missing networks as amyg_hipp_thal (amygdala, hippocampus, thalamus)
if (split_hipp_amyg){
  labels$network[str_detect(labels$subregion, regex("Hippocampus", ignore_case=TRUE))] <- "hippocampus"
  labels$network[str_detect(labels$subregion, regex("BLA", ignore_case=TRUE))] <- "amygdala"
  labels$network[str_detect(labels$subregion, regex("CMN", ignore_case=TRUE))] <- "amygdala"
} else {labels$network[is.na(labels$network)] <- "amyg_hipp_thal"}
#Merge labels with subcortical medusa data
# cortical_df <- cortical_df %>% inner_join(labels, by="atlas_value") #add labels to subcortical_df 
#Code the groupings
labels <- labels %>% mutate(group = as.factor(case_when(
  str_detect(subregion, regex("Anterior Putamen", ignore_case=TRUE)) ~ "Striatum",
  str_detect(subregion, regex("Caudate Tail and Lateral Putamen", ignore_case=TRUE)) ~ "Striatum",
  str_detect(subregion, regex("Caudate Head", ignore_case=TRUE)) ~ "Striatum",
  str_detect(subregion, regex("Ventral Striatum", ignore_case=TRUE)) ~ "Striatum",
  str_detect(subregion, regex("Posterior Putamen", ignore_case=TRUE)) ~ "Striatum",
  str_detect(subregion, regex("BLA", ignore_case=TRUE)) ~ "Amygdala",
  str_detect(subregion, regex("CMN", ignore_case=TRUE)) ~ "Amygdala",
  str_detect(subregion, regex("Hippocampus", ignore_case=TRUE)) ~ "Hippocampus",
  str_detect(subregion, regex("Thalamus", ignore_case=TRUE)) ~ "Thalamus"))) %>% 
  mutate(hippocampus_subregion = as.factor(case_when( #Code hippocampus subregions - anterior vs. posterior
    str_detect(subregion, regex("Anterior Putamen", ignore_case=TRUE)) ~ "Striatum",
    str_detect(subregion, regex("Caudate Tail and Lateral Putamen", ignore_case=TRUE)) ~ "Striatum",
    str_detect(subregion, regex("Caudate Head", ignore_case=TRUE)) ~ "Striatum",
    str_detect(subregion, regex("Ventral Striatum", ignore_case=TRUE)) ~ "Striatum",
    str_detect(subregion, regex("Posterior Putamen", ignore_case=TRUE)) ~ "Striatum",
    str_detect(subregion, regex("BLA", ignore_case=TRUE)) ~ "Amygdala",
    str_detect(subregion, regex("CMN", ignore_case=TRUE)) ~ "Amygdala",
    str_detect(subregion, regex("Thalamus", ignore_case=TRUE)) ~ "Thalamus",
    str_detect(subregion, regex("anterior", ignore_case=TRUE)) ~ "AH",
    str_detect(subregion, regex("posterior", ignore_case=TRUE)) ~ "PH")))
# add VTA label
labels[nrow(labels) + 1,] = rep(NA, 14)
labels$subregion[nrow(labels)] = "VTA"
labels$atlas_value[nrow(labels)] = 445
fwrite(labels, file = "Schaefer_445_region_labels.csv") # including VTA as 445
d <- d %>% inner_join(labels, by="atlas_value") #add labels to subcortical_df 

# cortical_df <- cortical_df %>% mutate(group = as.factor(case_when(
#   str_detect(subregion, regex("OFC1", ignore_case=TRUE)) ~ "mPFC",
#   str_detect(subregion, regex("OFC2", ignore_case=TRUE)) ~ "mPFC",
#   str_detect(subregion, regex("OFC3", ignore_case=TRUE)) ~ "mPFC",
#   str_detect(subregion, regex("PFC2", ignore_case=TRUE)) ~ "mPFC",
#   str_detect(subregion, regex("PFC4", ignore_case=TRUE)) ~ "mPFC",
#   str_detect(subregion, regex("PFC6", ignore_case=TRUE)) ~ "mPFC",
#   str_detect(subregion, regex("PFC7", ignore_case=TRUE)) ~ "mPFC",
#   str_detect(subregion, regex("PFCdPFCm1", ignore_case=TRUE)) ~ "mPFC",
#   str_detect(subregion, regex("PFCdPFCm2", ignore_case=TRUE)) ~ "mPFC",
#   str_detect(subregion, regex("PFCdPFCm4", ignore_case=TRUE)) ~ "mPFC",
#   str_detect(subregion, regex("PFCl1", ignore_case=TRUE)) ~ "mPFC",
#   str_detect(subregion, regex("PFCl2", ignore_case=TRUE)) ~ "mPFC")))



#Censor out data that bleeds into adjacent trials (added BY ANGELA 4/21/2023)
# AD: THIS SHOULD NOT BE NECESSARY IN MMCLOCK, DOUBLE-CHECK WITH MICHAEL
# if(align_to == "RT"){
#   Q$subcortical_decon[Q$evt_time > Q$iti_ideal] = NA #censors next trial
#   Q$subcortical_decon[Q$evt_time < -(Q$rt_csv + Q$iti_prev)] = NA # censors prior trial
# } else if(align_to == "clock"){
#   Q$subcortical_decon[Q$evt_time > Q$rt_csv + Q$iti_ideal] = NA # censors next trial
#   Q$subcortical_decon[Q$evt_time < -Q$iti_prev] = NA # censors prior trial (edited) 
# }

#To get basic demographics based on lethality group - before scaling age
#Q %>% filter(group_leth_c=="LL_Attempters") %>% select(id,sex,age) %>% unique() %>% summarise(group_n=n(), males_n=sum(sex=="M"),males_percent=100*(males_n/group_n))
#Q %>% filter(group_leth_c=="HL_Attempters") %>% select(id,sex,age) %>% unique() %>% summarise(group_n=n(), males_n=sum(sex=="M"),males_percent=100*(males_n/group_n))


# scale demographics stuff
#demo <- readRDS('/Users/angela/Documents/RESEARCH/ANALYSES/fMRI/Dec2022/extracted_values/12Dec2022/explore_n146.rds') %>% mutate(id = registration_redcapid)
#demos$id <- as.factor(demos$id)
#Q <- inner_join(Q,demos,by=c('id')) # takes a while

#Loop through all of the mdoels
for (model_to_run in models) {
  print(model_to_run)
  # left hand side of the GLM is outcome, which is specified in the mixed_by variable
  #9/22/2023 - just added extra controls to 5, 5.4, 5.5, 5.9 for complete sensitivity analysis
  #Add variant with all of the control stuff ()
  if(align_to == "RT"){ #need to align to v_max_wi_lead if response aligned
    if (model_to_run == 1) {decode_formula <- formula(~ v_max_lead_sc + rt_csv_sc + iti_sc + (1|id))} #model 1 - basic model, include RT and ISI to de-noise
    else if (model_to_run == 2) {decode_formula <- formula(~ v_max_lead_sc + rt_csv_sc + rt_csv_lag_sc + iti_sc + reward_rec + reward_lag_rec +  abs_pe_max_sc + (1|id))}
    else if (model_to_run == 3) {decode_formula <- formula(~ v_max_wi + rt_csv_sc + rt_csv_lag_sc + iti_sc + reward_rec + reward_lag_rec +  abs_pe_max_sc + (1|id))} 
    else if (model_to_run == 4) {decode_formula <- formula(~ v_max_wi + rt_csv_sc + rt_csv_lag_sc + iti_sc + reward_rec + reward_lag_rec +  abs_pe_max_sc + v_entropy_wi + v_entropy_wi_change + (1|id))} 
    else if (model_to_run == 5) {decode_formula <- formula(~ v_max_wi + rt_csv_sc + rt_csv_lag_sc + iti_sc + reward_rec + reward_lag_rec +  abs_pe_max_sc + v_entropy_wi + v_entropy_wi_change + 
                                                             (v_max_wi + rt_csv_sc + rt_csv_lag_sc + iti_sc + reward_rec + reward_lag_rec +  abs_pe_max_sc + v_entropy_wi + v_entropy_wi_change|id))} 
    
    #else if (model_to_run == 2) {decode_formula <- formula(~ rt_csv_sc + rt_csv_lag_sc + iti_sc + reward_rec + reward_lag_rec +  abs_pe_max_sc + (v_max_lead_sc + 1|id))}
    #else if (model_to_run == 3) {decode_formula <- formula(~ rt_csv_sc + rt_csv_lag_sc + iti_sc + reward_rec + reward_lag_rec +  abs_pe_max_sc + (v_max_wi + 1|id))} 
  } else if(align_to == "clock"){
    #   if (model_to_run == 1) {decode_formula <- formula(~ v_max_wi + rt_csv_sc + iti_sc + (1|id))} #model 1 - basic model, include RT and ISI to de-noise
    #   else if (model_to_run == 1.5) {decode_formula <- formula(~ v_max_wi + rt_csv_lag_sc + iti_prev_sc + (1|id))} #- use lagged rtcsv and iti variables
    #   #else if (model_to_run == 2) {decode_formula <- formula(~ v_max_wi + rt_csv_sc + rt_csv_lag_sc + iti_sc + reward_rec + reward_lag_rec +  abs_pe_max_sc + (1|id))}
    #   else if (model_to_run == 2) {decode_formula <- formula(~ v_max_wi + rt_csv_lag_sc + iti_prev_sc + reward_lag_rec +  abs_pe_max_lag_sc + (1|id))}
    #   else if (model_to_run == 3) {stop('REDUNDANT - THIS IS THE SAME AS MODEL 2 for CLOCK-ALIGNED DATA')}
  }
  #if (do_symmetry){ #does L and R of each region together, assuming they don't differ
  setwd(subcortical_output_dir)
  if (do_structure){
    splits = c('evt_time','group') #will do a separate MLM for each group in splits; evt_time is the TR; collapse across both hemispheres of subregions
    if(region_to_run == "subcortical") {dir <- paste0(subcortical_cache_dir,'/Structure_CombinedHemis')
    dir.create(file.path(subcortical_cache_dir, '/Structure_CombinedHemis'), showWarnings = FALSE)
    } else if(region_to_run == "cortical") {dir <- paste0(cortical_cache_dir,'/Structure_CombinedHemis')
    dir.create(file.path(cortical_cache_dir, '/Structure_CombinedHemis'), showWarnings = FALSE)}
    setwd(dir)
    df0 <- decode_formula   #df0 <- decode_formula[[i]]
    print(df0)
    ddf <- mixed_by(d, outcomes = "decon_mean", rhs_model_formulae = df0 , split_on = splits,
                    padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3, #ncores put 1 less than total cores on your machine or else it will slow down
                    tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE), #specify what output dataframe will have; will be looking mostly at fixed terms; lmer summary, will need all for that
                    calculate = "parameter_estimates_reml",
                    #Run post-hoc tests
                    emmeans_spec = list(
                      if(align_to == "RT"){ #need to align to v_max_lead_sc if response aligned
                        if (model_to_run == 3 | model_to_run == 5.3 | model_to_run == 5.8) {V = list(outcome="decon_mean", model_name=model_to_run,specs=c("v_max_wi"), at=list(v_max_wi=c(-1.5,1.5)))}
                        else {V = list(outcome="decon_mean", model_name=model_to_run,specs=c("v_max_lead_sc"), at=list(v_max_lead_sc=c(-1.5,1.5)))}
                      } else if(align_to == "clock"){
                        V = list(outcome="decon_mean", model_name=model_to_run,
                                 specs=c("v_max_wi"), at=list(v_max_wi=c(-1.5,1.5)))
                      }
                    )
    )
  } else if (do_subregion){
    splits = c('evt_time','subregion') #will do a separate MLM for each group in splits; evt_time is the TR; subregions - collapsed across hemispheres
    dir <- paste0(subcortical_output_dir,'/Subregion_CombinedHemis')
    dir.create(file.path(subcortical_output_dir, '/Subregion_CombinedHemis'), showWarnings = FALSE)
    setwd(dir)
    df0 <- decode_formula   #df0 <- decode_formula[[i]]
    print(df0)
    ddf <- mixed_by(d, outcomes = "decon_mean", rhs_model_formulae = df0 , split_on = splits,
                    padjust_by = c("term", "subregion"), padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3, #ncores put 1 less than total cores on your machine or else it will slow down
                    tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE), #specify what output dataframe will have; will be looking mostly at fixed terms; lmer summary, will need all for that
                    #Run post-hoc tests
                    emmeans_spec = list(
                      if(align_to == "RT"){ #need to align to v_max_lead_sc if response aligned
                        if (model_to_run == 3 | model_to_run == 5.3 | model_to_run == 5.8) {V = list(outcome="decon_mean", model_name=model_to_run,specs=c("v_max_lead_sc"), at=list(v_max_lead_sc=c(-1.5,1.5)))}
                        else {V = list(outcome="decon_mean", model_name=model_to_run,specs=c("v_max_lead_sc"), at=list(v_max_lead_sc=c(-1.5,1.5)))}
                      } else if(align_to == "clock"){
                        V = list(outcome="decon_mean", model_name=model_to_run,specs=c("v_max_wi"), at=list(v_max_wi=c(-1.5,1.5)))}
                    )
    )
  } else if (do_network){
    splits = c('evt_time','network') #will do a separate MLM for each group in splits; evt_time is the TR; collapse networks; but several are NA (change to amygdala_thalamus_hippocampus)
    if(region_to_run == "subcortical") {dir <- paste0(subcortical_cache_dir,'/Network_CombinedHemis')
    dir.create(file.path(subcortical_output_dir, '/Network_CombinedHemis'), showWarnings = FALSE)
    #AI added 5/16 to narrow down networks for subcortical models to just DMN, CONT, LIM, and Hipp/Amyg
    if (split_hipp_amyg){d <- d %>% filter(network %in% c("amygdala","hippocampus", "Cont","Default","Limbic"))
    } else {d <- d %>% filter(network %in% c("amyg_hipp_thal","Cont","Default","Limbic")) %>% filter(!group %in% c("Thalamus"))}
    } else if(region_to_run == "cortical") {dir <- paste0(cortical_cache_dir,'/Network_CombinedHemis')
    dir.create(file.path(subcortical_output_dir, '/Network_CombinedHemis'), showWarnings = FALSE)}
    setwd(dir)
    df0 <- decode_formula   #df0 <- decode_formula[[i]]
    print(df0)
    ddf <- mixed_by(Q, outcomes = "decon_mean", rhs_model_formulae = df0 , split_on = splits,
                    padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3, #ncores put 1 less than total cores on your machine or else it will slow down
                    tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE), #specify what output dataframe will have; will be looking mostly at fixed terms; lmer summary, will need all for that
                    #Run post-hoc tests
                    emmeans_spec = list(
                      if(align_to == "RT"){ #need to align to v_max_lead_sc if response aligned
                        if (model_to_run == 3 | model_to_run == 5.3 | model_to_run == 5.8) {V = list(outcome="decon_mean", model_name=model_to_run,specs=c("v_max_lead_sc"), at=list(v_max_lead_sc=c(-1.5,1.5)))}
                        else {V = list(outcome="decon_mean", model_name=model_to_run,specs=c("v_max_lead_sc"), at=list(v_max_lead_sc=c(-1.5,1.5)))}
                      } else if(align_to == "clock"){
                        V = list(outcome="decon_mean", model_name=model_to_run,
                                 specs=c("v_max_wi"), at=list(v_max_wi=c(-1.5,1.5)))
                      }
                    )
    )
  } else if (do_hippocampusAP) {
    splits = c('evt_time','hippocampus_subregion') #will do a separate MLM for each group in splits; evt_time is the TR; subregions - collapsed across hemispheres
    dir <- paste0(subcortical_output_dir,'/Structure_splitHippocampus_CombinedHemis')
    dir.create(file.path(subcortical_output_dir, '/Structure_splitHippocampus_CombinedHemis'), showWarnings = FALSE)
    setwd(dir)
    df0 <- decode_formula 
    print(df0)
    ddf <- mixed_by(d, outcomes = "decon_mean", rhs_model_formulae = df0 , split_on = splits,
                    padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3, #ncores put 1 less than total cores on your machine or else it will slow down
                    tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE), #specify what output dataframe will have; will be looking mostly at fixed terms; lmer summary, will need all for that
                    #Run post-hoc tests
                    emmeans_spec = list(
                      if(align_to == "RT"){ #need to align to v_max_lead_sc if response aligned
                        if (model_to_run == 3 | model_to_run == 5.3 | model_to_run == 5.8) {V = list(outcome="decon_mean", model_name=model_to_run,specs=c("v_max_wi"), at=list(v_max_lead_sc=c(-1.5,1.5)))}
                        else {V = list(outcome="decon_mean", model_name=model_to_run,specs=c("v_max_lead_sc"), at=list(v_max_lead_sc=c(-1.5,1.5)))}
                      } else if(align_to == "clock"){
                        V = list(outcome="decon_mean", model_name=model_to_run,specs=c("v_max_wi"), at=list(v_max_wi=c(-1.5,1.5)))}
                    )
    )
  } else {
    splits = c('evt_time','atlas_value') #will do a separate MLM for each group in splits; evt_time is the TR
    dir <- paste0(subcortical_cache_dir,'/AtlasValue_SeparateHemis')
    dir.create(file.path(subcortical_output_dir, '/AtlasValue_SeparateHemis'), showWarnings = FALSE)
    setwd(dir)
    df0 <- decode_formula   #df0 <- decode_formula[[i]]
    print(df0)
    ddf <- mixed_by(d, outcomes = "decon_mean", rhs_model_formulae = df0 , split_on = splits,
                    padjust_by = "term", padjust_method = "fdr", ncores = ncores, refit_on_nonconvergence = 3, #ncores put 1 less than total cores on your machine or else it will slow down
                    tidy_args = list(effects=c("fixed","ran_vals","ran_pars","ran_coefs"),conf.int=TRUE), #specify what output dataframe will have; will be looking mostly at fixed terms; lmer summary, will need all for that
                    #Run post-hoc tests
                    emmeans_spec = list(
                      if(align_to == "RT"){ #need to align to v_max_lead_sc if response aligned
                        if (model_to_run == 3 | model_to_run == 5.3 | model_to_run == 5.8) {V = list(outcome="decon_mean", model_name=model_to_run,specs=c("v_max_wi"), at=list(v_max_lead_sc=c(-1.5,1.5)))}
                        else {V = list(outcome="decon_mean", model_name=model_to_run,specs=c("v_max_lead_sc"), at=list(v_max_lead_sc=c(-1.5,1.5)))}
                      } else if(align_to == "clock"){
                        V = list(outcome="decon_mean", model_name=model_to_run,
                                 specs=c("v_max_wi"), at=list(v_max_wi=c(-1.5,1.5)))
                      }
                    )
    )
  }
  curr_date <- strftime(Sys.time(),format='%Y-%m-%d')
  #need to code in way to change the name of this file based on the model
  if(align_to == "RT"){
    if(region_to_run == "subcortical") {save(ddf,file=paste0(curr_date,'-subcortical-response-','model_',model_to_run,'.Rdata')) # removed i (i is for if you are running several models)
    } else if (region_to_run == "cortical") {save(ddf,file=paste0(curr_date,'-vmpfc-response-','model_',model_to_run,'.Rdata'))} # removed i (i is for if you are running several models)
  } else if(align_to == "clock"){
    if(region_to_run == "subcortical") {save(ddf,file=paste0(curr_date,'-subcortical-clock-','model_',model_to_run,'.Rdata')) # removed i (i is for if you are running several models)
    } else if (region_to_run == "cortical") {save(ddf,file=paste0(curr_date,'-vmpfc-clock-','model_',model_to_run,'.Rdata'))} # removed i (i is for if you are running several models)
  }
}

#Now run plot_mixed_by_AI.R