#the goal of this script is to run an entire fmri analysis for SCEPTIC data, including levels 1-3 in FSL

library(dependlab)
library(foreach)
library(parallel)
library(doParallel)
library(plyr)
library(tidyverse)

scripts_dir <- "/Volumes/bek/explore/clock_analysis/fmri/fsl_pipeline/"
box_dir<-"~/Box"
setwd(scripts_dir)

source(file.path(scripts_dir, "functions", "push_pipeline.R"))
source(file.path(scripts_dir, "functions", "finalize_pipeline_configuration.R"))
#Jun2017: further ICAs on the MMClock data suggest a short steady-state problem. Drop 2 volumes for good measure.

###
# SCEPTIC MMClock Y3
getMainDir<-function(){
  return("/Volumes/bek/explore")
}

#trial_df <- read.csv("/gpfs/group/mnh5174/default/temporal_instrumental_agent/clock_task/vba_fmri/vba_out/compiled_outputs/mmclock_fmri_decay_factorize_mfx_trial_statistics.csv.gz")
#trial_df <- read.csv("/gpfs/group/mnh5174/default/temporal_instrumental_agent/clock_task/vba_fmri/vba_out/compiled_outputs/mmclock_fmri_decay_mfx_trial_statistics.csv.gz")
subjectsdata_dir<-file.path(box_dir,'skinner','data','eprime','clock_reversal')
if(file.exists(file.path(subjectsdata_dir,"vba_out.rdata"))){
  load(file.path(subjectsdata_dir,"vba_out.rdata"))
} else {
  source('~/Documents/UPMC/RStation/temporal_instrumental_agent/clock_task/vba_fmri/parse_sceptic_outputs.R')
  vba_output<-parse_sceptic_outputs(outdir = "~/Box/skinner/projects_analyses/SCEPTIC/fMRI_paper/vba_output/",subjects_dir = subjectsdata_dir)
  save(vba_output,file = file.path(subjectsdata_dir,"vba_out.rdata"))
}
vba_output<-as_tibble(vba_output)
vba_output$trial<-as.double(vba_output$trial)
#factorized, selective maintenance, equal basis-generalization width
#trial_df <- read.csv("/gpfs/group/mnh5174/default/temporal_instrumental_agent/clock_task/vba_fmri/vba_out/compiled_outputs/mmclock_fmri_decay_factorize_selective_psequate_mfx_trial_statistics.csv.gz") %>%
trial_df <- vba_output %>%
  mutate(
    trial_rel=case_when(
    trial >= 1 & trial <= 120 ~ trial,
    trial >= 120 & trial <= 240 ~ trial - 120L, #dplyr/rlang has gotten awfully picky about data types!!
    TRUE ~ NA_real_),
    v_entropy_no5=if_else(trial_rel <= 5, NA_real_, v_entropy),
    d_auc_sqrt=if_else(d_auc > 0, NA_real_, sqrt(-1*d_auc)), #only compute the sqrt of d_auc for negative (i.e., reasonable) observations
    v_entropy_sqrt=sqrt(v_entropy),
    rew_om=if_else(score_vba > 0, 1, 0)
  ) %>% #for win/loss maps
  group_by(id, run) %>%
  dplyr::mutate(   #compute rt_swing within run and subject
    rt_vmax_lag = dplyr::lag(rt_vmax, 1, order_by=trial),
    rt_vmax_change = abs(rt_vmax - rt_vmax_lag),
    v_entropy_lag = dplyr::lag(v_entropy, 1, order_by=trial),
    v_entropy_change = v_entropy - v_entropy_lag, #change in entropy
    v_entropy_change_pos = v_entropy_change*(v_entropy_change > 0),
    v_entropy_change_neg = abs(v_entropy_change*(v_entropy_change < 0)),
    rt_swing = abs( c(NA, diff(rt_csv)))/1000,
    rt_swing_sqrt=sqrt(rt_swing)) %>%
  ungroup()

source("/Volumes/bek/explore/scripts/startup.R")
startup()
masterdemo<-bsrc.conredcap2(ptcs$masterdemo,online = T,output = T,batch_size = 1000)
explore_demo<-masterdemo$data[which(masterdemo$data$registration_ptcstat___explore=="1"),]
subject_df<-data.frame(redcapid=explore_demo$registration_redcapid,
           age=round(as.numeric((as.Date(explore_demo$reg_condate_explore) - as.Date(explore_demo$registration_dob))/365.25),2),
           female=explore_demo$registration_gender=="F",
           scandate=as.Date(explore_demo$reg_condate_explore),
           GroupATT=explore_demo$registration_group,
           GroupLH=NA)



subject_df<-subject_df %>%
  rename(id=redcapid, Age=age, Female=female, ScanDate=scandate) %>%
  mutate(mr_dir = paste0("/Volumes/bek/explore/MR_Proc/", id), #convert to Date, then reformat YYYYMMDD
    I_Age = -1000*1/Age,
    I_Age_c = I_Age - mean(I_Age, na.rm=TRUE),
    Age_c = Age - mean(Age, na.rm=TRUE),
    Q_Age = Age_c^2,
    Q_Age_c = Q_Age - mean(Q_Age, na.rm=TRUE)
  )

#from 2017:
##results from Mean SCEPTIC regressor correlation.pdf indicate that regressors for vchosen, ventropy_decay_matlab, dauc, and pemax are
##reasonably uncorrelated. The worst is dauc with vchosen (mean r = -0.31), which makes sense that as learning progresses, chosen values
##are higher and there is less residue to decay. These 4 regressors are also of greatest theoretical interest

#Setup the global configuration for the full FSL pipeline
fsl_model_arguments <- list(
  paraname="clockrev",
  #analysis_name="MMClock_aroma_preconvolve_fse",
  pipeline_home=scripts_dir,
  analysis_name="explore_clock",
  trial_statistics = trial_df,
  subject_covariates = subject_df,
  id_col = "id",
  fmri_dir = "/Volumes/bek/explore/MR_Proc",
  expectdir = "clockRev_proc", #subfolder name for processed data
  expectfile = "nfswudktm_clockrev[0-9]_7.nii.gz", #expected file name for processed clock data
  usepreconvolve=TRUE,
  ncpus=6,
  drop_volumes=0, #to handle steady state concerns
  tr=0.6, #seconds
  spikeregressors=FALSE, #don't include spike regressors in nuisance variables since we are using AROMA
  sceptic_run_variants=list(
#    c("clock", "feedback_bs")
#    c("clock_bs", "feedback")
#    c("clock", "feedback", "v_chosen", "v_entropy", "d_auc", "pe_max"), #all signals with entropy of weights
#    c("clock", "feedback", "v_chosen", "v_entropy_func", "d_auc", "pe_max"), #all signals with entropy of evaluated function
    c("clock", "feedback", "v_chosen"), #individual regressors
    c("clock", "feedback", "v_entropy"), #clock-aligned
#    c("clock", "feedback", "v_entropy_feedback"), #feedback-aligned
#    c("clock", "feedback", "v_entropy_func"),
#    c("clock", "feedback", "d_auc"), #feedback-aligned
#    c("clock", "feedback", "d_auc_clock"), #clock-aligned
    c("clock", "feedback", "pe_max"),
#    c("clock", "feedback", "v_entropy_no5"),
#    c("clock", "feedback", "v_auc"),
#    c("clock", "feedback", "d_auc_sqrt"),
#    c("clock", "feedback", "rt_swing"),
#    c("clock", "feedback", "rt_swing_sqrt"),
#    c("clock", "feedback", "v_max"),
#    c("clock", "feedback", "mean_kld"),
#    c("clock", "feedback", "intrinsic_discrepancy"),
#    c("clock", "feedback", "mean_kld_feedback"),
#    c("clock", "feedback", "intrinsic_discrepancy_feedback"),
#    c("clock", "feedback", "rt_vmax_change"),
    c("clock", "feedback", "v_entropy_change")
#    c("clock", "feedback", "v_entropy_change_pos"),
#    c("clock", "feedback", "v_entropy_change_neg")
#    c("clock", "feedback", "rew_om"),
#    c("clock", "feedback", "pe_max", "rew_om")
  ),
  group_model_variants=list(
    c("Intercept"),
    # c("Intercept", "Age"),
    c("Intercept", "Age", "Group")
#    c("Intercept", "I_Age"),
#    c("Intercept", "I_Age", "Female")
  ),
  execute_feat=FALSE, #passed through to fsl_sceptic_model to create fsf, but not run the model
  #model_suffix="_fse", #factorized, selective, equal generalization width
  model_suffix="_fse_groupfixed", #factorized, selective, equal generalization width
  root_workdir="/Volumes/bek/explore/clock_rev/",
  n_cluster_beta_cpus=8, #should be number of l2 contrasts, or lower
  badids = c("") #exclude ppl here
)

#validate and populate any other pipeline details before execution
fsl_model_arguments <- finalize_pipeline_configuration(fsl_model_arguments)

save(fsl_model_arguments, file=paste0("configuration_files/", fsl_model_arguments$analysis_name, ".RData"))

#this pushes the full analysis pipeline in parallel, where parallelism is across sceptic_run_variants
push_pipeline(fsl_model_arguments, ncpus=fsl_model_arguments$pipeline_cpus)


## Other pipelines go here
