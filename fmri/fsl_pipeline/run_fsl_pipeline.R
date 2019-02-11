#the goal of this script is to run an entire fmri analysis for SCEPTIC data, including levels 1-3 in FSL

library(dependlab)
library(foreach)
library(parallel)
library(doParallel)
library(plyr)
library(tidyverse)

scripts_dir <- "/gpfs/group/mnh5174/default/clock_analysis/fmri/fsl_pipeline"
setwd(scripts_dir)

source(file.path(scripts_dir, "functions", "push_pipeline.R"))
source(file.path(scripts_dir, "functions", "finalize_pipeline_configuration.R"))

#Jun2017: further ICAs on the MMClock data suggest a short steady-state problem. Drop 2 volumes for good measure.

###
# SCEPTIC MMClock Y3

#trial_df <- read.csv("/gpfs/group/mnh5174/default/temporal_instrumental_agent/clock_task/vba_fmri/vba_out/compiled_outputs/mmclock_fmri_decay_factorize_mfx_trial_statistics.csv.gz")
#trial_df <- read.csv("/gpfs/group/mnh5174/default/temporal_instrumental_agent/clock_task/vba_fmri/vba_out/compiled_outputs/mmclock_fmri_decay_mfx_trial_statistics.csv.gz")


#factorized, selective maintenance, equal basis-generalization width
#trial_df <- read.csv("/gpfs/group/mnh5174/default/temporal_instrumental_agent/clock_task/vba_fmri/vba_out/compiled_outputs/mmclock_fmri_decay_factorize_selective_psequate_mfx_trial_statistics.csv.gz") %>%
trial_df <- read.csv("/gpfs/group/mnh5174/default/temporal_instrumental_agent/clock_task/vba_fmri/vba_out/compiled_outputs/mmclock_fmri_decay_factorize_selective_psequate_fixedparams_ffx_trial_statistics.csv.gz") %>%
  mutate(trial_rel=case_when(
    trial >= 1 & trial <= 50 ~ trial,
    trial >= 51 & trial <= 100 ~ trial - 50L, #dplyr/rlang has gotten awfully picky about data types!!
    trial >= 101 & trial <= 150 ~ trial - 100L,
    trial >= 151 & trial <= 200 ~ trial - 150L,
    trial >= 201 & trial <= 250 ~ trial - 200L,
    trial >= 251 & trial <= 300 ~ trial - 250L,
    trial >= 301 & trial <= 350 ~ trial - 300L,
    trial >= 351 & trial <= 400 ~ trial - 350L,
    TRUE ~ NA_integer_
  ), v_entropy_no5=if_else(trial_rel <= 5, NA_real_, v_entropy),
  d_auc_sqrt=if_else(d_auc > 0, NA_real_, sqrt(-1*d_auc)), #only compute the sqrt of d_auc for negative (i.e., reasonable) observations
  v_entropy_sqrt=sqrt(v_entropy),
  rew_om=if_else(score_vba > 0, 1, 0)) %>% #for win/loss maps
  group_by(id, run) %>%  dplyr::mutate(rt_swing = abs( c(NA, diff(rt_csv)))/1000, rt_swing_sqrt=sqrt(rt_swing)) %>% ungroup() #compute rt_swing within run and subject

subject_df <- read.table("/gpfs/group/mnh5174/default/clock_analysis/fmri/data/mmy3_demographics.tsv", header=TRUE) %>%
  rename(ID=lunaid, Age=age, Female=female, ScanDate=scandate) %>%
  mutate(mr_dir = paste0("/gpfs/group/mnh5174/default/MMClock/MR_Proc/", ID, "_", format((as.Date(ScanDate, format="%Y-%m-%d")), "%Y%m%d")), #convert to Date, then reformat YYYYMMDD
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
  #analysis_name="MMClock_aroma_preconvolve_fse",
  analysis_name="MMClock_aroma_preconvolve_fse_groupfixed",
  trial_statistics = trial_df,
  subject_covariates = subject_df,
  fmri_dir = "/gpfs/group/mnh5174/default/MMClock/MR_Proc",
  expectdir = "mni_5mm_aroma", #subfolder name for processed data
  expectfile = "nfaswuktm_clock[0-9]_5.nii.gz", #expected file name for processed clock data
  usepreconvolve=TRUE,
  ncpus=20,
  drop_volumes=2, #to handle steady state concerns
  tr=1.0, #seconds
  spikeregressors=FALSE, #don't include spike regressors in nuisance variables since we are using AROMA
  idexpr=expression(subid), #how to match between the subject ID in trial_statistics and the folder structure on the file system
  idregex="([0-9]{5})_\\d+", #5 digit ID, followed by irrelevant date. A bit inelegant, but used in setup_feat_lvl2_inputs to infer the subject's ID from the MR folder structure
  sceptic_run_variants=list(
#    c("clock", "feedback", "v_chosen", "v_entropy", "d_auc", "pe_max"), #all signals with entropy of weights
#    c("clock", "feedback", "v_chosen", "v_entropy_func", "d_auc", "pe_max"), #all signals with entropy of evaluated function
    c("clock", "feedback", "v_chosen"), #individual regressors
    c("clock", "feedback", "v_entropy"), #clock-aligned
#    c("clock", "feedback", "v_entropy_feedback"), #feedback-aligned
#    c("clock", "feedback", "v_entropy_func"),
    c("clock", "feedback", "d_auc"), #feedback-aligned
#    c("clock", "feedback", "d_auc_clock"), #clock-aligned
    c("clock", "feedback", "pe_max"),
#    c("clock", "feedback", "v_entropy_no5"),
    c("clock", "feedback", "v_auc"),
#    c("clock", "feedback", "d_auc_sqrt"),
#    c("clock", "feedback", "rt_swing"),
#    c("clock", "feedback", "rt_swing_sqrt"),
    c("clock", "feedback", "v_max"),
    c("clock", "feedback", "mean_kld"),
    c("clock", "feedback", "intrinsic_discrepancy"),
    c("clock", "feedback", "mean_kld_feedback"),
    c("clock", "feedback", "intrinsic_discrepancy_feedback")
#    c("clock", "feedback", "rew_om"),
#    c("clock", "feedback", "pe_max", "rew_om"),
#    c("clock_bs", "feedback")
  ),
  group_model_variants=list(
    c("Intercept"),
    c("Intercept", "Age")
#    c("Intercept", "Age", "Female"),
#    c("Intercept", "I_Age"),
#    c("Intercept", "I_Age", "Female")
  ),    
  execute_feat=FALSE, #passed through to fsl_sceptic_model to create fsf, but not run the model
  #model_suffix="_fse", #factorized, selective, equal generalization width
  model_suffix="_fse_groupfixed", #factorized, selective, equal generalization width
  root_workdir="/gpfs/scratch/mnh5174/run_fsl_pipeline_qsub_tmp",
  n_cluster_beta_cpus=8 #should be number of l2 contrasts, or lower
)

#validate and populate any other pipeline details before execution
fsl_model_arguments <- finalize_pipeline_configuration(fsl_model_arguments)

save(fsl_model_arguments, file=paste0("configuration_files/", fsl_model_arguments$analysis_name, ".RData"))

#this pushes the full analysis pipeline in parallel, where parallelism is across sceptic_run_variants
push_pipeline(fsl_model_arguments, ncpus=fsl_model_arguments$pipeline_cpus)


## Other pipelines go here
