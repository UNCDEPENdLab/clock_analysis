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


#################
# FSL LEVEL 1

#N.B. in examining initial results from single subject analyses, it is clear that steady state magnetization is not achieved by the first volume acquired
#ICA analysis suggests that it takes up to 6 volumes to reach steady state, and the rel and mean uncertainty maps are being adversely affected by this problem
#because they also start high and decay... Mean uncertainty was consequently soaking up a huge amount of CSF in activation maps.
#Because the first presentation occurs at 8 seconds, it seems fine to drop 6 volumes (6s) 

#Jun2017: further ICAs on these data do not suggest a long steady-state problem. Drop 2 volumes for good measure

#abstract all necessary run specifications to this script
#write an RData object that contains all necessary arguments to call model_clock_fmri_lvl1.R

###
# SCEPTIC MMClock Y3

#trial_df <- read.csv("/gpfs/group/mnh5174/default/temporal_instrumental_agent/clock_task/vba_fmri/vba_out/compiled_outputs/mmclock_fmri_decay_factorize_mfx_trial_statistics.csv.gz")
#trial_df <- read.csv("/gpfs/group/mnh5174/default/temporal_instrumental_agent/clock_task/vba_fmri/vba_out/compiled_outputs/mmclock_fmri_decay_mfx_trial_statistics.csv.gz")


### LEFTOVERS FROM OLD LVL3

####
##SCEPTIC L3 ANALYSIS
##models <- c("sceptic_dauc_preconvolve", "sceptic_pemax_preconvolve", "sceptic_vmax_preconvolve", "sceptic_vchosen_preconvolve", "sceptic_ventropy_preconvolve")
## library(dplyr)
## options(dplyr.width=200)
## #preproc_dirname <- "mni_5mm_aroma"
## preproc_dirname <- "mni_5mm_3ddespike"
## fsl_dirname <- "FEAT_LVL2_runtrend.gfeat"
## ##models <- c("sceptic_vchosen_ventropy_decay_matlab_dauc_pemax_preconvolve")
## models <- c("sceptic_vchosen_ventropy_dauc_pemax_preconvolve")

## #dataset <- "MMY3"
## dataset <- "SPECC"

## if (dataset=="MMY3") {
##   subinfo <- read.table("/gpfs/group/mnh5174/default/clock_analysis/fmri/subinfo_db", header=TRUE)
##   subinfo <- subinfo %>% rename(ID=lunaid, Age=age, Female=female, ScanDate=scandate)
##   subinfo$mr_dir <- with(subinfo, paste0("/gpfs/group/mnh5174/default/MMClock/MR_Proc/", ID, "_", format((as.Date(ScanDate, format="%Y-%m-%d")), "%Y%m%d"))) #convert to Date, then reformat YYYYMMDD
##   outdir <- "/gpfs/group/mnh5174/default/MMClock/fsl_sceptic_group_jun2017"
## } else if (dataset=="SPECC") {
##   subinfo <- read.csv("/gpfs/group/mnh5174/default/SPECC/SPECC_Participant_Info.csv", header=TRUE)
##   subinfo <- subinfo %>% rowwise() %>% mutate(mr_dir=ifelse(LunaMRI==1,
##     paste0("/gpfs/group/mnh5174/default/MMClock/MR_Proc/", Luna_ID, "_", format((as.Date(ScanDate, format="%Y-%m-%d")), "%Y%m%d")), #convert to Date, then reformat YYYYMMDD
##     paste0("/gpfs/group/mnh5174/default/SPECC/MR_Proc/", tolower(SPECC_ID), "_", tolower(format((as.Date(ScanDate, format="%Y-%m-%d")), "%d%b%Y")))))
##   subinfo <- subinfo %>% rename(ID=SPECC_ID, Age=AgeAtScan)
##   outdir <- "/gpfs/group/mnh5174/default/SPECC/fsl_sceptic_group"  
## }


#factorized, selective maintenance, equal basis-generalization width
trial_df <- read.csv("/gpfs/group/mnh5174/default/temporal_instrumental_agent/clock_task/vba_fmri/vba_out/compiled_outputs/mmclock_fmri_decay_factorize_selective_psequate_mfx_trial_statistics.csv.gz")
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
  analysis_name="MMClock_aroma_preconvolve_fse",
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
    c("v_chosen", "v_entropy", "d_auc", "pe_max"), #all signals with entropy of weights
    c("v_chosen", "v_entropy_func", "d_auc", "pe_max"), #all signals with entropy of evaluated function
    c("v_chosen"), #individual regressors
    c("v_entropy"),
    c("v_entropy_func"),
    c("d_auc"),
    c("pe_max")
  ),
  group_model_variants=list(
    c("Intercept"),
    c("Intercept", "Age"),
    c("Intercept", "Age", "Female"),
    c("Intercept", "I_Age"),
    c("Intercept", "I_Age", "Female")
  ),    
  execute_feat=FALSE, #passed through to fsl_sceptic_model to create fsf, but not run the model
  model_suffix="_fse", #factorized, selective, equal generalization width
  root_workdir="/gpfs/scratch/mnh5174/run_fsl_pipeline_qsub_tmp"
)

#validate and populate any other pipeline details before execution
fsl_model_arguments <- finalize_pipeline_configuration(fsl_model_arguments)

save(fsl_model_arguments, file=paste0("configuration_files/", fsl_model_arguments$analysis_name, ".RData"))

#this pushes the full analysis pipeline in parallel, where parallelism is across sceptic_run_variants
push_pipeline(fsl_model_arguments, ncpus=fsl_model_arguments$pipeline_cpus)


## Other pipelines go here









### Some leftovers


## I have now converted all SPECC MR directory names to all lower case to allow for match on case-sensitive filesystem and to make the naming consistent
## idfile <- "/gpfs/group/mnh5174/default/SPECC/SPECC_Participant_Info.csv"
## idinfo <- read.csv(idfile)
## library(dplyr)
## options(dplyr.width=200)
## idinfo <- idinfo %>% rowwise() %>% mutate(mr_dir=ifelse(LunaMRI==1,
##   paste0("/gpfs/group/mnh5174/default/MMClock/MR_Proc/", Luna_ID, "_", format((as.Date(ScanDate, format="%Y-%m-%d")), "%Y%m%d")), #convert to Date, then reformat YYYYMMDD
##   paste0("/gpfs/group/mnh5174/default/SPECC/MR_Proc/", tolower(SPECC_ID), "_", tolower(format((as.Date(ScanDate, format="%Y-%m-%d")), "%d%b%Y"))))) %>% ungroup()

## ##verify that mr_dir is present as expected
## idinfo$dirfound <- file.exists(idinfo$mr_dir)
## subset(idinfo, dirfound==FALSE)

##subject CSVs in subjects/SPECC are names according to numeric SPECC_ID
##need to use idinfo data.frame to line up with MMClock, look in Luna dir as needed, etc.
#need SPECC here...
#trial_df <- read.csv("/gpfs/group/mnh5174/default/temporal_instrumental_agent/clock_task/vba_fmri/vba_out/compiled_outputs/mmclock_fmri_decay_factorize_mfx_trial_statistics.csv.gz")
#fit_all_fmri(trial_statistics=trial_df, iddf = idinfo, model="sceptic", usepreconvolve=TRUE, parmax1=TRUE, ncpus=1) #rescale to 1.0 max

###################################################
