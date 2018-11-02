#the goal of this script is to run an entire fmri analysis for SCEPTIC data, including levels 1-3 in FSL

library(dependlab)
library(foreach)
library(doSNOW)

setwd(file.path(getMainDir(), "clock_analysis", "fmri"))
source("fsl_sceptic_model.R")
source("glm_helper_functions.R")
##source("r_glm.R")
source("model_clock_fmri_lvl1.R")

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
#SCEPTIC MMClock Fit
#trial_df <- read.csv("/gpfs/group/mnh5174/default/temporal_instrumental_agent/clock_task/vba_fmri/vba_out/compiled_outputs/mmclock_fmri_decay_factorize_mfx_trial_statistics.csv.gz")
#trial_df <- read.csv("/gpfs/group/mnh5174/default/temporal_instrumental_agent/clock_task/vba_fmri/vba_out/compiled_outputs/mmclock_fmri_decay_mfx_trial_statistics.csv.gz")

#factorized, selective maintenance, equal basis-generalization width
trial_df <- read.csv("/gpfs/group/mnh5174/default/temporal_instrumental_agent/clock_task/vba_fmri/vba_out/compiled_outputs/mmclock_fmri_decay_factorize_selective_psequate_mfx_trial_statistics.csv.gz")

#from 2017:
##results from Mean SCEPTIC regressor correlation.pdf indicate that regressors for vchosen, ventropy_decay_matlab, dauc, and pemax are
##reasonably uncorrelated. The worst is dauc with vchosen (mean r = -0.31), which makes sense that as learning progresses, chosen values
##are higher and there is less residue to decay. These 4 regressors are also of greatest theoretical interest


#create a list of arguments, save it to an RData object, load this in a qsub, then use do.call() to populate arguments to fit_all_fmri
lvl1_model_arguments <- list(
  trial_statistics = trial_df,
  fmri_dir = "/gpfs/group/mnh5174/default/MMClock/MR_Proc",
  expectdir = "mni_5mm_aroma", #subfolder name for processed data
  expectfile = "nfaswuktm_clock[0-9]_5.nii.gz", #expected file name for processed clock data
  usepreconvolve=TRUE,
  runpar=TRUE,
  ncpus=20,
  drop_volumes=2, #to handle steady state concerns
  tr=1.0, #seconds
  spikeregressors=FALSE, #don't include spike regressors in nuisance variables since we are using AROMA
  idexpr=expression(subid), #how to match between the subject ID in trial_statistics and the folder structure on the file system
  sceptic_run_variants=list(
    c("v_chosen", "v_entropy", "d_auc", "pe_max"), #all signals with entropy of weights
    c("v_chosen", "v_entropy_func", "d_auc", "pe_max"), #all signals with entropy of evaluated function
    c("v_chosen"), #individual regressors
    c("v_entropy"),
    c("v_entropy_func"),
    c("d_auc")
  ),
  execute_feat=FALSE, #passed through to fsl_sceptic_model to create fsf, but not run the model
  model_suffix="_fse" #factorized, selective, equal generalization width
)

lvl1_model_arguments$outdir <- sapply(lvl1_model_arguments$sceptic_run_variants, function(x) {
  paste0("sceptic-", paste(x, collapse="-"),      #define output directory based on combination of signals requested
    ifelse(lvl1_model_arguments$usepreconvolve, "-preconvolve", ""),
    lvl1_model_arguments$model_suffix)
})

save(lvl1_model_arguments, file="fsl_pipeline_configurations/MMClock_aroma_preconvolve.RData")

#setup FSF files for all runs
setup_fsf_jobid <- qsub_file(script="execute_fsl_lvl1_pipeline.bash",
  pbs_args=c("-l nodes=1:ppn=20",
    "-l walltime=10:00:00",
    "-v fsl_lvl1_pipeline_file=fsl_pipeline_configurations/MMClock_aroma_preconvolve.RData"))








###


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
#fit_all_fmri(trial_statistics=trial_df, iddf = idinfo, model="sceptic", usepreconvolve=TRUE, parmax1=TRUE, runpar=FALSE, ncpus=1) #rescale to 1.0 max

###################################################
