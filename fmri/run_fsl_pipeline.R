#the goal of this script is to run an entire fmri analysis for SCEPTIC data, including levels 1-3 in FSL

library(dependlab)
library(foreach)
library(parallel)
library(doParallel)
library(plyr)
library(tidyverse)

scripts_dir <- "/gpfs/group/mnh5174/default/clock_analysis/fmri"
setwd(scripts_dir)
source("fsl_sceptic_model.R")
source("glm_helper_functions.R")
##source("r_glm.R")
source("model_clock_fmri_lvl1.R")
source(file.path(scripts_dir, "run_fsl_lvl1_sepqsub.R")) #executes FSF files in parallel batches

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
fsl_model_arguments <- list(
  analysis_name="MMClock_aroma_preconvolve_fse",
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
  idregex="([0-9]{5})_\\d+", #5 digit ID, followed by irrelevant date. A bit inelegant, but used in setup_fsl_lvl2_inputs to infer the subject's ID from the MR folder structure
  sceptic_run_variants=list(
    c("v_chosen", "v_entropy", "d_auc", "pe_max"), #all signals with entropy of weights
    c("v_chosen", "v_entropy_func", "d_auc", "pe_max"), #all signals with entropy of evaluated function
    c("v_chosen"), #individual regressors
    c("v_entropy"),
    c("v_entropy_func"),
    c("d_auc"),
    c("pe_max")
  ),
  execute_feat=FALSE, #passed through to fsl_sceptic_model to create fsf, but not run the model
  model_suffix="_fse", #factorized, selective, equal generalization width
  root_workdir="/gpfs/scratch/mnh5174/run_fsl_pipeline_qsub_tmp"
)

fsl_model_arguments$outdir <- sapply(fsl_model_arguments$sceptic_run_variants, function(x) {
  paste0("sceptic-", paste(x, collapse="-"), #define output directory based on combination of signals requested
    ifelse(fsl_model_arguments$usepreconvolve, "-preconvolve", ""),
    fsl_model_arguments$model_suffix)
})

fsl_model_arguments$n_l1_copes <- sapply(fsl_model_arguments$sceptic_run_variants, function(x) { length(x) + 2 }) #compute number of l1 copes for each variant
fsl_model_arguments$workdir <- file.path(fsl_model_arguments$root_workdir, fsl_model_arguments$outdir) #temp folder for each analysis variant

save(fsl_model_arguments, file=paste0("fsl_pipeline_configurations/", fsl_model_arguments$analysis_name, ".RData"))

#push each model variant through the full pipeline
#eventually, when we have l3 covariates, should probably allow l1 x l3 combinations

#LOOP OVER MODEL VARIANTS IN PARALLEL

push_pipeline <- function(fsl_model_arguments, ncpus=1, runpar=TRUE) { #shoudl this read from the fsl list, which also has runpar and ncpus??
  #this helper script walks through each variant of the level1 model and calls scripts to run the full analysis pipeline
  #this amounts to iterating over each of the sceptic_run_variants
  stopifnot(length(fsl_model_arguments$outdir) == length(fsl_model_arguments$sceptic_run_variants))

  require(parallel)
  require(doParallel)
  
  #setup parallel worker pool, if requested
  if (runpar) {
    cl <- makeCluster(ncpus)
    registerDoParallel(cl)
    
    on.exit(try(stopCluster(cl))) #cleanup pool upon exit of this function
  } else {
    registerDoSEQ() #formally register a sequential 'pool' so that dopar is okay
  }

  #iterate over lvl1/run variants
  nothing <- foreach(run_model_index=iter(1:length(fsl_model_arguments$outdir)), .export=c("qsub_file", "run_fsl_lvl1_sepqsub", "wait_for_job")) %dopar% {

    #setup FSF files for all runs in one qsub
    #then, call run_fsl_lvl1_sepqsub.R to execute all .fsf files using separate qsub instances
    ##THIS PRETTY MUCH JUST CALLS model_clock_fmri_lvl1
    
    setup_fsf_jobid <- qsub_file(script="qsub_fmri_r_script.bash",
      pbs_args=c("-l nodes=1:ppn=20", "-l walltime=10:00:00"),
      env_variables=c(R_SCRIPT="execute_fsl_lvl1_pipeline.R",
        run_model_index=run_model_index,
        fsl_pipeline_file=paste0("fsl_pipeline_configurations/", fsl_model_arguments$analysis_name, ".RData"))
    )

    #bring the LVL1 sep qsub out to here.
    #Create batched qsubs for all FSL LVL1 analyses after FSF setup completes
    ##Sys.setenv(TARGET=fsl_model_arguments$fmri_dir, MODEL_MATCH=fsl_model_arguments$outdir) #, WAIT_FOR=setup_fsf_jobid)

    #Okay, here's the qsub dilemma: We need for model_clock_fmri_lvl1 to complete before we try to run
    #run_fsl_lvl1_sepqsub. We could submit this as a job that depends on setup_fsf_jobid above
    #But! We need all of the subsidiary qsub jobs that are spawned by run_fsl_lvl1_sepqsub to be available/defined
    #when we call execute_fsl_lvl2_pipeline.R. So, if we queue a qsub job to run_fsl_lvl1_sepqsub that depends
    #on setup_fsf_jobid, the code will continue to the lvl2 pipeline before the multiple lvl1 pipeline jobids
    #are even defined. Thus, my provisional, if slightly hacky solution is to wait here in the code until the
    #model_clock_fmri_lvl1 job completes so that we can execute run_fsl_lvl1_sepqsub directly as a function
    #and recover the separate jobids that need to be passed to the level 2 script.

    wait_for_job(setup_fsf_jobid) #pause R script here until model_clock_fmri_lvl1 above completes

    #do not pass a wait_for here, as it will lead to an 'invalid job dependency' error during qsub (i.e., the parent job is already complete!).
    #this occurs because we now pause in the execution while the model_clock_fmri_lvl1 finishes in the previous step
    sep_lvl1_jobs <- run_fsl_lvl1_sepqsub(fsl_model_arguments, run_model_index, rerun=FALSE) #, wait_for=setup_fsf_jobid)

    ## lvl1_setup_jobid <- qsub_file(script="qsub_fmri_r_script.bash",
    ##   pbs_args=c("-l nodes=1:ppn=20", "-l walltime=10:00:00")
    ##   env_variables=c("R_SCRIPT=run_fsl_lvl1_sepqsub.R",
    ##     paste0("WAIT_FOR=", setup_fsf_jobid),
    ##     paste0("run_model_index=", run_model_index),
    ##     paste0("fsl_pipeline_file=fsl_pipeline_configurations/", fsl_model_arguments$analysis_name, ".RData"))
    ## )
    
    #NB. run_fsl_lvl1_sepqsub.R will create a file called sepqsub_lvl1_jobs.txt in the relevant temporary directory in the scratch folder
    
    # Run lvl2 analyses for the current model. The execute_fsl_lvl2_pipeline.R script will look for lvl1 jobs that are still running.
    # In terms of dependencies, because the qsub above fires immediately, the lvl2 needs to wait on:
    #    a) the job id from execute_fsl_lvl1_pipeline.bash
    #    b) all subsidiary job ids from each batch created by run_fsl_lvl1_sepqsub.R
    #
    # Notably, the subsidiary job ids for the separate qsubs should be available/spawned as part of the parent job (execute_fsl_lvl1_pipeline.bash),
    # so waiting on both sets for LVL 2 should be a safe approach.

    #the problem is that run_fsl_lvl1_sepqsub.R will be queued, but the subidiary jobs have not been created yet. Thus, the check below for separate qsub jobs will
    #probably be invalid.
    
    #this is now handled as a function return above
    #check in on the separate l1 jobs to wait on
    ## if (file.exists(l1_jobfile <- file.path(fsl_model_arguments$workdir[m], "sepqsub_lvl1_jobs.txt"))) {
    ##   sep_lvl1_jobs <- readLines(l1_jobfile)
    ## } else {
    ##   sep_lvl1_jobs <- NULL
    ## }

    if (!is.null(sep_lvl1_jobs)) {
      wait_string <- paste0("-W depend:afterok:", paste(sep_lvl1_jobs, collapse=":")) #setup_fsf_jobid, 
    } else { wait_string <- NULL }
      
    l2_execution_jobid <- qsub_file(script="qsub_fmri_r_script.bash",
      pbs_args=c("-l nodes=1:ppn=20", "-l walltime=10:00:00", wait_string),
      env_variables=c(
        R_SCRIPT="execute_fsl_lvl2_pipeline.R",
        run_model_index=run_model_index,
        fsl_pipeline_file=paste0("fsl_pipeline_configurations/", fsl_model_arguments$analysis_name, ".RData"))
    )
    
  }

}

push_pipeline(fsl_model_arguments, runpar=TRUE, ncpus=20)





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
