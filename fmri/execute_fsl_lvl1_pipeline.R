#load the specified RData object and call the LVL fitting function requested

to_run <- Sys.getenv("fsl_lvl1_pipeline_file")
if (nchar(to_run) == 0L) { stop("Cannot locate environment variable fsl_lvl1_pipeline_file") }
if (!file.exists(to_run)) { stop("Cannot locate configuration file", to_run) }

load(to_run)

library(dependlab)
library(foreach)
library(doSNOW)

scripts_dir <- "/gpfs/group/mnh5174/default/clock_analysis/fmri"

source(file.path(scripts_dir, "fsl_sceptic_model.R"))
source(file.path(scripts_dir, "glm_helper_functions.R"))
source(file.path(scripts_dir, "model_clock_fmri_lvl1.R"))

#use the list of arguments loaded from the configuration file to call the subject x run FEAT setup function
do.call(model_clock_fmri_lvl1, lvl1_model_arguments)

#Create batched qsubs for all FSL LVL1 analyses after FSF setup completes
for (v in 1:length(lvl1_model_arguments$outdir)) {
  Sys.setenv(TARGET=lvl1_model_arguments$fmri_dir, MODEL_MATCH=lvl1_model_arguments$outdir[v]) #, WAIT_FOR=setup_fsf_jobid)
  source(file.path(scripts_dir, "run_fsl_lvl1_sepqsub.R")) #executes FSF files in parallel batches
}
