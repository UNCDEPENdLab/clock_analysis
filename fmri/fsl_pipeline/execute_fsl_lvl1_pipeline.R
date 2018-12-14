#load the specified RData object and call the LVL fitting function requested

to_run <- Sys.getenv("fsl_pipeline_file")
run_model_index <- as.numeric(Sys.getenv("run_model_index")) #which variant to execute
if (nchar(to_run) == 0L) { stop("Cannot locate environment variable fsl_pipeline_file") }
if (!file.exists(to_run)) { stop("Cannot locate configuration file", to_run) }
if (is.na(run_model_index)) { stop("Couldn't identify usable run_model_index variable.") }

load(to_run)

library(dependlab)
library(foreach)

source(file.path(fsl_model_arguments$pipeline_home, "functions", "model_clock_fmri_lvl1.R"))
source(file.path(fsl_model_arguments$pipeline_home, "functions", "fsl_sceptic_model.R"))
source(file.path(fsl_model_arguments$pipeline_home, "functions", "glm_helper_functions.R"))


#specify the run model for the current batch, based on run_model_index (see run_fsl_pipeline.R)
fsl_model_arguments$sceptic_run_signals <- fsl_model_arguments$sceptic_run_variants[[run_model_index]]
fsl_model_arguments$sceptic_run_variants <- NULL #the broader set of model variants is not handled by subsidary scripts (i.e., all assume a single variant)
fsl_model_arguments$outdir <- fsl_model_arguments$outdir[run_model_index]
fsl_model_arguments$workdir <- NULL

#fsl_model_arguments$ncpus <- 1

#overuse of ... in subsidiary calls to model_clock_fmri_lvl1 will eventually blow up at the build_design_matrix step... need a better solution to make effective use of the ... argument
#this is my hack for now
lvl1_args <- fsl_model_arguments[c("trial_statistics", "fmri_dir", "idexpr",
  "drop_volumes", "ncpus", "expectdir", "expectfile",
  "sceptic_run_signals", "outdir", "usepreconvolve", "spikeregressors", "model_suffix",
  "execute_feat", "tr")]

#use the list of arguments loaded from the configuration file to call the subject x run FEAT setup function
do.call(model_clock_fmri_lvl1, lvl1_args)

