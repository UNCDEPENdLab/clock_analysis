to_run <- Sys.getenv("fsl_pipeline_file")

run_model_index <- as.numeric(Sys.getenv("run_model_index")) #which variant to execute
if (nchar(to_run) == 0L) { stop("Cannot locate environment variable fsl_pipeline_file") }
if (!file.exists(to_run)) { stop("Cannot locate configuration file", to_run) }
if (is.na(run_model_index)) { stop("Couldn't identify usable run_model_index variable.") }

load(to_run)

library(tidyverse)

scripts_dir <- "/gpfs/group/mnh5174/default/clock_analysis/fmri"

source(file.path(scripts_dir, "glm_helper_functions.R"))
source(file.path(scripts_dir, "setup_fsl_lvl2_inputs.R"))
source(file.path(scripts_dir, "run_feat_lvl2.R"))

#first, setup the inputs for LVL2 analysis (doesn't really benefit from parallel execution)
feat_l2_inputs_df <- setup_fsl_lvl2_inputs(fsl_model_arguments, run_model_index)

run_feat_lvl2(feat_l2_inputs_df, run=TRUE, force=FALSE, ncpus=20)
