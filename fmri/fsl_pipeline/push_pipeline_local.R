##Push pipeline (NON-CLUSTER VER):
##Maybe intergration in the future but not for now...
push_pipeline <- function(fsl_model_arguments, ncpus=1) {
stopifnot(length(fsl_model_arguments$outdir) == length(fsl_model_arguments$sceptic_run_variants))
stopifnot(is.numeric(ncpus) && ncpus >= 1)

require(parallel)
require(doParallel)
require(dependlab) #has qsub_file and wait_for_job
#require(Rniftilib)

source(file.path(fsl_model_arguments$pipeline_home, "functions", "glm_helper_functions.R"))
#source(file.path(fsl_model_arguments$pipeline_home, "functions", "run_feat_lvl1_sepqsub.R")) #executes FSF files in parallel batches


set_sysenvr<-function(env_variables=list()){
  do.call(Sys.setenv,env_variables)
  #source(R_SCRIPT)
}

setwd(fsl_model_arguments$pipeline_home) #to make qsub_file calls below happy with local paths

Sys.setenv(CLSTYLE="local")
#Let's be honest, without a cluster you won't be able to parallel models
#Generate all the fsf first;


for (run_model_index in 3:length(fsl_model_arguments$outdir)) {
  fsl_model_arguments$execute_feat<-FALSE
  set_sysenvr(env_variables=list(
    run_model_index=run_model_index,
    fsl_pipeline_file=file.path(fsl_model_arguments$pipeline_home, "configuration_files", paste0(fsl_model_arguments$analysis_name, ".RData")))
  )
  source("execute_fsl_lvl1_pipeline.R")

}
}
