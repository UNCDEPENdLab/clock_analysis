#This is a small helper function to validate the fsl_model_arguments list structure.
#It adds a few details such as the output directory to make it less burdensome for to setup a pipeline
finalize_pipeline_configuration <- function(fsl_model_arguments) {
  fsl_model_arguments$outdir <- sapply(fsl_model_arguments$sceptic_run_variants, function(x) {
    paste0("sceptic-", paste(x, collapse="-"), #define output directory based on combination of signals requested
      ifelse(fsl_model_arguments$usepreconvolve, "-preconvolve", ""),
      fsl_model_arguments$model_suffix)
  })

  fsl_model_arguments$n_l1_copes <- sapply(fsl_model_arguments$sceptic_run_variants, function(x) { length(x) + 2 }) #compute number of l1 copes for each variant
  fsl_model_arguments$workdir <- file.path(fsl_model_arguments$root_workdir, fsl_model_arguments$outdir) #temp folder for each analysis variant
  fsl_model_arguments$pipeline_cpus <- length(fsl_model_arguments$sceptic_run_variants) #number of workers to setup at the pipeline level (i.e., over run variants)
  if (is.null(fsl_model_arguments$l2_cpus)) { fsl_model_arguments$l2_cpus <- 20 } #number of cores to use in Feat LVL2 analyses (fixed effects combination of runs)
  if (is.null(fsl_model_arguments$pipeline_home)) { fsl_model_arguments$pipeline_home <- "/gpfs/group/mnh5174/default/clock_analysis/fmri/fsl_pipeline" }
  
  return(fsl_model_arguments)
}
