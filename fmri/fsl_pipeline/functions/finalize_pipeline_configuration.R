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
  if (is.null(fsl_model_arguments$group_output_dir)) { fsl_model_arguments$group_output_dir <- file.path(dirname(fsl_model_arguments$fmri_dir), "group_analyses", fsl_model_arguments$analysis_name) }
  if (is.null(fsl_model_arguments$center_l3_predictors)) { fsl_model_arguments$center_l3_predictors <- TRUE }
  if (is.null(fsl_model_arguments$l1_cope_names)) {
    fsl_model_arguments$l1_cope_names <- sapply(fsl_model_arguments$sceptic_run_variants, function(x) {
      signal_copes <- x
      names(signal_copes) <- paste0("cope", 2 + 1:length(x))
      cope_names <- c(cope1="clock_onset", cope2="feedback_onset", signal_copes)
      return(cope_names)
    })
  }

  if (is.null(fsl_model_arguments$zthresh)) { fsl_model_arguments$zthresh <- 3.09 }  #1-tailed p=.001 for z stat
  if (is.null(fsl_model_arguments$clustsize)) { fsl_model_arguments$clustsize <- 34 } #based on 3dClustSim using ACFs for first-level FEAT runs
  
  return(fsl_model_arguments)
}
