#This is a small helper function to validate the fsl_model_arguments list structure.
#It adds a few details such as the output directory to make it less burdensome for to setup a pipeline
finalize_pipeline_configuration <- function(fsl_model_arguments) {
  fsl_model_arguments$outdir <- sapply(fsl_model_arguments$sceptic_run_variants, function(x) {
    paste0("sceptic-", paste(x, collapse="-"), #define output directory based on combination of signals requested
      ifelse(fsl_model_arguments$usepreconvolve, "-preconvolve", ""),
      fsl_model_arguments$model_suffix)
  })

  fsl_model_arguments$n_l1_copes <- sapply(fsl_model_arguments$sceptic_run_variants, function(x) { length(x) }) #compute number of l1 copes for each variant
  fsl_model_arguments$workdir <- file.path(fsl_model_arguments$root_workdir, fsl_model_arguments$outdir) #temp folder for each analysis variant
  fsl_model_arguments$pipeline_cpus <- length(fsl_model_arguments$sceptic_run_variants) #number of workers to setup at the pipeline level (i.e., over run variants)
  if (is.null(fsl_model_arguments$l2_cpus)) { fsl_model_arguments$l2_cpus <- 20 } #number of cores to use in Feat LVL2 analyses (fixed effects combination of runs)
  if (is.null(fsl_model_arguments$pipeline_home)) { fsl_model_arguments$pipeline_home <- "/gpfs/group/mnh5174/default/clock_analysis/fmri/fsl_pipeline" }
  if (is.null(fsl_model_arguments$group_output_dir)) { fsl_model_arguments$group_output_dir <- file.path(dirname(fsl_model_arguments$fmri_dir), "group_analyses", fsl_model_arguments$analysis_name) }
  if (is.null(fsl_model_arguments$center_l3_predictors)) { fsl_model_arguments$center_l3_predictors <- TRUE }
  if (is.null(fsl_model_arguments$badids)) { fsl_model_arguments$badids <- c() }
  if (is.null(fsl_model_arguments$l1_cope_names)) {
    fsl_model_arguments$l1_cope_names <- lapply(fsl_model_arguments$sceptic_run_variants, function(x) {
      signal_copes <- x
      names(signal_copes) <- paste0("cope", 1:length(x))
      return(signal_copes)
    })
  }

  if (is.null(fsl_model_arguments$zthresh)) { fsl_model_arguments$zthresh <- 3.09 }  #1-tailed p=.001 for z stat
  if (is.null(fsl_model_arguments$clustsize)) { fsl_model_arguments$clustsize <- 34 } #based on 3dClustSim using ACFs for first-level FEAT runs

  #ensure that the user has specified some sort of clock event in the model
  for (v in fsl_model_arguments$sceptic_run_variants) {
    if (!any(c("clock", "clock_bs") %in% v)) {
      stop("No clock event is in the model: ", paste(v, collapse=","))
    }
    if (!any(c("feedback", "feedback_bs") %in% v)) {
      stop("No feedback event is in the model: ", paste(v, collapse=","))
    }
  }  

  #remove bad ids before running anything further
  if (!is.null(fsl_model_arguments$badids) && length(fsl_model_arguments$badids) > 0L) {
    fsl_model_arguments$subject_covariates <- fsl_model_arguments$subject_covariates %>% filter(!id %in% fsl_model_arguments$badids) #remove bad ids
  }
  
  return(fsl_model_arguments)
}
