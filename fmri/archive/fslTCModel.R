fslTCModel <- function(f, mrfiles, runlengths, mrrunnums, run=FALSE, force=FALSE, dropVolumes=0, f_value=NULL, outdir="fsl_hybrid_carryvalue") {
  require(Rniftilib)
  require(parallel)
  if (grepl('nomeanunc', outdir)) {
    fsfTemplate <- readLines(file.path(getMainDir(), "clock_analysis", "fmri", "feat_lvl1_clock_tc_template_nomeanunc.fsf"))  
  } else {
    fsfTemplate <- readLines(file.path(getMainDir(), "clock_analysis", "fmri", "feat_lvl1_clock_tc_template.fsf"))
  }
  
  
  fsldir <- file.path(normalizePath(file.path(dirname(mrfiles[1L]), "..")), outdir) #note: normalizePath will fail to evaluate properly if directory does not exist

  if (file.exists(fsldir) && force==FALSE) { message(fsldir, " exists. Skipping."); return(0) }
  cat("fsldir create: ", fsldir, "\n")
  dir.create(fsldir, showWarnings=FALSE) #one directory up from a given clock run
  timingdir <- file.path(fsldir, "run_timing_tc")

  d <- build_design_matrix(fitobj=f, regressors=c("clock", "mean_uncertainty", "rel_uncertainty", "ev", "feedback", "rpe_pos", "rpe_neg"), 
      event_onsets=c("clock_onset", "clock_onset", "clock_onset", "clock_onset", "feedback_onset", "feedback_onset", "feedback_onset"), 
      durations=c("clock_duration", "clock_duration", "clock_duration", "clock_duration", "feedback_duration", "feedback_duration", "feedback_duration"), 
      normalizations=c("durmax_1.0", "evtmax_1.0", "evtmax_1.0", "evtmax_1.0", "durmax_1.0", "evtmax_1.0", "evtmax_1.0"), 
      baselineCoefOrder=2, writeTimingFiles=c("FSL"), center_values=TRUE, convolve_wi_run=TRUE,      
      runVolumes=runlengths, runsToOutput=mrrunnums, output_directory=timingdir, dropVolumes=dropVolumes)
  
  #allow for hybrid TC + R-W model
  #if R-W fit is passed in, use estimates of EV, RPE+, and RPE- to generate FSL 3-column regressors
  #these will overwrite the TC-based estimates above and thus result in the desired hybrid model 
  if (!is.null(f_value)) {
    d_value <- build_design_matrix(fitobj=f_value, regressors=c("clock", "ev", "feedback", "rpe_neg", "rpe_pos"), 
        event_onsets=c("clock_onset", "clock_onset", "feedback_onset", "feedback_onset", "feedback_onset"), 
        durations=c("clock_duration", "clock_duration", "feedback_duration", "feedback_duration", "feedback_duration"),
        normalizations=c("durmax_1.0", "evtmax_1.0", "durmax_1.0", "evtmax_1.0", "evtmax_1.0"),
        baselineCoefOrder=2, writeTimingFiles=c("FSL"), center_values=TRUE, convolve_wi_run = TRUE,
        runVolumes=runlengths, runsToOutput=mrrunnums, output_directory=timingdir, dropVolumes=dropVolumes)
  }
  
  allFeatFiles <- list()
  
  #FSL computes first-level models on individual runs
  for (r in 1:length(mrfiles)) {
    stopifnot(file.exists(file.path(dirname(mrfiles[r]), "motion.par"))) #can't find motion parameters
    
    runnum <- sub("^.*/clock(\\d+)$", "\\1", dirname(mrfiles[r]), perl=TRUE)
    nvol <- nifti.image.read(mrfiles[r], read_data=0)$dim[4L]

    #just PCA motion on the current run
    mregressors <- pca_motion(mrfiles[r], runlengths[r], motion_parfile="motion.par", numpcs=3, dropVolumes=dropVolumes)$motion_pcs_concat
        
    #Add volumes to censor here. Use censor_intersection.mat, which flags fd > 0.9 and DVARS > 20
    censorfile <- file.path(dirname(mrfiles[r]), "motion_info", "censor_intersection.mat")
    if (file.exists(censorfile) && file.info(censorfile)$size > 0) {
      censor <- read.table(censorfile, header=FALSE)$V1
      censor <- censor[(1+dropVolumes):runlengths[r]]
      mregressors <- cbind(mregressors, censor)
    }
    
    #add CSF and WM regressors (with their derivatives)
    nuisancefile <- file.path(dirname(mrfiles[r]), "nuisance_regressors.txt")
    if (file.exists(nuisancefile)) {
      nuisance <- read.table(nuisancefile, header=FALSE)
      nuisance <- nuisance[(1+dropVolumes):runlengths[r],,drop=FALSE]
      nuisance <- as.data.frame(lapply(nuisance, function(col) { col - mean(col) })) #demean
      mregressors <- cbind(mregressors, nuisance)
    }
    
    motfile <- file.path(fsldir, paste0("run", runnum, "_confounds.txt"))
    write.table(mregressors, file=motfile, col.names=FALSE, row.names=FALSE)
    
    #search and replace within fsf file for appropriate sections
    ##.OUTPUTDIR. is the feat output location
    ##.NVOL. is the number of volumes in the run
    ##.FUNCTIONAL. is the fmri data to process (sans extension)
    ##.CONFOUNDS. is the confounds file for GLM
    ##.CLOCK_TIMES. is the three-column file for clock onset
    ##.FEEDBACK_TIMES. is the three-column file for feedback onset
    ##.RELUNC_TIMES. is the three-column file for relative uncertainty
    ##.MEANUNC_TIMES. is the three-column file for mean uncertainty
    ##.EV_TIMES. is the three-column file for expected value
    ##.RPEPOS_TIMES. is the three-column file for positive rpes
    ##.RPENEG_TIMES. is the three-column file for negative rpes
    
    thisTemplate <- fsfTemplate
    thisTemplate <- gsub(".OUTPUTDIR.", file.path(fsldir, paste0("FEAT_LVL1_run", runnum)), thisTemplate, fixed=TRUE)
    thisTemplate <- gsub(".NVOL.", nvol, thisTemplate, fixed=TRUE)
    thisTemplate <- gsub(".FUNCTIONAL.", gsub(".nii(.gz)*$", "", mrfiles[r]), thisTemplate, fixed=TRUE)
    thisTemplate <- gsub(".CONFOUNDS.", motfile, thisTemplate, fixed=TRUE)
    thisTemplate <- gsub(".CLOCK_TIMES.", file.path(timingdir, paste0("run", runnum, "_clock_FSL3col.txt")), thisTemplate, fixed=TRUE)
    thisTemplate <- gsub(".FEEDBACK_TIMES.", file.path(timingdir, paste0("run", runnum, "_feedback_FSL3col.txt")), thisTemplate, fixed=TRUE)
    thisTemplate <- gsub(".RELUNC_TIMES.", file.path(timingdir, paste0("run", runnum, "_rel_uncertainty_FSL3col.txt")), thisTemplate, fixed=TRUE)
    thisTemplate <- gsub(".MEANUNC_TIMES.", file.path(timingdir, paste0("run", runnum, "_mean_uncertainty_FSL3col.txt")), thisTemplate, fixed=TRUE) #for nomeanunc model, this should have no effect
    thisTemplate <- gsub(".EV_TIMES.", file.path(timingdir, paste0("run", runnum, "_ev_FSL3col.txt")), thisTemplate, fixed=TRUE)
    thisTemplate <- gsub(".RPEPOS_TIMES.", file.path(timingdir, paste0("run", runnum, "_rpe_pos_FSL3col.txt")), thisTemplate, fixed=TRUE)
    thisTemplate <- gsub(".RPENEG_TIMES.", file.path(timingdir, paste0("run", runnum, "_rpe_neg_FSL3col.txt")), thisTemplate, fixed=TRUE)
    
    featFile <- file.path(fsldir, paste0("FEAT_LVL1_run", runnum, ".fsf"))
    if (file.exists(featFile) && force==FALSE) { next } #skip re-creation of FSF and do not run below unless force==TRUE 
    cat(thisTemplate, file=featFile, sep="\n")
    
    allFeatFiles[[r]] <- featFile
  }
  
  if (run == TRUE) {    
    cl_fork <- makeForkCluster(nnodes=8)
    runfeat <- function(fsf) {
      runname <- basename(fsf)
      runFSLCommand(paste("feat", fsf), stdout=file.path(dirname(fsf), paste0("feat_stdout_", runname)), stderr=file.path(dirname(fsf), paste0("feat_stderr_", runname)))
    }
    clusterApply(cl_fork, allFeatFiles, runfeat)
    stopCluster(cl_fork)
  }
  
}
