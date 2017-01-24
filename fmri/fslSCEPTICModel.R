fslSCEPTICModel <- function(sceptic_signals, clockdata_subj, mrfiles, runlengths, mrrunnums, run=FALSE, force=FALSE, dropVolumes=0, outdir=NULL, usepreconvolve=FALSE, ...) {
  require(Rniftilib)
  require(parallel)

  if (is.null(outdir)) {
      outdir=paste0("sceptic_", paste(names(sceptic_signals)), collapse="_")
      if (usepreconvolve) { outdir=paste(outdir, "preconvolve", sep="_") }
  }
  
  univariate <- FALSE
  if (length(sceptic_signals) == 1L) {
    ##single model-based regressor
    if (usepreconvolve) {
        fsfTemplate <- readLines(file.path(getMainDir(), "clock_analysis", "fmri", "feat_lvl1_clock_sceptic_univariate_preconvolve_template.fsf"))
    } else {
        fsfTemplate <- readLines(file.path(getMainDir(), "clock_analysis", "fmri", "feat_lvl1_clock_sceptic_univariate_template.fsf"))
    }
    univariate <- TRUE
  } else if (length(sceptic_signals) == 4L) {
    if (usepreconvolve) {
      fsfTemplate <- readLines(file.path(getMainDir(), "clock_analysis", "fmri", "feat_lvl1_clock_sceptic_4param_preconvolve_template.fsf"))
    } else {
      stop("not implemented yet")
    }
  } else {
    stop("not implemented yet")
  }
  
  fsldir <- file.path(normalizePath(file.path(dirname(mrfiles[1L]), "..")), outdir) #note: normalizePath will fail to evaluate properly if directory does not exist

  if (file.exists(fsldir) && force==FALSE) { message(fsldir, " exists. Skipping."); return(0) }
  cat("fsldir create: ", fsldir, "\n")
  dir.create(fsldir, showWarnings=FALSE) #one directory up from a given clock run
  timingdir <- file.path(fsldir, "run_timing_sceptic")

  #create a fit object and populate with timings so that we can build a design matrix.
  f <- clock_fit()
  f$populate_fit(clockdata_subj)
  
  #add each sceptic_signal into the fit object with prefix sceptic_ so that design matrix code knows how to handle it
  signals_to_model <- rep(NA_character_, length(sceptic_signals))
  onsets <- rep(NA_character_, length(sceptic_signals))
  durations <- rep(NA_character_, length(sceptic_signals))
  normalizations <- rep(NA_character_, length(sceptic_signals))
  for (v in 1:length(sceptic_signals)) {
    thisName <- names(sceptic_signals)[v]
    signals_to_model[v] <- paste0("sceptic_", thisName) #name of regressor in build_design_matrix
    #clock_fit objects assume that each signal is a list of vectors where each element is a run (permits uneven runlengths)
    f$sceptic[[thisName]] <- split(sceptic_signals[[v]], row(sceptic_signals[[v]]))
    if (thisName %in% c("vauc", "vchosen", "ventropy", "vmax", "vsd", "ventropy_decay_matlab", "ventropy_fixed_matlab")) {
      onsets[v] <- "clock_onset"; durations[v] <- "clock_duration"; normalizations[v] <- "evtmax_1.0"
    } else if (thisName %in% c("pemax", "peauc", "dauc", "dsd")) {
      onsets[v] <- "feedback_onset"; durations[v] <- "feedback_duration"; normalizations[v] <- "evtmax_1.0"
    }
  }

  ##Use writingTimingFiles = "convolved" to receive HRF convolved parametric regressors
  d <- build_design_matrix(fitobj=f, regressors=c("clock", "feedback", signals_to_model), 
      event_onsets=c("clock_onset", "feedback_onset", onsets),
      durations=c("clock_duration", "feedback_duration", durations), 
      normalizations=c("durmax_1.0", "durmax_1.0", normalizations),
      baselineCoefOrder=2, writeTimingFiles=c("convolved", "FSL"), center_values=TRUE, convolve_wi_run=TRUE,
      runVolumes=runlengths, runsToOutput=mrrunnums, output_directory=timingdir, dropVolumes=dropVolumes, ...)

  save(d, f, file=file.path(fsldir, "designmatrix.RData"))
  
  allFeatFiles <- list()
  
  #FSL computes first-level models on individual runs
  for (r in 1:length(mrfiles)) {
    stopifnot(file.exists(file.path(dirname(mrfiles[r]), "motion.par"))) #can't find motion parameters
    
    runnum <- sub("^.*/clock(\\d+)$", "\\1", dirname(mrfiles[r]), perl=TRUE)
    nvol <- nifti.image.read(mrfiles[r], read_data=0)$dim[4L]

    ##just PCA motion on the current run
    ##mregressors <- pca_motion(mrfiles[r], runlengths[r], motion_parfile="motion.par", numpcs=3, dropVolumes=dropVolumes)$motion_pcs_concat

    ##Add volumes to censor here. Use censor_intersection.mat, which flags fd > 0.9 and DVARS > 20
    ##15Jun2016: Switch to FD > 0.9mm censoring in general (moving away from wavelet)
    ##If fd_0.9.mat doesn't exist, it means no spike regressors were generated at this threshold
    ##Thus, do not include in the nuisance set. Also do not include PCA motion regressors
    ##censorfile <- file.path(dirname(mrfiles[r]), "motion_info", "censor_intersection.mat")
    ##if (file.exists(censorfile) && file.info(censorfile)$size > 0) {
    ##  censor <- read.table(censorfile, header=FALSE)$V1
    ##  censor <- censor[(1+dropVolumes):runlengths[r]]
    ##  mregressors <- cbind(mregressors, censor)
    ##}

    mregressors <- NULL #start with NULL
    
    censorfile <- file.path(dirname(mrfiles[r]), "motion_info", "fd_0.9.mat")
    if (file.exists(censorfile) && file.info(censorfile)$size > 0) {
      censor <- read.table(censorfile, header=FALSE)
      censor <- censor[(1+dropVolumes):runlengths[r],,drop=FALSE] #need no drop here in case there is just a single volume to censor
      #if the spikes fall outside of the rows selected above, we will obtain an all-zero column. remove these
      censor <- censor[,sapply(censor, sum) > 0,drop=FALSE]
      if (ncol(censor) == 0L) { censor <- NULL } #no volumes to censor within valid timepoints
      mregressors <- censor
    }
    
    ##add CSF and WM regressors (with their derivatives)
    nuisancefile <- file.path(dirname(mrfiles[r]), "nuisance_regressors.txt")
    if (file.exists(nuisancefile)) {
      nuisance <- read.table(nuisancefile, header=FALSE)
      nuisance <- nuisance[(1+dropVolumes):runlengths[r],,drop=FALSE]
      nuisance <- as.data.frame(lapply(nuisance, function(col) { col - mean(col) })) #demean
      #cat("about to cbind with nuisance\n")
      #print(str(mregressors))
      #print(str(nuisance))
      if (!is.null(mregressors)) { mregressors <- cbind(mregressors, nuisance) #note that in R 3.3.0, cbind with NULL or c() is no problem...
      } else { mregressors <- nuisance }
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
    ##.VNAME. is the signal name in a univariate model
    ##.V_TIMES. is the three-column file for the signal
    ##.V_CON. is the contrast name for the signal
    
    thisTemplate <- fsfTemplate
    thisTemplate <- gsub(".OUTPUTDIR.", file.path(fsldir, paste0("FEAT_LVL1_run", runnum)), thisTemplate, fixed=TRUE)
    thisTemplate <- gsub(".NVOL.", nvol, thisTemplate, fixed=TRUE)
    thisTemplate <- gsub(".FUNCTIONAL.", gsub(".nii(.gz)*$", "", mrfiles[r]), thisTemplate, fixed=TRUE)
    thisTemplate <- gsub(".CONFOUNDS.", motfile, thisTemplate, fixed=TRUE)
    if (usepreconvolve) {
      thisTemplate <- gsub(".CLOCK_TIMES.", file.path(timingdir, paste0("run", runnum, "_clock.1D")), thisTemplate, fixed=TRUE)
      thisTemplate <- gsub(".FEEDBACK_TIMES.", file.path(timingdir, paste0("run", runnum, "_feedback.1D")), thisTemplate, fixed=TRUE)
    } else {
      thisTemplate <- gsub(".CLOCK_TIMES.", file.path(timingdir, paste0("run", runnum, "_clock_FSL3col.txt")), thisTemplate, fixed=TRUE)
      thisTemplate <- gsub(".FEEDBACK_TIMES.", file.path(timingdir, paste0("run", runnum, "_feedback_FSL3col.txt")), thisTemplate, fixed=TRUE)
    }
    
    for (s in 1:length(signals_to_model)) {
      if (usepreconvolve) {
        thisTemplate <- gsub(".V_TIMES.", file.path(timingdir, paste0("run", runnum, "_", signals_to_model[s], ".1D")), thisTemplate, fixed=TRUE)
      } else {
        thisTemplate <- gsub(".V_TIMES.", file.path(timingdir, paste0("run", runnum, "_", signals_to_model[s], "_FSL3col.txt")), thisTemplate, fixed=TRUE)
      }
      thisTemplate <- gsub(".VNAME.", sub("sceptic_", "", signals_to_model[s]), thisTemplate, fixed=TRUE) #remove sceptic_ prefix for brevity
      thisTemplate <- gsub(".V_CON.", sub("sceptic_", "", signals_to_model[s]), thisTemplate, fixed=TRUE)      
    }
    
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
