fslValueModel <- function(f_value, mrfiles, runlengths, mrrunnums, run=FALSE, force=FALSE) {
  #sloppily using f_value, mrfiles, runlengths, etc. from parent env
  require(Rniftilib)
  require(parallel)
  fsfTemplate <- readLines(file.path(getMainDir(), "clock_analysis", "fmri", "feat_lvl1_clock_template.fsf"))
  
  fsldir <- file.path(normalizePath(file.path(dirname(mrfiles[1L]), "..")), "fsl_value") #note: normalizePath will fail to evaluate properly if directory does not exist

  if (file.exists(fsldir) && force==FALSE) { message(fsldir, " exists. Skipping."); return(0) }
  cat("fsldir create: ", fsldir, "\n")
  dir.create(fsldir, showWarnings=FALSE) #one directory up from a given clock run
  timingdir <- file.path(fsldir, "run_timing_deltavalue")
  
  d_value <- f_value$build_design_matrix(regressors=c("clock", "feedback", "ev", "rpe_neg", "rpe_pos"), 
      event_onsets=c("clock_onset", "feedback_onset", "clock_onset", "feedback_onset", "feedback_onset"), 
      durations=c(0, 0, "clock_duration", "feedback_duration", "feedback_duration"), 
      baselineCoefOrder=2, writeTimingFiles=c("FSL"),
      runVolumes=runlengths, runsToOutput=mrrunnums, output_directory=timingdir)
  
  allFeatFiles <- list()
  
  #FSL computes first-level models on individual runs
  for (r in 1:length(mrfiles)) {
    stopifnot(file.exists(file.path(dirname(mrfiles[r]), "motion.par"))) #can't find motion parameters
    
    runnum <- sub("^.*/clock(\\d+)$", "\\1", dirname(mrfiles[r]), perl=TRUE)
    nvol <- nifti.image.read(mrfiles[r], read_data=0)$dim[4L]
    
    #for motion parameter regression, add temporal derivatives, then PCA and retain first 3 components (~95% of variance)
    mot <- read.table(file.path(dirname(mrfiles[r]), "motion.par"), col.names=c("r.x", "r.y", "r.z", "t.x", "t.y", "t.z"))
    motderiv <- as.data.frame(lapply(mot, function(col) { c(0, diff(col)) }))
    names(motderiv) <- paste0("d.", names(motderiv)) #add delta to names
    motall <- cbind(mot, motderiv)
    pc <- princomp(motall, scores=TRUE)
    cumvar <- cumsum(pc$sdev^2/sum(pc$sdev^2))
    #cat("first three motion principal components account for: ", plyr::round_any(cumvar[3], .001), "\n")
    cat("first two motion principal components account for: ", plyr::round_any(cumvar[2], .001), "\n")
    mregressors <- pc$scores[,1:2] #first two components (cf Churchill et al. 2012 PLoS ONE)
    
    #Add volumes to censor here. Use censor_intersection.mat, which flags fd > 0.9 and DVARS > 20
    censorfile <- file.path(dirname(mrfiles[r]), "motion_info", "censor_intersection.mat")
    if (file.exists(censorfile) && file.info(censorfile)$size > 0) {
      censor <- read.table(censorfile, header=FALSE)
      mregressors <- cbind(mregressors, censor)
    }
    
    motfile <- file.path(fsldir, paste0("run", runnum, "_confounds.txt"))
    write.table(mregressors[1:runlengths[r], ], file=motfile, col.names=FALSE, row.names=FALSE) #make sure confound output matches fMRI file in length
    
    #add deep ventricle time series as confound regressor?
    
    #search and replace within fsf file for appropriate sections
    #.OUTPUTDIR. is the feat output location
    #.NVOL. is the number of volumes in the run
    #.FUNCTIONAL. is the fmri data to process (sans extension)
    #.CONFOUNDS. is the confounds file for GLM
    #.CLOCK_TIMES. is the three-column file for clock onset
    #.FEEDBACK_TIMES. is the three-column file for feedback onset
    #.EV_TIMES. is the three-column file for expected value
    #.RPEPOS_TIMES. is the three-column file for positive rpes
    #.RPENEG_TIMES. is the three-column file for negative rpes
    
    thisTemplate <- fsfTemplate
    thisTemplate <- gsub(".OUTPUTDIR.", file.path(fsldir, paste0("FEAT_LVL1_run", runnum)), thisTemplate, fixed=TRUE)
    thisTemplate <- gsub(".NVOL.", nvol, thisTemplate, fixed=TRUE)
    thisTemplate <- gsub(".FUNCTIONAL.", gsub(".nii(.gz)*$", "", mrfiles[r]), thisTemplate, fixed=TRUE)
    thisTemplate <- gsub(".CONFOUNDS.", motfile, thisTemplate, fixed=TRUE)
    thisTemplate <- gsub(".CLOCK_TIMES.", file.path(timingdir, paste0("run", runnum, "_clock_FSL3col.txt")), thisTemplate, fixed=TRUE)
    thisTemplate <- gsub(".FEEDBACK_TIMES.", file.path(timingdir, paste0("run", runnum, "_feedback_FSL3col.txt")), thisTemplate, fixed=TRUE)
    thisTemplate <- gsub(".EV_TIMES.", file.path(timingdir, paste0("run", runnum, "_ev_FSL3col.txt")), thisTemplate, fixed=TRUE)
    thisTemplate <- gsub(".RPEPOS_TIMES.", file.path(timingdir, paste0("run", runnum, "_rpe_pos_FSL3col.txt")), thisTemplate, fixed=TRUE)
    thisTemplate <- gsub(".RPENEG_TIMES.", file.path(timingdir, paste0("run", runnum, "_rpe_neg_FSL3col.txt")), thisTemplate, fixed=TRUE)
    
    featFile <- file.path(fsldir, paste0("FEAT_LVL1_run", runnum, ".fsf"))
    if (file.exists(featFile) && force==FALSE) { next } #skip re-creation of FSF and do not run below unless force==TRUE 
    cat(thisTemplate, file=featFile, sep="\n")
    
    allFeatFiles[[r]] <- featFile
  }
  
  if (run == TRUE) {
    #N.B.: If FSLDIR is not setup properly, running feat using a system call hangs, whereas running it in bash -i is okay
    #cat("#!/bin/bash",
    #    paste0(file.path(Sys.getenv("FSLDIR"), "bin", "feat "), allFeatFiles, 
    #        " > ", sapply(allFeatFiles, dirname), "/feat_stdout_", sapply(allFeatFiles, "basename"),
    #        " 2> ", sapply(allFeatFiles, dirname), "/feat_stderr_", sapply(allFeatFiles, "basename"),
    #        " &"),
    #    "wait",
    #    file="runfeat.bash", sep="\n")
    #system("bash -i runfeat.bash")
    
    cl_fork <- makeForkCluster(nnodes=16)
    runfeat <- function(fsf) {
      runname <- basename(fsf)
      runFSLCommand(paste("feat", fsf), stdout=file.path(dirname(fsf), paste0("feat_stdout_", runname)), stderr=file.path(dirname(fsf), paste0("feat_stderr_", runname)))
    }
    clusterApply(cl_fork, allFeatFiles, runfeat)
    stopCluster(cl_fork)
  }
  
}
