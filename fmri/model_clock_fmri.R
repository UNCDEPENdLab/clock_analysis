# fit behavioral data for all participants who completed emo clock in scanner
# setup model-based fMRI GLM analysis for based on fitted data

library(fitclock)
library(Rniftilib) #has nice function for just reading header
source(file.path(getMainDir(), "Miscellaneous", "Global_Functions.R"))

fslValueModel <- function(mrfiles, fitobj, run=FALSE) {
  require(Rniftilib)
  require(parallel)
  fsfTemplate <- readLines(file.path(getMainDir(), "clock_analysis", "fmri", "feat_lvl1_clock_template.fsf"))
  
  fsldir <- file.path(normalizePath(file.path(dirname(mrfiles[1L]), "..")), "fsl_value") #note: normalizePath will fail to evaluate properly if directory does not exist
  cat("fsldir create: ", fsldir, "\n")
  dir.create(fsldir, showWarnings=FALSE) #one directory up from a given clock run
  timingdir <- file.path(fsldir, "run_timing_deltavalue")
  
  d_value <- fitobj$build_design_matrix(regressors=c("clock", "feedback", "ev", "rpe_neg", "rpe_pos"), 
      event_onsets=c("clock_onset", "feedback_onset", "clock_onset", "feedback_onset", "feedback_onset"), 
      durations=c(0, 0, "clock_duration", "feedback_duration", "feedback_duration"), 
      baselineCoefOrder=2, writeTimingFiles=c("FSL"),
      runVolumes=runlengths, output_directory=timingdir)
  
  allFeatFiles <- list()
  
  #FSL computes first-level models on individual runs
  for (r in mrfiles) {
    stopifnot(file.exists(file.path(dirname(r), "motion.par"))) #can't find motion parameters
    
    runnum <- sub("^.*/clock(\\d+)$", "\\1", dirname(r), perl=TRUE)
    nvol <- nifti.image.read(r, read_data=0)$dim[4L]
    
    #for motion parameter regression, add temporal derivatives, then PCA and retain first 3 components (~95% of variance)
    mot <- read.table(file.path(dirname(r), "motion.par"), col.names=c("r.x", "r.y", "r.z", "t.x", "t.y", "t.z"))
    motderiv <- as.data.frame(lapply(mot, function(col) { c(0, diff(col)) }))
    names(motderiv) <- paste0("d.", names(motderiv)) #add delta to names
    motall <- cbind(mot, motderiv)
    pc <- princomp(motall, scores=TRUE)
    cumvar <- cumsum(pc$sdev^2/sum(pc$sdev^2))
    cat("first three motion principal components account for: ", plyr::round_any(cumvar[3], .001), "\n")
    mregressors <- pc$scores[,1:3] #first three components
    
    #add volumes to censor here... with wavelet despiking, this becomes less clear
    motfile <- file.path(fsldir, paste0("run", runnum, "_mot_regress_3pc.txt"))
    write.table(mregressors, file=motfile, col.names=FALSE, row.names=FALSE)
    
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
    thisTemplate <- gsub(".FUNCTIONAL.", gsub(".nii(.gz)*$", "", r), thisTemplate, fixed=TRUE)
    thisTemplate <- gsub(".CONFOUNDS.", motfile, thisTemplate, fixed=TRUE)
    thisTemplate <- gsub(".CLOCK_TIMES.", file.path(timingdir, paste0("run", runnum, "_clock_FSL3col.txt")), thisTemplate, fixed=TRUE)
    thisTemplate <- gsub(".FEEDBACK_TIMES.", file.path(timingdir, paste0("run", runnum, "_feedback_FSL3col.txt")), thisTemplate, fixed=TRUE)
    thisTemplate <- gsub(".EV_TIMES.", file.path(timingdir, paste0("run", runnum, "_ev_FSL3col.txt")), thisTemplate, fixed=TRUE)
    thisTemplate <- gsub(".RPEPOS_TIMES.", file.path(timingdir, paste0("run", runnum, "_rpe_pos_FSL3col.txt")), thisTemplate, fixed=TRUE)
    thisTemplate <- gsub(".RPENEG_TIMES.", file.path(timingdir, paste0("run", runnum, "_rpe_neg_FSL3col.txt")), thisTemplate, fixed=TRUE)
    
    featFile <- file.path(fsldir, paste0("FEAT_LVL1_run", runnum, ".fsf"))
    cat(thisTemplate, file=featFile, sep="\n")
    
    allFeatFiles[[r]] <- featFile
  }
  
  if (run == TRUE) {
    #N.B.: for some reason, running feat using a system call hangs, whereas running it in bash is okay
    #I think it has something to do with trying to open the report in the browser at initialization, but I'm not sure...
    #looks like it has to do with system2 running a non-interactive shell
    #force interactive run using bash -i
    #for now, generate a bash script to run it.
#    cat("#!/bin/bash",
#        paste0(file.path(Sys.getenv("FSLDIR"), "bin", "feat "), allFeatFiles, 
#            " > ", sapply(allFeatFiles, dirname), "/feat_stdout_", sapply(allFeatFiles, "basename"),
#            " 2> ", sapply(allFeatFiles, dirname), "/feat_stderr_", sapply(allFeatFiles, "basename"),
#            " &"),
#        "wait",
#        file="runfeat.bash", sep="\n")

    #browser() 
    #system("bash -i runfeat.bash")
    #system2("bash", "runfeat.bash")
  
    cl_fork <- makeForkCluster(nnodes=8)
    runfeat <- function(fsf) {
      runname <- basename(fsf)
      #system2(file.path(Sys.getenv("FSLDIR"), "bin", "feat"), fsf, stdout=file.path(dirname(fsf), paste0("feat_stdout_", runname)), stderr=file.path(dirname(fsf), paste0("feat_stderr_", runname)))
      runFSLCommand(paste("feat", fsf), stdout=file.path(dirname(fsf), paste0("feat_stdout_", runname)), stderr=file.path(dirname(fsf), paste0("feat_stderr_", runname)))
    }
    #browser()
    clusterApply(cl_fork, allFeatFiles, runfeat)
    stopCluster(cl_fork)
  }
  
}

#function to setup and run afni value model
afniValueModel <- function(mrfiles, fitobj, run=FALSE) {
  #setup and run basic value model in AFNI
  afnidir <- file.path(normalizePath(file.path(dirname(mrfiles[1L]), "..")), "afni_value") #note: normalizePath will fail to evaluate properly if directory does not exist (e.g., dir not created yet)
  cat("afnidir create: ", afnidir, "\n")
  dir.create(afnidir, showWarnings=FALSE) #one directory up from a given clock run
  timingdir <- file.path(afnidir, "run_timing_deltavalue")
  
  d_value <- fitobj$build_design_matrix(regressors=c("clock", "feedback", "ev", "rpe_neg", "rpe_pos"), 
      event_onsets=c("clock_onset", "feedback_onset", "clock_onset", "feedback_onset", "feedback_onset"), 
      durations=c(0, 0, "clock_duration", "feedback_duration", "feedback_duration"), 
      baselineCoefOrder=2, writeTimingFiles=c("AFNI"),
      runVolumes=runlengths, output_directory=timingdir)
  
  #concat motion parameters
  motion_concat <- do.call(rbind, lapply(mrfiles, function(x)  {
            x <- read.table(file.path(dirname(x), "motion.par"), col.names=c("r.x", "r.y", "r.z", "t.x", "t.y", "t.z"))
            x
          }))
  
  write.table(motion_concat, file=file.path(afnidir, "motion_concat.par"), col.names=FALSE, row.names=FALSE)
  
  cat("correlation of motion parameters:\n\n")
  print(round(cor(motion_concat), 2))
  
  #concat motion censor
  censor_union_concat <- do.call(c, lapply(mrfiles, function(x)  {
            x <- read.table(file.path(dirname(x), "motion_info", "censor_union.1D"))$V1
            x
          }))
  
  write.table(censor_union_concat, file=file.path(afnidir, "censor_union_concat.1D"), col.names=FALSE, row.names=FALSE)
  
  afniScript <- paste(
      '#!/bin/bash',
      '3dDeconvolve',
      '-overwrite',
      '-input',
      paste(mrfiles, collapse=" \\\n"),
      '-rout -tout -allzero_OK -polort 3 -jobs 8',
      '-num_stimts 11', 
      paste('-stim_file 1', file.path(timingdir, "clock_concat.1D"), '-stim_label 1 clock_onset'),
      paste('-stim_file 2', file.path(timingdir, "feedback_concat.1D"), '-stim_label 2 feedback_onset'),
      paste('-stim_file 3', file.path(timingdir, "ev_concat.1D"), '-stim_label 3 ev'),
      paste('-stim_file 4', file.path(timingdir, "rpe_neg_concat.1D"), '-stim_label 4 rpe_neg'),
      paste('-stim_file 5', file.path(timingdir, "rpe_pos_concat.1D"), '-stim_label 5 rpe_pos'),
      '-stim_file 6 "motion_concat.par[0]" -stim_label 6 rx -stim_base 6',
      '-stim_file 7 "motion_concat.par[1]" -stim_label 7 ry -stim_base 7',
      '-stim_file 8 "motion_concat.par[2]" -stim_label 8 rz -stim_base 8',
      '-stim_file 9 "motion_concat.par[3]" -stim_label 9 tx -stim_base 9',
      '-stim_file 10 "motion_concat.par[4]" -stim_label 10 ty -stim_base 10',
      '-stim_file 11 "motion_concat.par[5]" -stim_label 11 tz -stim_base 11',
      '-censor censor_union_concat.1D',
      '-bucket glm_hrf_clock_preconvolve_valueModel_stats',
      '-fitts glm_hrf_clock_preconvolve_valueModel_fitts',
      '-errts glm_hrf_clock_preconvolve_valueModel_errts',
      '-x1D glm_hrf_clock_preconvolve_valueModel_x1D',
      '-cbucket glm_hrf_clock_preconvolve_valueModel_coefs',
      '-xjpeg glm_hrf_clock_preconvolve_valueModel_design.png',
      '-GOFORIT 10 2>&1', 
      sep=" \\\n")
  
  #add 3dREMLfit call
  afniScript <- c(afniScript,
      '',
      '#now rerun whole thing using REMLfit\n',
      'chmod +x glm_hrf_clock_preconvolve_valueModel_stats.REML_cmd\n',
      './glm_hrf_clock_preconvolve_valueModel_stats.REML_cmd 2>&1 | tee reml_glm_hrf_clock_preconvolve_valueModel_stats.$date.log')
  
  writeLines(afniScript, con=file.path(afnidir, "3ddecon_preconvolve_valueModel.bash"))
  
  if (run==TRUE) {
    oldwd <- getwd()
    setwd(afnidir)
    system2("bash", "3ddecon_preconvolve_valueModel.bash", stdout="3ddecon_stdout", stderr="3ddecon_stderr")
    setwd(oldwd)
  }
}

setwd(file.path(getMainDir(), "clock_analysis", "fmri"))
if (!file.exists("fmri_fits")) { dir.create("fmri_fits") }
setwd("fmri_fits")

#start with base Frank model
#force non-negative epsilon (no sticky choice)
posEps <- clock_model()
posEps$add_params(
    meanRT(max_value=4000),
    autocorrPrevRT(),
    goForGold(),
    go(),
    noGo(),
    meanSlowFast(),
    exploreBeta()
)

behavDir <- "/Volumes/rcn1/bea_res/Data/Tasks/EmoClockfMRI/Basic"
fmriDir <- "/Volumes/Serena/MMClock/MR_Raw"
behavFiles <- list.files(path=behavDir, pattern="*tcExport.csv", full.names=TRUE, recursive=TRUE)

for (b in behavFiles) {
  subid <- sub("^.*fMRIEmoClock_(\\d+)_tc_tcExport.csv$", "\\1", b, perl=TRUE)
  mrfiles <- c() #force clear of mr files over subjects to avoid potential persistence from one subject to the next
  
  ##identify corresponding fmri directory
  mrmatch <- grep(paste0(subid, "_\\d+"), list.files(fmriDir, full.names=TRUE), perl=TRUE, value=TRUE)
  if(length(mrmatch) != 1L) {
    warning("Unable to find fMRI directory for subid: ", subid)
    next
  }
  if (! file.exists(file.path(mrmatch, "MBclock_recon"))) {
    warning("Unable to find preprocessed data MBclock_recon for subid: ", subid)
    next
  }
  
  ##identify fmri run lengths (4th dimension)
  mrfiles <- list.files(mrmatch, pattern="nfswukdtm_clock.*_5.nii.gz", full.names=TRUE, recursive=TRUE)
  
  if (length(mrfiles) == 0L) {
    warning("Unable to find any preprocessed MB files in dir: ", mrmatch)
    next
  }
  
  runlengths <- unname(sapply(mrfiles, function(x) { nifti.image.read(x, read_data=0)$dim[4L] }))
  
  if (file.exists(paste0(subid, "_fitinfo.RData"))) { 
    cat("Fit data already present for: ", subid, "\n")
    load(paste0(subid, "_fitinfo.RData"))
  } else {
    s <- clockdata_subject(subject_ID=subid, dataset=b)
    
    cat("Fitting behavioral data for subject: ", subid, "\n")
    
    #set data for model fit
    posEps$set_data(s)
    
    incr_fit <- posEps$incremental_fit(njobs=7, plot=FALSE)
    
    png(file.path(paste0(subid, "_incrfit.png")), width=9, height=6, units="in", res=300)
    print(incr_fit$AICplot)
    dev.off()
    
    f <- posEps$fit(random_starts=5)
    
    #design matrix matching Badre et al. 2012 Neuron
    d <- f$build_design_matrix(regressors=c("mean_uncertainty", "rel_uncertainty", "rpe_pos", "rpe_neg", "rt"), 
        event_onsets=c("clock_onset", "clock_onset", "feedback_onset", "feedback_onset", "feedback_onset"), 
        durations=c("rt", "rt", "feedback_duration", "feedback_duration", 0), baselineCoefOrder=2, runVolumes=runlengths)
    
    #delta rule value model (simple)
    vm <- deltavalue_model(clock_data=s, alphaV=0.3, betaV=0.3) #N.B. This matches V matrix from full time-clock algorithm fit.
    f_value <- vm$fit() #estimate learning rate as a free parameter
    
    #Design with EV, clock onset, feedback_onset, PE+, PE-
    d_value <- f_value$build_design_matrix(regressors=c("clock", "feedback", "ev", "rpe_neg", "rpe_pos"), 
        event_onsets=c("clock_onset", "feedback_onset", "clock_onset", "feedback_onset", "feedback_onset"), 
        durations=c(0, 0, "clock_duration", "feedback_duration", "feedback_duration"), baselineCoefOrder=2, runVolumes=runlengths)
    
    save(f_value, d_value, f, d, s, incr_fit, mrfiles, file=paste0(subid, "_fitinfo.RData"))
  }
  
  #setup afni model
  afniValueModel(mrfiles, f_value, run=FALSE)
  fslValueModel(mrfiles, f_value, run=TRUE)
}