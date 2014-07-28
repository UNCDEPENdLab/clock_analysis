# fit behavioral data for all participants who completed emo clock in scanner
# setup model-based fMRI GLM analysis for based on fitted data

library(fitclock)

behavDir <- "/Volumes/bea_res/Data/Tasks/EmoClockfMRI/Basic"
fmriDir <- "/Volumes/Serena/MMClock/MR_Raw"
behavFiles <- list.files(path=behavDir, pattern="*tcExport.csv", full.names=TRUE, recursive=TRUE)

#function to setup and run afni value model
afniValueModel <- function(mrfiles, timingdir, run=FALSE) {
  #setup and run basic value model in AFNI
  afnidir <- suppressWarnings(normalizePath(file.path(dirname(mrfiles[1L]), "../afni_value")))
  dir.create(afnidir, showWarnings=FALSE) #one directory up from a given clock run
  
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

setwd("/Users/michael/CogEmoFaceReward/analysis")
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

library(Rniftilib) #has nice function for just reading header

for (b in behavFiles) {
  subid <- sub("^.*fMRIEmoClock_(\\d+)_tc_tcExport.csv$", "\\1", b, perl=TRUE)
  mrfiles <- c() #force clear of mr files over subjects to avoid potential persistence from one subject to the next
  
  if (file.exists(file.path(subid, paste0(subid, "_fitinfo.RData")))) { 
    cat("Fit data already present for: ", subid, "\n")
    load(file.path(subid, paste0(subid, "_fitinfo.RData")))
  } else {
    #identify corresponding fmri directory
    mrmatch <- grep(paste0(subid, "_\\d+"), list.files(fmriDir, full.names=TRUE), perl=TRUE, value=TRUE)
    if(length(mrmatch) != 1L) {
      warning("Unable to find fMRI directory for subid: ", subid)
      next
    }
    if (! file.exists(file.path(mrmatch, "MBclock_recon"))) {
      warning("Unable to find preprocessed data MBclock_recon for subid: ", subid)
      next
    }
    
    #identify fmri run lengths (4th dimension)
    mrfiles <- list.files(mrmatch, pattern="nfswuktmd_clock.*_5.nii.gz", full.names=TRUE, recursive=TRUE)
    
    if (length(mrfiles) == 0L) {
      warning("Unable to find any preprocessed MB files in dir: ", mrmatch)
      next
    }
    
    runlengths <- unname(sapply(mrfiles, function(x) { nifti.image.read(x, read_data=0)$dim[4L] }))
    
    dir.create(subid, showWarnings=FALSE)
    
    s <- clockdata_subject(subject_ID=subid, dataset=b)
    
    cat("Fitting behavioral data for subject: ", subid, "\n")
    
    #set data for model fit
    posEps$set_data(s)
    
    incr_fit <- posEps$incremental_fit(njobs=7)
    
    png(file.path(subid, paste0(subid, "_incrfit.png")), width=9, height=6, units="in", res=300)
    print(incr_fit$AICplot)
    dev.off()
    
    f <- posEps$fit(random_starts=5)
    
    #design matrix matching Badre et al. 2012 Neuron
    d <- f$build_design_matrix(regressors=c("mean_uncertainty", "rel_uncertainty", "rpe_pos", "rpe_neg", "rt"), 
        event_onsets=c("clock_onset", "clock_onset", "feedback_onset", "feedback_onset", "feedback_onset"), 
        durations=c("rt", "rt", "feedback_duration", "feedback_duration", 0), baselineCoefOrder=2, writeTimingFiles="AFNI",
        runVolumes=runlengths, output_directory=file.path(subid, "run_timing"))
    
    
    #delta rule value model (simple)
    vm <- deltavalue_model(clock_data=s, alphaV=0.3, betaV=0.3) #N.B. This matches V matrix from full time-clock algorithm fit.
    f_value <- vm$fit() #estimate learning rate as a free parameter
    
    #Design with EV, clock onset, feedback_onset, PE+, PE-
    d_value <- f_value$build_design_matrix(regressors=c("clock", "feedback", "ev", "rpe_neg", "rpe_pos"), 
        event_onsets=c("clock_onset", "feedback_onset", "feedback_onset", "feedback_onset", "feedback_onset"), 
        durations=c(0, 0, "feedback_duration", "feedback_duration", "feedback_duration"), baselineCoefOrder=2, writeTimingFiles="AFNI",
        runVolumes=runlengths, output_directory=file.path(subid, "run_timing_deltavalue"))
    
    save(f_value, d_value, f, d, s, incr_fit, mrfiles, file=file.path(subid, paste0(subid, "_fitinfo.RData")))
  }
  
  #setup afni model
  afniValueModel(mrfiles, timingdir=file.path(getwd(), subid, "run_timing_deltavalue"), run=FALSE)
  
}