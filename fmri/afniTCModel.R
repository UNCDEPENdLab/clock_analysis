#function to setup and run afni time clock model
afniTCModel <- function(f, mrfiles, runlengths, mrrunnums, run=FALSE, force=FALSE) {
  afnidir <- file.path(normalizePath(file.path(dirname(mrfiles[1L]), "..")), "afni_tc") #note: normalizePath will fail to evaluate properly if directory does not exist (e.g., dir not created yet)
  if (file.exists(afnidir) && force==FALSE) {
      message("afni_tc directory already exists. ", afnidir, ". Skipping subject")
      return(NULL)
  }
              
  cat("afnidir create: ", afnidir, "\n")
  dir.create(afnidir, showWarnings=FALSE) #one directory up from a given clock run
  timingdir <- file.path(afnidir, "run_timing_tc")

  d <- build_design_matrix(fitobj=f, regressors=c("clock", "mean_uncertainty", "rel_uncertainty", "rpe_pos", "rpe_neg"), 
                             event_onsets=c("clock_onset", "clock_onset", "clock_onset", "feedback_onset", "feedback_onset"), 
                             durations=c("rt", "rt", "rt", "feedback_duration", "feedback_duration"), 
                             baselineCoefOrder=2, writeTimingFiles="AFNI",
                             runVolumes=runlengths, runsToOutput=mrrunnums, output_directory=timingdir)

  ##concat motion parameters
  motion_concat <- do.call(rbind, lapply(1:length(mrfiles), function(i)  {
            mot <- read.table(file.path(dirname(mrfiles[i]), "motion.par"), col.names=c("r.x", "r.y", "r.z", "t.x", "t.y", "t.z"))
            mot <- mot[1:runlengths[i],]
            motderiv <- as.data.frame(lapply(mot, function(col) { c(0, diff(col)) }))
            names(motderiv) <- paste0("d.", names(mot)) #add delta to names
            cbind(mot, motderiv)
          }))
  
  pc <- princomp(motion_concat, scores=TRUE)
  cumvar <- cumsum(pc$sdev^2/sum(pc$sdev^2))
  
  cat("first three motion principal components account for: ", plyr::round_any(cumvar[3], .001), "\n")
  mregressors <- pc$scores[,1:3] #first three components (cf Churchill et al. 2012 PLoS ONE)
  
  cat("correlation of motion parameters:\n\n")
  print(round(cor(motion_concat), 2))
  
  write.table(mregressors, file=file.path(afnidir, 'motion_concat.par'), col.names=FALSE, row.names=FALSE)
  
  #concat motion censor
  #continue to use censor_union.1D, which will flag either FD > 0.9 OR DVARS > 20.
  censor_union_concat <- do.call(c, lapply(1:length(mrfiles), function(i)  {
            cen <- read.table(file.path(dirname(mrfiles[i]), "motion_info", "censor_union.1D"))$V1
            cen <- cen[1:runlengths[i]]
            cen
          }))
  
  write.table(censor_union_concat, file=file.path(afnidir, "censor_union_concat.1D"), col.names=FALSE, row.names=FALSE)
  
  #make between-session regressors based on emotion condition.
  fmriDir <- "/Volumes/Serena/MMClock/MR_Raw"
  fitDir <- file.path(getMainDir(), "clock_analysis", "fmri", "fmri_fits")
  #/Volumes/Serena/MMClock/MR_Raw/10997_20140308/MBclock_recon/clock1/nfswudktm_clock1_5_trunc282.nii.gz
  subid <- factor(sub(paste0(fmriDir, "/([0-9]{5})_\\d+/MBclock_recon/.*$"), "\\1", mrfiles[1L], perl=TRUE))
  
  loc <- local({load(file.path(fitDir, paste0(as.character(subid), "_fitinfo.RData"))); environment()})$f #time-clock fit object (load as local var)
  emocon <- data.frame(emotion=loc$run_condition[mrrunnums], contingency=loc$rew_function[mrrunnums]) #vector of emotion and contingency
  
  #generate interactions for ev, rpe_neg, and rpe_pos with emotion
  csum <- cumsum(runlengths)

  for (reg in c("mean_uncertainty", "rel_uncertainty", "rpe_neg", "rpe_pos")) {
    for (emo in c("fear", "scram", "happy")) {
      vec <- rep(0, sum(runlengths))
      rmatch <- which(emocon$emotion == emo)
      for (r in rmatch) {
        f <- read.table(file.path(timingdir, paste0("run", r, "_", reg, ".1D")))$V1
        start <- ifelse(r > 1, csum[r-1] + 1, 1)
        end <- csum[r] #if (r < nrow(emocon)) csum[r] else 
        vec[start:end] <- f
      }
      write.table(vec, file=file.path(timingdir, paste0(reg, "_", emo, "_concat.1D")), row.names=FALSE, col.names=FALSE)
    }
  }

  ##generate mask of mrfiles where temporal min is > 0 for all runs
  for (f in 1:length(mrfiles)) {
      runFSLCommand(paste0("fslmaths ", mrfiles[f], " -Tmin -bin ", afnidir, "/tmin", f))
  }

  #sum mins together over runs and threshold at number of runs
  runFSLCommand(paste0("fslmaths ", paste(paste0(afnidir, "/tmin", 1:length(mrfiles)), collapse=" -add "), " ", afnidir, "/tminsum"))
  runFSLCommand(paste0("fslmaths ", afnidir, "/tminsum -thr ", length(mrfiles), " -bin ", afnidir, "/runmask"))
  runFSLCommand(paste0("imrm ", afnidir, "/tmin*")) #cleanup 

  modelname <- "glm_hrf_clock_preconvolve_tcModel_emoint"
  afniScript <- c('#!/bin/bash',
      'source /Users/michael/.bashrc')
  
  afniScript <- c(afniScript, paste(
      '3dDeconvolve',
      '-overwrite',
      '-input',
      paste(mrfiles, collapse=" \\\n"),
      '-tout -allzero_OK -polort 3 -jobs 8',
      '-num_stimts 16 -mask runmask.nii.gz', 
      paste('-stim_file 1', file.path(timingdir, "clock_concat.1D"), '-stim_label 1 clock_onset'),
      paste('-stim_file 2', file.path(timingdir, "mean_uncertainty_concat.1D"), '-stim_label 2 mean_uncertainty'),
      paste('-stim_file 3', file.path(timingdir, "mean_uncertainty_happy_concat.1D"), '-stim_label 3 mean_uncertainty_happy'),
      paste('-stim_file 4', file.path(timingdir, "mean_uncertainty_fear_concat.1D"), '-stim_label 4 mean_uncertainty_fear'),
      paste('-stim_file 5', file.path(timingdir, "rel_uncertainty_concat.1D"), '-stim_label 5 rel_uncertainty'),
      paste('-stim_file 6', file.path(timingdir, "rel_uncertainty_happy_concat.1D"), '-stim_label 6 rel_uncertainty_happy'),
      paste('-stim_file 7', file.path(timingdir, "rel_uncertainty_fear_concat.1D"), '-stim_label 7 rel_uncertainty_fear'),
      paste('-stim_file 8', file.path(timingdir, "rpe_neg_concat.1D"), '-stim_label 8 rpe_neg'),
      paste('-stim_file 9', file.path(timingdir, "rpe_neg_happy_concat.1D"), '-stim_label 9 rpe_neg_happy'),
      paste('-stim_file 10', file.path(timingdir, "rpe_neg_fear_concat.1D"), '-stim_label 10 rpe_neg_fear'),
      paste('-stim_file 11', file.path(timingdir, "rpe_pos_concat.1D"), '-stim_label 11 rpe_pos'),
      paste('-stim_file 12', file.path(timingdir, "rpe_pos_happy_concat.1D"), '-stim_label 12 rpe_pos_happy'),
      paste('-stim_file 13', file.path(timingdir, "rpe_pos_fear_concat.1D"), '-stim_label 13 rpe_pos_fear'),
      '-stim_file 14 "motion_concat.par[0]" -stim_label 14 motpc1 -stim_base 14',
      '-stim_file 15 "motion_concat.par[1]" -stim_label 15 motpc2 -stim_base 15',
      '-stim_file 16 "motion_concat.par[2]" -stim_label 16 motpc3 -stim_base 16',      
      '-num_glt 28',
      "-gltsym 'SYM: +1*rel_uncertainty' -glt_label 1 m_rel_uncertainty_scram",
      "-gltsym 'SYM: +1*rel_uncertainty +1*rel_uncertainty_happy' -glt_label 2 m_rel_uncertainty_happy",
      "-gltsym 'SYM: +1*rel_uncertainty +1*rel_uncertainty_fear' -glt_label 3 m_rel_uncertainty_fear",
      "-gltsym 'SYM: +1*rel_uncertainty +0.25*rel_uncertainty_fear +0.25*rel_uncertainty_happy' -glt_label 4 m_rel_uncertainty_overall",
      "-gltsym 'SYM: +1*rel_uncertainty_fear' -glt_label 5 rel_uncertainty_fear_gt_scram",
      "-gltsym 'SYM: +1*rel_uncertainty_happy' -glt_label 6 rel_uncertainty_happy_gt_scram",
      "-gltsym 'SYM: +1*rel_uncertainty_fear -1*rel_uncertainty_happy' -glt_label 7 rel_uncertainty_fear_gt_happy",
      "-gltsym 'SYM: +1*mean_uncertainty' -glt_label 8 m_mean_uncertainty_scram",
      "-gltsym 'SYM: +1*mean_uncertainty +1*mean_uncertainty_happy' -glt_label 9 m_mean_uncertainty_happy",
      "-gltsym 'SYM: +1*mean_uncertainty +1*mean_uncertainty_fear' -glt_label 10 m_mean_uncertainty_fear",
      "-gltsym 'SYM: +1*mean_uncertainty +0.25*mean_uncertainty_fear +0.25*mean_uncertainty_happy' -glt_label 11 m_mean_uncertainty_overall",
      "-gltsym 'SYM: +1*mean_uncertainty_fear' -glt_label 12 mean_uncertainty_fear_gt_scram",
      "-gltsym 'SYM: +1*mean_uncertainty_happy' -glt_label 13 mean_uncertainty_happy_gt_scram",
      "-gltsym 'SYM: +1*mean_uncertainty_fear -1*mean_uncertainty_happy' -glt_label 14 mean_uncertainty_fear_gt_happy",
      "-gltsym 'SYM: +1*rpe_neg' -glt_label 15 m_rpe_neg_scram",
      "-gltsym 'SYM: +1*rpe_neg +1*rpe_neg_happy' -glt_label 16 m_rpe_neg_happy",
      "-gltsym 'SYM: +1*rpe_neg +1*rpe_neg_fear' -glt_label 17 m_rpe_neg_fear",
      "-gltsym 'SYM: +1*rpe_neg +0.25*rpe_neg_fear +0.25*rpe_neg_happy' -glt_label 18 m_rpe_neg_overall",
      "-gltsym 'SYM: +1*rpe_neg_fear' -glt_label 19 rpe_neg_fear_gt_scram",
      "-gltsym 'SYM: +1*rpe_neg_happy' -glt_label 20 rpe_neg_happy_gt_scram",
      "-gltsym 'SYM: +1*rpe_neg_fear -1*rpe_neg_happy' -glt_label 21 rpe_neg_fear_gt_happy",
      "-gltsym 'SYM: +1*rpe_pos' -glt_label 22 m_rpe_pos_scram",
      "-gltsym 'SYM: +1*rpe_pos +1*rpe_pos_happy' -glt_label 23 m_rpe_pos_happy",
      "-gltsym 'SYM: +1*rpe_pos +1*rpe_pos_fear' -glt_label 24 m_rpe_pos_fear",
      "-gltsym 'SYM: +1*rpe_pos +0.25*rpe_pos_fear +0.25*rpe_pos_happy' -glt_label 25 m_rpe_pos_overall",
      "-gltsym 'SYM: +1*rpe_pos_fear' -glt_label 26 rpe_pos_fear_gt_scram",
      "-gltsym 'SYM: +1*rpe_pos_happy' -glt_label 27 rpe_pos_happy_gt_scram",
      "-gltsym 'SYM: +1*rpe_pos_fear -1*rpe_pos_happy' -glt_label 28 rpe_pos_fear_gt_happy",
      '-censor censor_union_concat.1D',
      paste0('-bucket ', modelname, '_stats'),
      paste0('-fitts ', modelname, '_fitts'),
      paste0('-errts ', modelname, '_errts'),
      paste0('-x1D ', modelname, '_x1D'),
      paste0('-cbucket ', modelname, '_coefs'),
      paste0('-xjpeg ', modelname, '_design.png'),
      '-GOFORIT 10 2>&1', 
      sep=" \\\n"))
  
  #add 3dREMLfit call
  afniScript <- c(afniScript,
      '',
      '#now rerun whole thing using REMLfit\n',
      'chmod +x glm_hrf_clock_preconvolve_tcModel_emoint_stats.REML_cmd\n',
      './glm_hrf_clock_preconvolve_tcModel_emoint_stats.REML_cmd 2>&1 | tee reml_glm_hrf_clock_preconvolve_tcModel_emoint_stats.log')
  
  writeLines(afniScript, con=file.path(afnidir, "3ddecon_preconvolve_tcModel.bash"))

  if (run==TRUE) {
      oldwd <- getwd()
      setwd(afnidir)
      system2("bash", "3ddecon_preconvolve_tcModel.bash", stdout="3ddecon_stdout", stderr="3ddecon_stderr")
      setwd(oldwd)
  }
}
