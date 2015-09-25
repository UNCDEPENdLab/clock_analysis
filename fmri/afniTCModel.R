#function to setup and run afni time clock model
afniTCModel <- function(f, mrfiles, runlengths, mrrunnums, run=FALSE, force=FALSE, dropVolumes=0, outdir="afni_tc") {
  afnidir <- file.path(normalizePath(file.path(dirname(mrfiles[1L]), "..")), outdir) #note: normalizePath will fail to evaluate properly if directory does not exist (e.g., dir not created yet)
#  if (file.exists(afnidir) && force==FALSE) {
#    message(outdir, " directory already exists. ", afnidir, ". Skipping subject")
#    return(NULL)
#  }
  
  cat("afnidir create: ", afnidir, "\n")
  dir.create(afnidir, showWarnings=FALSE) #one directory up from a given clock run
  timingdir <- file.path(afnidir, "run_timing_tc")
  
  #matches model 8 from model_test_10895.R, which showed the most reliable and sensible activations
  #make sure to include .01 Hz high pass filter of design to match fmri preprocessing 
  d <- build_design_matrix(fitobj=f, regressors=c("clock", "mean_uncertainty", "rel_uncertainty", "ev", "feedback", "rpe_pos", "rpe_neg"), 
      event_onsets=c("clock_onset", "clock_onset", "clock_onset", "clock_onset", "feedback_onset", "feedback_onset", "feedback_onset"), 
      durations=c("clock_duration", "clock_duration", "clock_duration", "clock_duration", 0, 0, 0), 
      normalizations=c("durmax_1.0", "evtmax_1.0", "evtmax_1.0", "evtmax_1.0", "durmax_1.0", "evtmax_1.0", "evtmax_1.0"),
      center_values=TRUE, convolve_wi_run=TRUE, baselineCoefOrder=2,  runVolumes=runlengths, runsToOutput=mrrunnums,
      output_directory=timingdir, dropVolumes = dropVolumes, writeTimingFiles="AFNI", high_pass=.01)
  
  save(d, file=file.path(afnidir, "afni_tc_design.RData"))
  
  #concatenate motion parameters
  motpcs <- pca_motion(mrfiles, runlengths, motion_parfile="motion.par", numpcs=3, dropVolumes=dropVolumes)
  
  #add CSF and WM regressors (with their derivatives)
  nuisance_concat <- lapply(1:length(mrfiles), function(f) {
        nuisancefile <- file.path(dirname(mrfiles[f]), "nuisance_regressors.txt")
        stopifnot(file.exists(nuisancefile))
        nuisance <- read.table(nuisancefile, header=FALSE)
        nuisance <- nuisance[(1+dropVolumes):runlengths[f],,drop=FALSE]
        nuisance <- as.data.frame(lapply(nuisance, function(col) { col - mean(col) })) #demean
        nuisance
      })
  
  #apply .01 Hz bandpass filter to nuisance regressors, too
  allnuisance <- apply(cbind(motpcs$motion_pcs_concat, do.call(rbind, nuisance_concat)), 2, function(x) {
        fir1Bandpass(x, TR=1.0, low=.01, high=1/2, plotFilter=FALSE, forward_reverse=TRUE, padx=1, detrend=1)
      })
  
  write.table(allnuisance, file=file.path(afnidir, 'motion_pcs.txt'), col.names=FALSE, row.names=FALSE)    
  
  #concat motion censor
  #use censor_intersection.1D, which will flag FD > 0.9 AND DVARS > 20.
  #despiking should tamp down most of the problematic volumes.
  censor_intersection_concat <- do.call(c, lapply(1:length(mrfiles), function(i)  {
            cen <- read.table(file.path(dirname(mrfiles[i]), "motion_info", "censor_intersection.1D"))$V1
            cen <- cen[(1+dropVolumes):runlengths[i]]
            cen
          }))
  
  write.table(censor_intersection_concat, file=file.path(afnidir, "censor_intersection_concat.1D"), col.names=FALSE, row.names=FALSE)
  
  gen_emo_interaction_regressors(mrfiles[1L], regressors=c("clock", "mean_uncertainty", "rel_uncertainty", "ev", "feedback", "rpe_neg", "rpe_pos"), emotions=c("fear","scram","happy"),
      timingdir, mrrunnums, runlengths, dropVolumes)
  
  generateRunMask(mrfiles, outdir=afnidir, outfile="runmask")
  
  modelname <- "glm_hrf_clock_preconvolve_tcModel_emoint_normalized"
  afniScript <- c('#!/bin/bash',
                  'source /Users/michael/.bashrc',
                  'test -r /opt/ni_tools/ni_path.bash && . /opt/ni_tools/ni_path.bash')
  
  afniScript <- c(afniScript, paste(
          '3dDeconvolve',
          '-overwrite',
          '-input',
          paste(mrfiles, collapse=" \\\n"),
          '-tout -allzero_OK -polort 2 -jobs 8',
          '-num_stimts 28 -mask runmask.nii.gz', 
          paste('-stim_file 1', file.path(timingdir, "clock_concat.1D"), '-stim_label 1 clock'),
          paste('-stim_file 2', file.path(timingdir, "clock_happy_concat.1D"), '-stim_label 2 clock_happy'),
          paste('-stim_file 3', file.path(timingdir, "clock_fear_concat.1D"), '-stim_label 3 clock_fear'),
          paste('-stim_file 4', file.path(timingdir, "feedback_concat.1D"), '-stim_label 4 feedback'),
          paste('-stim_file 5', file.path(timingdir, "feedback_happy_concat.1D"), '-stim_label 5 feedback_happy'),
          paste('-stim_file 6', file.path(timingdir, "feedback_fear_concat.1D"), '-stim_label 6 feedback_fear'),      
          paste('-stim_file 7', file.path(timingdir, "mean_uncertainty_concat.1D"), '-stim_label 7 mean_uncertainty'),
          paste('-stim_file 8', file.path(timingdir, "mean_uncertainty_happy_concat.1D"), '-stim_label 8 mean_uncertainty_happy'),
          paste('-stim_file 9', file.path(timingdir, "mean_uncertainty_fear_concat.1D"), '-stim_label 9 mean_uncertainty_fear'),
          paste('-stim_file 10', file.path(timingdir, "rel_uncertainty_concat.1D"), '-stim_label 10 rel_uncertainty'),
          paste('-stim_file 11', file.path(timingdir, "rel_uncertainty_happy_concat.1D"), '-stim_label 11 rel_uncertainty_happy'),
          paste('-stim_file 12', file.path(timingdir, "rel_uncertainty_fear_concat.1D"), '-stim_label 12 rel_uncertainty_fear'),
          paste('-stim_file 13', file.path(timingdir, "ev_concat.1D"), '-stim_label 13 ev'),
          paste('-stim_file 14', file.path(timingdir, "ev_happy_concat.1D"), '-stim_label 14 ev_happy'),
          paste('-stim_file 15', file.path(timingdir, "ev_fear_concat.1D"), '-stim_label 15 ev_fear'),
          paste('-stim_file 16', file.path(timingdir, "rpe_neg_concat.1D"), '-stim_label 16 rpe_neg'),
          paste('-stim_file 17', file.path(timingdir, "rpe_neg_happy_concat.1D"), '-stim_label 17 rpe_neg_happy'),
          paste('-stim_file 18', file.path(timingdir, "rpe_neg_fear_concat.1D"), '-stim_label 18 rpe_neg_fear'),
          paste('-stim_file 19', file.path(timingdir, "rpe_pos_concat.1D"), '-stim_label 19 rpe_pos'),
          paste('-stim_file 20', file.path(timingdir, "rpe_pos_happy_concat.1D"), '-stim_label 20 rpe_pos_happy'),
          paste('-stim_file 21', file.path(timingdir, "rpe_pos_fear_concat.1D"), '-stim_label 21 rpe_pos_fear'),
          '-stim_file 22 "motion_pcs.txt[0]" -stim_label 22 motpc1 -stim_base 22',
          '-stim_file 23 "motion_pcs.txt[1]" -stim_label 23 motpc2 -stim_base 23',
          '-stim_file 24 "motion_pcs.txt[2]" -stim_label 24 motpc3 -stim_base 24',
          '-stim_file 25 "motion_pcs.txt[3]" -stim_label 25 csf -stim_base 25',
          '-stim_file 26 "motion_pcs.txt[4]" -stim_label 26 dcsf -stim_base 26',
          '-stim_file 27 "motion_pcs.txt[5]" -stim_label 27 dwm -stim_base 27',
          '-stim_file 28 "motion_pcs.txt[6]" -stim_label 28 wm -stim_base 28',
          '-num_glt 49',
          "-gltsym 'SYM: +1*clock' -glt_label 1 m_clock_scram",
          "-gltsym 'SYM: +1*clock +1*clock_happy' -glt_label 2 m_clock_happy",
          "-gltsym 'SYM: +1*clock +1*clock_fear' -glt_label 3 m_clock_fear",
          "-gltsym 'SYM: +1*clock +0.25*clock_fear +0.25*clock_happy' -glt_label 4 m_clock_overall",
          "-gltsym 'SYM: +1*clock_fear' -glt_label 5 clock_fear_gt_scram",
          "-gltsym 'SYM: +1*clock_happy' -glt_label 6 clock_happy_gt_scram",
          "-gltsym 'SYM: +1*clock_fear -1*clock_happy' -glt_label 7 clock_fear_gt_happy",
          "-gltsym 'SYM: +1*feedback' -glt_label 8 m_feedback_scram",
          "-gltsym 'SYM: +1*feedback +1*feedback_happy' -glt_label 9 m_feedback_happy",
          "-gltsym 'SYM: +1*feedback +1*feedback_fear' -glt_label 10 m_feedback_fear",
          "-gltsym 'SYM: +1*feedback +0.25*feedback_fear +0.25*feedback_happy' -glt_label 11 m_feedback_overall",
          "-gltsym 'SYM: +1*feedback_fear' -glt_label 12 feedback_fear_gt_scram",
          "-gltsym 'SYM: +1*feedback_happy' -glt_label 13 feedback_happy_gt_scram",
          "-gltsym 'SYM: +1*feedback_fear -1*feedback_happy' -glt_label 14 feedback_fear_gt_happy",
          "-gltsym 'SYM: +1*mean_uncertainty' -glt_label 15 m_mean_uncertainty_scram",
          "-gltsym 'SYM: +1*mean_uncertainty +1*mean_uncertainty_happy' -glt_label 16 m_mean_uncertainty_happy",
          "-gltsym 'SYM: +1*mean_uncertainty +1*mean_uncertainty_fear' -glt_label 17 m_mean_uncertainty_fear",
          "-gltsym 'SYM: +1*mean_uncertainty +0.25*mean_uncertainty_fear +0.25*mean_uncertainty_happy' -glt_label 18 m_mean_uncertainty_overall",
          "-gltsym 'SYM: +1*mean_uncertainty_fear' -glt_label 19 mean_uncertainty_fear_gt_scram",
          "-gltsym 'SYM: +1*mean_uncertainty_happy' -glt_label 20 mean_uncertainty_happy_gt_scram",
          "-gltsym 'SYM: +1*mean_uncertainty_fear -1*mean_uncertainty_happy' -glt_label 21 mean_uncertainty_fear_gt_happy",
          "-gltsym 'SYM: +1*rel_uncertainty' -glt_label 22 m_rel_uncertainty_scram",
          "-gltsym 'SYM: +1*rel_uncertainty +1*rel_uncertainty_happy' -glt_label 23 m_rel_uncertainty_happy",
          "-gltsym 'SYM: +1*rel_uncertainty +1*rel_uncertainty_fear' -glt_label 24 m_rel_uncertainty_fear",
          "-gltsym 'SYM: +1*rel_uncertainty +0.25*rel_uncertainty_fear +0.25*rel_uncertainty_happy' -glt_label 25 m_rel_uncertainty_overall",
          "-gltsym 'SYM: +1*rel_uncertainty_fear' -glt_label 26 rel_uncertainty_fear_gt_scram",
          "-gltsym 'SYM: +1*rel_uncertainty_happy' -glt_label 27 rel_uncertainty_happy_gt_scram",
          "-gltsym 'SYM: +1*rel_uncertainty_fear -1*rel_uncertainty_happy' -glt_label 28 rel_uncertainty_fear_gt_happy",
          "-gltsym 'SYM: +1*ev' -glt_label 29 m_ev_scram",
          "-gltsym 'SYM: +1*ev +1*ev_happy' -glt_label 30 m_ev_happy",
          "-gltsym 'SYM: +1*ev +1*ev_fear' -glt_label 31 m_ev_fear",
          "-gltsym 'SYM: +1*ev +0.25*ev_fear +0.25*ev_happy' -glt_label 32 m_ev_overall",
          "-gltsym 'SYM: +1*ev_fear' -glt_label 33 ev_fear_gt_scram",
          "-gltsym 'SYM: +1*ev_happy' -glt_label 34 ev_happy_gt_scram",
          "-gltsym 'SYM: +1*ev_fear -1*ev_happy' -glt_label 35 ev_fear_gt_happy",
          "-gltsym 'SYM: +1*rpe_neg' -glt_label 36 m_rpe_neg_scram",
          "-gltsym 'SYM: +1*rpe_neg +1*rpe_neg_happy' -glt_label 37 m_rpe_neg_happy",
          "-gltsym 'SYM: +1*rpe_neg +1*rpe_neg_fear' -glt_label 38 m_rpe_neg_fear",
          "-gltsym 'SYM: +1*rpe_neg +0.25*rpe_neg_fear +0.25*rpe_neg_happy' -glt_label 39 m_rpe_neg_overall",
          "-gltsym 'SYM: +1*rpe_neg_fear' -glt_label 40 rpe_neg_fear_gt_scram",
          "-gltsym 'SYM: +1*rpe_neg_happy' -glt_label 41 rpe_neg_happy_gt_scram",
          "-gltsym 'SYM: +1*rpe_neg_fear -1*rpe_neg_happy' -glt_label 42 rpe_neg_fear_gt_happy",
          "-gltsym 'SYM: +1*rpe_pos' -glt_label 43 m_rpe_pos_scram",
          "-gltsym 'SYM: +1*rpe_pos +1*rpe_pos_happy' -glt_label 44 m_rpe_pos_happy",
          "-gltsym 'SYM: +1*rpe_pos +1*rpe_pos_fear' -glt_label 45 m_rpe_pos_fear",
          "-gltsym 'SYM: +1*rpe_pos +0.25*rpe_pos_fear +0.25*rpe_pos_happy' -glt_label 46 m_rpe_pos_overall",
          "-gltsym 'SYM: +1*rpe_pos_fear' -glt_label 47 rpe_pos_fear_gt_scram",
          "-gltsym 'SYM: +1*rpe_pos_happy' -glt_label 48 rpe_pos_happy_gt_scram",
          "-gltsym 'SYM: +1*rpe_pos_fear -1*rpe_pos_happy' -glt_label 49 rpe_pos_fear_gt_happy",
          '-censor censor_intersection_concat.1D',
          paste0('-bucket ', modelname, '_stats'),
          ##paste0('-fitts ', modelname, '_fitts'),
          ##paste0('-errts ', modelname, '_errts'),
          paste0('-x1D ', modelname, '_x1D'),
          paste0('-cbucket ', modelname, '_coefs'),
          paste0('-xjpeg ', modelname, '_design.png'),
          '-GOFORIT 10 2>&1', 
          sep=" \\\n"))
  
  #add 3dREMLfit call
  afniScript <- c(afniScript,
      '',
      '#now rerun whole thing using REMLfit\n',
      paste0('chmod +x ', modelname, '_stats.REML_cmd\n'),
      paste0('./', modelname, '_stats.REML_cmd 2>&1 | tee reml_', modelname, '_stats.log'))
  
  writeLines(afniScript, con=file.path(afnidir, paste0("3ddecon_", modelname, ".bash")))
  
  if (run==TRUE && (!file.exists(file.path(afnidir, paste0("3ddecon_", modelname, "_stdout"))) || 
        (file.exists(file.path(afnidir, paste0("3ddecon_", modelname, "_stdout"))) &&  force == TRUE)) ) {
    oldwd <- getwd()
    setwd(afnidir)
    system2("bash", paste0("3ddecon_", modelname, ".bash"), stdout=paste0("3ddecon_", modelname, "_stdout"), stderr=paste0("3ddecon_", modelname, "_stderr"))
    setwd(oldwd)
  }
  
  ##model that excludes mean uncertainty given problems with collinearity with rel uncertainty and ev
  modelname <- "glm_hrf_clock_preconvolve_tcModel_emoint_nomeanunc_normalized"
  afniScript <- c('#!/bin/bash',
      'source /Users/michael/.bashrc',
      'test -r /opt/ni_tools/ni_path.bash && . /opt/ni_tools/ni_path.bash')
  
  afniScript <- c(afniScript, paste(
          '3dDeconvolve',
          '-overwrite',
          '-input',
          paste(mrfiles, collapse=" \\\n"),
          '-tout -allzero_OK -polort 2 -jobs 8',
          '-num_stimts 25 -mask runmask.nii.gz', 
          paste('-stim_file 1', file.path(timingdir, "clock_concat.1D"), '-stim_label 1 clock'),
          paste('-stim_file 2', file.path(timingdir, "clock_happy_concat.1D"), '-stim_label 2 clock_happy'),
          paste('-stim_file 3', file.path(timingdir, "clock_fear_concat.1D"), '-stim_label 3 clock_fear'),
          paste('-stim_file 4', file.path(timingdir, "feedback_concat.1D"), '-stim_label 4 feedback'),
          paste('-stim_file 5', file.path(timingdir, "feedback_happy_concat.1D"), '-stim_label 5 feedback_happy'),
          paste('-stim_file 6', file.path(timingdir, "feedback_fear_concat.1D"), '-stim_label 6 feedback_fear'),      
          paste('-stim_file 7', file.path(timingdir, "rel_uncertainty_concat.1D"), '-stim_label 7 rel_uncertainty'),
          paste('-stim_file 8', file.path(timingdir, "rel_uncertainty_happy_concat.1D"), '-stim_label 8 rel_uncertainty_happy'),
          paste('-stim_file 9', file.path(timingdir, "rel_uncertainty_fear_concat.1D"), '-stim_label 9 rel_uncertainty_fear'),
          paste('-stim_file 10', file.path(timingdir, "ev_concat.1D"), '-stim_label 10 ev'),
          paste('-stim_file 11', file.path(timingdir, "ev_happy_concat.1D"), '-stim_label 11 ev_happy'),
          paste('-stim_file 12', file.path(timingdir, "ev_fear_concat.1D"), '-stim_label 12 ev_fear'),
          paste('-stim_file 13', file.path(timingdir, "rpe_neg_concat.1D"), '-stim_label 13 rpe_neg'),
          paste('-stim_file 14', file.path(timingdir, "rpe_neg_happy_concat.1D"), '-stim_label 14 rpe_neg_happy'),
          paste('-stim_file 15', file.path(timingdir, "rpe_neg_fear_concat.1D"), '-stim_label 15 rpe_neg_fear'),
          paste('-stim_file 16', file.path(timingdir, "rpe_pos_concat.1D"), '-stim_label 16 rpe_pos'),
          paste('-stim_file 17', file.path(timingdir, "rpe_pos_happy_concat.1D"), '-stim_label 17 rpe_pos_happy'),
          paste('-stim_file 18', file.path(timingdir, "rpe_pos_fear_concat.1D"), '-stim_label 18 rpe_pos_fear'),
          '-stim_file 19 "motion_pcs.txt[0]" -stim_label 19 motpc1 -stim_base 19',
          '-stim_file 20 "motion_pcs.txt[1]" -stim_label 20 motpc2 -stim_base 20',
          '-stim_file 21 "motion_pcs.txt[2]" -stim_label 21 motpc3 -stim_base 21',
          '-stim_file 22 "motion_pcs.txt[3]" -stim_label 22 csf -stim_base 22',
          '-stim_file 23 "motion_pcs.txt[4]" -stim_label 23 dcsf -stim_base 23',
          '-stim_file 24 "motion_pcs.txt[5]" -stim_label 24 dwm -stim_base 24',
          '-stim_file 25 "motion_pcs.txt[6]" -stim_label 25 wm -stim_base 25',
          '-num_glt 42',
          "-gltsym 'SYM: +1*clock' -glt_label 1 m_clock_scram",
          "-gltsym 'SYM: +1*clock +1*clock_happy' -glt_label 2 m_clock_happy",
          "-gltsym 'SYM: +1*clock +1*clock_fear' -glt_label 3 m_clock_fear",
          "-gltsym 'SYM: +1*clock +0.25*clock_fear +0.25*clock_happy' -glt_label 4 m_clock_overall",
          "-gltsym 'SYM: +1*clock_fear' -glt_label 5 clock_fear_gt_scram",
          "-gltsym 'SYM: +1*clock_happy' -glt_label 6 clock_happy_gt_scram",
          "-gltsym 'SYM: +1*clock_fear -1*clock_happy' -glt_label 7 clock_fear_gt_happy",
          "-gltsym 'SYM: +1*feedback' -glt_label 8 m_feedback_scram",
          "-gltsym 'SYM: +1*feedback +1*feedback_happy' -glt_label 9 m_feedback_happy",
          "-gltsym 'SYM: +1*feedback +1*feedback_fear' -glt_label 10 m_feedback_fear",
          "-gltsym 'SYM: +1*feedback +0.25*feedback_fear +0.25*feedback_happy' -glt_label 11 m_feedback_overall",
          "-gltsym 'SYM: +1*feedback_fear' -glt_label 12 feedback_fear_gt_scram",
          "-gltsym 'SYM: +1*feedback_happy' -glt_label 13 feedback_happy_gt_scram",
          "-gltsym 'SYM: +1*feedback_fear -1*feedback_happy' -glt_label 14 feedback_fear_gt_happy",
          "-gltsym 'SYM: +1*rel_uncertainty' -glt_label 15 m_rel_uncertainty_scram",
          "-gltsym 'SYM: +1*rel_uncertainty +1*rel_uncertainty_happy' -glt_label 16 m_rel_uncertainty_happy",
          "-gltsym 'SYM: +1*rel_uncertainty +1*rel_uncertainty_fear' -glt_label 17 m_rel_uncertainty_fear",
          "-gltsym 'SYM: +1*rel_uncertainty +0.25*rel_uncertainty_fear +0.25*rel_uncertainty_happy' -glt_label 18 m_rel_uncertainty_overall",
          "-gltsym 'SYM: +1*rel_uncertainty_fear' -glt_label 19 rel_uncertainty_fear_gt_scram",
          "-gltsym 'SYM: +1*rel_uncertainty_happy' -glt_label 20 rel_uncertainty_happy_gt_scram",
          "-gltsym 'SYM: +1*rel_uncertainty_fear -1*rel_uncertainty_happy' -glt_label 21 rel_uncertainty_fear_gt_happy",
          "-gltsym 'SYM: +1*ev' -glt_label 22 m_ev_scram",
          "-gltsym 'SYM: +1*ev +1*ev_happy' -glt_label 23 m_ev_happy",
          "-gltsym 'SYM: +1*ev +1*ev_fear' -glt_label 24 m_ev_fear",
          "-gltsym 'SYM: +1*ev +0.25*ev_fear +0.25*ev_happy' -glt_label 25 m_ev_overall",
          "-gltsym 'SYM: +1*ev_fear' -glt_label 26 ev_fear_gt_scram",
          "-gltsym 'SYM: +1*ev_happy' -glt_label 27 ev_happy_gt_scram",
          "-gltsym 'SYM: +1*ev_fear -1*ev_happy' -glt_label 28 ev_fear_gt_happy",
          "-gltsym 'SYM: +1*rpe_neg' -glt_label 29 m_rpe_neg_scram",
          "-gltsym 'SYM: +1*rpe_neg +1*rpe_neg_happy' -glt_label 30 m_rpe_neg_happy",
          "-gltsym 'SYM: +1*rpe_neg +1*rpe_neg_fear' -glt_label 31 m_rpe_neg_fear",
          "-gltsym 'SYM: +1*rpe_neg +0.25*rpe_neg_fear +0.25*rpe_neg_happy' -glt_label 32 m_rpe_neg_overall",
          "-gltsym 'SYM: +1*rpe_neg_fear' -glt_label 33 rpe_neg_fear_gt_scram",
          "-gltsym 'SYM: +1*rpe_neg_happy' -glt_label 34 rpe_neg_happy_gt_scram",
          "-gltsym 'SYM: +1*rpe_neg_fear -1*rpe_neg_happy' -glt_label 35 rpe_neg_fear_gt_happy",
          "-gltsym 'SYM: +1*rpe_pos' -glt_label 36 m_rpe_pos_scram",
          "-gltsym 'SYM: +1*rpe_pos +1*rpe_pos_happy' -glt_label 37 m_rpe_pos_happy",
          "-gltsym 'SYM: +1*rpe_pos +1*rpe_pos_fear' -glt_label 38 m_rpe_pos_fear",
          "-gltsym 'SYM: +1*rpe_pos +0.25*rpe_pos_fear +0.25*rpe_pos_happy' -glt_label 39 m_rpe_pos_overall",
          "-gltsym 'SYM: +1*rpe_pos_fear' -glt_label 40 rpe_pos_fear_gt_scram",
          "-gltsym 'SYM: +1*rpe_pos_happy' -glt_label 41 rpe_pos_happy_gt_scram",
          "-gltsym 'SYM: +1*rpe_pos_fear -1*rpe_pos_happy' -glt_label 42 rpe_pos_fear_gt_happy",
          '-censor censor_intersection_concat.1D',
          paste0('-bucket ', modelname, '_stats'),
          ##paste0('-fitts ', modelname, '_fitts'),
          ##paste0('-errts ', modelname, '_errts'),
          paste0('-x1D ', modelname, '_x1D'),
          paste0('-cbucket ', modelname, '_coefs'),
          paste0('-xjpeg ', modelname, '_design.png'),
          '-GOFORIT 10 2>&1', 
          sep=" \\\n"))
  
  #add 3dREMLfit call
  afniScript <- c(afniScript,
      '',
      '#now rerun whole thing using REMLfit\n',
      paste0('chmod +x ', modelname, '_stats.REML_cmd\n'),
      paste0('./', modelname, '_stats.REML_cmd 2>&1 | tee reml_', modelname, '_stats.log'))
  
  writeLines(afniScript, con=file.path(afnidir, paste0("3ddecon_", modelname, ".bash")))
  
  if (run==TRUE && (!file.exists(file.path(afnidir, paste0("3ddecon_", modelname, "_stdout"))) || 
        (file.exists(file.path(afnidir, paste0("3ddecon_", modelname, "_stdout"))) &&  force == TRUE)) ) {
    oldwd <- getwd()
    setwd(afnidir)
    system2("bash", paste0("3ddecon_", modelname, ".bash"), stdout=paste0("3ddecon_", modelname, "_stdout"), stderr=paste0("3ddecon_", modelname, "_stderr"))
    setwd(oldwd)
  }
  
  ##models that include individual parametric regressors (including trial-ness regressors for clock and feedback)
  parametric <- c("mean_uncertainty", "rel_uncertainty", "ev", "rpe_pos", "rpe_neg")
  
  for (p in parametric) {
    modelname <- paste0("glm_hrf_clock_preconvolve_", p, "_only_emoint_normalized")
    afniScript <- c('#!/bin/bash',
                    'source /Users/michael/.bashrc',
                    'test -r /opt/ni_tools/ni_path.bash && . /opt/ni_tools/ni_path.bash')
    
    afniScript <- c(afniScript, paste(
            '3dDeconvolve',
            '-overwrite',
            '-input',
            paste(mrfiles, collapse=" \\\n"),
            '-tout -allzero_OK -polort 2 -jobs 8',
            '-num_stimts 16 -mask runmask.nii.gz', 
            paste('-stim_file 1', file.path(timingdir, "clock_concat.1D"), '-stim_label 1 clock'),
            paste('-stim_file 2', file.path(timingdir, "clock_happy_concat.1D"), '-stim_label 2 clock_happy'),
            paste('-stim_file 3', file.path(timingdir, "clock_fear_concat.1D"), '-stim_label 3 clock_fear'),
            paste('-stim_file 4', file.path(timingdir, "feedback_concat.1D"), '-stim_label 4 feedback'),
            paste('-stim_file 5', file.path(timingdir, "feedback_happy_concat.1D"), '-stim_label 5 feedback_happy'),
            paste('-stim_file 6', file.path(timingdir, "feedback_fear_concat.1D"), '-stim_label 6 feedback_fear'),      
            paste('-stim_file 7', file.path(timingdir, paste0(p, "_concat.1D")), '-stim_label 7', p),
            paste('-stim_file 8', file.path(timingdir, paste0(p, "_happy_concat.1D")), '-stim_label 8', paste0(p, "_happy")),
            paste('-stim_file 9', file.path(timingdir, paste0(p, "_fear_concat.1D")), '-stim_label 9', paste0(p, "_fear")),
            '-stim_file 10 "motion_pcs.txt[0]" -stim_label 10 motpc1 -stim_base 10',
            '-stim_file 11 "motion_pcs.txt[1]" -stim_label 11 motpc2 -stim_base 11',
            '-stim_file 12 "motion_pcs.txt[2]" -stim_label 12 motpc3 -stim_base 12',
            '-stim_file 13 "motion_pcs.txt[3]" -stim_label 13 csf -stim_base 13',
            '-stim_file 14 "motion_pcs.txt[4]" -stim_label 14 dcsf -stim_base 14',
            '-stim_file 15 "motion_pcs.txt[5]" -stim_label 15 dwm -stim_base 15',
            '-stim_file 16 "motion_pcs.txt[6]" -stim_label 16 wm -stim_base 16',
            '-num_glt 21',
            "-gltsym 'SYM: +1*clock' -glt_label 1 m_clock_scram",
            "-gltsym 'SYM: +1*clock +1*clock_happy' -glt_label 2 m_clock_happy",
            "-gltsym 'SYM: +1*clock +1*clock_fear' -glt_label 3 m_clock_fear",
            "-gltsym 'SYM: +1*clock +0.25*clock_fear +0.25*clock_happy' -glt_label 4 m_clock_overall",
            "-gltsym 'SYM: +1*clock_fear' -glt_label 5 clock_fear_gt_scram",
            "-gltsym 'SYM: +1*clock_happy' -glt_label 6 clock_happy_gt_scram",
            "-gltsym 'SYM: +1*clock_fear -1*clock_happy' -glt_label 7 clock_fear_gt_happy",
            "-gltsym 'SYM: +1*feedback' -glt_label 8 m_feedback_scram",
            "-gltsym 'SYM: +1*feedback +1*feedback_happy' -glt_label 9 m_feedback_happy",
            "-gltsym 'SYM: +1*feedback +1*feedback_fear' -glt_label 10 m_feedback_fear",
            "-gltsym 'SYM: +1*feedback +0.25*feedback_fear +0.25*feedback_happy' -glt_label 11 m_feedback_overall",
            "-gltsym 'SYM: +1*feedback_fear' -glt_label 12 feedback_fear_gt_scram",
            "-gltsym 'SYM: +1*feedback_happy' -glt_label 13 feedback_happy_gt_scram",
            "-gltsym 'SYM: +1*feedback_fear -1*feedback_happy' -glt_label 14 feedback_fear_gt_happy",
            paste0("-gltsym 'SYM: +1*", p, "' -glt_label 15 m_", p, "_scram"),
            paste0("-gltsym 'SYM: +1*", p, " +1*", p, "_happy' -glt_label 16 m_", p, "_happy"),
            paste0("-gltsym 'SYM: +1*", p, " +1*", p, "_fear' -glt_label 17 m_", p, "_fear"),
            paste0("-gltsym 'SYM: +1*", p, " +0.25*", p, "_fear +0.25*", p, "_happy' -glt_label 18 m_", p, "_overall"),
            paste0("-gltsym 'SYM: +1*", p, "_fear' -glt_label 19 ", p, "_fear_gt_scram"),
            paste0("-gltsym 'SYM: +1*", p, "_happy' -glt_label 20 ", p, "_happy_gt_scram"),
            paste0("-gltsym 'SYM: +1*", p, "_fear -1*", p, "_happy' -glt_label 21 ", p, "_fear_gt_happy"),
            '-censor censor_intersection_concat.1D',
            paste0('-bucket ', modelname, '_stats'),
            ##paste0('-fitts ', modelname, '_fitts'),
            ##paste0('-errts ', modelname, '_errts'),
            paste0('-x1D ', modelname, '_x1D'),
            paste0('-cbucket ', modelname, '_coefs'),
            paste0('-xjpeg ', modelname, '_design.png'),
            '-GOFORIT 10 2>&1', 
            sep=" \\\n"))
    
    #add 3dREMLfit call
    afniScript <- c(afniScript,
        '',
        '#now rerun whole thing using REMLfit\n',
        paste0('chmod +x ', modelname, '_stats.REML_cmd\n'),
        paste0('./', modelname, '_stats.REML_cmd 2>&1 | tee reml_', modelname, '_stats.log'))
    
    writeLines(afniScript, con=file.path(afnidir, paste0("3ddecon_", modelname, ".bash")))
    
    if (run==TRUE && (!file.exists(file.path(afnidir, paste0("3ddecon_", modelname, "_stdout"))) || 
          (file.exists(file.path(afnidir, paste0("3ddecon_", modelname, "_stdout"))) &&  force == TRUE)) ) {
      oldwd <- getwd()
      setwd(afnidir)
      system2("bash", paste0("3ddecon_", modelname, ".bash"), stdout=paste0("3ddecon_", modelname, "_stdout"), stderr=paste0("3ddecon_", modelname, "_stderr"))
      setwd(oldwd)
    }
    
  }
  
}
