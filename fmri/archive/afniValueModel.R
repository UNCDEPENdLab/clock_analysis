#function to setup and run afni value model
afniValueModel <- function(f_value, mrfiles, runlengths, mrrunnums, run=FALSE, force=FALSE, dropVolumes=0) {
  #setup and run basic value model in AFNI
  #afnidir <- file.path(normalizePath(file.path(dirname(mrfiles[1L]), "..")), "afni_value") #note: normalizePath will fail to evaluate properly if directory does not exist (e.g., dir not created yet)
  afnidir <- file.path(normalizePath(file.path(dirname(mrfiles[1L]), "..")), "afni_value_evtmax") #note: normalizePath will fail to evaluate properly if directory does not exist (e.g., dir not created yet)
  if (file.exists(afnidir) && force==FALSE) {
    message("afni_value directory already exists. ", afnidir, ". Skipping subject")
    return(NULL)
  }
  
  cat("afnidir create: ", afnidir, "\n")
  dir.create(afnidir, showWarnings=FALSE) #one directory up from a given clock run
  timingdir <- file.path(afnidir, "run_timing_deltavalue")
  
  #old version with 0 duration clock and feedback regressors
#  d_value <- build_design_matrix(fitobj=f_value, regressors=c("clock", "feedback", "ev", "rpe_neg", "rpe_pos"), 
#      event_onsets=c("clock_onset", "feedback_onset", "clock_onset", "feedback_onset", "feedback_onset"), 
#      durations=c(0, 0, "clock_duration", "feedback_duration", "feedback_duration"), 
#      baselineCoefOrder=2, writeTimingFiles=c("AFNI"),
#      runVolumes=runlengths, runsToOutput=mrrunnums, output_directory=timingdir)
  
  #new version with duration-convolved event regressors
  d_value <- build_design_matrix(fitobj=f_value, regressors=c("clock", "ev", "feedback", "rpe_neg", "rpe_pos"), 
                                 event_onsets=c("clock_onset", "clock_onset", "feedback_onset", "feedback_onset", "feedback_onset"), 
                                 durations=c("clock_duration", "clock_duration", "feedback_duration", "feedback_duration", "feedback_duration"),
                                 normalizations=c("durmax_1.0", "evtmax_1.0", "durmax_1.0", "evtmax_1.0", "evtmax_1.0"),
                                 baselineCoefOrder=2, writeTimingFiles=c("AFNI"), center_values=TRUE, convolve_wi_run = TRUE,
                                 runVolumes=runlengths, runsToOutput=mrrunnums, output_directory=timingdir, dropVolumes=dropVolumes)
  
  #another version with duration-convolved events and instant effects of parametric modulators at the time of feedback. 
#  d_value <- build_design_matrix(fitobj=f_value, regressors=c("clock", "ev", "feedback", "rpe_neg", "rpe_pos"), 
#      event_onsets=c("clock_onset", "feedback_onset", "feedback_onset", "feedback_onset", "feedback_onset"), 
#      durations=c("clock_duration", 0, "feedback_duration", 0, 0),
#      normalizations=c("durmax_1.0", "evtmax_1.0", "durmax_1.0", "evtmax_1.0", "evtmax_1.0"),
#      baselineCoefOrder=2, writeTimingFiles=c("AFNI"),
#      runVolumes=runlengths, runsToOutput=mrrunnums, output_directory=timingdir)
  
  motpcs <- pca_motion(mrfiles, runlengths, motion_parfile="motion.par", numpcs=3, dropVolumes=dropVolumes)
  write.table(motpcs$motion_pcs_concat, file=file.path(afnidir, 'motion_pcs.txt'), col.names=FALSE, row.names=FALSE)    
  
  #concat motion censor
  #use censor_intersection.1D, which will flag FD > 0.9 AND DVARS > 20.
  #despiking should tamp down most of the problematic volumes.
  censor_intersection_concat <- do.call(c, lapply(1:length(mrfiles), function(i)  {
            cen <- read.table(file.path(dirname(mrfiles[i]), "motion_info", "censor_intersection.1D"))$V1
            cen <- cen[(1+dropVolumes):runlengths[i]]
            cen
          }))
  
  write.table(censor_intersection_concat, file=file.path(afnidir, "censor_intersection_concat.1D"), col.names=FALSE, row.names=FALSE)
  
  gen_emo_interaction_regressors(mrfiles[1L], regressors=c("clock", "ev", "feedback", "rpe_neg", "rpe_pos"), emotions=c("fear","scram","happy"),
                                 timingdir, mrrunnums, runlengths, dropVolumes)

  generateRunMask(mrfiles, outdir=afnidir, outfile="runmask")
  
  #note that the 3ddeconvolve call is the same as earlier iterations
  #but we have significantly altered the parameterization by normalizing the HRF to be duration-convolved for event regressors
  #and 1.0-height convolved for parametric regressors. Regressors are also demeaned (without zeros) prior to convolution to reduce collinearity.
  
  #model ignoring emotion across runs
  modelname <- "glm_hrf_clock_preconvolve_valueModel_noemo_normalized"
  afniScript <- c('#!/bin/bash',
      'source /Users/michael/.bashrc')
  
  afniScript <- c(afniScript, paste(
          '3dDeconvolve',
          '-overwrite',
          '-input',
          paste(mrfiles, collapse=" \\\n"),
          '-tout -rout -allzero_OK -polort 2 -jobs 8',
          '-num_stimts 8 -mask runmask.nii.gz', 
          paste('-stim_file 1', file.path(timingdir, "clock_concat.1D"), '-stim_label 1 clock_onset'),
          paste('-stim_file 2', file.path(timingdir, "feedback_concat.1D"), '-stim_label 2 feedback_onset'),
          paste('-stim_file 3', file.path(timingdir, "ev_concat.1D"), '-stim_label 3 ev'),
          paste('-stim_file 4', file.path(timingdir, "rpe_neg_concat.1D"), '-stim_label 4 rpe_neg'),
          paste('-stim_file 5', file.path(timingdir, "rpe_pos_concat.1D"), '-stim_label 5 rpe_pos'),
          '-stim_file 6 "motion_pcs.txt[0]" -stim_label 6 motpc1 -stim_base 6',
          '-stim_file 7 "motion_pcs.txt[1]" -stim_label 7 motpc2 -stim_base 7',
          '-stim_file 8 "motion_pcs.txt[2]" -stim_label 8 motpc3 -stim_base 8',      
          '-censor censor_intersection_concat.1D',
          paste0('-bucket ', modelname, '_stats'),
          #paste0('-fitts ', modelname, '_fitts'),
          #paste0('-errts ', modelname, '_errts'),
          paste0('-x1D ', modelname, '_x1D'),
          paste0('-cbucket ', modelname, '_coefs'),
          paste0('-xjpeg ', modelname, '_design.png'),
          '2>&1', 
          sep=" \\\n"))
  
  #add 3dREMLfit call
  afniScript <- c(afniScript,
      '',
      '#now rerun whole thing using REMLfit\n',
      paste0('chmod +x ', modelname, '_stats.REML_cmd\n'),
      paste0('./', modelname, '_stats.REML_cmd 2>&1 | tee reml_', modelname, '_stats.log'))
  
  writeLines(afniScript, con=file.path(afnidir, paste0("3ddecon_", modelname, ".bash")))
  
  if (run==TRUE) {
    oldwd <- getwd()
    setwd(afnidir)
    system2("bash", paste0("3ddecon_", modelname, ".bash"), stdout=paste0("3ddecon_", modelname, "_stdout"), stderr=paste0("3ddecon_", modelname, "_stderr"))
    setwd(oldwd)
  }

  #model where interaction with emotion across runs is included
  modelname <- "glm_hrf_clock_preconvolve_valueModel_emoint_normalized"
  afniScript <- c('#!/bin/bash',
      'source /Users/michael/.bashrc')
  
  afniScript <- c(afniScript, paste(
      '3dDeconvolve',
      '-overwrite',
      '-input',
      paste(mrfiles, collapse=" \\\n"),
      '-tout -rout -allzero_OK -polort 3 -jobs 8',
      '-num_stimts 18 -mask runmask.nii.gz', 
      paste('-stim_file 1', file.path(timingdir, "clock_concat.1D"), '-stim_label 1 clock'),
      paste('-stim_file 2', file.path(timingdir, "clock_happy_concat.1D"), '-stim_label 2 clock_happy'),
      paste('-stim_file 3', file.path(timingdir, "clock_fear_concat.1D"), '-stim_label 3 clock_fear'),
      paste('-stim_file 4', file.path(timingdir, "feedback_concat.1D"), '-stim_label 4 feedback'),
      paste('-stim_file 5', file.path(timingdir, "feedback_happy_concat.1D"), '-stim_label 5 feedback_happy'),
      paste('-stim_file 6', file.path(timingdir, "feedback_fear_concat.1D"), '-stim_label 6 feedback_fear'),
      paste('-stim_file 7', file.path(timingdir, "ev_concat.1D"), '-stim_label 7 ev'),
      paste('-stim_file 8', file.path(timingdir, "ev_happy_concat.1D"), '-stim_label 8 ev_happy'),
      paste('-stim_file 9', file.path(timingdir, "ev_fear_concat.1D"), '-stim_label 9 ev_fear'),
      paste('-stim_file 10', file.path(timingdir, "rpe_neg_concat.1D"), '-stim_label 10 rpe_neg'),
      paste('-stim_file 11', file.path(timingdir, "rpe_neg_happy_concat.1D"), '-stim_label 11 rpe_neg_happy'),
      paste('-stim_file 12', file.path(timingdir, "rpe_neg_fear_concat.1D"), '-stim_label 12 rpe_neg_fear'),
      paste('-stim_file 13', file.path(timingdir, "rpe_pos_concat.1D"), '-stim_label 13 rpe_pos'),
      paste('-stim_file 14', file.path(timingdir, "rpe_pos_happy_concat.1D"), '-stim_label 14 rpe_pos_happy'),
      paste('-stim_file 15', file.path(timingdir, "rpe_pos_fear_concat.1D"), '-stim_label 15 rpe_pos_fear'),
      '-stim_file 16 "motion_pcs.txt[0]" -stim_label 16 motpc1 -stim_base 16',
      '-stim_file 17 "motion_pcs.txt[1]" -stim_label 17 motpc2 -stim_base 17',
      '-stim_file 18 "motion_pcs.txt[2]" -stim_label 18 motpc3 -stim_base 18',      
      '-num_glt 35',
      "-gltsym 'SYM: +1*clock' -glt_label 1 m_clock_scram",
      "-gltsym 'SYM: +1*clock +1*clock_happy' -glt_label 2 m_clock_happy",
      "-gltsym 'SYM: +1*clock +1*clock_fear' -glt_label 3 m_clock_fear",
      "-gltsym 'SYM: +1*clock +0.25*clock_fear +0.25*clock_happy' -glt_label 4 m_clock_overall",
      "-gltsym 'SYM: +1*clock_fear' -glt_label 5 clock_fear_gt_scram",
      "-gltsym 'SYM: +1*clock_happy' -glt_label 6 clock_happy_gt_scram",
      "-gltsym 'SYM: +1*clock_fear -1*ev_happy' -glt_label 7 clock_fear_gt_happy",
      "-gltsym 'SYM: +1*feedback' -glt_label 8 m_feedback_scram",
      "-gltsym 'SYM: +1*feedback +1*feedback_happy' -glt_label 9 m_feedback_happy",
      "-gltsym 'SYM: +1*feedback +1*feedback_fear' -glt_label 10 m_feedback_fear",
      "-gltsym 'SYM: +1*feedback +0.25*feedback_fear +0.25*feedback_happy' -glt_label 11 m_feedback_overall",
      "-gltsym 'SYM: +1*feedback_fear' -glt_label 12 feedback_fear_gt_scram",
      "-gltsym 'SYM: +1*feedback_happy' -glt_label 13 feedback_happy_gt_scram",
      "-gltsym 'SYM: +1*feedback_fear -1*ev_happy' -glt_label 14 feedback_fear_gt_happy",
      "-gltsym 'SYM: +1*ev' -glt_label 15 m_ev_scram",
      "-gltsym 'SYM: +1*ev +1*ev_happy' -glt_label 16 m_ev_happy",
      "-gltsym 'SYM: +1*ev +1*ev_fear' -glt_label 17 m_ev_fear",
      "-gltsym 'SYM: +1*ev +0.25*ev_fear +0.25*ev_happy' -glt_label 18 m_ev_overall",
      "-gltsym 'SYM: +1*ev_fear' -glt_label 19 ev_fear_gt_scram",
      "-gltsym 'SYM: +1*ev_happy' -glt_label 20 ev_happy_gt_scram",
      "-gltsym 'SYM: +1*ev_fear -1*ev_happy' -glt_label 21 ev_fear_gt_happy",
      "-gltsym 'SYM: +1*rpe_neg' -glt_label 22 m_rpe_neg_scram",
      "-gltsym 'SYM: +1*rpe_neg +1*rpe_neg_happy' -glt_label 23 m_rpe_neg_happy",
      "-gltsym 'SYM: +1*rpe_neg +1*rpe_neg_fear' -glt_label 24 m_rpe_neg_fear",
      "-gltsym 'SYM: +1*rpe_neg +0.25*rpe_neg_fear +0.25*rpe_neg_happy' -glt_label 25 m_rpe_neg_overall",
      "-gltsym 'SYM: +1*rpe_neg_fear' -glt_label 26 rpe_neg_fear_gt_scram",
      "-gltsym 'SYM: +1*rpe_neg_happy' -glt_label 27 rpe_neg_happy_gt_scram",
      "-gltsym 'SYM: +1*rpe_neg_fear -1*rpe_neg_happy' -glt_label 28 rpe_neg_fear_gt_happy",
      "-gltsym 'SYM: +1*rpe_pos' -glt_label 29 m_rpe_pos_scram",
      "-gltsym 'SYM: +1*rpe_pos +1*rpe_pos_happy' -glt_label 30 m_rpe_pos_happy",
      "-gltsym 'SYM: +1*rpe_pos +1*rpe_pos_fear' -glt_label 31 m_rpe_pos_fear",
      "-gltsym 'SYM: +1*rpe_pos +0.25*rpe_pos_fear +0.25*rpe_pos_happy' -glt_label 32 m_rpe_pos_overall",
      "-gltsym 'SYM: +1*rpe_pos_fear' -glt_label 33 rpe_pos_fear_gt_scram",
      "-gltsym 'SYM: +1*rpe_pos_happy' -glt_label 34 rpe_pos_happy_gt_scram",
      "-gltsym 'SYM: +1*rpe_pos_fear -1*rpe_pos_happy' -glt_label 35 rpe_pos_fear_gt_happy",
      '-censor censor_intersection_concat.1D',
      paste0('-bucket ', modelname, '_stats'),
      ##paste0('-fitts ', modelname, '_fitts'),
      ##paste0('-errts ', modelname, '_errts'),
      paste0('-x1D ', modelname, '_x1D'),
      paste0('-cbucket ', modelname, '_coefs'),
      paste0('-xjpeg ', modelname, '_design.png'),
      '2>&1', 
      sep=" \\\n"))
  
  #add 3dREMLfit call
  afniScript <- c(afniScript,
      '',
      '#now rerun whole thing using REMLfit\n',
      paste0('chmod +x ', modelname, '_stats.REML_cmd\n'),
      paste0('./', modelname, '_stats.REML_cmd 2>&1 | tee reml_', modelname, '_stats.log'))
  
  writeLines(afniScript, con=file.path(afnidir, paste0("3ddecon_", modelname, ".bash")))
  
  if (run==TRUE) {
    oldwd <- getwd()
    setwd(afnidir)
    system2("bash", paste0("3ddecon_", modelname, ".bash"), stdout=paste0("3ddecon_", modelname, "_stdout"), stderr=paste0("3ddecon_", modelname, "_stderr"))
    setwd(oldwd)
  }
  
  ##DMBLOCK MODEL
  #now using internal hrf_convolve_normalize to accomplish similar model with more flexibility
#    modelname <- "glm_hrf_clock_dmBLOCK_valueModel"
#    afniScript <- c('#!/bin/bash',
#            'source /Users/michael/.bashrc')
#    
#    afniScript <- c(afniScript, paste(
#                    '3dDeconvolve',
#                    '-overwrite',
#                    '-input',
#                    paste(mrfiles, collapse=" \\\n"),
#                    '-tout -allzero_OK -polort 3 -jobs 8',
#                    '-num_stimts 5 -mask runmask.nii.gz',
#                    paste('-stim_times_AM2 1', file.path(timingdir, "clock_ev_dmBLOCK.txt"), '\'dmUBLOCK(1)\' -stim_label 1 clock_ev'),
#                    paste('-stim_times_AM2 2', file.path(timingdir, "feedback_rpe_neg_rpe_pos_dmBLOCK.txt"), '\'dmUBLOCK(1)\' -stim_label 2 feedback_rpe'),
#                    '-stim_file 3 "motion_concat.par[0]" -stim_label 3 motpc1 -stim_base 3',
#                    '-stim_file 4 "motion_concat.par[1]" -stim_label 4 motpc2 -stim_base 4',
#                    '-stim_file 5 "motion_concat.par[2]" -stim_label 5 motpc3 -stim_base 5',      
#                    '-censor censor_intersection_concat.1D',
#                    paste0('-bucket ', modelname, '_stats'),
#                    #paste0('-fitts ', modelname, '_fitts'),
#                    #paste0('-errts ', modelname, '_errts'),
#                    paste0('-x1D ', modelname, '_x1D'),
#                    paste0('-cbucket ', modelname, '_coefs'),
#                    paste0('-xjpeg ', modelname, '_design.png'),
#                    '-GOFORIT 10 2>&1', 
#                    sep=" \\\n"))
#    
#    #add 3dREMLfit call
#    afniScript <- c(afniScript,
#            '',
#            '#now rerun whole thing using REMLfit\n',
#            paste0('chmod +x ', modelname, '_stats.REML_cmd\n'),
#            paste0('./', modelname, '_stats.REML_cmd 2>&1 | tee reml_', modelname, '_stats.log'))
#    
#    writeLines(afniScript, con=file.path(afnidir, "3ddecon_dmBLOCK_valueModel.bash"))
#    
#    if (run==TRUE) {
#        oldwd <- getwd()
#        setwd(afnidir)
#        system2("bash", "3ddecon_dmBLOCK_valueModel.bash", stdout="3ddecon_dmBLOCK_stdout", stderr="3ddecon_dmBLOCK_stderr")
#        setwd(oldwd)
#    }
  
}
