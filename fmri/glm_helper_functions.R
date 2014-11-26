#cutting down on redundancy across glm setup scripts

#wrapper for running an fsl command safely within R
#if FSL does not have its configuration setup properly, commands such as feat don't work, or hang strangely
runFSLCommand <- function(args, fsldir=NULL, stdout=NULL, stderr=NULL) {
  #look for FSLDIR in system environment if not passed in
  if (is.null(fsldir)) {
    #check for FSLDIR in sourced .bashrc
    bashrc_fsldir <- character(0)
    if (file.exists("~/.profile")) {
      bashrc_fsldir <- system("source ~/.profile && echo $FSLDIR", intern=TRUE)
    }
    
    #check for FSLDIR in current environment
    env <- system("env", intern=TRUE)
    if (length(fsldir <- grep("^FSLDIR=", env, value=TRUE)) > 0L) {
      fsldir <- sub("^FSLDIR=", "", fsldir)
    } else if (!identical(bashrc_fsldir, character(0))) {
      fsldir <- bashrc_fsldir      
    } else {
      warning("FSLDIR not found in environment. Defaulting to /usr/local/fsl.")
      fsldir <- "/usr/local/fsl"
    }
  }
  
  Sys.setenv(FSLDIR=fsldir) #export to R environment
  fslsetup=paste0("FSLDIR=", fsldir, "; PATH=${FSLDIR}/bin:${PATH}; . ${FSLDIR}/etc/fslconf/fsl.sh; ${FSLDIR}/bin/")
  fslcmd=paste0(fslsetup, args)
  if (!is.null(stdout)) { fslcmd=paste(fslcmd, ">", stdout) }
  if (!is.null(stderr)) { fslcmd=paste(fslcmd, "2>", stderr) }
  cat("FSL command: ", fslcmd, "\n")
  retcode <- system(fslcmd)
  return(retcode)
}

gen_emo_interaction_regressors <- function(examplefile, regressors, emotions=c("fear","scram","happy"), timingdir, mrrunnums, runlengths, dropVolumes) {
    ##make between-session regressors based on emotion condition.
    ##fmriDir <- "/Volumes/Serena/MMClock/MR_Raw"
    ##fmriDir <- "/Volumes/Serena/MMClock/MR_Proc"
    fmriDir <- "/Volumes/Serena/SPECC/MR_Proc"
    fitDir <- file.path(getMainDir(), "clock_analysis", "fmri", "fmri_fits")
    ##/Volumes/Serena/MMClock/MR_Raw/10997_20140308/MBclock_recon/clock1/nfswudktm_clock1_5_trunc282.nii.gz
    ##subid <- factor(sub(paste0(fmriDir, "/([0-9]{5})_\\d+/MBclock_recon/.*$"), "\\1", mrfiles[1L], perl=TRUE))
    ##subid <- factor(sub(paste0(fmriDir, "/([0-9]{5})_\\d+/mni_5mm_wavelet/.*$"), "\\1", examplefile, perl=TRUE))
    subid <- as.integer(sub(paste0(fmriDir, "/([0-9]{3})[A-z]{2}_.*/mni_5mm_wavelet/.*$"), "\\1", examplefile, perl=TRUE))
    
    
    loc <- local({load(file.path(fitDir, paste0(as.character(subid), "_fitinfo.RData"))); environment()})$f #time-clock fit object (load as local var)
    emocon <- data.frame(emotion=loc$run_condition[mrrunnums], contingency=loc$rew_function[mrrunnums]) #vector of emotion and contingency

    ##generate interactions for ev, rpe_neg, and rpe_pos with emotion
    csum <- cumsum(runlengths - dropVolumes)

    for (reg in regressors) {
        for (emo in emotions) {
            vec <- rep(0, max(csum))
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
}


truncateRuns <- function(s, mrfiles, mrrunnums, niftivols, dropVolumes=0) {
  ##Identify the last valid volume acquired in a given run.
  ##Subjects often exhibit head movement after run ends (MATLAB closes), but scan hasn't stopped
  ##This occurs because the MB raw transfer of the prior run is occurring, but does not finish before the current run
  ##Thus, truncate mr files to be 12 seconds after final feedback presentation, which is how the paradigm timing files are setup
  ##note that all of this would need to be reworked if TR were not 1.0 (i.e., 1 second = 1 volume)

  mrdf <- do.call(rbind, lapply(1:length(mrfiles), function(r) {
            iti_durations <- s$runs[[ mrrunnums[r] ]]$orig_data_frame$iti_ideal
            last_iti <- s$runs[[ mrrunnums[r] ]]$iti_onset[length(s$runs[[ mrrunnums[r] ]]$iti_onset)]
            last_vol_behavior <- floor(last_iti + iti_durations[length(iti_durations)]) #use floor to select last vol in the iti window
            first_vol <- dropVolumes #first volume to use for analysis 
            
            if (last_vol_behavior < niftivols[r]) {
              ##more vols were acquired than presented in paradigm. Thus, truncation may be needed
              ##check framewise displacement and truncate earlier than 12 second ITI if a big movement occurred
              fd <- read.table(file.path(dirname(mrfiles[r]), "motion_info", "fd.txt"))$V1
              badfd <- do.call(c, sapply(1:length(fd), function(x) { if (x >= last_iti && fd[x] > 0.9) x else NULL })) #flag volumes after last_iti with high FD
              if (length(badfd) == 0L) {
                ##no frames flagged in last volumes
                last_vol_analysis <- last_vol_behavior
              } else {
                ##use either the last volume of the task or the volume before the earliest bad movement 
                last_vol_analysis <- min(last_vol_behavior, (min(badfd) - 1))
              }
              
              #length of truncated file
              truncLength <- last_vol_analysis - first_vol
              
              #generate filename for truncated volume
              if (first_vol > 0) {
                truncfile <- sub("(^.*/[a-z]+_clock[0-9](?:_5)*)\\.nii\\.gz$", paste0("\\1_drop", dropVolumes, "_trunc", last_vol_analysis, ".nii.gz"), mrfiles[r], perl=TRUE)  
              } else {
                truncfile <- sub("(^.*/[a-z]+_clock[0-9](?:_5)*)\\.nii\\.gz$", paste0("\\1_trunc", last_vol_analysis, ".nii.gz"), mrfiles[r], perl=TRUE)
              }
              
              if (!file.exists(truncfile)) { runFSLCommand(paste("fslroi", mrfiles[r], truncfile, first_vol, truncLength)) } #create truncated volume
              mrfile_to_analyze <- truncfile
            } else {
              last_vol_analysis <- niftivols[r] 
              if (dropVolumes > 0) {
                truncLength <- niftivols[r] - dropVolumes
                truncfile <- sub("(^.*/[a-z]+_clock[0-9](?:_5)*)\\.nii\\.gz$", paste0("\\1_drop", dropVolumes, ".nii.gz"), mrfiles[r], perl=TRUE)
                if (!file.exists(truncfile)) { runFSLCommand(paste("fslroi", mrfiles[r], truncfile, first_vol, truncLength)) } #create truncated volume
                mrfile_to_analyze <- truncfile
              } else {
                mrfile_to_analyze <- mrfiles[r] #just use original file  
              }
              
            }
            #cat(paste0(paste(mrfiles[r], niftivols[r], floor(last_iti), truncLength, sep="\t"), "\n"), file="trunclog", append=TRUE)
            return(data.frame(last_vol_analysis, mrfile_to_analyze, stringsAsFactors=FALSE))
          }))
  
  mrdf
  
}

pca_motion <- function(mrfiles, runlengths, motion_parfile="motion.par", verbose=FALSE,  numpcs=3, dropVolumes=0) {
  #based on a vector of mr files to be analyzed, compute the PCA decomposition of motion parameters and their derivatives
  #do this for each run separately (e.g., for FSL or R glm), as well as concatenated files
  motion_runs <- lapply(1:length(mrfiles), function(i)  {
        mot <- read.table(file.path(dirname(mrfiles[i]), motion_parfile), col.names=c("r.x", "r.y", "r.z", "t.x", "t.y", "t.z"))
        mot <- mot[(1+dropVolumes):runlengths[i],]
        motderiv <- as.data.frame(lapply(mot, function(col) { c(0, diff(col)) }))
        names(motderiv) <- paste0("d.", names(mot)) #add delta to names
        cbind(mot, motderiv)
      })
  
  motion_pcs_runs <- lapply(1:length(motion_runs), function(r) {
        pc <- prcomp(motion_runs[[r]], retx=TRUE)
        cumvar <- cumsum(pc$sdev^2/sum(pc$sdev^2))
        if (verbose) message("first", numpcs, "motion principal components account for: ", round(cumvar[numpcs], 3))
        mregressors <- pc$x[,1:numpcs] #cf Churchill et al. 2012 PLoS ONE
      })
  
  motion_concat <- do.call(rbind, motion_runs)
  pc <- prcomp(motion_concat, retx=TRUE)
  cumvar <- cumsum(pc$sdev^2/sum(pc$sdev^2))
  
  if (verbose) message("first", numpcs, "motion principal components account for: ", round(cumvar[numpcs], 3))
  motion_pcs_concat <- pc$x[,1:numpcs] #cf Churchill et al. 2012 PLoS ONE
  
  if (verbose) {
    cat("correlation of motion parameters:\n\n")
    print(round(cor(motion_concat), 2))
  }
  
  return(list(motion_pcs_runs=motion_pcs_runs, motion_pcs_concat=motion_pcs_concat))
  
}

generateRunMask <- function(mrfiles, outdir=getwd(), outfile="runmask") {
  if (file.exists(file.path(outdir, paste0(outfile, ".nii.gz")))) { return(invisible(NULL)) }
  ##generate mask of mrfiles where temporal min is > 0 for all runs
  for (f in 1:length(mrfiles)) {
    runFSLCommand(paste0("fslmaths ", mrfiles[f], " -Tmin -bin ", outdir, "/tmin", f))#, fsldir="/usr/local/ni_tools/fsl")
  }
  
  ##sum mins together over runs and threshold at number of runs
  runFSLCommand(paste0("fslmaths ", paste(paste0(outdir, "/tmin", 1:length(mrfiles)), collapse=" -add "), " ", outdir, "/tminsum"))#, fsldir="/usr/local/ni_tools/fsl")
  runFSLCommand(paste0("fslmaths ", outdir, "/tminsum -thr ", length(mrfiles), " -bin ", outdir, "/", outfile))#, fsldir="/usr/local/ni_tools/fsl")
  runFSLCommand(paste0("imrm ", outdir, "/tmin*"))#, fsldir="/usr/local/ni_tools/fsl") #cleanup 
  
}
