# fit behavioral data for all participants who completed emo clock in scanner
# setup model-based fMRI GLM analysis for based on fitted data

library(fitclock)
library(Rniftilib) #has nice function for just reading header

#wrapper for running an fsl command safely within R
#if FSL does not have its configuration setup properly, commands such as feat don't work, or hang strangely
runFSLCommand <- function(args, fsldir=NULL, stdout=NULL, stderr=NULL) {
  #look for FSLDIR in system environment if not passed in
  if (is.null(fsldir)) {
    env <- system("env", intern=TRUE)
    if (length(fsldir <- grep("^FSLDIR=", env, value=TRUE)) > 0L) {
      fsldir <- sub("^FSLDIR=", "", fsldir)
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


setwd(file.path(getMainDir(), "clock_analysis", "fmri"))
source("afniValueModel.R")
source("afniTCModel.R")
source("fslValueModel.R")
if (!file.exists("fmri_fits")) { dir.create("fmri_fits") }
setwd("fmri_fits")

fit_all_fmri <- function(behavDir, fmriDir, idexpr) {
    ##start with base Frank model
    ##force non-negative epsilon (no sticky choice)
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

    behavFiles <- list.files(path=behavDir, pattern="*tcExport.csv", full.names=TRUE, recursive=TRUE)
    
    for (b in behavFiles) {
        #example location of file on bea_res, which contains scan date
        #/Volumes/bea_res/Data/Tasks/EmoClockfMRI/Basic/11229/20140521/Raw/fMRIEmoClock_11229_tc_tcExport.csv
        subid <- sub("^.*fMRIEmoClock_(\\d+)_tc_tcExport.csv$", "\\1", b, perl=TRUE)
        scandate <- sub("^.*/Basic/\\w+/(\\d+)/.*$", "\\1", b, perl=TRUE) 
        mrfiles <- c() #force clear of mr files over subjects to avoid potential persistence from one subject to the next
        
        ##identify corresponding fmri directory
        mrmatch <- grep(eval(idexpr), list.files(fmriDir, full.names=TRUE), perl=TRUE, value=TRUE) #MMClock Y3 format: 11273_20140610
        
        if(length(mrmatch) != 1L) {
            warning("Unable to find fMRI directory for subid: ", subid)
            next
        }
        if (! file.exists(file.path(mrmatch, "MBclock_recon"))) {
            warning("Unable to find preprocessed data MBclock_recon for subid: ", subid)
            next
        }
        
        ##identify fmri run lengths (4th dimension)
        mrfiles <- list.files(mrmatch, pattern="nfswudktm_clock\\d+_5.nii.gz", full.names=TRUE, recursive=TRUE)
        mrrunnums <- as.integer(sub(".*nfswudktm_clock(\\d+)_5.nii.gz$", "\\1", mrfiles, perl=TRUE))
        mrfiles <- mrfiles[order(mrrunnums)] #make absolutely sure that runs are ordered ascending
        
        if (length(mrfiles) == 0L) {
            warning("Unable to find any preprocessed MB files in dir: ", mrmatch)
            next
        }
        
        ##read number of volumes from NIfTI header
        runlengths <- unname(sapply(mrfiles, function(x) { nifti.image.read(x, read_data=0)$dim[4L] }))
        
        ##setup clock data subject object for fitting 
        s <- clockdata_subject(subject_ID=subid, dataset=b)
        
        ##Subjects often exhibit head movement after run ends (MATLAB closes), but scan hasn't stopped
        ##This occurs because the MB raw transfer of the prior run is occurring, but does not finish before the current run
        ##Thus, truncate mr files to be 12 seconds after final feedback presentation, which is how the paradigm timing files are setup
        ##note that all of this would need to be reworked if TR were not 1.0 (i.e., 1 second = 1 volume)

        mrdf <- do.call(rbind, lapply(1:length(mrfiles), function(r) {
            iti_durations <- s$runs[[ mrrunnums[r] ]]$orig_data_frame$iti_ideal
            last_iti <- s$runs[[ mrrunnums[r] ]]$iti_onset[length(s$runs[[ mrrunnums[r] ]]$iti_onset)]
            last_vol <- floor(last_iti + iti_durations[length(iti_durations)]) #use floor to select last vol in the iti window
            if (last_vol < runlengths[r]) {
                ##more vols were acquired than presented in paradigm. Thus, truncation may be needed
                ##check framewise displacement and truncate earlier than 12 second ITI if a big movement occurred...
                fd <- read.table(file.path(dirname(mrfiles[r]), "motion_info", "fd.txt"))$V1
                badfd <- do.call(c, sapply(1:length(fd), function(x) { if (x >= last_iti && fd[x] > 0.9) x else NULL })) #flag volumes after last_iti with high FD
                if (length(badfd) == 0L) {
                    ##no frames flagged in last volumes
                    truncLength <- last_vol
                } else {
                    ##use either the last volume of the task or the earliest bad movement 
                    truncLength <- min(last_vol, min(badfd))
                }
                truncfile <- sub("(^.*/nfswudktm_clock[0-9]_5)\\.nii\\.gz$", paste0("\\1_trunc", truncLength, ".nii.gz"), mrfiles[r], perl=TRUE)
                if (!file.exists(truncfile)) { runFSLCommand(paste("fslroi", mrfiles[r], truncfile, "0", truncLength)) }
                analyze <- truncfile
            } else {
                truncLength <- runlengths[r]
                analyze <- mrfiles[r] #just use original file
            }
            cat(paste0(paste(mrfiles[r], runlengths[r], floor(last_iti), truncLength, sep="\t"), "\n"), file="trunclog", append=TRUE) #very crude log file to track truncation
            return(data.frame(truncLength, analyze, stringsAsFactors=FALSE))
        }))
        
        mrfiles <- mrdf$analyze
        runlengths <- mrdf$truncLength
        
        if (file.exists(paste0(subid, "_fitinfo.RData"))) { 
            cat("Fit data already present for: ", subid, "\n")
            load(paste0(subid, "_fitinfo.RData"))
        } else {           
            cat("Fitting behavioral data for subject: ", subid, "\n")
            
            ##set data for model fit
            posEps$set_data(s)
            
            incr_fit <- posEps$incremental_fit(njobs=7, plot=FALSE)
            
            png(file.path(paste0(subid, "_incrfit.png")), width=9, height=6, units="in", res=300)
            print(incr_fit$AICplot)
            dev.off()
            
            f <- posEps$fit(random_starts=5)

            ##design matrix matching Badre et al. 2012 Neuron
            d <- f$build_design_matrix(regressors=c("mean_uncertainty", "rel_uncertainty", "rpe_pos", "rpe_neg", "rt"), 
                                       event_onsets=c("clock_onset", "clock_onset", "feedback_onset", "feedback_onset", "feedback_onset"), 
                                       durations=c("rt", "rt", "feedback_duration", "feedback_duration", 0), baselineCoefOrder=2, runVolumes=runlengths, runsToOutput=mrrunnums)
            
            ##delta rule value model (simple)
            vm <- deltavalue_model(clock_data=s, alphaV=0.3, betaV=0.3) #N.B. This matches V matrix from full time-clock algorithm fit.
            f_value <- vm$fit() #estimate learning rate as a free parameter
            
            ##Design with EV, clock onset, feedback_onset, PE+, PE-
            d_value <- f_value$build_design_matrix(regressors=c("clock", "feedback", "ev", "rpe_neg", "rpe_pos"), 
                                                   event_onsets=c("clock_onset", "feedback_onset", "clock_onset", "feedback_onset", "feedback_onset"), 
                                                   durations=c(0, 0, "clock_duration", "feedback_duration", "feedback_duration"), baselineCoefOrder=2, runVolumes=runlengths, runsToOutput=mrrunnums)

            save(f_value, d_value, f, d, s, incr_fit, file=paste0(subid, "_fitinfo.RData"))
        }
        
        ##setup afni and FSL models
        afniValueModel(f_value, mrfiles, runlengths, mrrunnums, run=TRUE)
        afniTCModel(f, mrfiles, runlengths, mrrunnums, run=TRUE)
        fslValueModel(f_value, mrfiles, runlengths, mrrunnums, run=TRUE, force=FALSE)
    }

}


#MMClock fit
fit_all_fmri(behavDir="/Volumes/bea_res/Data/Tasks/EmoClockfMRI/Basic",
             fmriDir="/Volumes/Serena/MMClock/MR_Raw",
             idexpr=expression(paste0(subid, "_", scandate))) ##MMClock/LunaID format: 10637_20140302

## fit_all_fmri(behavDir="/Users/michael/Dropbox/Hallquist_K01/Data/fMRI",
##              fmriDir="/Volumes/Serena/SPECC/MR_Raw",
##              idexpr=expression(paste0(sprintf("%03s", subid), "[A-z]{2}_\\d+"))) #SPECC format: 003aa_15Jul2014
