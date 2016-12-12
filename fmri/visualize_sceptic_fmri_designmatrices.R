#look at the design matrices for SCEPTIC data
library(fitclock)
setwd(file.path(getMainDir(), "clock_analysis", "fmri"))
source("glm_helper_functions.R")
if (!file.exists("fmri_fits")) { dir.create("fmri_fits") }
setwd("fmri_fits")

#updated version that used 24 basis functions and also 
if (file.exists("fmri_sceptic_signals_24basis.RData")) { sceptic <- local({load("fmri_sceptic_signals_24basis.RData"); as.list(environment())}) }

behavDir=file.path(getMainDir(), "temporal_instrumental_agent/clock_task/subjects")
behavFiles <- list.files(path=behavDir, pattern=".*tcExport.csv", full.names=TRUE, recursive=TRUE)

#for a moment...
sceptic[["ventropy_decay_matlab_3d"]] <- NULL
sceptic[["ventropy_fixedlrv_matlab_3d"]] <- NULL
sceptic[["rtumax"]] <- NULL

for (b in behavFiles) {
  
  subid <- sub("^.*fMRIEmoClock_(\\d+)_tc_tcExport.csv$", "\\1", b, perl=TRUE)
  mrfiles <- c() #force clear of mr files over subjects to avoid potential persistence from one subject to the next
  
  ##setup clock data subject object for fitting 
  s <- clockdata_subject(subject_ID=subid, dataset=b)
  
  mats3d <- sort(grep("ids", names(sceptic), value=TRUE, invert=TRUE))
  subj_sceptic <- lapply(sceptic[mats3d], function(mat) {
        mat[subid,,]
      })
  
  #create a fit object and populate with timings so that we can build a design matrix.
  f <- clock_fit()
  f$populate_fit(s)
  
  #loop over sceptic signals and plot them under various forms of convolution
  #1) no re-normalization, just mean center prior to convolution
  #2) duration max 1.0 normalization
  #3) event max 1.0 normalization
  #Also look at the effect of z-scoring prior to convolution
  #Plot the boxcars at their unit heights and parameter heights prior to convolution
  
  signals_to_model <- rep(NA_character_, length(subj_sceptic))
  onsets <- rep(NA_character_, length(subj_sceptic))
  durations <- rep(NA_character_, length(subj_sceptic))
  normalizations <- rep("none", length(subj_sceptic))
  for (v in 1:length(subj_sceptic)) {
    thisName <- names(subj_sceptic)[v] #for now, I have a silly divergence where the design matrix function expects a sceptic_ prefix, but the fit object should not have this...
    signals_to_model[v] <- paste0("sceptic_", thisName) #name of regressor in build_design_matrix
    f$sceptic[[thisName]] <- split(subj_sceptic[[v]], row(subj_sceptic[[v]]))
    
    if (thisName %in% c("vauc", "vchosen", "ventropy", "vmax", "vsd", "umax")) {
      onsets[v] <- "clock_onset"; durations[v] <- "clock_duration"
    } else if (thisName %in% c("pemax", "peauc", "dauc", "dsd")) {
      onsets[v] <- "feedback_onset"; durations[v] <- "feedback_duration"
    }
  }

  
  #visualize the design matrix as a series of boxcars
  
  d <- build_design_matrix(fitobj=f, regressors=c("clock", "feedback", signals_to_model), 
      event_onsets=c("clock_onset", "feedback_onset", onsets),
      durations=c("clock_duration", "feedback_duration", durations), 
      normalizations=c("none", "none", normalizations), #no re-normalization in case of no convolution
      #baselineCoefOrder=-1, center_values=FALSE, convolve=FALSE,
      baselineCoefOrder=-1, center_values=TRUE, convolve=TRUE,
      dropVolumes=0)
  
  pdf("run_designs_convolve_Dec2016.pdf", width=15, height=14)
  for (r in 1:length(d$design.convolve)) {
    g <- visualizeDesignMatrix(d$design.convolve[[r]], outfile=NULL, includeBaseline=FALSE)
    plot(g+ggtitle(paste("Run", r))+theme_bw(base_size=12))
    
    plot(corHeatmap(d$design.convolve[[r]], alphaSort=TRUE, tileTextSize=6))
    
  }
  dev.off()
  
  #now build up convolved variants...
  
  
#  d <- build_design_matrix(fitobj=f, regressors=c("clock", "feedback", signals_to_model), 
#      event_onsets=c("clock_onset", "feedback_onset", onsets),
#      durations=c("clock_duration", "feedback_duration", durations), 
#      normalizations=c("durmax_1.0", "durmax_1.0", normalizations),
#      baselineCoefOrder=2, writeTimingFiles=c("AFNI"), center_values=TRUE, convolve_wi_run=TRUE,
#      runVolumes=runlengths, runsToOutput=mrrunnums, output_directory=timingdir, dropVolumes=dropVolumes)
  
}
