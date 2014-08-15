#second-level model for BPD data
#generate group analysis (what FSL sometimes calls 3rd level) based on emotions of first level runs and multiple sessions per subject
#also drop runs where level of motion was unacceptable

fmriDir <- "/Volumes/Serena/SPECC/MR_Raw"
fitDir <- file.path(getMainDir(), "clock_analysis", "fmri", "fmri_fits")

featRuns <- list.files(path=fmriDir, pattern="FEAT_LVL1_run[0-9].feat", include.dirs=TRUE, recursive=TRUE, full.names=TRUE)

#flag runs with more than 10% volumes with FD 0.9mm or greater
#find fd.txt files corresponding to each FEAT run
fdFiles <- sub(paste0(fmriDir, "(.*)/fsl_value/FEAT_LVL1_run([0-9])\\.feat"), paste0(fmriDir, "\\1/clock\\2/motion_info/fd.txt"), featRuns, perl=TRUE)

library(plyr)
motexclude <- ldply(fdFiles, function(f) {
      fd <- read.table(f, header=FALSE)$V1
      propSpikes_0p9 <- sum(as.integer(fd > 0.9))/length(fd)
      spikeExclude <- if (propSpikes_0p9 > .25) 1 else 0
      maxFD <- max(fd)
      meanFD <- mean(fd)
      maxMotExclude <- if (maxFD > 15) 1 else 0
      data.frame(f, propSpikes_0p9, spikeExclude, meanFD, maxFD, maxMotExclude)
    })

motexclude$subid <- factor(sub(paste0(fmriDir, "/([^/]+)/MBclock_recon/.*$"), "\\1", featRuns, perl=TRUE))
motexclude <- ddply(motexclude, .(subid), function(subdf) {
      if (nrow(subset(subdf, maxMotExclude == 0 & spikeExclude == 0)) < 4) {
        subdf$lt4runs <- 1
      } else {
        subdf$lt4runs <- 0
      }
      subdf
    })

motexclude$anyExclude <- with(motexclude, as.integer(spikeExclude | maxMotExclude | lt4runs))
motexclude$featRun <- featRuns

featRuns_drops <- motexclude[-1*which(motexclude$anyExclude==1),]

inputs <- featRuns_drops$featRun



