#second-level model for MMY3 data
#generate group analysis (what FSL sometimes calls 3rd level) based on emotions of first level runs and multiple sessions per subject
#also drop runs where level of motion was unacceptable

setwd(file.path(getMainDir(), "clock_analysis", "fmri"))
source("glm_helper_functions.R")

##fmriDir <- "/Volumes/Serena/MMClock/MR_Proc/10873_20140918/mni_5mm_wavelet"
##fmriDir <- "/Volumes/Serena/MMClock/MR_Proc"
fmriDir <- "/storage/group/mnh5174_collab/MMClock/MR_Proc"

##not used: just look for designmatrix.RData in timing directory
##fitDir <- file.path(getMainDir(), "clock_analysis", "fmri", "fmri_fits")

featRuns_all <- list.files(path=fmriDir, pattern="FEAT_LVL1_run[0-9].feat", include.dirs=TRUE, recursive=TRUE, full.names=TRUE)
#featRuns <- featRuns_all[grepl("/fsl_tc_nomeanunc/", featRuns_all, fixed=TRUE)] #filter to only nomeanunc folders
#featRuns <- featRuns_all[grepl("/sceptic_vchosen/", featRuns_all, fixed=TRUE)]
featRuns <- featRuns_all[grepl("/sceptic_.*/", featRuns_all)]

cat("All runs identified:\n")
print(featRuns)

#flag runs with more than 15% volumes with FD 0.9mm or greater
#find fd.txt files corresponding to each FEAT run
fdFiles <- sub(paste0(fmriDir, "(.*)/[^/]+/FEAT_LVL1_run([0-9])\\.feat"), paste0(fmriDir, "\\1/clock\\2/motion_info/fd.txt"), featRuns, perl=TRUE)

#identify matching truncated 4d files (created by model_clock_fmri.R) since we are only concerned about movement within the run proper (not dead volumes)  
niTruncFiles <- Sys.glob(sub(paste0(fmriDir, "(.*)/[^/]+/FEAT_LVL1_run([0-9])\\.feat"), paste0(fmriDir, "\\1/clock\\2/nfsw*drop*.nii.gz"), featRuns, perl=TRUE))

#identify how many volumes were dropped at the beginning
dropVolumes <- as.integer(sub("^.*_drop(\\d+).*.nii.gz$", "\\1", niTruncFiles, perl=TRUE))

##some runs are not truncated (e.g., if they go all the way to 350 vols and the scanner stops)
##this results in NAs here
truncLengths <- as.integer(sub("^.*_trunc(\\d+).nii.gz$", "\\1", niTruncFiles, perl=TRUE))
nafiles <- which(is.na(truncLengths))

if (length(nafiles) > 0L) {
    library(Rniftilib)
    ##need to add back the dropped volumes since everything below assumes that the trunc length is the final volume in the original time series
    truncLengths[nafiles] <- unname(sapply(nafiles, function(x) { Rniftilib::nifti.image.read(niTruncFiles[x], read_data=0)$dim[4L] })) + dropVolumes[nafiles]
    detach("package:Rniftilib", unload=TRUE) #necessary to avoid dim() conflict with oro.nifti
}

library(plyr)
motexclude <- ldply(1:length(fdFiles), function(i) {
      fd <- read.table(fdFiles[i], header=FALSE)$V1
      fd <- fd[(dropVolumes[i]+1):truncLengths[i]] #only include volumes within run
      propSpikes_0p9 <- sum(as.integer(fd > 0.9))/length(fd)
      ##if (is.na(propSpikes_0p9[1])) browser()
      spikeExclude <- if (propSpikes_0p9 > .10) 1 else 0
      maxFD <- max(fd)
      meanFD <- mean(fd)
      maxMotExclude <- if (maxFD > 10) 1 else 0
      data.frame(f=fdFiles[i], propSpikes_0p9, spikeExclude, meanFD, maxFD, maxMotExclude)
    })

motexclude$subid <- factor(sub(paste0(fmriDir, "/([0-9]{5})_\\d+/mni_5mm_wavelet/.*$"), "\\1", featRuns, perl=TRUE))
##motexclude$subid <- factor(paste0("10873", sub(".*/mni_5mm_wavelet/(fsl_.*)/.*$", "\\1", featRuns)))
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

#10637 has pretty bad head movement in runs 5-8... in runs 7 and 8, it falls just below 10% FD > 0.9mm, so exclude subject altogether
motexclude[which(motexclude$subid == "10637"),"anyExclude"] <- 1

nrow(motexclude[which(motexclude$anyExclude == 1),])

motexclude[which(motexclude$subid == "10711"),] #runs 6, 7, 8 are bad
motexclude[which(motexclude$subid == "11324"),] #a lot of movement in runs 3 and 4, but otherwise very still...
motexclude[which(motexclude$subid == "11336"),] #run 1 has a 14.5 mm FD, run 4 has an 8.5mm movement, but otherwise still


#generate data frame of runs to analyze
featL1Df <- motexclude[which(motexclude$anyExclude==0),] #only retain good runs
featL1Df$runnums <- as.integer(sub("^.*/[^/]+/FEAT_LVL1_run([0-9])\\.feat", "\\1", featL1Df$featRun, perl=TRUE))

#figure out emotion and rew contingency for all runs
run_conditions <- do.call(rbind, lapply(1:nrow(featL1Df), function(i) {
    ##loc <- local({load(file.path(fitDir, paste0(as.character(featL1Df$subid[i]), "_fitinfo.RData"))); environment()})$f #time-clock fit object (load as local var)
    ##loc <- local({load(file.path(fitDir, paste0(as.character("10873"), "_fitinfo.RData"))); environment()})$f #time-clock fit object (load as local var)
    designmat <- file.path(dirname(featL1Df[i,"featRun"]), "designmatrix.RData")
    loc <- local({load(designmat); environment()})$f #load designmatrix.RData    
    data.frame(emotion=loc$run_condition[featL1Df$runnums[i]], contingency=loc$rew_function[featL1Df$runnums[i]]) #vector of emotion and contingency
}))

#build design matrix
featL1Df <- cbind(featL1Df, run_conditions)
featL1Df$emotion <- relevel(featL1Df$emotion, ref="scram")
##featL1Df$model <- sub(paste0(fmriDir, "/.*/mni_5mm_wavelet/fsl_([^/]+)/FEAT.*$"), "\\1", featL1Df$featRun, perl=TRUE)
featL1Df$model <- sub(paste0(fmriDir, "/.*/mni_5mm_wavelet/([^/]+)/FEAT.*$"), "\\1", featL1Df$featRun, perl=TRUE)

save(featL1Df, file="Feat_runinfo_sceptic.RData")

