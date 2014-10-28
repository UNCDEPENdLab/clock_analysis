#cutting down on redundancy across glm setup scripts

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


setup_mrfiles <- function(s, mrfiles) {
    #identify 
    
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
                                ##use either the last volume of the task or the volume before the earliest bad movement 
                                truncLength <- min(last_vol, (min(badfd) - 1))
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
    
}

pca_motion <- function(mrfiles, runlengths, motion_parfile="motion.par") {
    #based on a vector of mr files to be analyzed, look 
    #concat motion parameters
    motion_concat <- do.call(rbind, lapply(1:length(mrfiles), function(i)  {
                        mot <- read.table(file.path(dirname(mrfiles[i]), motion_parfile), col.names=c("r.x", "r.y", "r.z", "t.x", "t.y", "t.z"))
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
    
    
}