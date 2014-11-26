#library(fmri)
#library(fitclock)
#load a subject
#setwd(file.path(getMainDir(), "clock_analysis", "fmri", "fmri_fits"))


r_valueModel <- function(f_value, mrfiles, runlengths, mrrunnums, force=FALSE, njobs=8) {
  require(abind)
  require(fmri)
  require(oro.nifti)
  require(pracma)
  
  rglmdir <- file.path(normalizePath(file.path(dirname(mrfiles[1L]), "..")), "rglm_feedback_0dur") #note: normalizePath will fail to evaluate properly if directory does not exist (e.g., dir not created yet)
  if (file.exists(rglmdir) && force==FALSE) {
    message("rglm_feedback_0dur directory already exists. ", rglmdir, ". Skipping subject")
    return(NULL)
  }
  
  cat("rglmdir create: ", rglmdir, "\n")
  dir.create(rglmdir, showWarnings=FALSE) #one directory up from a given clock run
  ## regressornames <- c("clock", "ev", "feedback", "rpe_pos", "rpe_neg")
  ## d_value <- build_design_matrix(fitobj=f_value, regressors=regressornames, 
  ##     event_onsets=c("clock_onset", "clock_onset", "feedback_onset", "feedback_onset", "feedback_onset"), 
  ##     durations=c("clock_duration", "clock_duration", 0, 0, 0),
  ##     normalizations=c("durmax_1.0", "evtmax_1.0", "evtmax_1.0", "evtmax_1.0", "evtmax_1.0"), baselineCoefOrder=3,
  ##     runVolumes=runlengths, runsToOutput=mrrunnums)
  regressornames <- c("clock", "feedback", "rel_uncertainty", "mean_uncertainty", "ev", "rpe_pos", "rpe_neg")
  
  d_value <- build_design_matrix(fitobj=f_value, regressors=regressornames, 
      event_onsets=c("clock_onset", "feedback_onset", "clock_onset", "clock_onset", "feedback_onset", "feedback_onset", "feedback_onset"), 
      durations=c("clock_duration", "feedback_duration", "clock_duration", "clock_duration", "feedback_duration", 0, 0),
      normalizations=c("durmax_1.0", "durmax_1.0", "evtmax_1.0", "evtmax_1.0", "evtmax_1.0", "evtmax_1.0", "evtmax_1.0"), baselineCoefOrder=3,
      runVolumes=runlengths, runsToOutput=mrrunnums)
  
  ##read in run mask!
  ##because glms are estimated on single runs, do not assume that a Tmin across runs makes any sense. Haven't coregistered runs to each other.
  #generateRunMask(mrfiles, outdir=rglmdir, outfile="runmask")
  #brainmask <- readNIfTI(file.path(rglmdir, "runmask.nii.gz"), reorient=FALSE)
  #maskIndices <- which(brainmask == 0.0, arr.ind=TRUE)
  
  #motion parameters
  motpcs <- pca_motion(mrfiles, runlengths, motion_parfile="motion.par", numpcs=3)
  
  contrastmat <- c()
  
  for (r in 1:length(regressornames)) {
    thisCon <- rep(0.0, length(regressornames))
    thisCon[r] <- 1.0
    contrastmat <- rbind(contrastmat,thisCon)
    rownames(contrastmat)[r] <- paste("me", regressornames[r], sep=".") 
  }
  
  require(foreach)
  require(doSNOW)
  
  setDefaultClusterOptions(master="localhost", port=10290)
  clusterobj <- makeSOCKcluster(njobs)
  registerDoSNOW(clusterobj)
  
  on.exit(stopCluster(clusterobj))
  
  res <- foreach(f=1:length(mrfiles), .packages=c("oro.nifti", "fmri", "abind"), .export="runFSLCommand") %dopar% {
    #for (f in 1:length(mrfiles)) {
    
    #mrdat <- read.NIFTI(mrfiles[f], setmask=FALSE)
    mrdat <- readNIfTI(mrfiles[f], reorient=FALSE)
    
    #copy NIfTI header information for stats output (fmri write.NIFTI giving bad grid)
    mrout <- mrdat
    mrout@.Data <- array(0) #delete data and preserve header alone
    mrout@dim_[5] <- 5 #stats output includes b, se, t, 1-p, seg
    mrdat <- oro2fmri(mrdat, setmask=FALSE) #convert to fmri object for analysis
    
    ##use tmin=0 as mask
    runFSLCommand(paste0("fslmaths ", mrfiles[f], " -Tmin -bin ", rglmdir, "/runmask", f))#, fsldir="/usr/local/ni_tools/fsl")
    brainmask <- readNIfTI(file.path(rglmdir, paste0("runmask", f, ".nii.gz")), reorient=FALSE)    
    mrdat$mask <- as.logical(brainmask@.Data) #set mask directly
    
    #apply run mask
    #this generates a massive array (162994378 x 4)... probably more effective to loop and assign
    #mask4d <- cbind(repmat(maskIndices, dim(mrdat)[4], 1), rep(1:dim(mrdat)[4], each=nrow(maskIndices))) 
    #mrdat[mask4d] <- NA_real_
    
    #masking the functional image one volume at a time is very slow
    #but using one big 4d array leads to a huge RAM overhead
    #middle ground is to divide into blocks of 25 volumes to be masked
#    blocksize <- 25
#    nvols <- dim(mrdat)[4L]
#    blocks <- split(1:nvols, ceiling((1:nvols)/blocksize))
#    
#    tim <- system.time(lapply(blocks, function(b) {
#              mask4d <- cbind(repmat(maskIndices, length(b), 1), rep(b, each=nrow(maskIndices)))
#              mrdat[mask4d] <- NA_real_
#            }))
    
#    time <- system.time(for (vol in 1:dim(mrdat)[4L]) {
#      mrdat[cbind(maskIndices, vol)] <- NA_real_
#    })
    
    for (v in 1:nrow(contrastmat)) {
      dmat <- as.matrix(cbind(d_value$design.convolve[[f]], motpcs$motion_pcs_runs[[f]]))
      result <- fmri.lm(mrdat, dmat, actype="smooth", contrast=c(contrastmat[v,], rep(0, ncol(dmat) - length(contrastmat[v,]))), keep="all")
      
      #adaptive smoothing
      sm <- fmri.smooth(result, adaptation="fullaws")
      
      #summary.fmridata(sm)
      #local p-values based on adaptive smoothing
      pv <- fmri.pvalue(sm, mode="local")
      #plot(pv, anatomic=anat, maxpvalue=.05)
      
      #structural segmentation smoothing approach
      #the $segm field has values 0 (null), 1 (positive activation), and -1 (negative activation) for each voxel basd on segmentation.
      sm_seg <- tryCatch(fmri.smooth(result, adaptation="segment", alpha=.05), error=function(e) { warning("segmentation algorithm failed: ", mrfiles[f]); return(array(-5, dim=dim(sm$cbeta))) } )
      
      ##old fmri output code
      #b_se_t_p_seg <- abind(sm$cbeta, sm$var, sm$cbeta/sqrt(sm$var), 1 - pv$pvalue, sm_seg$segm, along=4)
      #niftiHead <- sm$header
      #niftiHead$slicecode <- "" #slight mishap in converting oro2fmri
      #niftiHead$xyztunits <- "\n"
      #niftiHead$dimension[5] <- 5 #b, se, t, 1-p, seg
      #write.NIFTI(b_se_t_p_seg, niftiHead, filename=file.path(rglmdir, paste0("run", mrrunnums[f], "_", rownames(contrastmat)[v])))
      
      ##new code based on oro.nifti
      mrout@.Data <- abind(sm$cbeta, sm$var, sm$cbeta/sqrt(sm$var), 1 - pv$pvalue, sm_seg$segm, along=4)
      mrout@cal_min <- min(mrout) #need to have min and max vals in header that match data for oro.nifti to succeed
      mrout@cal_max <- max(mrout)
      writeNIfTI(mrout, filename=file.path(rglmdir, paste0("run", mrrunnums[f], "_", rownames(contrastmat)[v])))
    }
    
    rm(mrdat)
    return("done")
  }
  
  ##      stopCluster(clusterobj)
#need to get an emotion-modulated matrix at some point, although this is between runs
  
}
