# fit behavioral data for all participants who completed emo clock in scanner
# setup model-based fMRI GLM analysis for based on fitted data

library(fitclock)
library(foreach)
library(doSNOW)

setwd(file.path(getMainDir(), "clock_analysis", "fmri"))
source("afniValueModel.R")
source("afniTCModel.R")
source("fslValueModel.R")
source("fslTCModel.R")
source("fslSCEPTICModel.R")
source("glm_helper_functions.R")
source("r_glm.R")
if (!file.exists("fmri_fits")) { dir.create("fmri_fits") }
setwd("fmri_fits")

#this version has rtvmax etc., as well as specific variants of entropy
if (file.exists("fmri_sceptic_signals_24basis.RData")) { sceptic <- local({load("fmri_sceptic_signals_24basis.RData"); as.list(environment())}) }

#N.B. in examining initial results from single subject analyses, it is clear that steady state magnetization is not achieved by the first volume acquired
#ICA analysis suggests that it takes up to 6 volumes to reach steady state, and the rel and mean uncertainty maps are being adversely affected by this problem
#because they also start high and decay... Mean uncertainty was consequently soaking up a huge amount of CSF in activation maps.
#Because the first presentation occurs at 8 seconds, it seems fine to drop 6 volumes (6s) 
fit_all_fmri <- function(behavDir, fmriDir, idexpr, dropVolumes=6, usenative=FALSE, model="sceptic", runpar=FALSE, ...) {
  
  behavFiles <- list.files(path=behavDir, pattern=".*tcExport.csv", full.names=TRUE, recursive=FALSE)
  #behavFiles <- grep("11317", behavFiles, value=TRUE) #temporarily here for running a single subject

  if (runpar) {
    require(doSNOW)
    setDefaultClusterOptions(master="localhost") #move away from 10187 to avoid collisions
    clusterobj <- makeSOCKcluster(10)
    registerDoSNOW(clusterobj)
    
    on.exit(try(stopCluster(clusterobj)))
  }

  ll <- foreach(b = iter(behavFiles), .inorder=FALSE, .packages=c("fitclock"), .export=c("truncateRuns", "r_valueModel", "afniTCModel", "fslSCEPTICModel", "runFSLCommand", "sceptic") ) %dopar% {
  #for (b in behavFiles) {
    #example location of file on bea_res, which contains scan date
    #/Volumes/bea_res/Data/Tasks/EmoClockfMRI/Basic/11229/20140521/Raw/fMRIEmoClock_11229_tc_tcExport.csv
    subid <- sub("^.*fMRIEmoClock_(\\d+)_tc_tcExport.csv$", "\\1", b, perl=TRUE)
    #scandate <- sub("^.*/Basic/\\w+/(\\d+)/.*$", "\\1", b, perl=TRUE)
    scandate <- sub("^.*/Basic/\\w+/(\\d+)/.*$", "\\1", b, perl=TRUE) 
    mrfiles <- c() #force clear of mr files over subjects to avoid potential persistence from one subject to the next
    
    ##identify corresponding fmri directory
    mrmatch <- grep(eval(idexpr), list.files(fmriDir, full.names=TRUE), perl=TRUE, value=TRUE)

    if(length(mrmatch) != 1L) {
      warning("Unable to find fMRI directory for subid: ", subid)
      ##next
      return(NULL)
    }
    
    if (usenative==TRUE) {
      expectdir <- "native_nosmooth"
      ##expectfile <- "nfudktm_clock(\\d+).nii.gz"
      expectfile <- "nfudktm_clock[0-9].nii.gz"
    } else {
        ##expectdir <- "mni_5mm_wavelet"
        expectdir <- "mni_5mm_3ddespike"
      ##expectfile <- "nfswudktm_clock(\\d+)_5.nii.gz"
      expectfile <- "nfswudktm_clock[0-9]_5.nii.gz"
    }
    
    if (! file.exists(file.path(mrmatch, expectdir))) {
      warning("Unable to find preprocessed data ", expectdir, " for subid: ", subid)
      ##next
      return(NULL)
    }
    
    ##identify fmri run lengths (4th dimension)
    ##mrfiles <- list.files(mrmatch, pattern=expectfile, full.names=TRUE, recursive=TRUE)
    ##cat(paste0("command: find ", mrmatch, " -iname '", expectfile, "' -ipath '*", expectdir, "*' -type f\n"))
    mrfiles <- system(paste0("find ", mrmatch, " -iname '", expectfile, "' -ipath '*", expectdir, "*' -type f | sort -n"), intern=TRUE)
    mrfiles <- mrfiles[!grepl("(exclude|bbr_noref|old)", mrfiles, ignore.case=TRUE)] #if exclude is in path/filename, then skip
    ##mrrunnums <- as.integer(sub(paste0(".*", expectfile, "$"), "\\1", mrfiles, perl=TRUE))
    mrrunnums <- as.integer(sub(paste0(".*clock(\\d+)_.*$"), "\\1", mrfiles, perl=TRUE))

    ##NB. If we reorder the mrfiles, then the run numbers diverge unless we sort(mrrunnums). Remove for now for testing
    ##mrfiles <- mrfiles[order(mrrunnums)] #make absolutely sure that runs are ordered ascending

    if (length(mrfiles) == 0L) {
      warning("Unable to find any preprocessed MB files in dir: ", mrmatch)
      ##next
      return(NULL)
    }

    ##read number of volumes from NIfTI header
    library(Rniftilib)
    runlengths <- unname(sapply(mrfiles, function(x) { Rniftilib::nifti.image.read(x, read_data=0)$dim[4L] }))
    detach("package:Rniftilib", unload=TRUE) #necessary to avoid dim() conflict with oro.nifti
    
    ##setup clock data subject object for fitting 
    s <- clockdata_subject(subject_ID=subid, dataset=b)
    
    ##create truncated run files to end analysis 12s after last ITI (or big head movement)
    ##also handle removal of 6 volumes from the beginning of each run due to steady state magnetization
    mrdf <- truncateRuns(s, mrfiles, mrrunnums, runlengths, dropVolumes=dropVolumes)
    mrfiles <- mrdf$mrfile_to_analyze
    runlengths <- mrdf$last_vol_analysis
    
    if (model == "tc") {
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
        
        f <- posEps$fit(random_starts=20)
        
        #refit without value carryover
        posEps$carryover_value <- FALSE
        f_nocarryover <- posEps$fit(random_starts=20)
        
      }
    } else if (model == "value") {
      ##delta rule value model (simple)
      vm <- deltavalue_model(clock_data=s, alphaV=0.3, betaV=0.3) #N.B. This matches V matrix from full time-clock algorithm fit.
      f_value <- vm$fit() #estimate learning rate as a free parameter
      
      vm$carryover_value <- FALSE
      f_value_nocarryover <- vm$fit()
    } else if (model == "sceptic") {
      #everything but the ID vector is a 3d matrix subjects x runs x trials
      mats3d <- sort(grep("ids", names(sceptic), value=TRUE, invert=TRUE))
      subj_sceptic <- lapply(sceptic[mats3d], function(mat) {
            mat[subid,,]
          })
      #now have a list where each element is a runs x trials matrix and the elements are the various sceptic signals available
    }
    
    ##setup afni and/or FSL models
    message("About to analyze the following files:")
    print(mrfiles)
    
    if (model=="tc") {
      if (usenative) {
        r_valueModel(f, mrfiles, runlengths, mrrunnums, force=FALSE, dropVolumes=dropVolumes, outdir="rglm_tc_carry")      
      } else {
        #afniValueModel(f_value, mrfiles, runlengths, mrrunnums, run=TRUE, dropVolumes=dropVolumes)
        afniTCModel(f, mrfiles, runlengths, mrrunnums, run=TRUE, dropVolumes=dropVolumes, outdir="afni_tc")
        ##fslValueModel(f_value, mrfiles, runlengths, mrrunnums, run=TRUE, force=FALSE, dropVolumes=dropVolumes)
        ##fslTCModel(f_nocarryover, mrfiles, runlengths, mrrunnums, run=TRUE, force=FALSE, dropVolumes=dropVolumes, outdir="fsl_tc_nocarry") #, f_value=f_value) #hybrid model with EV, RPE+, and RPE- from R-W (carryover value)
        ##fslTCModel(f, mrfiles, runlengths, mrrunnums, run=TRUE, force=FALSE, dropVolumes=dropVolumes, outdir="fsl_tc_nomeanunc")
      }
      
    } else if (model=="sceptic") {
      #for now, trying out a handful of univariate model-based regressors
      #fslSCEPTICModel(subj_sceptic["vmax"], s, mrfiles, runlengths, mrrunnums, run=FALSE, dropVolumes=dropVolumes, ...)
      #fslSCEPTICModel(subj_sceptic["pemax"], s, mrfiles, runlengths, mrrunnums, run=FALSE, dropVolumes=dropVolumes, ...)
      #fslSCEPTICModel(subj_sceptic["vchosen"], s, mrfiles, runlengths, mrrunnums, run=FALSE, dropVolumes=dropVolumes, ...)
      #fslSCEPTICModel(subj_sceptic["ventropy"], s, mrfiles, runlengths, mrrunnums, run=FALSE, dropVolumes=dropVolumes, ...)
      #fslSCEPTICModel(subj_sceptic["vsd"], s, mrfiles, runlengths, mrrunnums, run=FALSE, dropVolumes=dropVolumes, ...)
      #fslSCEPTICModel(subj_sceptic["dauc"], s, mrfiles, runlengths, mrrunnums, run=FALSE, dropVolumes=dropVolumes, ...)
      #fslSCEPTICModel(subj_sceptic["dsd"], s, mrfiles, runlengths, mrrunnums, run=FALSE, dropVolumes=dropVolumes, ...)
  
      ##results from Mean SCEPTIC regressor correlation.pdf indicate that regressors for vchosen, ventropy_decay_matlab, dauc, and pemax are
      ##reasonably uncorrelated. The worst is dauc with vchosen (mean r = -0.31), which makes sense that as learning progresses, chosen values
      ##are higher and there is less residue to decay. These 4 regressors are also of greatest theoretical interest
      subj_sceptic[["dauc"]] <- -1*subj_sceptic[["dauc"]] #invert decay such that higher values indicate greater decay
      fslSCEPTICModel(subj_sceptic[c("vchosen", "ventropy_decay_matlab", "dauc", "pemax")], s, 
                      mrfiles, runlengths, mrrunnums, run=FALSE, dropVolumes=dropVolumes, ...)
    }
  }  
}

#SCEPTIC MMClock Fit
fit_all_fmri(behavDir="/storage/group/mnh5174_collab/temporal_instrumental_agent/clock_task/subjects",
    fmriDir="/storage/group/mnh5174_collab/MMClock/MR_Proc",
    idexpr=expression(subid), ##MMClock/LunaID format: 10637_20140302
    model="sceptic", usepreconvolve=TRUE, parmax1=TRUE, runpar=TRUE) #rescale to 1.0 max

#MMClock fit
# fit_all_fmri(behavDir="/Volumes/bea_res/Data/Tasks/EmoClockfMRI/Basic",
#     fmriDir="/Volumes/Serena/MMClock/MR_Proc",
#     idexpr=expression(paste0(subid, "_", scandate))) ##MMClock/LunaID format: 10637_20140302

## fit_all_fmri(behavDir="/Volumes/bea_res/Data/Tasks/EmoClockfMRI/Basic",
##    fmriDir="/Volumes/Serena/MMClock/MR_Proc",
##    idexpr=expression(paste0(subid, "_", scandate)), usenative=TRUE) ##MMClock/LunaID format: 10637_20140302

## fit_all_fmri(behavDir="/Users/michael/Dropbox/Hallquist_K01/Data/fMRI",
##              fmriDir="/Volumes/Serena/SPECC/MR_Proc",
##              idexpr=expression(paste0(sprintf("%03s", subid), "[A-z]{2}_\\d+"))) #SPECC format: 003aa_15Jul2014
