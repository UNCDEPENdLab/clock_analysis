# fit behavioral data for all participants who completed emo clock in scanner
# setup model-based fMRI GLM analysis for based on fitted data

library(fitclock)

setwd(file.path(getMainDir(), "clock_analysis", "fmri"))
source("afniValueModel.R")
source("afniTCModel.R")
source("fslValueModel.R")
source("fslTCModel.R")
source("glm_helper_functions.R")
source("r_glm.R")
if (!file.exists("fmri_fits")) { dir.create("fmri_fits") }
setwd("fmri_fits")

#N.B. in examining initial results from single subject analyses, it is clear that steady state magnetization is not achieved by the first volume acquired
#ICA analysis suggests that it takes up to 6 volumes to reach steady state, and the rel and mean uncertainty maps are being adversely afected by this problem
#because they also start high and decay... Mean uncertainty was consequently soaking up a huge amount of CSF in activation maps.
#Because the first presentation occurs at 8 seconds, it seems fine to drop 6 volumes (6s) 
fit_all_fmri <- function(behavDir, fmriDir, idexpr, dropVolumes=6) {
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
        mrmatch <- grep(eval(idexpr), list.files(fmriDir, full.names=TRUE), perl=TRUE, value=TRUE)
        
        if(length(mrmatch) != 1L) {
            warning("Unable to find fMRI directory for subid: ", subid)
            next
        }
        #if (! file.exists(file.path(mrmatch, "MBclock_recon"))) {
        if (! file.exists(file.path(mrmatch, "mni_5mm_wavelet"))) {
        ##if (! file.exists(file.path(mrmatch, "native_nosmooth"))) {
            #warning("Unable to find preprocessed data mni_5mm_wavelet for subid: ", subid)
            warning("Unable to find preprocessed data native_nosmooth for subid: ", subid)
            next
        }

        ##identify fmri run lengths (4th dimension)
        mrfiles <- list.files(mrmatch, pattern="nfswudktm_clock\\d+_5.nii.gz", full.names=TRUE, recursive=TRUE)
        ##mrfiles <- list.files(mrmatch, pattern="nfudktm_clock\\d+.nii.gz", full.names=TRUE, recursive=TRUE)
        mrrunnums <- as.integer(sub(".*nfswudktm_clock(\\d+)_5.nii.gz$", "\\1", mrfiles, perl=TRUE))
        ##mrrunnums <- as.integer(sub(".*nfudktm_clock(\\d+).nii.gz$", "\\1", mrfiles, perl=TRUE))
        mrfiles <- mrfiles[order(mrrunnums)] #make absolutely sure that runs are ordered ascending
        
        if (length(mrfiles) == 0L) {
            warning("Unable to find any preprocessed MB files in dir: ", mrmatch)
            next
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
            
            #refit without value carryover
            posEps$carryover_value <- FALSE
            f_nocarryover <- posEps$fit(random_starts=5)
            
            ##delta rule value model (simple)
            vm <- deltavalue_model(clock_data=s, alphaV=0.3, betaV=0.3) #N.B. This matches V matrix from full time-clock algorithm fit.
            f_value <- vm$fit() #estimate learning rate as a free parameter

            vm$carryover_value <- FALSE
            f_value_nocarryover <- vm$fit()
            
            save(f_value, f_value_nocarryover, f, f_nocarryover, s, incr_fit, file=paste0(subid, "_fitinfo.RData"))
        }
        
        ##setup afni and FSL models
        message("About to analyze the following files:")
        print(mrfiles)

        #afniValueModel(f_value, mrfiles, runlengths, mrrunnums, run=TRUE, dropVolumes=dropVolumes)
        afniTCModel(f, mrfiles, runlengths, mrrunnums, run=TRUE, dropVolumes=dropVolumes, outdir="afni_tc_carry_filt")
        ##fslValueModel(f_value, mrfiles, runlengths, mrrunnums, run=TRUE, force=FALSE, dropVolumes=dropVolumes)
        ##fslTCModel(f_nocarryover, mrfiles, runlengths, mrrunnums, run=TRUE, force=FALSE, dropVolumes=dropVolumes, outdir="fsl_tc_nocarry") #, f_value=f_value) #hybrid model with EV, RPE+, and RPE- from R-W (carryover value)
        
        ##r_valueModel(f, mrfiles, runlengths, mrrunnums, force=FALSE)
    }
    
}


#MMClock fit
##fit_all_fmri(behavDir="/Volumes/bea_res/Data/Tasks/EmoClockfMRI/Basic",
##        fmriDir="/Volumes/Serena/MMClock/MR_Proc",
##        idexpr=expression(paste0(subid, "_", scandate))) ##MMClock/LunaID format: 10637_20140302

 fit_all_fmri(behavDir="/Users/michael/Dropbox/Hallquist_K01/Data/fMRI",
              fmriDir="/Volumes/Serena/SPECC/MR_Proc",
              idexpr=expression(paste0(sprintf("%03s", subid), "[A-z]{2}_\\d+"))) #SPECC format: 003aa_15Jul2014
