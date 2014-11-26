#second-level model for MMY3 data
#generate group analysis (what FSL sometimes calls 3rd level) based on emotions of first level runs and multiple sessions per subject
#also drop runs where level of motion was unacceptable

setwd(file.path(getMainDir(), "clock_analysis", "fmri"))
source("glm_helper_functions.R")

##fmriDir <- "/Volumes/Serena/MMClock/MR_Proc/10873_20140918/mni_5mm_wavelet"
fmriDir <- "/Volumes/Serena/MMClock/MR_Proc"
fitDir <- file.path(getMainDir(), "clock_analysis", "fmri", "fmri_fits")

featRuns <- list.files(path=fmriDir, pattern="FEAT_LVL1_run[0-9].feat", include.dirs=TRUE, recursive=TRUE, full.names=TRUE)

#flag runs with more than 10% volumes with FD 0.9mm or greater
#find fd.txt files corresponding to each FEAT run
fdFiles <- sub(paste0(fmriDir, "(.*)/fsl_.*/FEAT_LVL1_run([0-9])\\.feat"), paste0(fmriDir, "\\1/clock\\2/motion_info/fd.txt"), featRuns, perl=TRUE)

#identify matching truncated 4d files (created by model_clock_fmri.R) since we are only concerned about movement within the run proper (not dead volumes)  
niTruncFiles <- Sys.glob(sub(paste0(fmriDir, "(.*)/fsl_.*/FEAT_LVL1_run([0-9])\\.feat"), paste0(fmriDir, "\\1/clock\\2/nfsw*drop*.nii.gz"), featRuns, perl=TRUE))

#identify how many volumes were dropped at the beginning
dropVolumes <- as.integer(sub("^.*_drop(\\d+).*.nii.gz$", "\\1", niTruncFiles, perl=TRUE))

##some runs are not truncated (e.g., if they go all the way to 350 vols and the scanner stops)
##this results in NAs here
truncLengths <- as.integer(sub("^.*_trunc(\\d+).nii.gz$", "\\1", niTruncFiles, perl=TRUE))
nafiles <- which(is.na(truncLengths))

library(Rniftilib)
#need to add back the dropped volumes since everything below assumes that the trunc length is the final volume in the original time series
truncLengths[nafiles] <- unname(sapply(nafiles, function(x) { Rniftilib::nifti.image.read(niTruncFiles[x], read_data=0)$dim[4L] })) + dropVolumes[nafiles]
detach("package:Rniftilib", unload=TRUE) #necessary to avoid dim() conflict with oro.nifti

library(plyr)
motexclude <- ldply(1:length(fdFiles), function(i) {
      fd <- read.table(fdFiles[i], header=FALSE)$V1
      fd <- fd[(dropVolumes[i]+1):truncLengths[i]] #only include volumes within run
      propSpikes_0p9 <- sum(as.integer(fd > 0.9))/length(fd)
      ##if (is.na(propSpikes_0p9[1])) browser()
      spikeExclude <- if (propSpikes_0p9 > .15) 1 else 0
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

motexclude[which(motexclude$anyExclude == 1),]

#generate data frame of runs to analyze
featL1Df <- motexclude[which(motexclude$anyExclude==0),] #only retain good runs
featL1Df$runnums <- as.integer(sub("^.*/fsl_.*/FEAT_LVL1_run([0-9])\\.feat", "\\1", featL1Df$featRun, perl=TRUE))

#figure out emotion and rew contingency for all runs
run_conditions <- do.call(rbind, lapply(1:nrow(featL1Df), function(i) {
    loc <- local({load(file.path(fitDir, paste0(as.character(featL1Df$subid[i]), "_fitinfo.RData"))); environment()})$f #time-clock fit object (load as local var)
    ##loc <- local({load(file.path(fitDir, paste0(as.character("10873"), "_fitinfo.RData"))); environment()})$f #time-clock fit object (load as local var)
    data.frame(emotion=loc$run_condition[featL1Df$runnums[i]], contingency=loc$rew_function[featL1Df$runnums[i]]) #vector of emotion and contingency
}))

#build design matrix
featL1Df <- cbind(featL1Df, run_conditions)
featL1Df$emotion <- relevel(featL1Df$emotion, ref="scram")
featL1Df$model <- sub(paste0(fmriDir, "/.*/mni_5mm_wavelet/fsl_([^/]+)/FEAT.*$"), "\\1", featL1Df$featRun, perl=TRUE)

save(featL1Df, file="Feat_runinfo.RData")

#generate dummy DV to get model.matrix from lm
featL1Df$dummy <- rnorm(nrow(featL1Df), 0, 1)

#need to pull in age here

## designmat <- lm(dummy ~ emotion + subid, data=featL1Df)

## mm <- model.matrix(designmat)
## write.table(mm, file="fsl_LVL3_emo_design.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)
## cat(featL1Df$featRun, sep="\n", file="feat_runlist.txt")
## print(dimnames(mm)[2]) #column names (for FSL)

## library(lsmeans)

##generate per-subject second-level FE analyses to get contrasts of interest for group analysis
#saved file does not include relevel for ref of scram (hence copy from above -- redundant)



run_feat_lvl2 <- function(featL1Df, run=TRUE, force=FALSE) {
    allFeatRuns <- list()
    require(plyr)
    require(parallel)
    d_ply(featL1Df, .(subid, model), function(subdf) {
        subdf <- subdf[order(subdf$runnums),] #verify that we have ascending runs
        dummy <- lm(runnums ~ emotion, subdf)
        mm <- model.matrix(dummy)
#      library(lsmeans)
#      v <- lsmeans(dummy, list(pairwise ~ emotion))
#      v[[1]]@linfct #condition means
#      v[[2]]@linfct #emo diffs (pairwise)
#      lsm <- lsmeans(dummy, "emotion")
#      contrast(lsm, method="eff")@linfct #cell versus gm
#      
#      #grand mean contrast
#      #unbalanced design, so need to set weights based on relative frequency
        props <- prop.table(table(subdf$emotion))
        gm_coef <- c(1, props[2:length(props)]) #1 for intercept/reference, proportions for other cells
#      
#      summary(glht(dummy, linfct=matrix(gm_coef, 1)))
        
        #generate and run lvl2 for this subject
        fsfTemplate <- readLines(file.path(getMainDir(), "clock_analysis", "fmri", "feat_lvl2_clock_template.fsf"))

        #depending on lower-level model (e.g., TC versus value, will have different number of copes to compute

                         
        
        if (subdf$model[1] == "value") {
            ##value model has 5 copes: clock onset, feedback onset, ev, rpe+, rpe-
            fsfTemplate <- c(fsfTemplate,
                             "# Number of lower-level copes feeding into higher-level analysis",
                             "set fmri(ncopeinputs) 5",
                             "# Use lower-level cope 1 for higher-level analysis",
                             "set fmri(copeinput.1) 1",
                             "# Use lower-level cope 2 for higher-level analysis",
                             "set fmri(copeinput.2) 1",
                             "# Use lower-level cope 3 for higher-level analysis",
                             "set fmri(copeinput.3) 1",
                             "# Use lower-level cope 4 for higher-level analysis",
                             "set fmri(copeinput.4) 1",
                             "# Use lower-level cope 5 for higher-level analysis",
                             "set fmri(copeinput.5) 1"
                             )
        } else if (subdf$model[1] =="tc_nocarry") {
            ##TC model has 7 copes: clock onset, feedback onset, ev, rpe+, rpe-, mean_unc, rel_unc
            fsfTemplate <- c(fsfTemplate,
                             "# Number of lower-level copes feeding into higher-level analysis",
                             "set fmri(ncopeinputs) 7",
                             "",
                             "# Use lower-level cope 1 for higher-level analysis",
                             "set fmri(copeinput.1) 1",
                             "",
                             "# Use lower-level cope 2 for higher-level analysis",
                             "set fmri(copeinput.2) 1",
                             "",
                             "# Use lower-level cope 3 for higher-level analysis",
                             "set fmri(copeinput.3) 1",
                             "",
                             "# Use lower-level cope 4 for higher-level analysis",
                             "set fmri(copeinput.4) 1",
                             "",
                             "# Use lower-level cope 5 for higher-level analysis",
                             "set fmri(copeinput.5) 1",
                             "# Use lower-level cope 6 for higher-level analysis",
                             "set fmri(copeinput.6) 1",
                             "",
                             "# Use lower-level cope 7 for higher-level analysis",
                             "set fmri(copeinput.7) 1"
                             )
        } else { warning("unable to match model: ", subdf$model[1]); return(NULL) }
        
        if (nrow(subdf) != 8) { warning("can't handle less than 8 runs at the moment! ", subdf$subid[1]); return(NULL) }
        #search and replace within fsf file for appropriate sections
        #.OUTPUTDIR. is the feat output location
        
        thisTemplate <- fsfTemplate
        thisTemplate <- gsub(".OUTPUTDIR.", file.path(dirname(subdf$featRun[1L]), "FEAT_LVL2"), thisTemplate, fixed=TRUE)
        for (i in 1:nrow(subdf)) {
          thisTemplate <- gsub(paste0(".INPUT", i, "."), subdf$featRun[i], thisTemplate, fixed=TRUE)
        }
        
        #EV1 is intercept (scrambled is reference)
        #EV2 is emofear
        #EV3 is emohappy
        for (i in 1:nrow(mm)) {
          for (j in 1:ncol(mm)) {
            thisTemplate <- gsub(paste0(".I", i, "EV", j, "."), mm[i,j], thisTemplate, fixed=TRUE)
          }
        }
        
        #grand mean contrast (depends on number of runs of each emotion)
        thisTemplate <- gsub(".GMCOL1.", gm_coef[1], thisTemplate, fixed=TRUE)
        thisTemplate <- gsub(".GMCOL2.", gm_coef[2], thisTemplate, fixed=TRUE)
        thisTemplate <- gsub(".GMCOL3.", gm_coef[3], thisTemplate, fixed=TRUE)
        
        featFile <- file.path(dirname(subdf$featRun[1L]), "FEAT_LVL2.fsf")
        if (file.exists(featFile) && force==FALSE) { return(NULL) } #skip re-creation of FSF and do not run below unless force==TRUE 
        cat(thisTemplate, file=featFile, sep="\n")      
        
        allFeatRuns[[featFile]] <<- featFile
      })

  if (run == TRUE) {
    cl_fork <- makeForkCluster(nnodes=8)
    runfeat <- function(fsf) {
      runname <- basename(fsf)
      runFSLCommand(paste("feat", fsf), stdout=file.path(dirname(fsf), paste0("feat_stdout_", runname)), stderr=file.path(dirname(fsf), paste0("feat_stderr_", runname)))
    }
    clusterApply(cl_fork, allFeatRuns, runfeat)
    stopCluster(cl_fork)
  } 
  
}

run_feat_lvl2(featL1Df, run=TRUE, force=TRUE)
    
    
    
    
    
    
    
    
#soooo.... fsl is very slow and I don't think it handles the multiple sessions per subject very well (since we have dummy regressors for subject)
#what about using weighted lmer approach ala David Paulsen?

#loop over copes and varcopes and load into 4-d array (x y z subject)
library(oro.nifti)
library(doSNOW)
library(foreach)
library(iterators)

load(file="Feat_runinfo.RData")

#about 1.29 GB in RAM per array
exampCope <- readNIfTI(file.path(featL1Df$featRun[1], "stats", "cope1.nii.gz"), reorient=FALSE)@.Data
mask <- readNIfTI("/Users/michael/standard/mni_icbm152_nlin_asym_09c/mni_icbm152_t1_tal_nlin_asym_09c_mask_2.3mm.nii", reorient=FALSE)@.Data

#generate a group cope and varcope image for all subjects
for (copenum in 1:5) {
  cope_array <- array(NA_real_, c(dim(exampCope),nrow(featL1Df)))
  varcope_array <- array(NA_real_, c(dim(exampCope),nrow(featL1Df)))
  
  for (i in 1:nrow(featL1Df)) {
    cope <- readNIfTI(file.path(featL1Df$featRun[i], "stats", paste0("cope", copenum, ".nii.gz")), reorient=FALSE)@.Data
    varcope <- readNIfTI(file.path(featL1Df$featRun[i], "stats", paste0("varcope", copenum, ".nii.gz")), reorient=FALSE)@.Data
    
    #apply mask
    cope[which(mask==0)] <- NA_real_
    varcope[which(mask==0)] <- NA_real_
    
    cope_array[,,,i] <- cope
    varcope_array[,,,i] <- varcope
    
  }
  
  save(cope_array, varcope_array, file=paste0("allsubj_cope", copenum, ".RData"))
}


#now have 4d files for copes and varcopes in same order as featL1Df
#loop over voxels that are within the mask

#only consider unmasked voxels
vindices <- which(mask==1, arr.ind=TRUE)

copenum <- 1
load(paste0("allsubj_cope", copenum, ".RData"))

#each 2-D array is ~283MB
cope_voxels <- apply(cope_array, 4, '[', vindices) #obtain 154985 x 240 matrix (voxels x sessions)
varcope_voxels <- apply(varcope_array, 4, '[', vindices) #obtain 154985 x 240 matrix (voxels x sessions)

rm(cope_array)
rm(varcope_array)

#maskData <- maskData[1:1000,] #for testing

#setDefaultClusterOptions(master="localhost", port=10290) #move away from 10187 to avoid collisions
#clusterobj <- makeSOCKcluster(8)
#registerDoSNOW(clusterobj)

#as a memory-saving hack, cbind cope and varcope matrices together so that we can iterate by row (and avoid exporting full data to each worker)
vanalyze <- cbind(cope_voxels, varcope_voxels)
rownames(vanalyze) <- 1:nrow(vanalyze) #allows for lookup of mask indices inside dopar loop

#parallel loop over mask indices
fourDMat <- foreach(vox=iter(vanalyze, by="row"), .inorder=TRUE, .packages=c("lme4", "pbkrtest"), .combine=rbind, 
        .noexport=c("cope_array", "varcope_array", "cope_voxes", "varcope_voxels")) %do% {
      
      #each v is 1 x nsessions*2 matrix
      #first half of columns are betas, second half are variances
      vdf <- cbind(featL1Df[,c("subid", "emotion", "contingency", "runnums")], data.frame(b=vox[1, 1:(ncol(vox)/2)], v=vox[1,(ncol(vox)/2 + 1):ncol(vox)]))
      
      #for some edge voxels, some subjects to not have any values. In this case, make a note of the voxel coordinates and analyze those with data present... 
      if (any(vdf$v==0.0)) {
        cat("Vmiss: ", paste(vindices[rownames(vox),], collapse=" "), "indices: ", paste(which(vdf$v==0.0), collapse=","), "\n", file="misslog", append=TRUE)
        vdf <- subset(vdf, v > 0)
      }
      
      browser()
      require(pbkrtest)
      
      #normalize reciprocal of variances to have sd = 1.0 (and mean of ~0.5-0.7)
      #vdf$weight <- scale(1/vdf$v, center=FALSE)
      
      #ldiv <- lmer(b ~ 1 + emotion + (1 | subid), data=vdf, weights=1/(v/1000000))
      l <- lmer(b ~ 1 + emotion + (1 | subid), data=vdf, weights=weight)
      #l <- lmer(b ~ 1 + emotion + (1 | subid), data=vdf, weights=1/(v/10000))
      #l <- lmer(b ~ 1 + emotion + (1 | subid), data=vdf)
      
      #Kenward-Roger approximation to df
      #note: this only produces one value, not a per-parameter estimate
      #df.KR <- get_ddf_Lb(l, fixef(l))
      # get p-values from the t-distribution using the t-values and approximated
      # degrees of freedom
      #coefs$p.KR <- 2 * (1 - pt(abs(coefs$t.value), df.KR))
      #coefs
      
      ltest <- summary(lmerTest::lmer(b ~ 1 + emotion + (1 | subid), data=vdf))$coefficients[,"df"]
      l2 <- lmerTest::lmer(b ~ 1 + emotion + (1 | subid), data=vdf, weights=1/v)
      #mod <- lmer(b ~ 1 + emotion + (1 | subid), data=vdf, weights=1/v)
      rl <- rlmer(b ~ 1 + emotion + (1 | subid), data=vdf)
      wt_huber <- getME(rl, "w_e")
      
      #nlme needs ~v, not 1/v, to get variance weighting correct
      test <- nlme::lme(b ~ 1 + emotion, random=~1|subid, data=vdf, weights=~v)
      ranef(test)
      
      #ltest <- lmer(b ~ 1 + emotion + (1 | subid), data=vdf, weights=wt_huber)
      
      af <- afex::mixed(b ~ 1 + emotion + (1 | subid), data=vdf, method="KR", test.intercept=TRUE, per.parameter=".")
      
      af <- mixed(b ~ 1 + emotion + (1 | subid), data=vdf, method="KR", test.intercept=TRUE, per.parameter=".")
      
      library(varComp)
      
      test2 <- varComp::varComp(b ~ 1 + emotion, vdf, ~subid)
      
      vout <- ltest$coefficients["(Intercept)", c("Estimate", "t value", "df", "Pr(>|t|)")]
      return(vout)      
      #afex::mixed breaks down for intercept only model
      #voxDf$dummy <- 0
      #ltest <- mixed(v ~ 1 + dummy + (1 | LunaID), voxDf, method="KR", test.intercept=TRUE, per.parameter=".")      
    }


colnames(fourDMat) <- c("b", "t", "df", "p")

#allocate memory
#nifti.image.alloc.data(out_hdr)

#recycle each row of the mask indices 4 times
#this is ugly
#maskIndicesMod <- cbind(mi[rep(1:nrow(mi), each=ncol(fourDMat)),], dim4=1:ncol(fourDMat))
#fourDUnmat <- gdata::unmatrix(fourDMat, byrow=TRUE) #convert into a vector by row (b, t, df, p) to match indices above

#insert results into 4d array
#system.time(out_hdr[maskIndicesMod] <- fourDUnmat) #as.matrix(fourDMat))

#aha! Use abind to build 4-d stats mat
#lapply over columns of matrix (take advantage of data.frame inheriting list)
out_hdr@.Data <- abind(lapply(data.frame(fourDMat), function(col) {
          m <- array(0, dim(mask)[1:3]) #empty 3d matrix matching dims of images
          m[mi] <- col
          return(m)
        }), along=4)



#from here: http://mindingthebrain.blogspot.com/2014/02/three-ways-to-get-parameter-specific-p.html
require(pbkrtest)
# get the KR-approximated degrees of freedom
df.KR <- get_ddf_Lb(m.sem, fixef(m.sem))
# get p-values from the t-distribution using the t-values and approximated
# degrees of freedom
coefs$p.KR <- 2 * (1 - pt(abs(coefs$t.value), df.KR))
coefs
