setwd(file.path(getMainDir(), "clock_analysis", "fmri"))
source("glm_helper_functions.R")
##load("Feat_runinfo_sceptic_mmclock.RData")
load("Feat_runinfo_sceptic_mmclock.RData")
##load("Feat_runinfo_sceptic_specc.RData")

#generate dummy DV to get model.matrix from lm
featL1Df$dummy <- rnorm(nrow(featL1Df), 0, 1)
#designmat <- lm(dummy ~ emotion + subid, data=featL1Df)

## mm <- model.matrix(designmat)
## write.table(mm, file="fsl_LVL3_emo_design.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)
## cat(featL1Df$featRun, sep="\n", file="feat_runlist.txt")
## print(dimnames(mm)[2]) #column names (for FSL)

## library(lsmeans)

##generate per-subject second-level FE analyses to get contrasts of interest for group analysis
#saved file does not include relevel for ref of scram (hence copy from above -- redundant)

run_feat_lvl2 <- function(featL1Df, run=TRUE, force=FALSE, ncpus=8) {
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
    #fsfTemplate <- readLines(file.path(getMainDir(), "clock_analysis", "fmri", "fsf_templates", "feat_lvl2_clock_template.fsf"))
    fsfTemplate <- readLines(file.path(getMainDir(), "clock_analysis", "fmri", "fsf_templates", "feat_lvl2_clock_template_runtrend.fsf"))
    #depending on lower-level model (e.g., TC versus value, will have different number of copes to compute

    
    if (subdf$model[1] == "value") {
      ##value model has 5 copes: clock onset, feedback onset, ev, rpe+, rpe-
      ncopes <- 5      
    } else if (subdf$model[1] == "tc_nocarry") {
      ##TC model has 7 copes: clock onset, feedback onset, ev, rpe+, rpe-, mean_unc, rel_unc
      ncopes <- 7
    } else if (subdf$model[1] =="tc_nomeanunc") {
      ##TC no mean unc model has 6 copes: clock onset, feedback onset, ev, rpe+, rpe-, rel_unc
      ncopes <- 6
    } else if (subdf$model[1] %in% c("sceptic_vchosen_ventropy_decay_matlab_dauc_pemax_preconvolve", "sceptic_vchosen_ventropy_dauc_pemax_preconvolve")) {
      ncopes <- 6
    } else if (subdf$model[1] == "sceptic_vchosen_ventropy_dauc_pemax_vtime_preconvolve") {
      ncopes <- 7
    } else {
      ##assuming single single parametric regressor model (clock onset, feedback onset, parameter)
      ##warning("unable to match model: ", subdf$model[1]); return(NULL)
      ##message("Assuming single parametric modulator")
      ncopes <- 3
    }

    fsfTemplate <- c(fsfTemplate,
      "# Number of lower-level copes feeding into higher-level analysis",
      paste0("set fmri(ncopeinputs) ", ncopes))
    
    for (n in 1:ncopes) {
      fsfTemplate <- c(fsfTemplate,
      paste0("# Use lower-level cope ", n, " for higher-level analysis"),
      paste0("set fmri(copeinput.", n, ") 1"))      
    }
    
    if (nrow(subdf) != 8) { warning("can't handle less than 8 runs at the moment! ", subdf$subid[1]); return(NULL) }
    #search and replace within fsf file for appropriate sections
    #.OUTPUTDIR. is the feat output location
    
    thisTemplate <- fsfTemplate
    ##thisTemplate <- gsub(".OUTPUTDIR.", file.path(dirname(subdf$featRun[1L]), "FEAT_LVL2"), thisTemplate, fixed=TRUE)
    thisTemplate <- gsub(".OUTPUTDIR.", file.path(dirname(subdf$featRun[1L]), "FEAT_LVL2_runtrend"), thisTemplate, fixed=TRUE)
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

    featOutDir <- file.path(dirname(subdf$featRun[1L]), "FEAT_LVL2_runtrend.gfeat")
    featFile <- file.path(dirname(subdf$featRun[1L]), "FEAT_LVL2_runtrend.fsf")
    if (file.exists(featOutDir) && force==FALSE) { return(NULL) } #skip re-creation of FSF and do not run below unless force==TRUE 
    cat(thisTemplate, file=featFile, sep="\n")      
    
    allFeatRuns[[featFile]] <<- featFile
  })

  print(allFeatRuns)
  if (run == TRUE) {
    cl_fork <- makeForkCluster(nnodes=ncpus)
    runfeat <- function(fsf) {
      runname <- basename(fsf)
      runFSLCommand(paste("feat", fsf), stdout=file.path(dirname(fsf), paste0("feat_stdout_", runname)), stderr=file.path(dirname(fsf), paste0("feat_stderr_", runname)))
    }
    clusterApply(cl_fork, allFeatRuns, runfeat)
    stopCluster(cl_fork)
  } 
  
}

run_feat_lvl2(featL1Df, run=TRUE, force=FALSE, ncpus=20)
print(warnings())
