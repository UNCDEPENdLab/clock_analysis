#look at the design matrices for SCEPTIC data
library(fitclock)
library(abind)
library(ggplot2)
library(dplyr)
setwd(file.path(getMainDir(), "clock_analysis", "fmri"))
source("glm_helper_functions.R")
if (!file.exists("fmri_fits")) { dir.create("fmri_fits") }
setwd("fmri_fits")

source(file.path(getMainDir(), "Miscellaneous", "Global_Functions.R"))

#updated version that used 24 basis functions and also has additional signals (rtvmax etc.)
#if (file.exists("fmri_sceptic_signals_24basis.RData")) { sceptic <- local({load("fmri_sceptic_signals_24basis.RData"); as.list(environment())}) }
if (file.exists("fmri_sceptic_signals_24basis_mmclock_Jun2017.RData")) { sceptic <- local({load("fmri_sceptic_signals_24basis_mmclock_Jun2017.RData"); as.list(environment())}) }

behavDir=file.path(getMainDir(), "temporal_instrumental_agent/clock_task/subjects")
behavFiles <- list.files(path=behavDir, pattern=".*tcExport.csv", full.names=TRUE, recursive=FALSE)

#11282 was on floor of RTs for all runs (also leads to low magnitude of clock regressor)
#11246 had huge and unforgivable head movements (centimeters). RTs were also on the floor and had poor variability.
behavFiles <- grep("11282|11246", behavFiles, value=TRUE, invert=TRUE)
ids <- grep("11282|11246", sceptic$ids, value=TRUE, invert=TRUE)

#for a moment...
#sceptic[["ventropy_decay_matlab"]] <- NULL
#sceptic[["ventropy_fixedlrv_matlab"]] <- NULL
#sceptic[["rtumax"]] <- NULL
#pdf("run_designs_Jun2017.pdf", width=15, height=14)

require(doSNOW)
setDefaultClusterOptions(master="localhost") #move away from 10187 to avoid collisions
clusterobj <- makeSOCKcluster(8)
registerDoSNOW(clusterobj)
#registerDoSEQ()
#for (b in behavFiles[1:2]) {
dmats <- foreach(b=iter(behavFiles), .packages=c("fitclock", "ggplot2"), .export=c("corHeatmap")) %dopar% {
  
  subid <- sub("^.*fMRIEmoClock_(\\d+)_tc_tcExport.csv$", "\\1", b, perl=TRUE)
  mrfiles <- c() #force clear of mr files over subjects to avoid potential persistence from one subject to the next
  
  ##setup clock data subject object for fitting 
  s <- clockdata_subject(subject_ID=subid, dataset=b)
  
  #mats3d <- sort(grep("ids", names(sceptic), value=TRUE, invert=TRUE))
  
  #rtumax and umax are wonky because they can be quite constant in some cases, which blows up various regressors
  #"rtumax", "umax",
  
  #also, for collinearity diagnostics, stick with the ventropy_decay_matlab used in the behavioral paper
  #"ventropy", "ventropy_fixedlrv_matlab",
  
  #mats3d <- sort(c("pemax", "vchosen", "vmax", "rtvmax", "ventropy_decay_matlab", "dauc")) #manual specification
  mats3d <- sort(c("pemax", "vchosen", "ventropy", "dauc")) #manual specification "vmax", "rtvmax", 
  
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
  #normalizations <- rep("none", length(subj_sceptic))
  normalizations <- rep("evtmax_1.0", length(subj_sceptic)) #for now, normalize all to have height 1.0 per event
  for (v in 1:length(subj_sceptic)) {
    thisName <- names(subj_sceptic)[v] #for now, I have a silly divergence where the design matrix function expects a sceptic_ prefix, but the fit object should not have this...
    signals_to_model[v] <- paste0("sceptic_", thisName) #name of regressor in build_design_matrix
    f$sceptic[[thisName]] <- split(subj_sceptic[[v]], row(subj_sceptic[[v]]))
    
    if (thisName %in% c("vauc", "vchosen", "ventropy", "ventropy_fixedlrv_matlab", "ventropy_decay_matlab", "vmax", "rtvmax", "vsd", "umax", "rtumax")) {
      onsets[v] <- "clock_onset"; durations[v] <- "clock_duration"
    } else if (thisName %in% c("pemax", "peauc", "dauc", "dsd")) {
      onsets[v] <- "feedback_onset"; durations[v] <- "feedback_duration"
    }
  }
  
  f$sceptic[["vtime"]] <- sceptic[["vtime_list"]][subid,,]
  signals_to_model[v+1] <- "custom_vtime"
  onsets[v+1] <- NA #not used for custom (but expects same lengths)
  durations[v+1] <- NA #not used for custom
  normalizations[v+1] <- "none"
  
  #visualize the design matrix as a series of boxcars  
  dnocon <- build_design_matrix(fitobj=f, regressors=c("clock", "feedback", signals_to_model), 
      event_onsets=c("clock_onset", "feedback_onset", onsets),
      durations=c("clock_duration", "feedback_duration", durations), 
      normalizations=rep("none", length(signals_to_model) + 3), #no re-normalization in case of no convolution
      baselineCoefOrder=-1, center_values=TRUE, convolve=FALSE,
      dropVolumes=0)
  
#  pdf("run_designs_noconvolve_Dec2016.pdf", width=15, height=14)
#  for (r in 1:length(dnocon$design.convolve)) {
#    g <- visualizeDesignMatrix(dnocon$design.convolve[[r]], outfile=NULL, includeBaseline=FALSE)
#    plot(g+ggtitle(paste("Subj", subid, "Run", r))+theme_bw(base_size=12))
#    
#    plot(corHeatmap(dnocon$design.convolve[[r]], alphaSort=TRUE, tileTextSize=6))
#    
#  }
  #dev.off()
  
  #now build up convolved variants...
  
  dcon <- build_design_matrix(fitobj=f, regressors=c("clock", "feedback", signals_to_model), 
      event_onsets=c("clock_onset", "feedback_onset", onsets),
      durations=c("clock_duration", "feedback_duration", durations), 
      normalizations=c("durmax_1.0", "durmax_1.0", normalizations),
      baselineCoefOrder=-1, center_values=TRUE, convolve=TRUE, convolve_wi_run=TRUE, dropVolumes=0)
  
#  for (r in 1:length(dcon$design.convolve)) {
#    g <- visualizeDesignMatrix(dcon$design.convolve[[r]], outfile=NULL, includeBaseline=FALSE)
#    plot(g+ggtitle(paste("Subj", subid, "Run", r))+theme_bw(base_size=12))
#    plot(corHeatmap(dcon$design.convolve[[r]], alphaSort=TRUE, tileTextSize=6))
#  }
  
#  dcon <- build_design_matrix(fitobj=f, regressors=c("clock", "feedback", signals_to_model), 
#      event_onsets=c("clock_onset", "feedback_onset", onsets),
#      durations=c("clock_duration", "feedback_duration", durations), 
#      normalizations=c("durmax_1.0", "durmax_1.0", normalizations),
#      baselineCoefOrder=2, writeTimingFiles=c("AFNI"), center_values=TRUE, convolve_wi_run=TRUE,
#      runVolumes=runlengths, runsToOutput=mrrunnums, output_directory=timingdir, dropVolumes=dropVolumes)
  
  return(list(file=b, id=subid, dcon=dcon, dnocon=dnocon))
}
#dev.off()

try(stopCluster(clusterobj))

#save(file="n76_sceptic_design_matrices_Jan2017.RData", dmats)
#save(file="n76_sceptic_design_matrices_vtime_Jun2017.RData", dmats)
#save(file="n76_sceptic_design_matrices_vtime_novmax_Jun2017.RData", dmats)
save(file="n76_sceptic_design_matrices_novtime_novmax_Jun2017.RData", dmats)

#aggregate unconvolved and convolved design matrices by: 1) runs within subject, 2) across subjects
rawcollin <- lapply(dmats, function(d) { lapply(d$dcon$collin.raw, "[[", "r") })
rawcollin <- lapply(dmats, function(d) { lapply(d$dcon$collin.raw, "[[", "vif") })

sink("convolved vif sceptic.txt")
convcollin <- lapply(dmats, function(d) { lapply(d$dcon$collin.convolve, "[[", "vif") })
convcollin <- do.call(rbind, lapply(convcollin, function(subj) {
      colMeans(do.call(rbind, subj), na.rm=TRUE)
    }))
rownames(convcollin) <- ids
print(convcollin)
sink()

dmats[[1]]$dnocon$collin.raw$run1 #this is the correlations among decision signals independent of event timing (on trial grid, not on MRI time grid)
dmats[[1]]$dcon$collin.raw$run1 #identical to above

dmats[[1]]$dnocon$collin.convolve$run1 #this is the unconvolved signals laid down onto the MRI time grid, mean centered, and correlated
dmats[[1]]$dcon$collin.convolve$run1 #this is the correlations among *convolved* signals on the MRI time grid

#aggregate the convolved design matrices by subject
alldcon_agg <- lapply(dmats, function(d) {
      allruns <- do.call(abind, list(lapply(d$dcon$collin.convolve, function(run) { run$r }), along=0))
      allaggm <- apply(allruns, c(2,3), mean, na.rm=TRUE)
      allaggsd <- apply(allruns, c(2,3), sd, na.rm=TRUE)
      allaggmin <- apply(allruns, c(2,3), function(cell) { cell[which.min(abs(cell))] }) #least extreme value
      allaggmax <- apply(allruns, c(2,3), function(cell) { cell[which.max(abs(cell))] }) #most extreme value
      allaggmiss <- apply(allruns, c(2,3), function(cell) { sum(is.na(cell)) }) #missing values occur when a regressor is constant. This occurs when the unconvolve parametric values are also constant. So far, rtvmax
      summaries <- abind(allaggm, allaggsd, allaggmin, allaggmax, allaggmiss, along=3, new.names=list(NULL, NULL, statistic=c("mean", "sd", "min", "max", "miss")))
      return(summaries)
    })

#obtain similar summaries for unconvolved designs (still on MRI time grid, just prior to convolution)
alldnocon_agg <- lapply(dmats, function(d) {
      allruns <- do.call(abind, list(lapply(d$dnocon$collin.convolve, function(run) { run$r }), along=0))
      allaggm <- apply(allruns, c(2,3), mean, na.rm=TRUE)
      allaggsd <- apply(allruns, c(2,3), sd, na.rm=TRUE)
      allaggmin <- apply(allruns, c(2,3), function(cell) { cell[which.min(abs(cell))] }) #least extreme value
      allaggmax <- apply(allruns, c(2,3), function(cell) { cell[which.max(abs(cell))] }) #most extreme value
      allaggmiss <- apply(allruns, c(2,3), function(cell) { sum(is.na(cell)) }) #missing values occur when a regressor is constant. This occurs when the unconvolve parametric values are also constant. So far, rtvmax
      summaries <- abind(allaggm, allaggsd, allaggmin, allaggmax, allaggmiss, along=3, new.names=list(NULL, NULL, statistic=c("mean", "sd", "min", "max", "miss")))
      return(summaries)
    })

#aggregate mean matrices across subjects for display
alldcon_groupmean <- do.call(abind, list(lapply(alldcon_agg, function(d) {
      d[,,"mean"]
    }), along=0))

#melt for display
pdf("Mean SCEPTIC regressor correlation vtime.pdf", width=17, height=13)
gmelt <- reshape2::melt(alldcon_groupmean, varnames=c("ID", "V1", "V2"))
gmeans <- gmelt %>% filter(V1 != V2) %>% group_by(V1, V2) %>% summarize(value=mean(value))
ggplot(filter(gmelt, V1 != V2), aes(x=value)) + geom_density() + facet_grid(V1 ~ V2) + theme_bw(base_size=16) +
    geom_vline(xintercept=0) + geom_vline(data=gmeans, mapping=aes(xintercept=value), color="blue") +
    geom_text(data=gmeans, mapping=aes(x=value + 0.2, y=8, label=round(value, 2)))
dev.off()

alldcon_groupmax <- do.call(abind, list(lapply(alldcon_agg, function(d) {
              d[,,"max"]
            }), along=0))

gmelt <- reshape2::melt(alldcon_groupmax, varnames=c("ID", "V1", "V2"))
ggplot(gmelt, aes(x=value)) + stat_density() + facet_grid(V1 ~ V2)


#also look at regressor scaling differences
alldcon_regheights <- lapply(dmats, function(d) {
      #return a runs x signals x min/max matrix
      mall <- do.call(abind, list(lapply(d$dcon$design.convolve, function(run) { t(sapply(run, range)) }),
              along=0, new.names=list(NULL, NULL, stat=c("min", "max"))))
      
      maggruns <- apply(mall, c(2,3), mean)
      return(maggruns)
    })

#now combine over subjects to look at distribution
alldcon_regheights_group <- do.call(abind, list(alldcon_regheights, along=0))
alldcon_regheights_regrange <- apply(alldcon_regheights_group, c(1,2), function(dd) { dd[2] - dd[1]})
#try transform of dauc
alldcon_regheights_regrange[,c("sceptic_dauc")] <- sqrt(alldcon_regheights_regrange[,c("sceptic_dauc")])
gmelt <- reshape2::melt(alldcon_regheights_regrange, varnames=c("ID", "regressor"))

pdf("SCEPTIC Regressor Scaling Differences.pdf", width=12, height=8)
ggplot(gmelt, aes(x=value)) + geom_histogram(bins=10) + facet_wrap(~regressor, scales="free")
dev.off()

behavFiles[which(alldcon_regheights_regrange[,"clock"] < 0.6)] #outlier on clock magnitude. Yes, this is due to bad subject 11282 (floor of RTs)
behavFiles[which(alldcon_regheights_regrange[,"sceptic_pemax"] > 75)] #high outliers on PE magnitudes: 10717, 11328, 11366
behavFiles[which(alldcon_regheights_regrange[,"sceptic_ventropy_decay_matlab"] > 2.5)] #high outliers on sceptic entropy: 11314, 11329

#z score convolved regressors
alldcon_regheights_normalized <- lapply(dmats, function(d) {
      #return a runs x signals x min/max matrix
      mall <- do.call(abind, list(lapply(d$dcon$design.convolve, function(run) { xx <- data.frame(lapply(run, scale)); t(sapply(xx, range)) }), #z score each regressor in each run
              along=0, new.names=list(NULL, NULL, stat=c("min", "max"))))
      
      maggruns <- apply(mall, c(2,3), mean)
      
      rescaled <- lapply(d$dcon$design.convolve, function(run) {
            df <- data.frame(lapply(1:ncol(run), function(reg) {
                      run[,reg]/sd(run[,reg])
                    }))
            names(df) <- names(run)
            df
          })
      
      return(list(rescaled_design=rescaled, maggruns=maggruns))
    })

alldcon_regheights_regrange_normalized <- apply(do.call(abind, list(alldcon_regheights_normalized, along=0)), c(1,2), function(dd) { dd[2] - dd[1]})
gmelt <- reshape2::melt(alldcon_regheights_regrange_normalized, varnames=c("ID", "regressor"))

pdf("SCEPTIC Regressor Scaling Differences Normalized.pdf", width=12, height=8)
ggplot(gmelt, aes(x=value)) + geom_histogram(bins=15) + facet_wrap(~regressor, scales="free")
dev.off()


#z score convolved regressors at the subject level, not run level
#if indeed parameter scaling differences are to blame for differences in regressor magnitude, this should correct the problem
#without potentially undoing meaningful differences in scaling between runs (e.g., due to contingency)
alldcon_regheights_subj_normalized <- lapply(dmats, function(d) {
      #a bit of a hack, but compute SD in each run, then divide by the mean SD
      #should probably be some weighted mean by run length
      runsds <- sapply(d$dcon$design.convolve, function(run) { sapply(run, sd) })
      avgsd <- apply(runsds, 1, mean)
      
      rescaled <- lapply(d$dcon$design.convolve, function(run) {
            df <- data.frame(lapply(1:length(avgsd), function(reg) {
                      run[,reg]/avgsd[reg]
                    }))
            names(df) <- names(run)
            df
          })
      
      #lapply(rescaled, function(run) { lapply(run, function(col) { sd(col) })}) #check looks basically right (unit variance on average)
#
#      
#      dconcat <- do.call(rbind, d$dcon$design.convolve)
#      colm <- colMeans(dconcat)
#      colsd <- apply(dconcat, 2, sd)
#      browser()
#      rescaled <- sapply(1:ncol(dconcat), function(col) { (dconcat[,col] - colm[col])/colsd[col]})
#      apply(rescaled, 2, mean)
#      apply(rescaled, 2, sd)
      
      mall <- do.call(abind, list(lapply(rescaled, function(run) { t(sapply(run, range)) }),
              along=0, new.names=list(NULL, NULL, stat=c("min", "max"))))
      
      maggruns <- apply(mall, c(2,3), mean)
      return(list(rescaled_design=rescaled, maggruns=maggruns))
    })

alldcon_regheights_regrange_subj_normalized <- apply(do.call(abind, list(lapply(alldcon_regheights_subj_normalized, "[[", "maggruns"), along=0)), c(1,2), function(dd) { dd[2] - dd[1]})
gmelt <- reshape2::melt(alldcon_regheights_regrange_subj_normalized, varnames=c("ID", "regressor"))

pdf("SCEPTIC Regressor Scaling Differences Subject-Level Normalized.pdf", width=12, height=8)
ggplot(gmelt, aes(x=value)) + geom_histogram(bins=15) + facet_wrap(~regressor, scales="free")
dev.off()

str(dmats)
orig <- lapply(dmats, function(d) { d$dcon$design.convolve })
str(alldcon_regheights_normalized)
alldcon_regheights_subj_normalized


#between subjects rank stability in ordering of regressor ranges
#looks like the z scoring really shuffles up the rank ordering across subjects of parametric regressors
sapply(1:ncol(alldcon_regheights_regrange), function(col) { cor(alldcon_regheights_regrange[,col], alldcon_regheights_regrange_normalized[,col], use="pairwise.complete.obs")})
#[1] -0.67899198  0.30186849 -0.46698857  0.28872240 -0.36188949 -0.17916680 -0.07682001 -0.23570318

sapply(1:ncol(alldcon_regheights_regrange), function(col) { cor(alldcon_regheights_regrange[,col], alldcon_regheights_regrange_subj_normalized[,col], use="pairwise.complete.obs")})
#[1] -0.6854808  0.2837072 -0.5149371  0.5415190 -0.1321394 -0.1552493 -0.2681861 -0.2552397

#run-level versus subject-level normalization is pretty similar
sapply(1:ncol(alldcon_regheights_regrange), function(col) { cor(alldcon_regheights_regrange_normalized[,col], alldcon_regheights_regrange_subj_normalized[,col], use="pairwise.complete.obs")})
#[1] 0.9833575 0.9978105 0.9474152 0.9195217 0.9194900 0.9806524 0.8813005 0.9603892


#what about transforming the DAUC variable before convolution to normalize while not changing ranks?

