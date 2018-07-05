#plot the fitted data from VBA
library(R.matlab)
library(ggplot2)
library(plyr)
library(dplyr)
library(reshape2)
library(entropy)
options(dplyr.print_max = 1000)
options(max.print = 1000)
# setwd("~/Data_Analysis/temporal_instrumental_agent/clock_task/vba_fmri")
setwd("~/code/temporal_instrumental_agent/clock_task/vba_fmri")

#fit <- readMat("/Users/michael/Data_Analysis/temporal_instrumental_agent/clock_task/vba_fmri/posterior_states_decay_nomultisession.mat")
#fit <- readMat("/Users/michael/Data_Analysis/temporal_instrumental_agent/clock_task/vba_fmri/posterior_states_decay_nomultisession_constrain0p025.mat")
#fit <- readMat("/Users/michael/Data_Analysis/temporal_instrumental_agent/clock_task/vba_fmri/posterior_states_decay_nomultisession_constrain0p0125_niv.mat")
# fit <- readMat("~/Data_Analysis/temporal_instrumental_agent/clock_task/vba_fmri/posterior_states_decay_nomultisession_psfixed0p0125_k24.mat") #24 basis pre-niv (PLoS Comp Bio submission)
fit <- readMat("~/code/temporal_instrumental_agent/clock_task/vba_fmri/posterior_states_decay_nomultisession_psfixed0p0125_k24.mat") #24 basis pre-niv (PLoS Comp Bio submission)

#fit <- readMat("~/Data_Analysis/temporal_instrumental_agent/clock_task/vba_fmri/posterior_states_decay_nomultisession_specc_decay_psfixed0p0125_k24.mat") #specc n=94 dataset

#basis <- readMat("/Users/michael/Data_Analysis/temporal_instrumental_agent/clock_task/vba_fmri/sceptic_fmri_basis_setup.mat")
# basis <- readMat("~/Data_Analysis/temporal_instrumental_agent/clock_task/sceptic_fmri_basis_setup_k24_p0125.mat")
basis <- readMat("~/code/temporal_instrumental_agent/clock_task/sceptic_fmri_basis_setup_k24_p0125.mat")

source("clock_functions.R")

#pull the uncertainties from kalman_uv_sum. Make sure this is not predictive of RT swings
#udata <- readMat("/Users/michael/Google Drive/skinner/projects_analyses/SCEPTIC/subject_fitting/uncertainty_results/sigma_matrix.mat")
#udata <- readMat(file.path(GoogleDriveDir(), "skinner/projects_analyses/SCEPTIC/subject_fitting/uncertainty_results/multisession_uncertainty_fixed_uv_kalman_uv_sum.mat"))

#this one should be correct for resetting U at run boundaries
#udata <- readMat(file.path(GoogleDriveDir(), "skinner/projects_analyses/SCEPTIC/subject_fitting/uncertainty_results/updated_multisession_u_matrix.mat"))

udata <- readMat("~/Box Sync/skinner/projects_analyses/SCEPTIC/subject_fitting/uncertainty_results/multisession_uncertainty_fixed_uv_kalman_uv_sum.mat")

#sigtrials <- udata[[1]]["kalman.uv.sum",,][[1]]["sigma.all.trials",,][[1]] #hideous syntax, but that's how we get it from the .mat!
#sigtrials <- udata[[1]]["fixed.uv",,][[1]]["sigma.all.trials",,][[1]] #hideous syntax, but that's how we get it from the .mat!
#the other data in the fit objects tend to be subjects x runs x basis functions x trials
#sigtrials <- aperm(sigtrials, c(3,2,1))
#sigrestruct <- array(aperm(sigtrials, c(3,2,1)), dim=c(76, 8, 24, 50)) #this doesn't get the ordering quite right...

#m <- melt(sigtrials, varnames=c("basis", "trial", "rowID")) %>% arrange(rowID, trial, basis) %>% select(-trial)#rename(trial_abs=trial)
#m$run <- rep(1:8, each=50*24) #1:8 numbering
#m$trial <- rep(1:50, each=24) #1:50 numbering as elsewhere

#sigrestruct <- acast(m, rowID ~ run ~ basis ~ trial, value.var="value") #tada! 

#sigrestruct[1,1,1:24, 1:10] #all basis functions for 10 trials
#sigtrials[1:24,1:10, 1]

ids <- gsub(".*CORRECT_(\\d+)_fixed_decay.*", "\\1", unlist(fit$fitfiles), perl=TRUE)

#plot RTs for clock data
#allData <- getClockGroupData(path="~/Data_Analysis/temporal_instrumental_agent/clock_task/subjects/SPECC") #for SPECC
# allData <- getClockGroupData(path="~/Data_Analysis/temporal_instrumental_agent/clock_task/subjects") #for MMY3
allData <- getClockGroupData(path="~/code/temporal_instrumental_agent/clock_task/subjects") #for MMY3

allData$timestep <- plyr::round_any(allData$rt, 100)/100 #1:40 coding


#these .mat structures are all MMClock Y3, no equivalent for SPECC at the moment.
#newEntropy <- readMat(file.path(GoogleDriveDir(), "skinner/projects_analyses/SCEPTIC/subject_fitting/entropy_analysis/val_based_shannon_H.mat"))
#newEntropy <- readMat(file.path(GoogleDriveDir(), "skinner/projects_analyses/SCEPTIC/subject_fitting/entropy_analysis/H_values_decay_and_fixed_multisession.mat"))
#newEntropy <- readMat(file.path(GoogleDriveDir(), "skinner/projects_analyses/SCEPTIC/subject_fitting/entropy_analysis/unisession_decay_and_fixed_rand_priors.mat"))
#newEntropy <- readMat(file.path(GoogleDriveDir(), "skinner/projects_analyses/SCEPTIC/subject_fitting/entropy_analysis/unisession_decay_and_fixed_rand_priors_full.mat"))

#value distribution entropy under decay and fixed LR V 
#ventropy_decay_matlab <- melt(newEntropy$decay.H, varnames=c("rowID", "trial"), value.name="entropyH")
#ventropy_fixedlrv_matlab <- melt(newEntropy$fixed.H, varnames=c("rowID", "trial"), value.name="entropyFixed")

#ventropy_decay_matlab <- ventropy_decay_matlab %>% arrange(rowID, trial)
#ventropy_decay_matlab$run <- rep(1:8, each=50)
#ventropy_decay_matlab$trial <- 1:50

#ventropy_fixedlrv_matlab <- ventropy_fixedlrv_matlab %>% arrange(rowID, trial)
#ventropy_fixedlrv_matlab$run <- rep(1:8, each=50)
#ventropy_fixedlrv_matlab$trial <- 1:50

#with(allData, table(ID, rewFunc, emotion))
#allData$rewarded <- factor(allData$score > 0, levels=c(TRUE, FALSE), labels=c("Reward", "Omission"))
#pdf("ClockfMRI_AllSubjRTs_withMax.pdf", width=14, height=8)
#for (s in split(allData, allData$ID)) {
#  sm <- reshape2::melt(s[,c("trial", "rt", "emotion", "rewFunc", "bestRewRT", "bestEVRT", "run", "rewarded")], id.vars=c("trial", "emotion", "rewFunc", "run", "rewarded"))
#  g <- ggplot(sm, aes(x=trial, y=value, color=variable)) + geom_line() + geom_point(data=subset(sm, variable=="rt"), aes(shape=rewarded), size=2) + facet_grid(emotion ~ rewFunc) + ggtitle(s$ID[1L]) + ylab("RT") + xlab("Trial") +
#      geom_text(aes(x=45, y=3700, label=run, color=NULL), color="black", family="Helvetica", size=6) + theme_bw(base_size=21) + scale_shape("Outcome") + scale_color_discrete("") + coord_cartesian(ylim=c(0,4000))
#  print(g)
#}
#dev.off()


#plot of V

#trialsToPlot <- c(1, 5, 10, 20, 30, 40, 49)
#Vsub <- fit$V[,,,trialsToPlot]
#
##PE and D are right-shifted by 1 such that PE(1) falls in the second column since it is used to compute V(2) = V(1) + alpha*(r(1) - V(1)) 
#PEsub <- fit$PE[,,,trialsToPlot + 1]
#Dsub <- fit$D[,,,trialsToPlot + 1]
#
#toplot <- list(V=Vsub, PE=PEsub, D=Dsub)
#
#for (v in 1:length(toplot)) {
#
#  #multiply against basis to get distribution/representation
#  Vexpand <- aaply(toplot[[v]], c(1,2,4), function(b) {
#        v = outer(b, rep(1, dim(basis$tvec)[2])) * basis$gaussmat #use vector outer product to replicate weight vector
#        vfunc <- colSums(v)
#        return(vfunc)
#      })
#  
#  pdf(paste(names(toplot)[v], "unfolding.pdf"), width=12, height=10)
#  for (s in 1:dim(Vexpand)[1]) {
#    subj <- Vexpand[s,,,]
#    m <- melt(subj, varnames=c("run", "trial", "timestep"))
#    m$ID <- factor(ids[s])
#    m$trialFac <- factor(m$trial, levels=1:length(trialsToPlot), labels=paste0("i = ", trialsToPlot))
#    m$trial <- plyr::mapvalues(m$trial, 1:length(trialsToPlot), trialsToPlot)
#    rtdf <- subset(allData, ID == ids[s] & trial %in% trialsToPlot, select=c(run, trial, rt, score)) #[,c("ID", "run", "trial", "rt", "score")], by=c("ID", "run", "trial"))
#    rtdf$timestep <- plyr::round_any(rtdf$rt, 100)/100
#    rtdf$trialFac <- factor(rtdf$trial, levels=trialsToPlot, labels=paste0("i = ", trialsToPlot))
#    g <- ggplot(m, aes(x=timestep, y=value)) + geom_hline(yintercept=0, color="blue") + geom_line() + geom_text(data=rtdf, aes(label=score, y=0)) + facet_grid(run ~ trialFac, scales="free_y") + ggtitle(ids[s])
#    plot(g)
#  }
#  dev.off()
#  
#}

#note: the problems with "positive decay" and other weird effects in decay and PE were due to an error in the use of an incorrect refspread
#and the max sigma of the eligibility of 1.0 (instead of this being truly free). I have fixed these 8Jun2016 and plots now look correct.
#hence, code below only applies to previous problem.

#diagnose a few weird things
#positive decay:
#SUB17: id 11178, run 7, trial 20
#V7 <-  fit$V[17,7,,]
#
#V7ex <- apply(V7, c(2), function(b) {
#      v = outer(b, rep(1, dim(basis$tvec)[2])) * basis$gaussmat #use vector outer product to replicate weight vector
#      vfunc <- colSums(v)
#      return(vfunc)
#    })
#
#D7 <-  fit$D[17,7,,]
#
#D7ex <- apply(D7, c(2), function(b) {
#      v = outer(b, rep(1, dim(basis$tvec)[2])) * basis$gaussmat #use vector outer product to replicate weight vector
#      vfunc <- colSums(v)
#      return(vfunc)
#    })
#
#PE7 <-  fit$PE[17,7,,]
#
#PE7ex <- apply(PE7, c(2), function(b) {
#      v = outer(b, rep(1, dim(basis$tvec)[2])) * basis$gaussmat #use vector outer product to replicate weight vector
#      vfunc <- colSums(v)
#      return(vfunc)
#    })
#
#par(mfrow=c(3,1))
#plot(1:40, V7ex[,20], type="l")
#plot(1:40, D7ex[,21], type="l")
#plot(1:40, PE7ex[,21], type="l")


#generate regressors for PE, value, and decay at each trial
tofmri <- list(V=fit$V, PE=fit$PE, D=fit$D)#, U=sigrestruct) #dropping U for SPECC
#tofmri <- list(V=fit$V, PE=fit$PE, D=fit$D, U=sigrestruct) #MMY3
fmriexpanded <- list()
for (v in 1:length(tofmri)) {  
  #multiply against basis to get distribution/representation
  fmriexpanded[[ names(tofmri)[v] ]] <- aaply(tofmri[[v]], c(1,2,4), function(b) {
        v = outer(b, rep(1, dim(basis$tvec)[2])) * basis$gaussmat #use vector outer product to replicate weight vector
        vfunc <- colSums(v)
        return(vfunc)
      })
  
  dimnames(fmriexpanded[[ names(tofmri)[v] ]]) <- list(ID=ids, run=1:8, trial=1:50, timestep=1:40) #list(ID=ids, run=1:8, trial=1:50, timestep=1:40)
}

#fmriexpanded now has arrays that are subjects x runs x trials x timesteps
#make regressors of these

#1) Maximum value
options(error=recover)
vmax <- apply(fmriexpanded[["V"]], c(1,2,3), function(trial) {
      if (all(is.na(trial)) || #this happens when the subject didn't complete a given trial at all 
          sd(trial, na.rm=TRUE) < .01) { #essentially no variability, so hard to say there's a max
        return(NA)
      } else {
        return(max(trial, na.rm=TRUE))
      } 
    })

#1.5) RT of max value
rtvmax <- apply(fmriexpanded[["V"]], c(1,2,3), function(trial) {
      if (all(is.na(trial)) || sd(trial) < .01) { #essentially no variability, so hard to say there's a max 
        return(NA) 
      } else {
        return(which.max(trial))
      } 
    })

#2) AUC of value function
vauc <- apply(fmriexpanded[["V"]], c(1,2,3), function(trial) {
      if (all(is.na(trial)) || sd(trial) < .01) { #essentially no variability, so hard to say there's a max 
        return(NA) 
      } else {
        return(sum(trial))
      } 
    })

#3) Entropy of value function
#essentially looking at this as a problem of "how much risk to I incur choosing one of these 40 options" (time bins)
#when all are equal, risk is maximal. Here, we arbitrarily carve up the value dimension into 20 bins so that there isn't
#too much sparsity. Divide by log(20) since this is the maximum entropy -- rescales entropy to be 0..1

library(parallel)
#cl <- makeCluster(getOption("cl.cores", 4))

ventropy <- apply(fmriexpanded[["V"]], c(1,2,3), function(trial) {
#ventropy <- parApply(cl=cl, X=fmriexpanded[["V"]], c(1,2,3), function(trial) {
      if (all(is.na(trial)) || sd(trial) < .01) { #essentially no variability, so hard to say there's a max 
        return(NA) 
      } else {
        #normalize entropy by discretizing into 20 bins
        dd <- discretize(trial, numBins=20)
        entropy.empirical(dd)/log(20) #ML-based Shannon entropy
      } 
    })

#plot ventropy
#mventropy = melt(ventropy, varnames=c("subject", "run", "trial"), value.name="entropy")
#
#mvsplit = split(mventropy, mventropy$subject)
#pdf("ventropy_unfolding.pdf", width=8, height=11)
#for (vdf in mvsplit) {
#  g = ggplot(vdf, aes(x=trial, y=entropy)) + geom_line(size=1.3) + ggtitle(paste0("Subject: ", vdf$subject[1])) + theme_bw(base_size=12) + facet_wrap(~run, ncol=1)
#  plot(g)
#}
#dev.off()

#4) standard deviation of value function
vsd <- apply(fmriexpanded[["V"]], c(1,2,3), function(trial) {
      if (all(is.na(trial)) || sd(trial) < .01) { #essentially no variability, so hard to say there's a max 
        return(NA) 
      } else {
        return(sd(trial))
      } 
    })

#5) value of chosen option
mv <- melt(fmriexpanded[["V"]], varnames=c("ID", "run", "trial", "timestep"))

##NB. this line was implemented before ID was added directly (and in the correct order) to fmriexpanded elements
##also, it causes serious problems (terrible idea) when ids are out of order, as in SPECC data (MMY3 appears okay)
##mv$ID <- factor(mv$subj, levels=sort(unique(mv$subj)), labels=ids)

#um <- melt(fmriexpanded[["U"]], value.name="U")
#um$subject <- factor(um$ID, levels=sort(unique(um$ID)), labels=1:76)
#dcast(ID ~ run + trial + timestep, value.name="U")

#6) max uncertainty
#problem is that subjects can't really respond fast enough to choose very early, so U will always go there.
#pull back to say min of 10?.. this gets pretty hacky, but going with 6 and 6 for now
# rtumax <- apply(fmriexpanded[["U"]], c(1,2,3), function(trial) {
#       if (all(is.na(trial)) || sd(trial) < 1) { #essentially no variability, so hard to say there's a max 
#         return(NA) 
#       } else {
#         minoffset <- 6 #first five time steps (of 40) deemed ineligible
#         maxoffset <- 6 #last three time steps also ineligible
#         #browser()
#         trial <- trial[minoffset:(length(trial)-maxoffset)]
#         return(which.max(trial) + minoffset - 1)
#       } 
#     })
# 
# umax <- apply(fmriexpanded[["U"]], c(1,2,3), function(trial) {
#       if (all(is.na(trial)) || sd(trial) < 1) { #essentially no variability, so hard to say there's a max 
#         return(NA) 
#       } else {
#         minoffset <- 4 #first three time steps (of 40) deemed ineligible
#         maxoffset <- 4 #last three time steps also ineligible
#         trial <- trial[minoffset:(length(trial)-maxoffset)]
#         return(max(trial))
#       } 
#     })


#there are two RTs less than 50ms that round down to zero, leading to their being missing in the merge
#this may be right -- these may not be valid rts... but subjects do get feedback
#there are also 119 RTs above 4s, representing non-response. Again, subjects get feedback (always zero) but is a bit different
#for now:
allData$timestep[which(allData$timestep==0)] <- 1
allData$timestep[which(allData$timestep==41)] <- 40 #overflow on RT

#30400 rows = 76 subjects x 8 runs x 50 trials
#mv now has one row per subject and trial representing the value at the chosen option.

#problem with a handful of subjects. Look at subject 18, run 8
#dplyr::anti_join(mv, allData, by=c("ID", "run", "trial", "timestep"))

#from SPECC
#tt1 <- filter(mv, ID=="18")
#tt2 <- filter(allData, ID=="18") %>% arrange(run, trial)
#mtt <- merge(tt1, tt2, by=c("run", "trial", "timestep"), all.y=TRUE) %>% arrange(run, trial)

#this check out... no concerns here
# > with(tt1, value[run==1 & trial==1 & timestep==8])
# [1] 0
# > with(tt1, value[run==1 & trial==2 & timestep==5])
# [1] 7.680692
# > with(tt1, value[run==1 & trial==3 & timestep==14])
# [1] 0.4172741
# > with(tt1, value[run==1 & trial==4 & timestep==9])
# [1] 16.08931
# with(tt1, value[run==1 & trial==5 & timestep==20])
# [1] 5.716574e-07
# 
# > head(mtt, n=5)[,c("run", "trial", "timestep", "value")]
# run trial timestep        value
# 1   1     1        8 0.000000e+00
# 2   1     2        5 7.680692e+00
# 3   1     3       14 4.172741e-01
# 4   1     4        9 1.608931e+01
# 5   1     5       20 5.716574e-07


mv_chosen <- merge(mv, allData, by=c("ID", "run", "trial", "timestep"))
with(mv_chosen, table(ID, run))

#mv_chosen <- subset(mv_chosen, !ID==74)
#recast into 3d array to match other regressors.
vchosen_suspicious <- acast(mv_chosen, ID ~ run ~ trial, value.var="value")

#I have no idea what's happening in the recast... let's just try a brute force rearrangement
#but it is generating very weird missingness (264 missing values interspersed in strange places)
#in addition, when I look at the value estimates themselves, they seem out of order in vchosen, but good in the brute force setup

#on further analysis, it appears that there are 264 missing values (legitimately)
#but in the original there were 465 missing values. As best I can tell, this links back to the bad mv$ID line above (that inadvertently messed up IDs)
ids_arranged <- as.character(sort(as.numeric(ids)))
vchosen <- array(NA, dim=c(length(ids), 8, 50), dimnames=list(ID=ids_arranged, run=1:8, trial=1:50))
for (s in 1:length(ids)) {
  for (r in 1:8) {
    df <- subset(mv_chosen, ID==ids_arranged[s] & run==r)
    df <- df %>% arrange(trial) #ensure that the data are ordered by trial
    if (nrow(df) > 0) { vchosen[s,r,1:nrow(df)] <- df$value
    } else { message("no data for id:", ids_arranged[s], ", run: ", r) }
  }
}

#further checks reveal no substantial problems in the acast version
all.equal(vchosen, vchosen_suspicious)
#[1] "Attributes: < Component “dimnames”: names for target but not for current >"
sum(vchosen-vchosen_suspicious, na.rm=T)
#[1] 0


#compute an evolving value signal
mmdf <- melt(fmriexpanded[["V"]])
mmdf <- merge(mmdf, allData %>% dplyr::rename(timestep_chosen=timestep), by=c("ID", "run", "trial")) %>% arrange(ID, run, trial, timestep)

#so now we have the timestep chosen merged in. Seems we would just need to filter to observations before chosen value
mmdf <- mmdf %>% filter(timestep <= timestep_chosen)

#if the EV is held for 100ms, then we add timestep - 1 to the trial onset
mmdf$onset <- mmdf$clock_onset + (mmdf$timestep-1)/10
mmdf$duration = 0.1

mmdf <- select(mmdf, ID, run, trial, onset, value, duration) #timestep, timestep_chosen,


#prototype: checks out

#vtime <- split(mmdf, list(mmdf$ID, mmdf$run))
#xx <- vtime[[1]]
#xr1 <- filter(xx, run==1)
#
#ev_convolve <- fmri.stimulus(scans=3000, values=xr1$value, times=xr1$onset, durations=xr1$duration, rt=0.1)
#ev_noconvolve <- fmri.stimulus(scans=3000, values=xr1$value, times=xr1$onset, durations=xr1$duration, rt=0.1, convolve=FALSE)
#par(mfrow=c(2,1))
#plot(ev_convolve, type="l")
#plot(ev_noconvolve, type="l")


#creates a nested list structure... but not multidimensional
#mylist <- dlply(mmdf, .(ID), dlply, .(run), dlply, .(trial), identity)
subjs <- sort(unique(mmdf$ID))
runs <- sort(unique(mmdf$run))
trials <- sort(unique(mmdf$trial))
vtime_list <- array(list(), dim=c(length(subjs), length(runs), length(trials)), dimnames=list(ID=subjs, run=runs, trial=trials))

#brute force
for (s in as.character(subjs)) {
  for (r in as.character(runs)) {
    for (t in as.character(trials)) {
      vtime_list[[s,r,t]] <- filter(mmdf, ID==s & run==r & trial==t)
    }
  }
}





#this should give us what we need for convolution
#just need to reshape to subjects x runs x trials array (and to have a 4th dimension of timestep)
#or maybe just break with this convention and allow a custom pass through for full manual specification


#now we need to lay these down onto the MRI time grid...
#it seems like essentially we just need to support passing in a parametric regressor that can have many values per trial
#and have the convolution lay down all of those timesteps wrt trial onset

#vdf <- dcast(mmdf, ID + run + trial ~ timestep)



#convert to 4000ms time series that can be convolved
#conv <- vdf %>% group_by(ID, run) %>%
#    do({
#          vmat <- as.matrix(.[,as.character(1:40)]) # trials x timesteps matrix of evaluated value function
#          vmat <- vmat - mean(vmat) #subtract off run-wise mean before convolution
#          .[,as.character(1:40)] <- vmat
#          return(.)
#        }) %>% rowwise() %>%
#    do({
#          v <- unlist(.[as.character(1:40)])
#          #this is basically correct
#          vconv <- fmri.stimulus(scans=50, times=seq(0, 3900, 100)/1000, durations=0.1, values=v, rt=0.1, convolve = TRUE)
#          browser()
#        })


#but really, we want run-wise convolved regressors, which means concatenating trials

#For PE- AND D-based statistics, these have been right-shifted by 1 position to line up the inputs and outputs in SCEPTIC VBA fitting
#Consequently, the PE of trial 1 is actually in position 2, etc. And, for now, the PE and D of the last trial (50) is not collected/estimated
#Thus, flip the first trial (value of 0) to the end (currently unspecified/not estimated)

rightshift1 <- function(mat3d) {
  #need aaply to keep dimensions correct (otherwise collapses runs and trials)
  plyr::aaply(mat3d, 1, function(subj) {
        v <- gdata::unmatrix(subj, byrow=TRUE)
        vlag <- Hmisc::Lag(v, -1)
        matrix(vlag, nrow=nrow(subj), ncol=ncol(subj), byrow=TRUE)
      })
}

#PE (need to right shift to get PE(i)
pemax <- apply(fmriexpanded[["PE"]], c(1,2,3), function(trial) {
      #note that some trials have essentially no PE
      if (all(is.na(trial))) {
        return(NA)
      } else if (min(trial) >= 0){
        return(max(trial))
      } else {
        return(min(trial))
      }
    })

pemax <- rightshift1(pemax)

#AUC of PE (seems like it will be highly correlated with max PE, so probably not informative)... yeah, looks like r = .99 for the first run!
peauc <- apply(fmriexpanded[["PE"]], c(1,2,3), sum)
peauc <- rightshift1(peauc)

#AUC of decay
dauc <- apply(fmriexpanded[["D"]], c(1,2,3), sum)
dauc <- rightshift1(dauc)

#SD of decay, which essentially captures how much the amount of decay varies by response option
#Should have a high inverse relationship with AUC, I think.
dsd <- apply(fmriexpanded[["D"]], c(1,2,3), sd)
dsd <- rightshift1(dsd)

#copy across luna ids
#ventropy_decay_matlab_3d <- acast(ventropy_decay_matlab, rowID ~ run ~ trial, value.var="entropyH")
#ventropy_fixedlrv_matlab_3d <- acast(ventropy_fixedlrv_matlab, rowID ~ run ~ trial, value.var="entropyFixed")

#dimnames(ventropy_decay_matlab_3d) <- dimnames(ventropy_fixedlrv_matlab_3d) <- dimnames(rtvmax)

save(vmax, rtvmax, vauc, vsd, ventropy, vchosen, vtime_list, pemax, peauc, dauc, dsd, ids, file="fmri_sceptic_signals_24basis_mmclock_Jun2017.RData")
#save(vmax, rtvmax, vauc, vsd, ventropy, vchosen, pemax, peauc, dauc, dsd, umax, rtumax, ventropy_decay_matlab, ventropy_fixedlrv_matlab, ids, file="fmri_sceptic_signals_24basis_specc.RData")
#save(vmax, rtvmax, vauc, vsd, ventropy, vchosen, vchosen_suspicious, pemax, peauc, dauc, dsd, ids, file="fmri_sceptic_signals_24basis_specc_correctedvchosen.RData")

#detach("package:plyr", unload=TRUE)
library(dplyr)

#try out model of entropy exploration
#need RT, RT vmax, trialwise entropy, and reward/omission
edf <- melt(ventropy, varnames=c("ID", "run", "trial"), value.name="entropy")
rtvmaxdf <- melt(rtvmax, varnames=c("ID", "run", "trial"), value.name="rtvmax")
vchosendf <- melt(vchosen, varnames=c("ID", "run", "trial"), value.name="vchosen")
vmaxdf <- melt(vmax, varnames=c("ID", "run", "trial"), value.name="vmax")
#rtumaxdf <- melt(rtumax, varnames=c("ID", "run", "trial"), value.name="rtumax")
#umaxdf <- melt(umax, varnames=c("ID", "run", "trial"), value.name="umax")
pemaxdf <- melt(pemax, varnames=c("ID", "run", "trial"), value.name="pemax")

bdf <- merge(edf, rtvmaxdf, by=c("ID", "run", "trial"))
bdf <- merge(bdf, vchosendf, by=c("ID", "run", "trial"))
bdf <- merge(bdf, vmaxdf, by=c("ID", "run", "trial"))
#bdf <- merge(bdf, rtumaxdf, by=c("ID", "run", "trial"))
#bdf <- merge(bdf, umaxdf, by=c("ID", "run", "trial"))
bdf <- merge(bdf, pemaxdf, by=c("ID", "run", "trial"))
bdf <- merge(bdf, allData, by=c("ID", "run", "trial"))

str(bdf)

bdf <- bdf %>% arrange(ID, run, trial)

bdf$rowID <- rep(1:76, each=400) #Jon says decayH sorted ascending by subject
bdf <- merge(bdf, ventropy_decay_matlab, by=c("rowID", "run", "trial"))

bdf <- merge(bdf, ventropy_fixedlrv_matlab, by=c("rowID", "run", "trial"))


#bdf <- bdf %>% group_by(ID, run) %>% arrange(ID, run, trial) %>% mutate(timesteplag = Hmisc::Lag(timestep)) #do({.$tlag = Hmisc::Lag(.$timestep); return(.) }) %>% ungroup() # %>% ungroup()
#bdf <- bdf %>% group_by(ID, run) %>% arrange(ID, run, trial) %>% mutate(timesteplag = c(NA, diff(timestep))) #do({.$tlag = Hmisc::Lag(.$timestep); return(.) }) %>% ungroup() # %>% ungroup()

#need to be explicit about mutate because it was picking up mutate from plyr, which doesn't work properly here
bdf <- bdf %>% group_by(ID, run) %>%
    dplyr::mutate(
        rtlag=lag(rt, order_by=trial),
        rtlag2=lag(rt, order_by=trial, n=2),
        rtlag3=lag(rt, order_by=trial, n=3),
        pemaxlag=lag(pemax, order_by=trial),
        abspe=abs(pemax),
        abspelag=lag(abspe, order_by=trial),
        ppe=abs((pemax > 0)*pemax),
        ppelag=lag(ppe, order_by=trial),
        npe=abs((pemax <= 0)*pemax), #absolute value so that bigger numbers represent larger PEs
        npelag=lag(npe, order_by=trial),
        timesteplag = lag(timestep, order_by=trial), 
        timestepchange = c(NA, diff(timestep)),
        timestepchangelag = lag(timestepchange, n=1, order_by=trial),
        abstschange = abs(timestepchange),
        abstschangelag = lag(abstschange, n=1, order_by=trial),
        abstschangelag2 = lag(abstschange, n=2, order_by=trial),
        abstschangelag3 = lag(abstschange, n=3, order_by=trial),
        abstschangelag4 = lag(abstschange, n=4, order_by=trial),
        vdev = timestep - rtvmax,
        vdevlag = lag(vdev, order_by=trial),
        rtvmaxlag = lag(rtvmax, order_by=trial),
        #udev = timestep - rtumax,
        #udevlag = lag(udev, order_by=trial),
        #rtumaxlag = lag(rtumax, order_by=trial),
        absvdevlag = abs(vdevlag), #just absolute distance from max value (unsigned)
        evdev = vmax - vchosen,
        evdevlag = lag(evdev, order_by=trial),
        omission = factor(as.numeric(score > 0), levels=c(0,1), labels=c("Omission", "Reward")),
        omissionlag = lag(omission, order_by=trial),
        entropylag = lag(entropy, order_by=trial),
        #entropyHlag = lag(entropyH, order_by=trial),
        #entropyFlag = lag(entropyFixed, order_by=trial),
        wizentropy = as.vector(scale(entropy)),
        wizabstschange = as.vector(scale(abstschange)),
        distfromedge = sapply(timestep, function(x) { min(40-x,x) }),
        distfromedgelag = lag(distfromedge, order_by=trial)
    ) %>% ungroup() %>% arrange(ID, run, trial)

# counting process dataframe for survival analysis
library(tidyr)
library(tibble)
mmdf <- as.tibble(mmdf)
bdf$trial <- as.numeric(bdf$trial)
sdf <- as.tibble(merge(mmdf,bdf, by = c("ID", "run", "trial")))
sdf <- sdf %>% arrange(ID, run, trial)
sdf <- sdf %>% group_by(ID,run,trial) %>% 
  dplyr::mutate(
  t2 = seq(n())/10,
  t1 = t2-.1
  ) %>% ungroup() %>% arrange(ID, run, trial)

# sanity check: plot value by condition.  Checks out.
p <- ggplot(sdf,aes(t2,value, color = rewFunc)) + geom_jitter() + facet_wrap(~ID)
ggsave("binned_value_for_coxme_by_subject.pdf", p, width = 20, height = 20)
p <- ggplot(sdf,aes(t2,value, color = rewFunc)) + geom_jitter(size = 0.25) + facet_wrap(rewFunc~trial>10)
ggsave("binned_value_for_coxme_early_late.pdf", p, width = 10, height = 10)

# save data for coxme analyses
setwd("~/code/clock_analysis/coxme")
save(file="clock_for_coxme_value_only_070518.RData", list = c('sdf', 'bdf'))


bdf_toplot <- bdf %>% select(ID, run, trial, rtvmax, rtumax, timestep) %>% gather(key=signal, value=value, rtvmax, rtumax, timestep)

subids <- 1:length(unique(bdf_toplot$ID))

#pdf("signals.pdf", width=10, height=15)
#for (s in subids) {
#  g <- ggplot(filter(bdf_toplot, ID==s), aes(x=trial, y=value, color=signal, shape=signal)) + facet_wrap(~run, ncol=1) + geom_line() + geom_point(position=position_jitter(width=0.5)) + ggtitle(s)
#  plot(g)
#}
#dev.off()

#pdf("ufunc.pdf", width=10, height=15)
#for (s in subids) {
#  g <- ggplot(filter(um, ID==s), aes(x=trial, y=timestep, color=U, fill=U)) + facet_wrap(~run, ncol=1) +  geom_raster(interpolate = TRUE, alpha=0.7) + ggtitle(s) #geom_tile() 
#  plot(g)
#}

#pdf("signals_withudist.pdf", width=10, height=15)
#for (s in subids) {
#  #color=signal, 
##  g <- ggplot(filter(bdf_toplot, ID==s), aes(x=trial, y=value, shape=signal)) + facet_wrap(~run, ncol=1) + geom_line(position=position_dodge(width=0.5)) + geom_point(position=position_dodge(width=0.5)) + ggtitle(s) +
##      theme_bw(base_size=15) + geom_raster(data=filter(um, ID==s), aes(x=trial, y=timestep, fill = U, color=U, shape=NULL), interpolate = TRUE) +
##      scale_x_continuous("Trial") + scale_y_continuous("Timestep") 
##  plot(g)
#  g <- ggplot(data=filter(um, ID==s), aes(x=trial, y=timestep, fill = U, color=NULL, shape=NULL)) + theme_bw(base_size=15) + 
#      geom_raster(interpolate = TRUE) + facet_wrap(~run, ncol=1) +
#      geom_point(data=filter(bdf_toplot, ID==s), aes(x=trial, y=value, shape=signal, color=signal, fill=NULL), position=position_dodge(width=0.5), size=1.5) + #, fill=NULL
#      geom_line(data=filter(bdf_toplot, ID==s), aes(x=trial, y=value, color=signal, fill=NULL), position=position_dodge(width=0.5), size=1.5) + #, fill=NULL #shape=signal,  
#      ggtitle(s) + scale_x_continuous("Trial") + scale_y_continuous("Timestep") 
#  
#  plot(g)
#  
#}
#dev.off()


#g <- ggplot() +  geom_raster(data=filter(um, ID==s), aes(x=trial, y=timestep, fill = U, color=U, shape=NULL), interpolate = TRUE) + facet_wrap(~run, ncol=1) + #+ geom_line(position=position_dodge(width=0.5)) + geom_point(position=position_dodge(width=0.5)) + ggtitle(s) +
#    theme_bw(base_size=15) + geom_line(data)
#plot(g)


#looks correct
bdf %>% select(ID, run, trial, rt, timestep, timesteplag, timestepchange, rtvmax, vdev, vdevlag, score, omission) %>% arrange(ID, run, trial) %>% print(n=101)

bdf %>% select(ID, run, trial, rt, timestep, score, omission, entropy, wizentropy) %>% arrange(ID, run, trial) %>% print(n=101)


#save(file="dataframe_for_entropy_analysis_Mar2017.RData", bdf)
save(file="dataframe_for_entropy_analysis_specc_Mar2017.RData", bdf)

#save(file="dataframe_for_entropy_analysis_SPECC_Dec2016.RData", bdf)
#write.csv(file="sceptic_fmri_behavior_Nov2016.csv", bdf, row.names=FALSE)

#now fit basic model: timestep change is predicted by prior deviation from value max and entropy 
library(lme4)




#summary(test_cur)
#anova(test_past, test_cur)

#example of entropy
mixing = seq(0, 1, .01)
mixing <- seq(1, 10)
u = runif(3000, min=3, max=6)
y = rnorm(3000, mean=4)
#y = c(rep(100, 1000), rep(0,2000))
m = rep(NA, length(mixing))
s = rep(NA, length(mixing))
e = rep(NA, length(mixing))
library(entropy)
library(np)
dd <- matrix(NA_real_, nrow=length(mixing), ncol=length(u))


enorm = rep(NA, length(mixing))
esoft = rep(NA, length(mixing))
options(error=recover, warn=0)
for (i in 1:length(mixing)) {
  #d = mixing[i]*u + (1-mixing[i])*y
  d <- y/mixing[i] 
  dd[i,] <- d
  #hist(d)
  #dd = discretize(d, numBins=100)
  #maximum entropy for a 100-bit vector is log(100) = 4.605
  #m[i] = entropy(dd)
  #s[i] = sd(d)
  #e[i] = npsymtest(d,boot.num=200)
  dnorm <- d/sum(d)
  dsoft <- exp(d)/sum(exp(d))
  enorm[i] <- -sum(dnorm*log(dnorm))
  esoft[i] <- -sum(dsoft*log(dsoft))
}

matplot(dd)
dm <- melt(dd)
ggplot(dm, aes(x=value)) + geom_histogram() + facet_wrap(~Var1)

cor(enorm, esoft)

scatter3d(1:10, enorm, esoft, point.col = "blue", surface=FALSE)
library(plotly)

plot(1:101, enorm, main="normalize")
plot(1:101, esoft, main="softmax")

par(mfrow=c(3,1))
plot(mixing, m/log(100), type="l")
plot(mixing, s, type="l")
plot(mixing, e, type="l")

cor(cbind(m, s, e))
