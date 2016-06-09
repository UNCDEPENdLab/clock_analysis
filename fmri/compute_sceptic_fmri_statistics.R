#plot the fitted data from VBA
library(R.matlab)
library(ggplot2)
library(plyr)
library(reshape2)
library(entropy)
setwd("/Users/michael/Data_Analysis/temporal_instrumental_agent/clock_task/vba_fmri")
#fit <- readMat("/Users/michael/Data_Analysis/temporal_instrumental_agent/clock_task/vba_fmri/posterior_states_decay_nomultisession.mat")
fit <- readMat("/Users/michael/Data_Analysis/temporal_instrumental_agent/clock_task/vba_fmri/posterior_states_decay_nomultisession_constrain0p025.mat")
basis <- readMat("/Users/michael/Data_Analysis/temporal_instrumental_agent/clock_task/vba_fmri/sceptic_fmri_basis_setup.mat")
source("clock_functions.R")

ids <- gsub(".*CORRECT_(\\d+)_fixed_decay.*", "\\1", unlist(fit$fitfiles), perl=TRUE)

#plot RTs for clock data
allData <- getClockGroupData(path="/Users/michael/Data_Analysis/temporal_instrumental_agent/clock_task/subjects")

with(allData, table(LunaID, rewFunc, emotion))
allData$rewarded <- factor(allData$score > 0, levels=c(TRUE, FALSE), labels=c("Reward", "Omission"))
pdf("ClockfMRI_AllSubjRTs_withMax.pdf", width=14, height=8)
for (s in split(allData, allData$LunaID)) {
  sm <- reshape2::melt(s[,c("trial", "rt", "emotion", "rewFunc", "bestRewRT", "bestEVRT", "run", "rewarded")], id.vars=c("trial", "emotion", "rewFunc", "run", "rewarded"))
  g <- ggplot(sm, aes(x=trial, y=value, color=variable)) + geom_line() + geom_point(data=subset(sm, variable=="rt"), aes(shape=rewarded), size=2) + facet_grid(emotion ~ rewFunc) + ggtitle(s$LunaID[1L]) + ylab("RT") + xlab("Trial") +
      geom_text(aes(x=45, y=3700, label=run, color=NULL), color="black", family="Helvetica", size=6) + theme_bw(base_size=21) + scale_shape("Outcome") + scale_color_discrete("") + coord_cartesian(ylim=c(0,4000))
  print(g)
}
dev.off()


#plot of V

trialsToPlot <- c(1, 5, 10, 20, 30, 40, 49)
Vsub <- fit$V[,,,trialsToPlot]

#PE and D are right-shifted by 1 such that PE(1) falls in the second column since it is used to compute V(2) = V(1) + alpha*(r(1) - V(1)) 
PEsub <- fit$PE[,,,trialsToPlot + 1]
Dsub <- fit$D[,,,trialsToPlot + 1]

toplot <- list(V=Vsub, PE=PEsub, D=Dsub)

for (v in 1:length(toplot)) {

  #multiply against basis to get distribution/representation
  Vexpand <- aaply(toplot[[v]], c(1,2,4), function(b) {
        v = outer(b, rep(1, dim(basis$tvec)[2])) * basis$gaussmat #use vector outer product to replicate weight vector
        vfunc <- colSums(v)
        return(vfunc)
      })
  
  pdf(paste(names(toplot)[v], "unfolding.pdf"), width=12, height=10)
  for (s in 1:dim(Vexpand)[1]) {
    subj <- Vexpand[s,,,]
    m <- melt(subj, varnames=c("run", "trial", "timestep"))
    m$LunaID <- factor(ids[s])
    m$trialFac <- factor(m$trial, levels=1:length(trialsToPlot), labels=paste0("i = ", trialsToPlot))
    m$trial <- plyr::mapvalues(m$trial, 1:length(trialsToPlot), trialsToPlot)
    rtdf <- subset(allData, LunaID == ids[s] & trial %in% trialsToPlot, select=c(run, trial, rt, score)) #[,c("LunaID", "run", "trial", "rt", "score")], by=c("LunaID", "run", "trial"))
    rtdf$timestep <- plyr::round_any(rtdf$rt, 100)/100
    rtdf$trialFac <- factor(rtdf$trial, levels=trialsToPlot, labels=paste0("i = ", trialsToPlot))
    g <- ggplot(m, aes(x=timestep, y=value)) + geom_hline(yintercept=0, color="blue") + geom_line() + geom_text(data=rtdf, aes(label=score, y=0)) + facet_grid(run ~ trialFac, scales="free_y") + ggtitle(ids[s])
    plot(g)
  }
  dev.off()
  
}

#note: the problems with "positive decay" and other weird effects in decay and PE were due to an error in the use of an incorrect refspread
#and the max sigma of the eligibility of 1.0 (instead of this being truly free). I have fixed these 8Jun2016 and plots now look correct.
#hence, code below only applies to previous problem.

#diagnose a few weird things
#positive decay:
# SUB17: id 11178, run 7, trial 20
V7 <-  fit$V[17,7,,]

V7ex <- apply(V7, c(2), function(b) {
      v = outer(b, rep(1, dim(basis$tvec)[2])) * basis$gaussmat #use vector outer product to replicate weight vector
      vfunc <- colSums(v)
      return(vfunc)
    })

D7 <-  fit$D[17,7,,]

D7ex <- apply(D7, c(2), function(b) {
      v = outer(b, rep(1, dim(basis$tvec)[2])) * basis$gaussmat #use vector outer product to replicate weight vector
      vfunc <- colSums(v)
      return(vfunc)
    })

PE7 <-  fit$PE[17,7,,]

PE7ex <- apply(PE7, c(2), function(b) {
      v = outer(b, rep(1, dim(basis$tvec)[2])) * basis$gaussmat #use vector outer product to replicate weight vector
      vfunc <- colSums(v)
      return(vfunc)
    })

par(mfrow=c(3,1))
plot(1:40, V7ex[,20], type="l")
plot(1:40, D7ex[,21], type="l")
plot(1:40, PE7ex[,21], type="l")


#generate regressors for PE, value, and decay at each trial
tofmri <- list(V=fit$V, PE=fit$PE, D=fit$D)
fmriexpanded <- list()
for (v in 1:length(tofmri)) {
  
  #multiply against basis to get distribution/representation
  fmriexpanded[[ names(tofmri)[v] ]] <- aaply(tofmri[[v]], c(1,2,4), function(b) {
        v = outer(b, rep(1, dim(basis$tvec)[2])) * basis$gaussmat #use vector outer product to replicate weight vector
        vfunc <- colSums(v)
        return(vfunc)
      })
  
  dimnames(fmriexpanded[[ names(tofmri)[v] ]]) <- list(LunaID=ids, run=1:8, trial=1:50, timestep=1:40)
}

#fmriexpanded now has arrays that are subjects x runs x trials x timesteps
#make regressors of these

#1) Maximum value
vmax <- apply(fmriexpanded[["V"]], c(1,2,3), function(trial) {
      if (sd(trial) < .01) { #essentially no variability, so hard to say there's a max 
        return(NA) 
      } else {
        return(max(trial))
      } 
    })
  
#2) AUC of value function
vauc <- apply(fmriexpanded[["V"]], c(1,2,3), function(trial) {
      if (sd(trial) < .01) { #essentially no variability, so hard to say there's a max 
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
      if (sd(trial) < .01) { #essentially no variability, so hard to say there's a max 
        return(NA) 
      } else {
        #normalize entropy by discretizing into 20 bins
        dd <- discretize(trial, numBins=20)
        entropy(dd)/log(20)
      } 
    })

#4) standard deviation of value function
vsd <- apply(fmriexpanded[["V"]], c(1,2,3), function(trial) {
      if (sd(trial) < .01) { #essentially no variability, so hard to say there's a max 
        return(NA) 
      } else {
        return(sd(trial))
      } 
    })

#5) value of chosen option
mv <- melt(fmriexpanded[["V"]], varnames=c("subj", "run", "trial", "timestep"))
mv$LunaID <- factor(mv$subj, levels=sort(unique(mv$subj)), labels=ids)
allData$timestep <- plyr::round_any(allData$rt, 100)/100

#there are two RTs less than 50ms that round down to zero, leading to their being missing in the merge
#this may be right -- these may not be valid rts... but subjects do get feedback
#there are also 119 RTs above 4s, representing non-response. Again, subjects get feedback (always zero) but is a bit different
#for now:
allData$timestep[which(allData$timestep==0)] <- 1

#30400 rows = 76 subjects x 8 runs x 50 trials
#mv now has one row per subject and trial representing the value at the chosen option.
mv <- merge(mv, allData, by=c("LunaID", "run", "trial", "timestep"))
with(mv, table(LunaID, run))

#recast into 3d array to match other regressors.
vchosen <- acast(mv, LunaID ~ run ~ trial, value.var="value")

#PE (need to right shift to get PE(i)
pemax <- apply(fmriexpanded[["PE"]], c(1,2,3), function(trial) {
      #note that some trials have essentially no PE
      if (min(trial) >= 0){
        return(max(trial))
      } else {
        return(min(trial))
      }
    })

#AUC of PE (seems like it will be highly correlated with max PE, so probably not informative)... yeah, looks like r = .99 for the first run!
peauc <- apply(fmriexpanded[["PE"]], c(1,2,3), sum)

#AUC of decay
dauc <- apply(fmriexpanded[["D"]], c(1,2,3), sum)

#SD of decay, which essentially captures how much the amount of decay varies by response option
#Should have a high inverse relationship with AUC, I think.
dsd <- apply(fmriexpanded[["D"]], c(1,2,3), sd)

save(vmax, vauc, vsd, ventropy, vchosen, pemax, peauc, dauc, dsd, ids, file="fmri_sceptic_signals.RData")


#example of entropy
mixing = seq(0, 1, .01)
u = runif(3000, min=-3, max=3)
y = rnorm(3000)
#y = c(rep(100, 1000), rep(0,2000))
m = rep(NA, length(mixing))
s = rep(NA, length(mixing))
e = rep(NA, length(mixing))
library(entropy)
library(np)

for (i in 1:length(mixing)) {
  d = mixing[i]*u + (1-mixing[i])*y
  #hist(d)
  dd = discretize(d, numBins=100)
  #maximum entropy for a 100-bit vector is log(100) = 4.605
  m[i] = entropy(dd)
  s[i] = sd(d)
  #e[i] = npsymtest(d,boot.num=200)
}

par(mfrow=c(3,1))
plot(mixing, m/log(100), type="l")
plot(mixing, s, type="l")
plot(mixing, e, type="l")

cor(cbind(m, s, e))