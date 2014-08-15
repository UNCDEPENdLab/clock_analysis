#simple comparisons of parameters between BPD and control groups
setwd(file.path(getMainDir(), "clock_analysis", "fmri", "fmri_fits"))
fitobjs <- list.files(pattern=".*fitinfo.RData")
#cw behavior looks pretty ugly... increasing AIC for more parameters, no go for gold or beyond
fitobjs <- fitobjs[-1*which(fitobjs=="15_fitinfo.RData")]
subids <- as.integer(sub("^(\\d+)_fitinfo.RData", "\\1", fitobjs, perl=TRUE))
bpdsub <- as.integer(subids < 10000)


allfits <- c()
for (i in 1:length(fitobjs)) {
  loc <- local({load(fitobjs[i]); environment()})$f #time-clock fit object
  emo <- loc$run_condition
  rew <- loc$rew_function
  pars <- loc$theta[,"cur_value"]
  
  rw <- local({load(fitobjs[i]); environment()})$f_value #rescorla-wagner fit object
  rewhappy <- sum(rw$Reward[which(rw$run_condition == "happy"),])
  rewfear <- sum(rw$Reward[which(rw$run_condition == "fear"),])
  rewscram <- sum(rw$Reward[which(rw$run_condition == "scram"),])
  
  evhappy <- sum(rw$ev[which(rw$run_condition == "happy"),])
  evfear <- sum(rw$ev[which(rw$run_condition == "fear"),])
  evscram <- sum(rw$ev[which(rw$run_condition == "scram"),])
  
  pars <- c(pars, rw_alphaV=rw$theta["alphaV", "cur"], rw_betaV=rw$theta["betaV", "cur"], avg_ev=mean(rw$ev), totreward=sum(rw$Reward),
      rewhappy=rewhappy, rewfear=rewfear, rewscram=rewscram,
      evhappy=evhappy, evfear=evfear, evscram=evscram)
  allfits <- rbind(allfits, c(id=subids[i], bpd=bpdsub[i], pars))
}

allfits <- data.frame(allfits)

t.test(K ~ bpd, allfits)
t.test(lambda ~ bpd, allfits)
t.test(scale ~ bpd, allfits)
t.test(alphaG ~ bpd, allfits)
t.test(alphaN ~ bpd, allfits)
t.test(rho ~ bpd, allfits)
t.test(epsilonBeta ~ bpd, allfits)
t.test(epsilonBeta > 0 ~ bpd, allfits)
wilcox.test(epsilonBeta ~ bpd, allfits)

t.test(rw_alphaV ~ bpd, allfits)
t.test(rw_betaV ~ bpd, allfits)
t.test(avg_ev ~ bpd, allfits)
t.test(totreward ~ bpd, allfits)
t.test(rewhappy ~ bpd, allfits)
t.test(rewfear ~ bpd, allfits)
t.test(rewscram ~ bpd, allfits)


t.test(evhappy ~ bpd, allfits)
t.test(evfear ~ bpd, allfits)
t.test(evscram ~ bpd, allfits)
