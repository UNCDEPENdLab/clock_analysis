#simple comparisons of parameters between BPD and control groups
setwd(file.path(getMainDir(), "clock_analysis", "fmri", "fmri_fits"))
fitobjs <- list.files(pattern=".*fitinfo.RData")
#cw behavior looks pretty ugly... increasing AIC for more parameters, no go for gold or beyond
#fitobjs <- fitobjs[!which(fitobjs=="15_fitinfo.RData")]
subids <- as.integer(sub("^(\\d+)_fitinfo.RData", "\\1", fitobjs, perl=TRUE))
bpdsub <- as.integer(subids < 10000)
fitobjs <- fitobjs[subids > 10000] #only keep controls from MMY3


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

#add age, narrow to sample under consideration for fMRI (n = 37)
subjage <- local({load(file.path(getMainDir(), "clock_analysis", "fmri", "subjdata_8Sep2014.RData")); environment()})$x
subjage <- subset(subjage, subj != 11258) #something wrong with fmri stats, excluded
subjage <- plyr::rename(subjage, c(subj="id"))

allfits <- merge(allfits, subjage, all.y=TRUE, by="id")
table(allfits$epsilonBeta > 0)

cor.test(~ K + age, allfits)
cor.test(~ lambda + age, allfits)
cor.test(~ scale + age, allfits)
cor.test(~ alphaG + age, allfits)
cor.test(~ alphaG + age, subset(allfits, alphaG>0.01))
cor.test(as.integer(allfits$alphaG > .01), allfits$age)
cor.test(~ alphaN + age, allfits)
cor.test(~ rho + age, allfits)
cor.test(~ epsilonBeta + age, allfits)
cor.test(as.integer(allfits$epsilonBeta > 0), allfits$age)
cor.test(as.integer(allfits$epsilonBeta > 0), allfits$totreward)

cor.test(~ K + adult, allfits)
cor.test(~ lambda + adult, allfits)
cor.test(~ scale + adult, allfits)
cor.test(~ alphaG + adult, allfits)
cor.test(~ alphaG + adult, subset(allfits, alphaG>0.01))
cor.test(as.integer(allfits$alphaG > .01), allfits$adult)
cor.test(~ alphaN + adult, allfits)
cor.test(~ rho + adult, allfits)
cor.test(~ epsilonBeta + adult, allfits)
cor.test(as.integer(allfits$epsilonBeta > 0), allfits$adult)
cor.test(as.integer(allfits$epsilonBeta > 0), allfits$totreward)



library(lattice)
xyplot(alphaG ~ age, allfits)
xyplot(alphaN ~ age, allfits)
xyplot(epsilonBeta ~ age, allfits)
xyplot(totreward ~ age, allfits)
xyplot(totreward ~ epsilonBeta, allfits)
xyplot(avg_ev ~ age, allfits)

allfits$explorer <- as.integer(allfits$epsilonBeta > 0)
allfits$female.c <- allfits$female - mean(allfits$female)
allfits$age.c <- allfits$age - mean(allfits$age)

allfits$exp_age.c <- with(allfits, age.c*explorer)

#use mean-centered female to treat as nuisance confound (roughly +.5/-.5)
write.table(allfits[,c("id", "female.c", "age.c", "explorer", "exp_age.c")], file="memacov_ageexplore.txt", row.names=FALSE, quote=FALSE, sep="\t")

cor.test(~ rw_alphaV + age, allfits)
cor.test(~ rw_betaV + age, allfits)
cor.test(~ avg_ev + age, allfits)
cor.test(~ totreward + age, allfits)
cor.test(~ rewhappy + age, allfits)
cor.test(~ rewfear + age, allfits)
cor.test(~ rewscram + age, allfits)





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
