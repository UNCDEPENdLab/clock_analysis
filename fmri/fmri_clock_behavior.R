setwd(file.path(getMainDir(), "clock_analysis", "fmri"))
subinfo <- read.table("subinfo_db", header=TRUE)
subinfo <- subset(subinfo, lunaid != 10637) #high movement exclude subject

allfits <- merge(allfits, subinfo, by="lunaid")

table(allfits$epsilonBeta > 0)

prop.table(table(allfits$epsilonBeta > 0))

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
xyplot(alphaG - alphaN ~ age, allfits)
xyplot(epsilonBeta ~ age, allfits)
xyplot(totreward ~ age, allfits)
cor.test(~ totreward + age, allfits)
xyplot(totreward ~ epsilonBeta, allfits)
cor.test(~ totreward + epsilonBeta, allfits)
cor.test(~ totreward + as.integer(epsilonBeta > 0), allfits)
xyplot(avg_ev ~ age, allfits)
cor.test(~ avg_ev + age, allfits)

summary(lm(totreward ~age*as.integer(epsilonBeta>0), allfits))


allfits$explorer <- as.integer(allfits$epsilonBeta > 0)
allfits$female.c <- allfits$female - mean(allfits$female)
allfits$age.c <- allfits$age - mean(allfits$age)

allfits$exp_age.c <- with(allfits, age.c*explorer)

hist(allfits$alphaG)
hist(allfits$alphaN)

allfits$alpha_diff <- allfits$alphaG - allfits$alphaN #difference in learning rate is approximately normal
allfits$alpha_diff.c <- allfits$alpha_diff - mean(allfits$alpha_diff) #difference in learning rate is approximately normal

allfits$avg_ev.c <- allfits$avg_ev - mean(allfits$avg_ev)

#use mean-centered female to treat as nuisance confound (roughly +.5/-.5)
write.table(allfits[,c("lunaid", "adult", "female.c", "age.c", "explorer", "exp_age.c", "alpha_diff.c", "avg_ev.c")], file="memacov_ageexplore_center.txt", row.names=FALSE, quote=FALSE, sep="\t")
write.table(allfits[,c("lunaid", "adult", "female", "age", "explorer", "alpha_diff")], file="memacov_ageexplore.txt", row.names=FALSE, quote=FALSE, sep="\t")


cor.test(~ rw_alphaV + age, allfits)
cor.test(~ rw_betaV + age, allfits)
cor.test(~ avg_ev + age, allfits)
cor.test(~ totreward + age, allfits)
cor.test(~ rewhappy + age, allfits)
cor.test(~ rewfear + age, allfits)
cor.test(~ rewscram + age, allfits)


###APR2015: Basic behavioral analyses of RTs and rewards








##NOV 2014: comparisons of behavioral data between controls and BPD participants
##this is the n=15 BPD sample with 15 age and sex-matched controls.
##see match_bpd_controls.R for approach
setwd(file.path(getMainDir(), "clock_analysis", "fmri"))
subjinfo <- read.table("subjinfo_bpd15_control15.txt", header=TRUE)
subjinfo$numid <- as.integer(sub("(\\d{3})[A-z]{2}", "\\1", subjinfo$lunaid, perl=TRUE))

fitobjs <- file.path(getMainDir(), "clock_analysis", "fmri", "fmri_fits", paste0(subjinfo$numid, "_fitinfo.RData"))

n30fits <- get_fit_array(fitobjs)
n30fits <- merge(n30fits, subjinfo[,c("numid", "female", "age")], by.x="lunaid", by.y="numid") #put in female and age

t.test(K ~ bpd, n30fits)
t.test(lambda ~ bpd, n30fits)
t.test(scale ~ bpd, n30fits)
t.test(alphaG ~ bpd, n30fits)
t.test(alphaN ~ bpd, n30fits)
t.test(alphaG - alphaN ~ bpd, n30fits)
t.test(rho ~ bpd, n30fits)
t.test(epsilonBeta ~ bpd, n30fits)
t.test(log(epsilonBeta + 0.5) ~ bpd, n30fits)
t.test(epsilonBeta > 0 ~ bpd, n30fits)
wilcox.test(epsilonBeta ~ bpd, n30fits)

n30fits$exploreGt0 <- as.integer(n30fits$epsilonBeta > 0)
library(exact2x2)

exact2x2(n30fits$exploreGt0, n30fits$bpd)



library(lattice)
histogram(~epsilonBeta | bpd, n30fits)
histogram(~alphaG-alphaN| bpd, n30fits)

xyplot(epsilonBeta ~ age | bpd, n30fits)

t.test(rw_alphaV ~ bpd, n30fits)
t.test(rw_betaV ~ bpd, n30fits)
t.test(avg_ev ~ bpd, n30fits)
t.test(totreward ~ bpd, n30fits)
t.test(rewhappy ~ bpd, n30fits)
t.test(rewfear ~ bpd, n30fits)
t.test(rewscram ~ bpd, n30fits)


t.test(evhappy ~ bpd, n30fits)
t.test(evfear ~ bpd, n30fits)
t.test(evscram ~ bpd, n30fits)


fitobjs <- file.path(getMainDir(), "clock_analysis", "fmri", "fmri_fits", paste0(subjinfo$numid, "_fitinfo.RData"))



alldf$bpd <- factor(as.numeric(sub("(\\d+)_fitinfo.RData", "\\1", alldf$id, perl=TRUE)) < 10000, levels=c(FALSE,TRUE), labels=c("Control", "BPD"))

options(width=140)
head(alldf, n=50)

library(lme4)
alldf$rpes_lag_pos <- factor(as.numeric(alldf$rpes_lag > 0), levels=c(0,1), labels=c("Prior RPE-", "Prior RPE+"))
alldf$rpes_lag_abs <- abs(alldf$rpes_lag)

test <- lmer(rtchange ~ rpes_lag*rpes_lag_pos*bpd*cond + (1|id), alldf)

cm <- lmerCellMeans(test, cont.pts=list(rpes_lag=c(-100, -75, -50, -25, 0, 25, 50, 75, 100)))
cm_mod <- rbind(subset(cm, rpes_lag_pos=="Prior RPE-" & rpes_lag <= 0), subset(cm, rpes_lag_pos=="Prior RPE+" & rpes_lag >= 0))

#just look at emotion x contingency a prior RPE direction effect on delta RT
pdf("RT change by RPE.pdf", width=12, height=6)
ggplot(cm_mod, aes(x=rpes_lag, y=rtchange, color=bpd)) + geom_line(size=2) + 
    geom_point(size=5) + facet_wrap(~contingency) + 
    geom_errorbar(aes(ymin=rtchange-se, ymax=rtchange+se), size=2, width=5) + geom_vline(xintercept=0) + geom_hline(yintercept=0) + 
    ylab("Change in RT (ms)") + xlab("Magnitude of Reward Prediction Error (points)") + theme_bw(base_size=24) + scale_color_brewer("Group", palette="Set1")

dev.off()

#test <- lmer(rtchange ~ rpes_lag*rpes_lag_pos*bpd*contingency + (1|id), alldf)
#test <- lmer(rtchange ~ rpes_lag_pos*cond*contingency*bpd + (1|id), alldf)

levels(alldf$cond) <- c("Fear", "Happy", "Scrambled")

#test <- lmer(rtchange ~ rpes_lag_pos*bpd*cond + (1|id), alldf)

car::Anova(test)
library(ez)
#ezANOVA(data=alldf, dv=rtchange, within=.(cond, rpes_lag_pos, rpes_lag), wid=id, between=bpd)
#(rtchange ~ rpes_lag*rpes_lag_pos*bpd*cond + (1|id), alldf)

source(file.path(getMainDir(), "Miscellaneous", "Global_Functions.R"))
cm <- lmerCellMeans(test, cont.pts=list(rpes_lag=c(-100, -75, -50, -25, 0, 25, 50, 75, 100)))

getDescriptivesByGroup(alldf, "rtchange", "rpes_lag_pos", zscale=FALSE, tscale=FALSE)

#just look at emotion x contingency a prior RPE direction effect on delta RT
pdf("RT change by PE emo and contingency.pdf", width=12, height=6)
ggplot(cm, aes(x=contingency, y=rtchange, color=cond)) +  
    geom_point(size=5, position=position_dodge(width=0.7)) + facet_grid(bpd~rpes_lag_pos) + 
    geom_errorbar(aes(ymin=rtchange-se, ymax=rtchange+se), size=2, width=0.7, position="dodge") + geom_vline(xintercept=0) + geom_hline(yintercept=0) + 
    ylab("Change in RT (ms)") + xlab("Contingency") + theme_bw(base_size=24) + scale_color_brewer("Emotion", palette="Set1")

dev.off()


cm_mod <- rbind(subset(cm, rpes_lag_pos=="Prior RPE-" & rpes_lag <= 0), subset(cm, rpes_lag_pos=="Prior RPE+" & rpes_lag >= 0))

cm_mod$rpes_lag_jitter <- NA_real_
for (i in 1:nrow(cm_mod)) {
  if (cm_mod[i,"rpes_lag_pos"]=="Prior RPE-") {
    if (cm_mod[i,"bpd"] == "Control") {
      v <- cm_mod[i,"rpes_lag"] - 5
    } else {
      v<- cm_mod[i,"rpes_lag"]- 10
    }
  } else if (cm_mod[i,"rpes_lag_pos"]=="Prior RPE+") {
    if (cm_mod[i,"bpd"] == "Control") {
      v<- cm_mod[i,"rpes_lag"] +  10
    } else {
      v<- cm_mod[i,"rpes_lag"] + 5
    }
  }
  cm_mod$rpes_lag_jitter[i] <- v
}

cm_mod$grouping <- paste0(cm_mod$bpd, cm_mod$rpes_lag_pos)
pdf("RT change by RPE.pdf", width=12, height=6)
ggplot(cm_mod, aes(x=rpes_lag_jitter, y=rtchange, color=bpd, group=grouping)) + geom_line(size=2) + 
    geom_point(size=5) + facet_wrap(~cond) + 
    geom_errorbar(aes(ymin=rtchange-se, ymax=rtchange+se), size=2, width=5) + geom_vline(xintercept=0) + geom_hline(yintercept=0) + 
    ylab("Change in RT (ms)") + xlab("Magnitude of Reward Prediction Error (points)") + theme_bw(base_size=24) + scale_color_brewer("Group", palette="Set1")

dev.off()
#ggplot(cm, aes(x=rpes_lag_pos, y=rtchange, color=bpd)) + 
#    geom_point(size=5) + facet_wrap(~cond) + 
#    geom_errorbar(aes(ymin=rtchange-se, ymax=rtchange+se))


test <- lmer(rtchange ~ bpd*cond*rpes_lag_abs + (1|id), subset(alldf, rpes_lag > 0)) #change after positive rpes only
test <- lmer(rtchange ~ bpd*cond*rpes_lag_abs + (1|id), subset(alldf, rpes_lag < 0)) #change after positive rpes only

str(n30fits)


#look at slow down in RT after RPE-
afterNeg <- subset(alldf, rpes_lag < 0)
afterPos <- subset(alldf, rpes_lag > 0)

cor(afterNeg$rtchange, afterNeg$rpes_lag)
cor(afterPos$rtchange, afterPos$rpes_lag)

#simplest descriptives...
mean(afterNeg$rtchange)
mean(afterPos$rtchange)






m1 <- lmer(epsilonBeta ~ emotion + (1|subid), allf)

library(lme4)
car::Anova(m1 <- lmer(epsilonBeta ~ reward*emotion + (1|subid), allf))
car::Anova(m2 <- lmer(scale ~ reward*emotion + (1|subid), allf))
car::Anova(m3 <- lmer(alphaG ~ reward*emotion + (1|subid), allf))
car::Anova(m4 <- lmer(alphaN ~ reward*emotion + (1|subid), allf))
allf$alpha_diff <- allf$alphaG -allf$alphaN 
car::Anova(m5 <- lmer(alpha_diff ~ reward*emotion + (1|subid), allf))
car::Anova(m6 <- lmer(K ~ reward*emotion + (1|subid), allf))
car::Anova(m7 <- lmer(lambda ~ reward*emotion + (1|subid), allf))
car::Anova(m8 <- lmer(rho ~ reward*emotion + (1|subid), allf))

subinfo <- read.table("../../subinfo_db", header=TRUE)



###ANALYSIS OF FRANK LEARNING PARAMETERS
#does exploration vary by emotion?
load("fits_epsilon_emo_controls_5May2015.RData")

theta_params <- do.call(rbind, lapply(epsemo_controls, function(x) {
          x$theta[,"cur_value"]      
        }))

str(theta_params)

eps_emo <- melt(theta_params[,c("epsilonBeta/run_condition:fear", "epsilonBeta/run_condition:scram", "epsilonBeta/run_condition:happy")], varnames=c("ID", "Condition"), value.name="Epsilon")
eps_emo$ID <- factor(eps_emo$ID)
#significant by ANOVA: p = .006
ezANOVA(data=eps_emo, dv=Epsilon, wid=ID, within=Condition)

#descriptives
ezStats(data=eps_emo, dv=Epsilon, wid=ID, within=Condition)

histogram(~Epsilon | Condition, eps_emo)

#highly non-normal distributions. look at nonparametric test
#yes, significant by Friedman: p = .017
friedman_test(Epsilon ~ Condition | ID, eps_emo)

library(doMC);options(cores=4);registerDoMC()
#significant by permutation test: p = .005
perm_results <- ezPerm(data=eps_emo, dv=Epsilon, wid=ID, within=Condition, perms=10000, parallel=TRUE)
print(perm_results)

#basic pairwise tests
friedman_test(Epsilon ~ Condition | ID, gdata::drop.levels(subset(eps_emo, Condition != "epsilonBeta/run_condition:fear"))) #p = .02 happy > scram
friedman_test(Epsilon ~ Condition | ID, gdata::drop.levels(subset(eps_emo, Condition != "epsilonBeta/run_condition:scram"))) #p = .80, fear vs. happy
friedman_test(Epsilon ~ Condition | ID, gdata::drop.levels(subset(eps_emo, Condition != "epsilonBeta/run_condition:happy"))) #p = .009, fear > scram

t.test(Epsilon ~ Condition, gdata::drop.levels(subset(eps_emo, Condition != "epsilonBeta/run_condition:scram")), paired=TRUE)
