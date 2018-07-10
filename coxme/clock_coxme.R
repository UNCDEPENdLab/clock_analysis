setwd("~/code/clock_analysis/coxme")
library(readr)
library(ggplot2)
library(tidyr)
library(tibble)
library(reshape2)
library(emmeans)
# library(factoextra)
# library(ggfortify)
# library(RColorBrewer)
# library(MASS)
# library(readr)
# library(VIM)
# library(mice)
# library(multcompView)
# library(stargazer)
library(dplyr)
library(lme4)
library(survival)
library(coxme)
library(survminer)
# library(OIsurv)
library(ggpubr)

load(file="clock_for_coxme_value_only_070518.RData")

# mark events (responses)
sdf$response <- round(sdf$rt/1000, digits = 1)==sdf$t2


# inspect piecewize hazard functions
library(muhaz)
# piecewise <- pehaz(ddf$latency, delta=ddf$quit, width=NA, min.time=0, max.time=20.1)
# plot(piecewise, xlab="Time, seconds", ylab="Quit Rate")
pwI <- pehaz(na.omit(bdf$rt[bdf$rewFunc=='IEV']), width=200, min.time=0, max.time=4000)
pwD <- pehaz(na.omit(bdf$rt[bdf$rewFunc=='DEV']), width=200, min.time=0, max.time=4000)
pwR <- pehaz(na.omit(bdf$rt[bdf$rewFunc=='CEVR']), width=200, min.time=0, max.time=4000)
pwC <- pehaz(na.omit(bdf$rt[bdf$rewFunc=='CEV']), width=200, min.time=0, max.time=4000)

h <- c(pwI$Hazard,pwD$Hazard,pwC$Hazard,pwR$Hazard)
times <- c(pwI$Cuts[2:length(pwI$Cuts)],pwD$Cuts[2:length(pwD$Cuts)],pwC$Cuts[2:length(pwC$Cuts)],pwR$Cuts[2:length(pwR$Cuts)])
condition <- c(rep('IEV',length(pwI$Hazard)), rep('DEV',length(pwD$Hazard)),rep('CEV',length(pwC$Hazard)),rep('CEVR',length(pwR$Hazard)))
H <- as.tibble(cbind(h,times,condition))
H$h <- as.numeric(H$h)
H$logh <- log(H$h)
H$times <- as.numeric(H$times)
# pdf("piecewise_hazard_fx_all_groups.pdf", width = 12, height = 4)
ggplot(H,aes(times, h, color = condition)) + geom_line()
ggplot(H,aes(times, logh, color = condition)) + geom_line()
# not what I expected: shared underlying hazard across contingencies

badfit <- survfit(Surv(t2) ~ rewFunc, type = 'right', origin = .1, data=sdf)
badfit1 <- survfit(Surv(rt) ~ rewFunc,  data=bdf)


plot(badfit1, mark.time=FALSE, lty=1:4,
     xlab="ms", ylab="Proportion still waiting")
legend(3000, .85, c("CEV", "CEVR", "DEV", "IEV"),
       lty=1:4, bty='n')

# most basic model
c00 <- coxme(Surv(rt) ~ rtlag + rewFunc + (1|ID/run), bdf)
c1 <- coxme(Surv(rt) ~ rtvmaxlag + (1|ID/run), bdf)
c2 <- coxme(Surv(rt) ~ rtlag + rtvmaxlag + entropylag + distfromedgelag + (1|ID/run), bdf)

cf0 <- coxph(Surv(t2) ~ 1, sdf)
# plot(cf0)
c0 <- coxme(Surv(t1,t2,response) ~ (1|ID), sdf)
summary(c1 <- coxme(Surv(t1,t2,response) ~ rtlag + rewFunc + (1|ID), sdf))
summary(c1a <- coxme(Surv(t1,t2,response) ~ rtlag + rewFunc + (1|ID/run), sdf))

summary(c2 <- coxme(Surv(t1,t2,response) ~ rtlag + rtvmaxlag + entropylag + (1|ID), sdf))
summary(c3 <- coxme(Surv(t1,t2,response) ~ rtlag + value + entropylag +  (1|ID), sdf))
summary(c4 <- coxme(Surv(t1,t2,response) ~ rtlag + value + uncertainty + entropylag +  (1|ID), sdf))
summary(c5 <- coxme(Surv(t1,t2,response) ~ rtlag + value + uncertainty + entropylag + distfromedgelag +  (1|ID), sdf))

# limit analysis to middle 3s to eliminate speed constraints and avoidance of interval end: U-aversion holds
summary(c6 <- coxme(Surv(t1,t2,response) ~ rtlag + value + uncertainty + entropylag + distfromedgelag +  (1|ID/run), sdf[sdf$t1>.5 & sdf$t1<3.5,]))

# limit analysis to middle 2s to eliminate speed constraints and avoidance of interval end: U-aversion holds
summary(c6a <- coxme(Surv(t1,t2,response) ~ rtlag + value + uncertainty + entropylag + distfromedgelag +  (1|ID/run), sdf[sdf$t1>1 & sdf$t1<3,]))

# limit to IEV to r/o uncertainty/value tradeoff explanation: U-aversion holds
summary(c7 <- coxme(Surv(t1,t2,response) ~ rtlag + value + uncertainty + entropylag + distfromedgelag +  (1|ID/run), sdf[sdf$rewFunc=='IEV',]))

# limit to CEVR to r/o probability/magnitude tradeoff explanation: U-aversion holds
summary(c8 <- coxme(Surv(t1,t2,response) ~ rtlag + value + uncertainty + entropylag + distfromedgelag +  (1|ID/run), sdf[sdf$rewFunc=='CEVR',]))


anova(c2,c3,c4, c5)
