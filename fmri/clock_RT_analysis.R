#analyses of RTs in clock task during fMRI
#these are more basic analyses of RTs (relatively model free, not learning-based per se)
#fmri_clock_behavior.R contains analyses of model-based parameters etc.

library(lme4)
library(ggplot2)
library(scales)
library(grid)
library(abind)
library(reshape2)
library(lattice)
library(ez)
library(coin)
library(lsmeans)
source(file.path(getMainDir(), "Miscellaneous", "Global_Functions.R"))
setwd(file.path(getMainDir(), "clock_analysis", "fmri"))

load("fmri_fits/allRTs_withLaggedVars_1Jun2015.RData")

subinfo <- read.table("subinfo_db", header=TRUE)
subinfo <- subset(subinfo, lunaid != 10637) #high movement exclude subject

setdiff(alldf$lunaid, subinfo$lunaid)
 
#allfits <- merge(allfits, subinfo, all.y=TRUE, by="lunaid")
alldf <- merge(alldf, subinfo, all.y=TRUE, by="lunaid")

alldf <- alldf[order(alldf$lunaid, alldf$trial),]
alldf$age.c <- alldf$age - mean(alldf$age, na.rm=TRUE) 
alldf$runningrtdev_lag.c <- alldf$runningrtdev_lag - mean(alldf$runningrtdev_lag, na.rm=TRUE)
alldf$vdevlag.c <- alldf$vdevlag - mean(alldf$vdevlag, na.rm=TRUE)
alldf$adult <- factor(alldf$adult, levels=c(0,1), labels=c("adolescent", "adult"))
#alldf$lunaid <- factor(alldf$lunaid)

hist(alldf$udevlag)
library(lattice)
histogram(~udevlag | contingency, alldf)
hist(alldf$vdevlag)
hist(alldf$maxu)
histogram(~maxu | contingency, alldf)


test <- lmer(rts ~ rt_lag + vdevlag + (1|lunaid), alldf)#, na.action=na.exclude)
summary(test)

test2 <- lmer(rts ~ rt_lag * vdevlag + (1|lunaid), alldf)
summary(test2)
anova(test, test2)

#0) look at rt change as a function of prior reward outcome (omitted/received)
test <- lmer(rtchange ~ omission_lag + (1|lunaid), alldf)
summary(test)

cm <- lmerCellMeans(test, n.cont=10)
pdf("figures/Effect of prior rew outcome on RT change.pdf", width=3.5, height=3.5)
ggplot(cm, aes(x=omission_lag, y=rtchange, ymin=rtchange-se, ymax=rtchange+se)) + geom_line() + geom_bar(stat="identity", fill="darkblue") + geom_errorbar(width=0.4) + geom_hline(yintercept=0) +
    ylab("Delta RT (ms)") + xlab("Outcome on Prior Trial") + theme_bw(base_size=16)
dev.off()



#1) look at rt change as a function of deviation from value on prior trial
#yes, if prior RT was far from the value max, then people tend to shift back toward the value max
#interpretation of vdevlag, RT(i) - Vmax(i):
# - Large positive values indicate that prior RT was much *slower* than Vmax (e.g., RT was 3000, Vmax was 1000)
# - Large negative values indicate that prior RT was much *faster* than Vmax (e.g., RT was 1000, Vmax was 3000)
test <- lmer(rtchange ~ vdevlag + (1|lunaid), alldf)
summary(test)

#huge effect: if prior RT was ~2000ms faster than Vmax, RTchange will be about +1000ms (slow down)
#							if prior RT was ~2000ms slower than Vmax, RTchange will be about -1000ms (speed up)
cm <- lmerCellMeans(test, n.cont=10)
ggplot(cm, aes(x=vdevlag, y=rtchange, ymin=rtchange-se, ymax=rtchange+se)) + geom_line() + geom_pointrange()

#2) If the prior RT was far from value max AND a reward omission occurred, does this enhance shift toward value max?
levels(alldf$omission_lag) <- c("Reward", "Omission")
test <- lmer(rtchange ~ omission_lag*vdevlag + (1|lunaid), alldf)
summary(test)

#asymmetric cross-over interaction: if prior RT was much slower than Vmax (person waited), then omissions lead to much bigger speed ups than reward (left side of X)
# if prior RT was faster than Vmax (person responded more quickly), then omissions lead to even longer waits than reward (right side of X)
# This is intuitive: If you waited longer than Vmax and got an omission, you move even harder toward Vmax (fast) 
# When you respond more quickly than indicated by learning/Vmax and don't get a reward, you shift harder toward Vmax by waiting longer
# These are suggestive of trialwise learning/updating of value
cm <- lmerCellMeans(test, n.cont=10)
pdf("figures/Effect of prior vdev and omission on RT change.pdf", width=8, height=6)
ggplot(cm, aes(x=vdevlag, y=rtchange, color=omission_lag, ymin=rtchange-se, ymax=rtchange+se)) + geom_line(size=2.5) + theme_bw(base_size=24) + #geom_linerange(size=1.0, width=250) + 
    ylab("Mean trialwise RT change (ms)") + xlab("RT deviation from maximum value\non prior trial (ms)") + scale_color_brewer("Prior Outcome", palette="Set2") + geom_vline(xintercept=0) + geom_hline(yintercept=0)
dev.off()

#3) Ala Moustafa, do we see a change in RT depending on whether the prior RT was faster or slower than (running) average?
test <- lmer(rtchange ~ priorrtfast_running + (1|lunaid), alldf)
summary(test)

# Yes, strong effect: Speed up after slower than average, slow down after faster than average
cm <- lmerCellMeans(test, n.cont=10)
pdf("figures/Effect of prior RT compared to average on RT change.pdf", width=3.5, height=3.5)
ggplot(cm, aes(x=priorrtfast_running, y=rtchange, ymin=rtchange-se, ymax=rtchange+se)) + geom_bar(stat="identity", fill="darkblue") + geom_errorbar(width=0.4) + geom_hline(yintercept=0) +
    ylab("Delta RT (ms)") + xlab("Speed of Prior RT\ncompared to running average") + theme_bw(base_size=14)
dev.off()


#quantitative version of the same thing (how far was RT from average)
#yes, an even bigger effect here... values run from about -3000 -- 3000
#big negative values indicate that prior RT was 
test <- lmer(rtchange ~ runningrtdev_lag + (1|lunaid), alldf)
summary(test)

cm <- lmerCellMeans(test, n.cont=10)
ggplot(cm, aes(x=runningrtdev_lag, y=rtchange, ymin=rtchange-se, ymax=rtchange+se)) + geom_line() + geom_pointrange(size=2)

#look at joint effects of prior RT versus average and prior reward outcome
test <- lmer(rtchange ~ runningrtdev_lag*omission_lag + (1|lunaid), alldf)
car::Anova(test) #yes, all whopping effects...
summary(test)

cm <- lmerCellMeans(test, n.cont=10)
pdf("figures/Effects of prior outcome and prior RT speed on RT change.pdf", width=6.5, height=3.5)
ggplot(cm, aes(x=runningrtdev_lag, y=rtchange, ymin=rtchange-se, ymax=rtchange+se, color=omission_lag)) + geom_line() + geom_linerange(size=2) +
    ylab("Delta RT (ms)") + xlab("Speed of Prior RT\ncompared to running average") + theme_bw(base_size=14)  + geom_hline(yintercept=0)
dev.off()


#put together this with #1 (vdevlag) and #2 (prior omission)
test <- lmer(rtchange ~ runningrtdev_lag*omission_lag*vdevlag + (1|lunaid), alldf)
car::Anova(test) #yes, all whopping effects...
summary(test)

#runningrtdev_lag M=75, SD=821
#looks like there is NOT a vdevlag * omission_lag effect for rewarded trials (parallel lines in priorRew panel)
#for omission trials, looks like there is a big effect of prior RT deviation from average
#for reward trials, there is a small effect of prior RT deviation from average, but a much more substantial effect of deviation from Vmax

cm <- lmerCellMeans(test, n.cont=10, divide="runningrtdev_lag")
pdf("figures/Effect of prior outcome, prior RT rel avg, and prior deviance from Vmax on RT change.pdf", width=7, height=4)
ggplot(cm, aes(x=vdevlag, y=rtchange, ymin=rtchange-se, ymax=rtchange+se, color=runningrtdev_lag)) + geom_line() + geom_pointrange(size=0.5) + 
    facet_wrap(~omission_lag) + ylab("Delta RT (ms)") + xlab("Deviation of prior RT from Vmax (ms)") + theme_bw(base_size=14)  + geom_hline(yintercept=0)
dev.off()



#two-way vdevlag*runningrtdev_lag interaction for omissions versus rewards
summary(lmer(rtchange ~ runningrtdev_lag*vdevlag + (1|lunaid), subset(alldf, omission_lag=="PriorOmission"))) #yes, two-way interaction
summary(mrew <- lmer(rtchange ~ runningrtdev_lag.c*vdevlag.c + (1|lunaid), subset(alldf, omission_lag=="PriorRew"))) #just barely significant interaction here...

cm_sep <- lmerCellMeans(mrew, divide="runningrtdev_lag", n.cont=10)
head(cm_sep)
head(subset(cm, omission_lag=="PriorRew"))
dev.new()
ggplot(cm, aes(x=vdevlag, y=rtchange, color=runningrtdev_lag)) + geom_line()


#results for prior rewards only
#Estimate Std. Error t value
#(Intercept)               1.421e+02  1.128e+01   12.60
#runningrtdev_lag         -3.167e-01  8.556e-03  -37.02
#vdevlag                  -1.685e-01  9.651e-03  -17.46
#runningrtdev_lag:vdevlag -1.147e-05  5.711e-06   -2.01

#tests of vdevlag x runningrtdev_lag separately for omission and reward 
cmat <- rbind(
    "vdevlag x rtdevlag for prior omissions"=c(
        0, #intercept
        0, #runningrtdev_lag
        0, #omission_lag (prior omission)
        0, #vdevlag
        0, #runningrtdev_lag x omission_lag
        1, #runningrtdev_lag x vdevlag
        0, #omission_lag x vdevlag
        1  #runningrtdev_lag x omission_lag x vdevlag
    ),
    "vdevlag x rtdevlag for prior rewards"=c(
        0, #intercept
        0, #runningrtdev_lag
        0, #omission_lag (prior omission)
        0, #vdevlag
        0, #runningrtdev_lag x omission_lag
        1, #runningrtdev_lag x vdevlag
        0, #omission_lag x vdevlag
        0  #runningrtdev_lag x omission_lag x vdevlag
    )
)

#test whether ME of prior RT deviation is larger for omissions that rewards (looks like it by eye)
cmat <- rbind(
    "rtdevlag effect after omissions"=c(
        0, #intercept
        1, #runningrtdev_lag
        0, #omission_lag (prior omission)
        0, #vdevlag
        1, #runningrtdev_lag x omission_lag
        0, #runningrtdev_lag x vdevlag
        0, #omission_lag x vdevlag
        0  #runningrtdev_lag x omission_lag x vdevlag
    ),
    "rtdevlag effect after rewards"=c(
        0, #intercept
        1, #runningrtdev_lag
        0, #omission_lag (prior omission)
        0, #vdevlag
        0, #runningrtdev_lag x omission_lag
        0, #runningrtdev_lag x vdevlag
        0, #omission_lag x vdevlag
        0  #runningrtdev_lag x omission_lag x vdevlag
    )
)

#model without vdevlag
#summary(lmer(rtchange ~ runningrtdev_lag*omission_lag + (1|lunaid), alldf))

library(multcomp)
summary(glht(test, linfct=cmat))

#confirms that there is a vdevlag x runningrtdevlag interaction for omissions, but not rewards

lsmeans(test, "runningrtdev_lag")

lsmeans(test, ~ vdevlag + runningrtdev_lag | omission_lag)

#l <- lsmeans(test, ~ vdevlag:runningrtdev_lag | omission_lag)
l <- lsmeans(test, ~ runningrtdev_lag | omission_lag)
summary(l, infer=TRUE, level=.90, adjust="bon")


#interpretation:
#  - 




#similar model with original RTs, not change. But harder to interpret because of adding lagged Rt effect alone (too many predictors to wrap my head around)
#test <- lmer(rts ~ rt_lag + runningrtdev_lag*omission_lag*vdevlag + (1|lunaid), alldf)
#car::Anova(test) #yes, all whopping effects...

#how related are prior RT deviation from average versus prior RT deviation from Vmax?
cor.test(~runningrtdev_lag + vdevlag, alldf) #r = .70
summary(lmer(runningrtdev_lag ~ vdevlag + (1|lunaid), alldf))

#do these effects depend on emotion?
test1 <- lmer(rtchange ~ runningrtdev_lag*vdevlag*omission_lag + contingency + cond + (1|lunaid), alldf)
summary(test1)
car::Anova(test1)

test2 <- lmer(rtchange ~ runningrtdev_lag*vdevlag*omission_lag + contingency*cond + (1|lunaid), alldf)
car::Anova(test2)

test3 <- lmer(rtchange ~ runningrtdev_lag*vdevlag*omission_lag*contingency + (1|lunaid), alldf)
test4 <- lmer(rtchange ~ runningrtdev_lag*vdevlag*omission_lag*cond + (1|lunaid), alldf)
car::Anova(test4)
test5 <- lmer(rtchange ~ runningrtdev_lag*vdevlag*omission_lag*cond + runningrtdev_lag*vdevlag*omission_lag*contingency + (1|lunaid), alldf)
car::Anova(test5)
test6 <- lmer(rtchange ~ runningrtdev_lag*vdevlag*omission_lag*contingency + cond + (1|lunaid), alldf)

anova(test1, test2, test3, test4, test5, test6)

car::Anova(test3)

#test5 best
pdf("figures/Effects of vdevlag omissionlag runningrtdevlag contingency and emotion on RT change.pdf", width=8, height=8)
cm <- lmerCellMeans(test5, divide="runningrtdev_lag", n.cont=10)
cm <- subset(cm, !(contingency %in% c("CEV", "CEVR") & cond %in% c("happy", "fear")))
ggplot(cm, aes(x=vdevlag, y=rtchange, color=runningrtdev_lag, ymin=rtchange-se, ymax=rtchange+se, shape=omission_lag)) + geom_line() + 
    geom_pointrange() + facet_grid(contingency~cond)
dev.off()

#look only at IEV and DEV (make a complete factorial design)
test5 <- lmer(rtchange ~ runningrtdev_lag*vdevlag*omission_lag*cond + runningrtdev_lag*vdevlag*omission_lag*contingency + (1|lunaid), subset(alldf, contingency %in% c("IEV", "DEV")))
car::Anova(test5)

pdf("figures/Effects of vdevlag omissionlag runningrtdevlag contingency and emotion on RT change IEV DEV ONLY.pdf", width=8, height=8)
cm <- lmerCellMeans(test5, divide="runningrtdev_lag", n.cont=10)
ggplot(cm, aes(x=vdevlag, y=rtchange, color=runningrtdev_lag, ymin=rtchange-se, ymax=rtchange+se, shape=omission_lag)) + geom_line() + 
    geom_pointrange() + facet_grid(contingency~cond)
dev.off()


#age effects...
test <- lmer(rtchange ~ runningrtdev_lag*omission_lag*vdevlag*age + (1|lunaid), alldf)
test <- lmer(rtchange ~ runningrtdev_lag*omission_lag*vdevlag*age.c + (1|lunaid), alldf)
car::Anova(test) #looks like prior RT deviation*age and prior value deviation*age effects


#remove runningrtdev_lag since I don't knows its interpretation for now (esp. when included with vdevlag
test <- lmer(rtchange ~ vdevlag*omission_lag*age.c + (1|lunaid), alldf)
summary(test)
car::Anova(test)

test <- lmer(rtchange ~ vdevlag*omission_lag + (1|lunaid), alldf)
summary(test)

pieceage <- lmer(rtchange ~ vdevlag*omission_lag*adult*age + (1|lunaid), alldf)
summary(test)
car::Anova(test)


#no 3-way interactions with age are significant (nor is 4-way)
afex::mixed(rtchange ~ runningrtdev_lag.c*omission_lag*vdevlag.c*age.c + (1|lunaid), alldf)

simpler <- lmer(rtchange ~ runningrtdev_lag*omission_lag*vdevlag + age.c + runningrtdev_lag:age.c + omission_lag:age.c + vdevlag:age.c + (1|lunaid), alldf)
simpler <- lmer(rtchange ~ runningrtdev_lag*omission_lag*vdevlag + age + runningrtdev_lag:age + omission_lag:age + vdevlag:age + (1|lunaid), alldf)
anova(simpler, test) #confirms lack of 3-way

pdf("Age effects on RT change.pdf", width=8, height=6)
cm <- lmerCellMeans(simpler, divide=c("vdevlag", "runningrtdev_lag"), n.cont=10)
ggplot(cm, aes(x=age, y=rtchange, color=vdevlag, shape=runningrtdev_lag)) + geom_point(size=5) + geom_line(size=1.5) + facet_wrap(~omission_lag)
dev.off()


#adolescent versus adult version
test <- lmer(rtchange ~ runningrtdev_lag*omission_lag*vdevlag*adult + (1|lunaid), alldf)
test <- lmer(rtchange ~ runningrtdev_lag.c*omission_lag*vdevlag.c*adult + (1|lunaid), alldf)
summary(test)
car::Anova(test)

afex::mixed(rtchange ~ runningrtdev_lag.c*omission_lag*vdevlag.c*adult.c + (1|lunaid), alldf)

#piecewise adolescent/adult age model
test3 <- lmer(rtchange ~ runningrtdev_lag*omission_lag*vdevlag*adult*age + (1|lunaid), alldf)
summary(test)
car::Anova(test)

#model comparison
linage <- lmer(rtchange ~ runningrtdev_lag*omission_lag*vdevlag*age + (1|lunaid), alldf)
binage <- lmer(rtchange ~ runningrtdev_lag*omission_lag*vdevlag*adult + (1|lunaid), alldf)
linbin_noint <- lmer(rtchange ~ runningrtdev_lag*omission_lag*vdevlag*age + runningrtdev_lag*omission_lag*vdevlag*adult + (1|lunaid), alldf)
linbin_int <- lmer(rtchange ~ runningrtdev_lag*omission_lag*vdevlag*adult*age + (1|lunaid), alldf)
splineage <- lmer(rtchange ~ runningrtdev_lag*omission_lag*vdevlag*ns(age,5) + (1|lunaid), alldf)

anova(linage, binage, linbin_noint, linbin_int, splineage)

#ugh, 5-way interaction?!
car::Anova(test4)

pdf("Age effects on RT change.pdf", width=8, height=6)
cm_teen <- lmerCellMeans(pieceage, divide=c("vdevlag"), cont.pts=list(age=c(14, 16, 18)))
cm_teen <- subset(cm_teen, adult=="adolescent")
cm_adult <- lmerCellMeans(pieceage, divide=c("vdevlag"), cont.pts=list(age=c(18, 20, 22, 24, 26, 28, 30, 32)))
cm_adult <- subset(cm_adult, adult=="adult")
#cm <- subset(cm, (adult=="adult" & age >= 18) | (adult=="adolescent" & age < 18))
ggplot(cm_adult, aes(x=age, y=rtchange, color=vdevlag)) + #shape=runningrtdev_lag) 
    geom_line(data=cm_teen, aes(x=age, y=rtchange, color=vdevlag)) + #, shape=runningrtdev_lag
    geom_point(size=5) + geom_line(size=1.5) + facet_wrap(~omission_lag)
dev.off()


library(splines)
str(ns(alldf$age, 5))

test1 <- lmer(rtchange ~ runningrtdev_lag*omission_lag*vdevlag*ns(age,5) + (1|lunaid), alldf)
summary(test1)
car::Anova(test1)

#try to get predicted values off of spline model
newdata <- data.frame()

newdata <- with(cbpp, expand.grid(period=unique(period), herd=unique(herd)))
str(p2 <- predict(gm1,newdata))    # new data, all RE
str(p3 <- predict(gm1,newdata,REform=NA)) # new data, level-0
str(p4 <- predict(gm1,newdata,REform=~(1|herd))) # explicitly specify RE


ggplot(cm_adult, aes(x=age, y=rtchange, color=vdevlag, shape=runningrtdev_lag)) + 
    geom_line(data=cm_teen, aes(x=age, y=rtchange, color=vdevlag, shape=runningrtdev_lag)) +
    geom_point(size=5) + geom_line(size=1.5) + facet_wrap(~omission_lag)

str()



cm <- lmerCellMeans(test2, divide="rt_lag", n.cont=10)
ggplot(cm, aes(x=vdevlag, y=rts, color=rt_lag, ymin=rts-se, ymax=rts+se)) + geom_line() + geom_pointrange()


#work in progress
#test3 <- lmer(rts ~ rt_lag + vdevlag*trial*contingency + udevlag*trial*contingency + trial*contingency + (1|lunaid), alldf) #, na.action=na.exclude)
#test3 <- lmer(rts ~ rt_lag + vdevlag + udevlag + (1|lunaid), alldf) #, na.action=na.exclude)
#test3 <- lmer(rts ~ rt_lag + trial*contingency + (1|lunaid), alldf) #, na.action=na.exclude)
#summary(test3)
#car::Anova(test3)
#
#cm <- lmerCellMeans(test3, divide=c("rt_lag", "contingency", "trial"))
#
#ggplot(cm, aes(x=udevlag, y=rts, ymax=rts+se, ymin=rts-se, color=rt_lag)) + facet_wrap(~vdevlag) + geom_pointrange(size=1.6) + geom_line()

#f <- fitted(test)
#alldf$f <- fitted(test)
#alldf$res <- residuals(test)
#
#subj1 <- subset(alldf, lunaid==10637)
#subj1 <- subj1[order(subj1$trial),]
#subj1 <- subset(subj1, trial <= 50)
#m <- melt(subj1[,c("trial", "rts", "f", "res")], id.vars="trial")
#ggplot(m, aes(x=trial, y=value, color=variable)) + geom_line()
#
##c, #divide="vdevlag", n.divide=5)
#ggplot(cm, aes(x=rt_lag, y=rts, color=vdevlag, ymin=rts-se, ymax=rts+se)) + geom_pointrange() + geom_line()


cm <- lmerCellMeans(test, cont.pts=list(rpes_lag=c(-100, -75, -50, -25, 0, 25, 50, 75, 100)))
cm_mod <- rbind(subset(cm, rpes_lag_pos=="Prior RPE-" & rpes_lag <= 0), subset(cm, rpes_lag_pos=="Prior RPE+" & rpes_lag >= 0))

#just look at emotion x contingency a prior RPE direction effect on delta RT
pdf("RT change by RPE.pdf", width=12, height=6)
ggplot(cm_mod, aes(x=rpes_lag, y=rtchange, color=cond)) + geom_line(size=2) + 
    geom_point(size=5) + facet_wrap(~rpes_lag_pos) + 
    geom_errorbar(aes(ymin=rtchange-se, ymax=rtchange+se), size=2, width=5) + geom_vline(xintercept=0) + geom_hline(yintercept=0) + 
    ylab("Change in RT (ms)") + xlab("Magnitude of Reward Prediction Error (points)") + theme_bw(base_size=24) + scale_color_brewer("Group", palette="Set1")
dev.off()


test <- lmer(rtchange ~ omission_lag*cond + (1|lunaid), alldf)
summary(test)

#after reward omissions, people speed up (gradient for happy > fear > scram)
pdf("Effect of emotional cue on RT change.pdf", width=6, height=5)
cm <- lmerCellMeans(test)
ggplot(cm, aes(x=omission_lag, y=rtchange, ymin=rtchange-se, ymax=rtchange+se, color=cond)) + geom_pointrange(size=1.2, position=position_dodge(width=0.2)) +
    xlab("Prior Outcome") + ylab("RT change (ms)") + scale_color_brewer("Condition", palette="Set2")
dev.off()

test <- lmer(rtchange ~ priorrtfast*cond + (1|lunaid), alldf)
summary(test)

#yes, if prior RT was fast, people slow down and vice versa (Moustafa 2008)
cm <- lmerCellMeans(test)
ggplot(cm, aes(x=priorrtfast, y=rtchange, ymin=rtchange-se, ymax=rtchange+se, color=cond)) + geom_pointrange(size=2, position=position_dodge(width=0.2))

#looking at running average RT (computed as average of RTs in run so far)
test <- lmer(rtchange ~ priorrtfast_running*cond + (1|lunaid), alldf)
summary(test)

#do both processes occurs (omission-related swings and prior RT fast/slow swings?)
#yes: people speed up much more if a) they were previously slow and b) a reward omission occurred
test <- lmer(rtchange ~ omission_lag*priorrtfast*cond + (1|lunaid), alldf)
summary(test)

car::Anova(test)

cm <- lmerCellMeans(test)
ggplot(cm, aes(x=omission_lag, y=rtchange, ymin=rtchange-se, ymax=rtchange+se, color=priorrtfast)) + geom_pointrange(size=2, position=position_dodge(width=0.2)) + facet_wrap(~cond)

test <- lmer(rtchange ~ omission_lag*priorrtfast_running*cond + (1|lunaid), alldf)
summary(test)

cm <- lmerCellMeans(test)
ggplot(cm, aes(x=omission_lag, y=rtchange, ymin=rtchange-se, ymax=rtchange+se, color=priorrtfast_running)) + geom_pointrange(size=2, position=position_dodge(width=0.2)) + facet_wrap(~cond)


withoutcond <- lmer(rtchange ~ omission_lag*priorrtfast_running + (1|lunaid), alldf)
anova(test, withoutcond) #condition is important to capture cond x omission effects

#look at previous max value
test <- lmer(rtchange ~ omission_lag*priorrtfast*cond*vdevlag+ (1|lunaid), alldf)
summary(test)
car::Anova(test)


#does it vary by contingency?
addconting <- lmer(rtchange ~ omission_lag*priorrtfast_running*cond*contingency + (1|lunaid), alldf)
summary(addconting)

car::Anova(addconting)

pdf("Effect of Contingency, Emotion, Prior Omission, and Prior Fast-Slow on RT change.pdf", width=12, height=10)
cm <- lmerCellMeans(addconting)

#incomplete nesting such that CEV and CEVR were only completed for scrambled
#           cond
#contingency fear happy scram
#CEV            0     0    50
#CEVR           0     0    50
#DEV           50    50    50
#IEV           50    50    50
cm <- subset(cm, !(cond %in% c("fear", "happy") & contingency %in% c("CEV", "CEVR")))

ggplot(cm, aes(x=omission_lag, y=rtchange, ymin=rtchange-se, ymax=rtchange+se, color=priorrtfast_running)) + geom_pointrange(size=0.8, position=position_dodge(width=0.2)) + 
    facet_grid(contingency~cond) + geom_hline(yintercept=0) + ylab("RT Change") + xlab("Prior Trial Outcome") + theme(panel.margin = unit(2, "lines"))
dev.off()


afex::mixed(rtchange ~ omission_lag*priorrtfast*cond*contingency + (1|lunaid), alldf)

#age effect
test <- lmer(rtchange ~ omission_lag*priorrtfast*age + (1|lunaid), alldf)
summary(test)
car::Anova(test)

cm <- lmerCellMeans(test, n.cont=10)
ggplot(cm, aes(x=age, color=omission_lag, y=rtchange, ymin=rtchange-se, ymax=rtchange+se, shape=priorrtfast)) + geom_pointrange(size=1.5, position=position_dodge(width=0.2)) + geom_line() + 
    geom_hline(yintercept=0) + ylab("RT Change") + xlab("Age")

test <- lmer(rtchange ~ omission_lag*priorrtfast*age*cond + (1|lunaid), alldf)
summary(test)
car::Anova(test)

cm <- lmerCellMeans(test, n.cont=10)
ggplot(cm, aes(x=age, color=omission_lag, y=rtchange, ymin=rtchange-se, ymax=rtchange+se, shape=priorrtfast)) + geom_pointrange(size=1.5, position=position_dodge(width=0.2)) + geom_line() + 
    geom_hline(yintercept=0) + ylab("RT Change") + xlab("Age") + facet_wrap(~cond)

library(plyr)
aggby <- ddply(alldf, .(lunaid, omission_lag, priorrtfast, cond), summarise, rtchange=mean(rtchange, na.rm=TRUE), age=age[1])
aggby <- subset(aggby, !is.na(omission_lag)) #remove rows where lagged omission was missing (since undefined on first trial)

#table(firstsubj$cond, firstsubj$omission_lag, firstsubj$priorrtfast)
ggplot(aggby, aes(x=age, y=rtchange, color=omission_lag, shape=priorrtfast)) + geom_point() + stat_smooth(size=4) + facet_wrap(~cond)


##New variant in Sep2016: entropy
library(R.matlab)
library(reshape2)
vars <- readMat("/Users/michael/Data_Analysis/temporal_instrumental_agent/entropy_analysis.mat")
str(vars)

rt <- melt(vars$rt.decay, varnames=c("subject", "trial"), value.name="rt") #76 x 400
vmax <- melt(vars$value.max.decay, varnames=c("subject", "trial"), value.name="vmax") #76 x 400
entropy <- melt(vars$S.decay, varnames=c("subject", "trial"), value.name="entropy") #76 x 400

library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
df <- full_join(rt, vmax, by=c("subject", "trial"))
df <- full_join(df, entropy, by=c("subject", "trial")) %>% arrange(subject, trial)

mdf <- df %>% gather(key=signal, value=value, -subject, -trial)

plotsubj <- function(data) {
#  browser()
  g <- ggplot(data, aes(x=trial, y=value)) + geom_line() + facet_wrap(~signal, ncol=1, scales="free_y") #+ ggtitle(paste0("ID: ", data$subject[1]))
  plot(g)
  return(g)
}

pdf("entropy diagnostics.pdf", width=11, height=8)
xx = mdf %>% group_by(subject) %>% nest() %>% mutate(obj = map(data, plotsubj))

#ddply(   { ggplot(., aes(x=trial, y=value)) + geom_line() + facet_wrap(~signal, ncol=1, scales="free_y") })
dev.off()
#diagnostic plots


