#kai check
sig <- read.csv("/Users/mnh5174/Downloads/db_testcsv.csv")
library(dplyr)
library(lme4)
library(ggplot2)
filter(sig, Pow < -100)
str(sig)
mm <- lmer(Pow ~ Faces  + Age + Faces*Age + Trial + Pe + Pe*Faces + (0 + Trial | Subject) + (1 | Subject/Run), sig)

mm2 <- lmer(Pow ~ Faces  + Age + Faces*Age + Trial + Pe + Pe*Faces + (0 + Trial | Subject) + (1 | Subject:Run), sig)
mm3 <- lmer(Pow ~ Faces  + Age + Faces*Age + Trial + Pe + Pe*Faces + (1 | Subject/Run), sig)
summary(mm3)

mm4 <- lmer(Pow ~ 1 + (1 | Subject/Run), sig %>% filter (Pow > -300), REML=FALSE)
summary(mm4)

mm9 <- lmer(Pow ~ 1 + Trial + (0 + Trial | Subject) + (1 | Subject/Run), sig %>% filter (Pow > -300))
summary(mm9)


mm4_all <- lmer(Pow ~ 1 + (1 | Subject/Run), sig)
summary(mm4_all)


mm5 <- lmer(Pow ~ 1 + (1 | Run), sig %>% filter (Pow > -300))
summary(mm5)

mm6 <- lmer(Pow ~ 1 + (1 | Subject:Run), sig %>% filter (Pow > -300))
summary(mm6)

anova(mm4, mm5, mm6)

anova(mm,mm2)
logLik(mm)
logLik(mm2)


xtabs(~Run+Subject, sig)
runmeans <- sig %>% group_by(Subject, Run) %>% summarize(pow_runmean=mean(Pow)) %>% 
  group_by(Subject) %>% mutate(pow_personmean=mean(pow_runmean)) %>% ungroup()

ggplot(runmeans, aes(x=pow_runmean)) + geom_histogram()


pmeans <- sig %>% group_by(Subject) %>% summarize(pow_personmean=mean(Pow)) %>% ungroup()
ggplot(pmeans, aes(x=pow_personmean)) + geom_histogram()
