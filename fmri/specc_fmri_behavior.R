#mlms of SPECC final fMRI data

setwd("/Users/michael/Data_Analysis/temporal_instrumental_agent/clock_task/vba_fmri")
#load(file="dataframe_for_entropy_analysis_SPECC_Dec2016.RData")
load(file="dataframe_for_entropy_analysis_specc_Mar2017.RData")
library(lme4)
library(ggplot2)
library(dplyr)
library(tidyr)
library(psych)
library(gdata)
source(file.path(getMainDir(), "Miscellaneous", "Global_Functions.R"))
library(R.matlab)

specc_info <- gdata::read.xls("/Users/michael/Box_Sync/DEPENd/Projects/SPECC/ID Management/SPECC_Participant_Info.xlsx")
specc_info$ID <- as.numeric(substr(specc_info$SPECC_ID, 1, 3)) #for merge with bdf
specc_info$group <- factor(specc_info$BPD, levels=c(0,1), labels=c("Control", "BPD"))
specc_info$sex <- factor(specc_info$Female, levels=c(0,1), labels=c("Male", "Female"))
bdf <- left_join(bdf, specc_info, by="ID")

bdf$abspe <- abs(bdf$pemax)
bdf$emotion <- relevel(bdf$emotion, ref="scram")
bdf <- filter(bdf, !ID==23) #exclude participant who lied about age
bdf$iage_c <- (1/bdf$AgeAtScan - mean(1/bdf$AgeAtScan))*100

#try to build a simple model of good choice

#here's the old emotion model
summary(memo <- lmer(rt ~ rtlag + rtlag2 + rtlag3 + emotion + rtvmaxlag + abspelag + rewFunc + 
            emotion*omissionlag + emotion*abspelag + omissionlag*abspelag + (1 | ID) + (1 | run) + (1 | run:ID), filter(bdf, rewFunc %in% c("IEV", "DEV") & !ID==11282), REML=FALSE))


#simpler updated model with bpd in the mix
summary(memo <- lmer(rt ~ rtlag + rtlag2 + emotion*group + rtvmaxlag + abspelag + rewFunc + 
            emotion*omissionlag + emotion*abspelag*group + (1 | ID) + (1 | run) + (1 | run:ID), filter(bdf, rewFunc %in% c("IEV", "DEV") & !ID==11282), REML=FALSE))

car::Anova(memo)


summary(memo <- lmer(rt ~ rtlag + rtlag2 + emotion*group + rtvmaxlag*group*emotion + rewFunc*group + 
            + emotion*ppelag*group + emotion*npelag*group + (1 | ID/run), bdf, REML=FALSE))

car::Anova(memo)

#look under the hood a bit at value reversion increase in BPD
summary(msimple <- lmer(rt ~ rtlag + rtlag2 + emotion*group + rtvmaxlag*group*emotion + rewFunc*group + entropy*group + (1 | ID/run), bdf, REML=FALSE))
car::Anova(msimple)

lstrends(msimple, ~ group , var="rtvmaxlag")
lstrends(msimple, ~ group | emotion, var="rtvmaxlag") #just barely an interaction p = .05 in the Anova when rewFunc 3-way entered
pairs(lstrends(msimple, ~ group | emotion, var="rtvmaxlag"))

pairs(lstrends(msimple, ~ group , var="entropy"))

#huh, is entropy different in BPD? no, it appears not, so it's the relationship of entropy to choice
summary(mmm <- lmer(entropy ~ emotion*group + rewFunc*group + (1 | ID/run), bdf, REML=FALSE))
car::Anova(mmm)


##what about aggregate performance

bdfruns = bdf %>% group_by(SPECC_ID, run) %>% filter(entropy > 0) %>% arrange(SPECC_ID, run, trial) %>% dplyr::summarize(runreward=sum(score), 
        earlyentropy=mean(entropy[trial > 1 & trial < 10])/mean(entropy), lateentropy=mean(entropy[trial > 40 & trial <= 50])/mean(entropy),
        allentropy=mean(entropy), midentropy=mean(entropy[trial > 10 & trial <=40]/mean(entropy)),
        #earlyentropyF=mean(entropyFixed[trial > 1 & trial < 10])/mean(entropyFixed), lateentropyF=mean(entropyFixed[trial > 40 & trial <= 50])/mean(entropyFixed),
        #allentropyF=mean(entropyFixed), midentropyF=mean(entropyFixed[trial > 10 & trial <=40])/mean(entropyFixed),
        elratio=earlyentropy/lateentropy, rewFunc=head(rewFunc, n=1), emotion=head(emotion, n=1), group=head(group, n=1), AgeAtScan=head(AgeAtScan, n=1)) %>%
    group_by(SPECC_ID) %>% #elratioF=earlyentropyF/lateentropyF, 
    mutate(elratioLag = lag(elratio, order_by=run)) %>% ungroup() #elratioLagF = lag(elratioF, order_by=run

bdfruns$AgeAtScan.c <- bdfruns$AgeAtScan - mean(bdfruns$AgeAtScan)

summary(mm <- lmer(runreward ~ AgeAtScan.c*group + run + rewFunc*group*emotion + (1 | SPECC_ID), bdfruns))
car::Anova(mm)

#runlevel entropy?
summary(mm <- lmer(runreward ~ AgeAtScan.c*group + run + group*emotion + group*rewFunc + earlyentropy*emotion*group + (1 | SPECC_ID), bdfruns))
car::Anova(mm)

pairs(lsmeans(mm, ~ group | emotion)) #close to significant!
lsmeans(mm, ~ group | emotion)

#not much purchase on group effects here, but what if we revert back to quality of choice in trial data?
summary(mmtrial <- lmer(ev ~ emotion*group + rewFunc*group + AgeAtScan.c + run + omissionlag*group*emotion + abspelag*group*emotion + (1 | SPECC_ID/run), bdf))
car::Anova(mmtrial)

#some emerging evidence of better? quality decisions in BPD under specific emotion x omission and emotion x PE effects
#pairs(lsmeans(mmtrial, ~ emotion*omissionlag )) 

#unpack group x emotion x omission effect
#so: after an omission in fear, controls choose a higher value optoin compared to BPDs
#trend/weak: after reward in fear, controls tend to choose lower value option
#so maybe more plasticity in fear conditions for controls than bpd
pairs(lsmeans(mmtrial, ~ group | emotion + omissionlag )) 
lsmeans(mmtrial, ~ emotion | group + omissionlag) 

pairs(lsmeans(mmtrial, ~ group | emotion + abspelag )) 

#conclusion: after large PE, shift toward good choices is larger in controls than BPD for fear only
#this seems to dovetail with omission effect above
pairs(lstrends(mmtrial, ~ group | emotion, var="abspelag"))

#so both wrt PEs and omissions, BPD tends to show poorer adaptation in fear
#what if we go back to signed PE alone...
summary(mmtrial2 <- lmer(ev ~ emotion*group + rewFunc*group + AgeAtScan.c + run + pemaxlag*group*emotion + (1 | SPECC_ID/run), bdf))
car::Anova(mmtrial2)

#separate into positive and negative PEs
summary(mmtrial2 <- lmer(ev ~ emotion*group + rewFunc*group + AgeAtScan.c + run + ppelag*group*emotion + (1 | SPECC_ID/run), bdf))
car::Anova(mmtrial2)

summary(mmtrial2 <- lmer(ev ~ emotion*group + rewFunc*group*ppelag + rewFunc*group*npelag + AgeAtScan.c + run + ppelag*group*emotion + npelag*group*emotion + (1 | SPECC_ID/run), bdf))
car::Anova(mmtrial2)

pairs(lstrends(mmtrial2, ~ group | emotion, var="ppelag"))
pairs(lstrends(mmtrial2, ~ group | emotion, var="npelag"))

#looks like there are also emotion condition differences within group
pairs(lstrends(mmtrial2, ~ emotion | group, var="ppelag"))
pairs(lstrends(mmtrial2, ~ emotion | group, var="npelag"))



plot(lstrends(mmtrial2, ~ group | emotion, var="ppelag", infer=c(TRUE, TRUE)))
plot(lstrends(mmtrial2, ~ group | emotion, var="npelag", infer=c(TRUE, TRUE)))

ppeeff <- summary(lstrends(mmtrial2, ~ group | emotion, var="ppelag", infer=c(TRUE, TRUE)))
npeeff <- summary(lstrends(mmtrial2, ~ group | emotion, var="npelag", infer=c(TRUE, TRUE)))

both <- rbind(ppeeff %>% rename(trend=ppelag.trend) %>% mutate(eff="Positive PE"), npeeff %>% rename(trend=npelag.trend) %>% mutate(eff="Negative PE"))

pdf("PE shifts toward high EV BPD diff.pdf", width=11, height=8)
ggplot(both, aes(x=emotion, y=trend, ymin=trend-SE, ymax=trend+SE, color=group)) + facet_wrap(~eff, scales="free_y") + 
    geom_pointrange(position=position_dodge(width=0.5), size=1.6) +
    ylab("Shift toward high-value option after PE") + xlab("Emotion Stimulus") + scale_color_brewer("Group", palette="Set1") + scale_x_discrete(labels=c("Scrambled", "Fear", "Happy")) +
    theme_bw(base_size=24) + theme(axis.title.x=element_text(margin=margin(t=20, r=0,l=0, b=0)), axis.title.y=element_text(margin=margin(t=0, r=15, l=0, b=0)))
dev.off()

lstrends(mmtrial2, ~ rewFunc | group, var="ppelag")

#shouldn't we just run separate models for trials following NPEs versus PPEs? Seems way simpler.... though the above is effective




bdfagg = bdf %>% group_by(SPECC_ID) %>% filter(entropy > 0) %>% #zero entropy indicative of very early learning (not valid)
    arrange(SPECC_ID, run, trial) %>% dplyr::summarize(runreward=sum(score), 
        earlyentropy=mean(entropy[trial > 1 & trial < 10])/mean(entropy), lateentropy=mean(entropy[trial > 40 & trial <= 50])/mean(entropy),
        allentropy=mean(entropy), midentropy=mean(entropy[trial > 10 & trial <=40]/mean(entropy)),
        #earlyentropyF=mean(entropyFixed[trial > 1 & trial < 10])/mean(entropyFixed), lateentropyF=mean(entropyFixed[trial > 40 & trial <= 50])/mean(entropyFixed),
        #allentropyF=mean(entropyFixed), midentropyF=mean(entropyFixed[trial > 10 & trial <=40])/mean(entropyFixed),
        elratio=earlyentropy/lateentropy, group=head(group, n=1), AgeAtScan=head(AgeAtScan, n=1)) #elratioF=earlyentropyF/lateentropyF, 

bdfagg$AgeAtScan.c <- bdfagg$AgeAtScan - mean(bdfagg$AgeAtScan)

summary(mval <- lm(runreward ~ earlyentropy + lateentropy, bdfagg))
car::Anova(mval)

bdfagg$elratio_wins <- winsor(bdfagg$elratio, trim=0.02) #trim a few outliers
#bdfagg$elratioF_wins <- winsor(bdfagg$elratioF, trim=0.02) #trim a few outliers
bdfagg$runreward_wins <- winsor(bdfagg$runreward, trim=0.02) #trim a few outliers

cor.test(~ runreward + elratio, bdfagg) 
cor.test(~ runreward_wins + elratio_wins, bdfagg)
 
summary(lm(runreward_wins ~ elratio_wins*group, bdfagg))

summary(lm(runreward_wins ~ AgeAtScan.c, bdfagg))
summary(lm(runreward_wins ~ elratio_wins + AgeAtScan.c, bdfagg))

#does elratio act as mediator of age effect?
#nope
library(lavaan)
m1 <- '
    runreward_wins ~ c*AgeAtScan.c + b*elratio_wins 
    elratio_wins ~ a*AgeAtScan.c 
    
    ide := a*b
    total := c + a*b
    '

fsem <- sem(m1, data=bdfagg, bootstrap=1000, se="bootstrap")
summary(fsem, fit.measures=TRUE, standardize=TRUE, rsquare=TRUE)



cor(select(bdfagg, AgeAtScan, elratio))
##leftovers

m1 <- lmer(rt ~ 1 + rtlag + rtlag2 + emotion*abspe + omissionlag + rewFunc + rtvmaxlag + #condition_num + emotion*rewardlag + abspe*rewardlag + rtvmaxlag + rtlag + rtlag2 +rtlag3 + emotion*abspe +
        (1|ID) + (1|run), bdf)

summary(m1)

m2 <- lmer(rt ~ 1 + rtlag + rtlag2 + rtlag3 + emotion*group*omissionlag*abspe*iage_c + rewFunc + rtvmaxlag + #condition_num + emotion*rewardlag + abspe*rewardlag + rtvmaxlag + rtlag + rtlag2 +rtlag3 + emotion*abspe +
        (1 + omissionlag |ID) + (1|run), filter(bdf, rewFunc %in% c("IEV", "DEV")))

summary(m2)

car::Anova(m2)