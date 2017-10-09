setwd(file.path(getMainDir(), "temporal_instrumental_agent/clock_task/vba_fmri"))
#load(file="dataframe_for_entropy_analysis_Oct2016.RData")
#this contains data with 24 basis functions and post-Niv learning rule
#load(file="dataframe_for_entropy_analysis_Nov2016.RData")
load(file="dataframe_for_entropy_analysis_Mar2017.RData") #has the random priors entropy
library(lme4)
library(ggplot2)
library(dplyr)
library(tidyr)
library(psych)
library(gdata)
source(file.path(getMainDir(), "Miscellaneous", "Global_Functions.R"))
library(R.matlab)

#bring in Luna IDs
idlist <- read.xls("id list.xlsx", header=FALSE, col.names="LunaID")
idlist$subject <- 1:nrow(idlist) #numeric ID

bdf <- bdf %>% rename(subject=rowID)
bdf <- left_join(bdf, idlist, by="subject")

#compute total rewards and cumulative rewards by subject
bdf = bdf %>% group_by(subject) %>% arrange(subject, run, trial) %>% mutate(totreward=sum(score), cumreward=cumsum(score)) %>% ungroup() %>%
		mutate(medreward=median(totreward), #between subjects
				msplit=factor(as.numeric(totreward > medreward), levels=c(0,1), labels=c("< Median", ">= Median")))


##11282 is at the floor of responding on all runs -- invalid, exclude
bdf <- filter(bdf, LunaID != 11282)

#plot average entropy for fixed and decay models
edescriptives <- bdf %>% 
  dplyr::select(LunaID, run, trial, entropyH, entropyFixed) %>% 
  gather(key=Model, value=entropy, entropyFixed, entropyH) %>% 
  mutate(Model=recode(Model, entropyFixed="Fixed LR V", entropyH="Fixed LR V Sel. Maint."))

edesc_aggruns <- edescriptives %>% filter(run > 1) %>%
  group_by(LunaID, Model, trial) %>% dplyr::summarize(entropy=mean(entropy))

ggplot(edesc_aggruns, aes(x=trial, y=entropy, color=Model)) + 
  stat_smooth(size=2) + theme_bw(base_size=25)

summary((mdiff = lmer(entropy ~ Model + trial + (1|LunaID) + (1|run), edescriptives)))
car::Anova(mdiff)

#just for merging
mdf <- bdf %>% group_by(LunaID) %>% summarize(msplit=dplyr::first(msplit))

#This just gets the 95% bootstrapped CIs. Can't put a multi-return call in summarize
edesc_aggtrialsruns <- edescriptives %>% filter(run > 1) %>% 
  group_by(Model, trial) %>% do({data.frame(rbind(Hmisc::smean.cl.boot(.$entropy))) } )  %>% 
  ungroup() #%>% left_join(mdf, by="LunaID")
    #summarize(m=mean(value), se=plotrix::std.error(value))
    #, b=list(Hmisc::smean.cl.boot(value))) %>% unnest()

edesc_aggtrialsruns_msplit <- edescriptives %>% filter(run > 1) %>% left_join(mdf, by="LunaID") %>% 
    group_by(Model, msplit, trial) %>% do({data.frame(rbind(Hmisc::smean.cl.boot(.$entropy))) } )  %>% 
    ungroup() %>% mutate(msplit=recode(msplit, "< Median" = "Total~earnings<median", ">= Median"="Total~earnings>=median"))

pdf("entropy curves raw means.pdf", width=5, height=3.5)
#ggplot(filter(edesc_aggtrialsruns), aes(x=trial, y=m, ymin=m-se, ymax=m+se, color=Model)) +
ggplot(edesc_aggtrialsruns, aes(x=trial, y=Mean, ymin=Lower, ymax=Upper, color=Model)) +
    theme_bw(base_size=16) + xlab("Trial") + ylab("Entropy") +
    geom_line(size=1.5) + geom_ribbon(aes(fill=Model, color=NULL), alpha=0.3) +
  #geom_ribbon(aes(color=NULL, fill=Model), alpha=0.3)
    scale_color_brewer("Model", palette="Dark2") + scale_fill_brewer("Model", palette="Dark2") +
    theme(legend.position = c(0.75, 0.5), legend.background = element_rect(color = "grey40", 
            fill = "grey98", size = 0.5, linetype = "solid")) +
    theme(axis.title.y=element_text(margin=margin(r=15)), axis.title.x=element_text(margin=margin(t=10)))
dev.off()

pdf("entropy curves raw means inset.pdf", width=3.5, height=2.5)
#ggplot(filter(edesc_aggtrialsruns), aes(x=trial, y=m, ymin=m-se, ymax=m+se, color=Model)) +
ggplot(edesc_aggtrialsruns, aes(x=trial, y=Mean, ymin=Lower, ymax=Upper, color=Model)) +
    theme_bw(base_size=13) + xlab("Trial") + ylab("Entropy") +
    geom_line(size=1.5, show.legend=FALSE) + geom_ribbon(aes(color=NULL, fill=Model), alpha=0.3, show.legend=FALSE) +
    scale_color_brewer("Model", palette="Dark2") + scale_fill_brewer("Model", palette="Dark2") +
    theme(legend.position = c(0.75, 0.5), legend.background = element_rect(color = "grey40", 
            fill = "grey98", size = 0.5, linetype = "solid")) +
    theme(axis.title.y=element_text(margin=margin(r=15)), axis.title.x=element_text(margin=margin(t=10))) +
    ggtitle("Averaged over runs (2-8)")
dev.off()

pdf("entropy curves raw means inset msplit.pdf", width=7.9, height=2.5)
#ggplot(filter(edesc_aggtrialsruns), aes(x=trial, y=m, ymin=m-se, ymax=m+se, color=Model)) +
ggplot(edesc_aggtrialsruns_msplit, aes(x=trial, y=Mean, ymin=Lower, ymax=Upper, color=Model)) +
    theme_bw(base_size=16) + xlab("Trial within run") + ylab("Entropy") +
    geom_line(size=1.5, show.legend=FALSE) + geom_ribbon(aes(color=NULL, fill=Model), alpha=0.3, show.legend=FALSE) +
    scale_color_brewer("Model", palette="Dark2") + scale_fill_brewer("Model", palette="Dark2") +
    theme(legend.position = c(0.75, 0.5), legend.background = element_rect(color = "grey40", 
            fill = "grey98", size = 0.5, linetype = "solid")) +
    theme(axis.title.y=element_text(margin=margin(r=15)), axis.title.x=element_text(margin=margin(t=10))) +
    facet_wrap(~msplit, labeller=label_parsed) +
    #ggtitle("Averaged over runs (2-8)") 
    theme(strip.background = element_blank(), 
        strip.text = element_text(face = "bold", size=12, margin=margin(t=0, b=0)),
        plot.margin = margin(l=7, r=7, t=0, b=6))
#strip.placement = "outside", 
dev.off()


#what about plotting all trials? (No run aggregation and keep run 1)
#we should probably drop true zeros... doesn't seem to have much effect: filter(entropy > 0) %>%
edesc_aggtrials <- edescriptives %>% group_by(Model, run, trial) %>% do({data.frame(rbind(Hmisc::smean.cl.boot(.$entropy))) } )  %>% ungroup() %>%
    arrange(Model, run, trial) %>% mutate(trialabs=rep(1:400,2))

pdf("entropy curves raw means alltrials.pdf", width=10.25, height=5)
ggplot(filter(edesc_aggtrials, trialabs > 1), aes(x=trialabs, y=Mean, ymin=Lower, ymax=Upper, color=Model)) +
#ggplot(edesc_aggtrials, aes(x=trialabs, y=Mean, ymin=Lower, ymax=Upper, color=Model)) +
    theme_bw(base_size=16) + xlab("Trial") + ylab("Entropy") +
    geom_line(size=1.5) + geom_ribbon(aes(color=NULL, fill=Model), alpha=0.3) +
    scale_color_brewer("Model", palette="Dark2") + scale_fill_brewer("Model", palette="Dark2") +
    #theme(legend.position = c(0.75, 0.3), legend.background = element_rect(color = "grey40", 
    #        fill = "grey98", size = 0.5, linetype = "solid")) +
    theme(axis.title.y=element_text(margin=margin(r=15)), axis.title.x=element_text(margin=margin(t=10))) +
    geom_vline(xintercept=seq(50, 350, 50), color="gray60")
dev.off()

subset(edesc_aggtrials, trialabs < 10)



pdf("entropy curves with spaghetti.pdf", width=8, height=6)
ggplot(filter(edescriptives, run > 1), aes(x=trial, y=entropy, color=Model)) + stat_smooth(size=2, alpha=0.2, method="gam", formula = y ~ s(x), size = 1) + theme_bw(base_size=20) + xlab("Trial") + ylab("Entropy") + 
    #geom_line(data=subset(edesc_aggruns), aes(x=trial, y=entropy, color=Model, group=interaction(LunaID, Model)), size=0.4, alpha=0.3)
    stat_smooth(data=edesc_aggruns, aes(x=trial, y=entropy, color=Model, group=interaction(LunaID, Model)), size=0.4, alpha=0.1, se=FALSE, method="loess")
    #geom_jitter(data=edesc_aggruns, size=0.4, alpha=0.2, width=0.5) + 
#    layer(mapping=NULL,
#        data=edesc_aggruns,
#        geom=geom_point(alpha=0.1),
#        stat="identity", position="identity")
dev.off()

pdf("entropy curves gam.pdf", width=8, height=6)
ggplot(filter(edescriptives, run > 1), aes(x=trial, y=entropy, color=Model)) + theme_bw(base_size=20) + xlab("Trial") + ylab("Entropy") +
    stat_smooth(size=2, alpha=0.2, method="gam", formula = y ~ s(x), method.args=list(method="REML"), size = 1)
    #stat_smooth(size=2, alpha=0.2, method="loess")
    #geom_line(data=subset(edesc_aggruns), aes(x=trial, y=entropy, color=Model, group=interaction(LunaID, Model)), size=0.4, alpha=0.3)
#    stat_smooth(data=edesc_aggruns, aes(x=trial, y=entropy, color=Model, group=interaction(LunaID, Model)), size=0.4, alpha=0.1, se=FALSE, method="loess")
#geom_jitter(data=edesc_aggruns, size=0.4, alpha=0.2, width=0.5) + 
#    layer(mapping=NULL,
#        data=edesc_aggruns,
#        geom=geom_point(alpha=0.1),
#        stat="identity", position="identity")
dev.off()


ggplot(subset(edesc_aggruns, LunaID == 10637), aes(x=trial, y=entropy, color=Model, group=Model)) + geom_line() + geom_point()

wicenter <- c("entropyHlag", "distfromedgelag", "evdevlag", "vdevlag", "abstschangelag", "ev", "rtumax", "rtumaxlag", "rtvmaxlag") #predictors to center within subject and run
predictors <- c(wicenter, "trial", "omissionlag", "ev") #full set of predictors (GM center)

#compute entropy of trials 2-10 and 41-50 versus the run-average mean (normalize to remove between run and between subjects effects)
bdfruns = bdf %>% group_by(subject, run) %>% filter(entropyH > 0) %>% arrange(subject, run, trial) %>% dplyr::summarize(runreward=sum(score), 
    earlyentropy=mean(entropyH[trial > 1 & trial < 10])/mean(entropyH), lateentropy=mean(entropyH[trial > 40 & trial <= 50])/mean(entropyH),
    allentropy=mean(entropyH), midentropy=mean(entropyH[trial > 10 & trial <=40]/mean(entropyH)),
    earlyentropyF=mean(entropyFixed[trial > 1 & trial < 10])/mean(entropyFixed), lateentropyF=mean(entropyFixed[trial > 40 & trial <= 50])/mean(entropyFixed),
    allentropyF=mean(entropyFixed), midentropyF=mean(entropyFixed[trial > 10 & trial <=40])/mean(entropyFixed),
    elratioF=earlyentropyF/lateentropyF, elratio=earlyentropy/lateentropy, rewFunc=head(rewFunc, n=1), emotion=head(emotion, n=1)) %>% group_by(subject) %>%
    mutate(elratioLag = lag(elratio, order_by=run), elratioLagF = lag(elratioF, order_by=run)) %>% ungroup()

#ggplot(bdfruns, aes(x=factor(run), y=earlyentropy)) + geom_boxplot()
ggplot(bdfruns, aes(x=interaction(rewFunc, emotion, run), y=earlyentropy)) + geom_boxplot()
xtabs(~rewFunc + run, bdfruns) #run order by rewFunc

df1 <- bdfruns %>% dplyr::select(subject, run, elratio) %>%
    mutate(run=paste0("r", run)) %>%
    spread(key=run, value=elratio, convert=TRUE) %>%
    left_join(dplyr::select(bdfruns, subject, run, runreward) %>% group_by(subject) %>% dplyr::summarize(totreward=sum(runreward)), by=c("subject"))  

corstarsl(df1)
summary(lm(totreward ~ r1, df1))
summary(lm(totreward ~ r2, df1))
summary(lm(totreward ~ r3, df1))
summary(lm(totreward ~ r4, df1))
summary(lm(totreward ~ r5, df1))
summary(lm(totreward ~ r6, df1))
summary(lm(totreward ~ r7, df1))
summary(lm(totreward ~ r8, df1))

# %>% filter(run > 1) 
#bdfagg = bdf %>% group_by(subject) %>% filter(entropyH > 0) %>% #zero entropy indicative of very early learning (not valid) 
bdfagg = bdf %>% group_by(LunaID) %>% filter(entropyH > 0) %>% #zero entropy indicative of very early learning (not valid)
    arrange(subject, run, trial) %>% dplyr::summarize(runreward=sum(score), 
    earlyentropy=mean(entropyH[trial > 1 & trial < 10])/mean(entropyH), lateentropy=mean(entropyH[trial > 40 & trial <= 50])/mean(entropyH),
    allentropy=mean(entropyH), midentropy=mean(entropyH[trial > 10 & trial <=40]/mean(entropyH)),
    earlyentropyF=mean(entropyFixed[trial > 1 & trial < 10])/mean(entropyFixed), lateentropyF=mean(entropyFixed[trial > 40 & trial <= 50])/mean(entropyFixed),
    allentropyF=mean(entropyFixed), midentropyF=mean(entropyFixed[trial > 10 & trial <=40])/mean(entropyFixed),
    elratioF=earlyentropyF/lateentropyF, elratio=earlyentropy/lateentropy)

#median instead of mean (no big deal)
#bdfagg = bdf %>% filter(run > 1) %>% group_by(subject) %>% arrange(subject, run, trial) %>% dplyr::summarize(runreward=sum(score), 
#    earlyentropy=median(entropyH[trial > 1 & trial < 10])/median(entropyH), lateentropy=median(entropyH[trial > 40 & trial <= 50])/median(entropyH),
#    allentropy=median(entropyH), midentropy=median(entropyH[trial > 10 & trial <=40]/median(entropyH)),
#    earlyentropyF=median(entropyFixed[trial > 1 & trial < 10])/median(entropyFixed), lateentropyF=median(entropyFixed[trial > 40 & trial <= 50])/median(entropyFixed),
#    allentropyF=median(entropyFixed), midentropyF=median(entropyFixed[trial > 10 & trial <=40])/median(entropyFixed),
#    elratioF=earlyentropyF/lateentropyF, elratio=earlyentropy/lateentropy)

    
    
summary(mval <- lm(runreward ~ earlyentropy + lateentropy, bdfagg))
car::Anova(mval)


bdfagg$elratio_wins <- winsor(bdfagg$elratio, trim=0.02) #trim a few outliers
bdfagg$elratioF_wins <- winsor(bdfagg$elratioF, trim=0.02) #trim a few outliers
bdfagg$runreward_wins <- winsor(bdfagg$runreward, trim=0.02) #trim a few outliers

ggplot(bdfagg, aes(x=runreward, y=elratio)) + geom_point()
ggplot(bdfagg, aes(x=runreward_wins, y=runreward_wins)) + geom_point()
ggplot(bdfagg, aes(x=runreward_wins, y=elratioF)) + geom_point()

bothmodels <- bdfagg %>% dplyr::select(runreward_wins, runreward_wins, elratio_wins, elratioF_wins) %>% 
    gather(key="model", value="elratio", elratioF_wins, elratio_wins) %>% mutate(model=recode(model, elratioF_wins="Fixed LR V", elratio_wins="Fixed LR V Sel. Maint."))

pdf("Beauty.pdf", width=10, height=5)
ggplot(bothmodels, aes(x=runreward_wins, y=elratio)) + geom_point() + facet_wrap(~model) + stat_smooth(se=FALSE, method="lm") +
    geom_hline(yintercept=1, color="gray50") +
    theme_bw(base_size=20) +
    ylab("Early:Late Entropy Ratio") + xlab("Total points earned in task") +
    theme(axis.title.y=element_text(margin=margin(r=15)), axis.title.x=element_text(margin=margin(t=15))) +
    geom_text(data=data.frame(runreward_wins=13800, elratio=1.38, model=c("Fixed LR V", "Fixed LR V Sel. Maint."), text=c("italic(r) == -.16", "italic(r) == .56")), 
        aes(label=text), parse=TRUE, size=7)
dev.off()


library(robust)
covRob(dplyr::select(bdfagg, runreward, elratio), corr=TRUE)
covRob(dplyr::select(bdfagg, runreward, elratioF), corr=TRUE)

cor.test(~ runreward + elratio, bdfagg) 
cor.test(~ runreward_wins + elratio_wins, bdfagg)
cor.test(~ runreward + elratioF, bdfagg)
cor.test(~ runreward_wins + elratioF_wins, bdfagg)


cor.test(~ runreward + earlyentropy, bdfagg)
cor.test(~ runreward + lateentropy, bdfagg)
cor.test(~ runreward + earlyentropyF, bdfagg)
cor.test(~ runreward + lateentropyF, bdfagg)

#pull in age for a brief jaunt
#ran the top chunk of sceptic_external_correlates.R (should really just cache, but being lazy)
bdfagg <- left_join(bdfagg, df, by=c(LunaID="lunaid"))

cor.test(~runreward_wins + age, bdfagg)
plot(~runreward_wins + age, bdfagg)
ggplot(bdfagg, aes(x=age, y=runreward_wins)) + geom_point()

cor.test(~runreward_wins + fmri_gamma_t, bdfagg) #uhh, why didn't we see this before?!
cor.test(~runreward_wins + fmri_beta_t, bdfagg)
cor.test(~runreward_wins + fmri_alpha_t, bdfagg)
cor.test(~fmri_gamma_t + fmri_beta_t, bdfagg) #higher decay, lower beta

plot(~runreward_wins + PerformanceTScore, bdfagg)
cor.test(~fmri_gamma_t + PerformanceTScore, bdfagg)

#correlation of entropy and IQ
cor.test(~ PerformanceTScore + elratio, bdfagg)
cor.test(~ PerformanceTScore + elratio_wins, bdfagg)
cor.test(~ PerformanceTScore + elratioF, bdfagg)
cor.test(~ PerformanceTScore + elratioF_wins, bdfagg)

cor.test(~ PerformanceTScore + elratio, bdfagg)
cor.test(~ PerformanceTScore + elratio_wins, bdfagg)
cor.test(~ PerformanceTScore + elratioF, bdfagg)
cor.test(~ PerformanceTScore + elratioF_wins, bdfagg)

bdfagg$PerformanceTScore_wins <- winsor(bdfagg$PerformanceTScore, trim=0.02) #trim a few outliers

cor.test(~ PerformanceTScore_wins + elratio_wins, bdfagg)

ggplot(bdfagg, aes(x=PerformanceTScore_wins, y=elratio_wins)) + geom_point() + stat_smooth()

#conclusion: no relationship between IQ and EL ratio

#but what about mediation?
# elratio -> total rewards is mediated by gamma?  elratio -> gamma -> reward

cor.test(~runreward_wins + elratio_wins, bdfagg) #X -> Y: .55

cor.test(~runreward_wins + fmri_gamma_t, bdfagg) #M -> Y: .37

cor.test(~elratio_wins + fmri_gamma_t, bdfagg) #X -> M: .52

#what about relationships under fixed only?
cor.test(~fixed_alpha + runreward_wins, bdfagg)
cor.test(~fixed_beta + runreward_wins, bdfagg)

cor.test(~fixed_decay_alpha + runreward_wins, bdfagg)
cor.test(~fixed_decay_gamma + runreward_wins, bdfagg)
cor.test(~fixed_decay_beta + runreward_wins, bdfagg)

corstarsl(dplyr::select(bdfagg, runreward_wins, fixed_decay_alpha, fixed_decay_gamma, fixed_decay_beta, fixed_alpha, fixed_beta))
corstarsl(dplyr::select(bdfagg, runreward_wins, elratio_wins, fixed_decay_alpha, fixed_decay_gamma, fixed_decay_beta))
corstarsl(dplyr::select(bdfagg, runreward_wins, elratio_wins, fixed_decay_alpha, fixed_decay_gamma, fixed_decay_beta, fixed_alpha, fixed_beta))

model_withm <- lm(runreward_wins ~ elratio_wins + fmri_gamma_t, data=bdfagg) # + fmri_beta_t
model_withoutm_elratio <- lm(runreward_wins ~ elratio_wins, data=bdfagg) # + fmri_beta_t
model_withoutm_gamma <- lm(runreward_wins ~ fmri_gamma_t, data=bdfagg) # + fmri_beta_t

library(mediation)

#technically, I'm supposed to specify means for control and treatment on elratio_wins -- haven't done this ... not sure how to since it seems arbitrary
med <- mediate(model_withoutm_elratio, model_withm, sims=1000, treat="elratio_wins", mediator="fmri_gamma_t", control.value=.95, treat.value=1.2)
summary(med)

med <- mediate(model_withoutm_gamma, model_withm, sims=1000, mediator="elratio_wins", treat="fmri_gamma_t", treat.value=0.5, control.value=.05)
summary(med)


library(MBESS)
round(with(bdfagg, mediation(x=elratio_wins, mediator=fmri_gamma_t, dv=runreward_wins, bootstrap=TRUE, B=1000)), 3)
round(with(bdfagg, mediation(x=elratio_wins, mediator=fmri_gamma_t, dv=runreward_wins)), 3)

cor.test(~runreward_wins + age, bdfagg)
cor.test(~fmri_gamma_t + age, bdfagg)
cor.test(~fmri_alpha_t + age, bdfagg)
cor.test(~fmri_beta_t + age, bdfagg)
cor.test(~elratio_wins + age, bdfagg)
cor.test(~elratio_wins + age, bdfagg)

corwithtarget(as.data.frame(bdfagg), target="age", omit=c("CompletionDate", "Notes", "lunaid", "scandate"))

round(medmodel <- with(bdfagg, mediation(mediator=elratio_wins, x=fmri_gamma_t, dv=runreward_wins)), 3) #isn't gamma really X since it controls things?
m2 <- with(bdfagg, mediation(mediator=elratio_wins, x=fmri_gamma_t, dv=runreward_wins, bootstrap=TRUE, B=1000, complete.set=TRUE)) #isn't gamma really X since it controls things?
upsilon(bdfsem, x="fmri_gamma_t", m="elratio_wins", y="runreward_wins")
upsilon(bdfsem, m="fmri_gamma_t", x="elratio_wins", y="runreward_wins")

mediation.effect.bar.plot(x=bdfagg$elratio_wins, mediator=bdfagg$fmri_gamma_t, dv=bdfagg$runreward_wins)
mediation.effect.plot(x=bdfagg$elratio_wins, mediator=bdfagg$fmri_gamma_t, dv=bdfagg$runreward_wins)




bdfsem <- bdfagg
bdfsem$runreward_wins <- bdfsem$runreward_wins/1000 #normalize variances 

#lavaan (of course!)
library(lavaan)
m1 <- '
runreward_wins ~ c*elratio_wins + b*fmri_gamma_t
fmri_gamma_t ~ a*elratio_wins

ide := a*b
total := c + a*b
'

fsem <- sem(m1, data=bdfsem, bootstrap=1000, se="bootstrap")
summary(fsem, fit.measures=TRUE, standardize=TRUE, rsquare=TRUE)

#elratio is mediator
m2 <- '
    runreward_wins ~ c*fmri_gamma_t + b*elratio_wins 
    elratio_wins ~ a*fmri_gamma_t
    
    ide := a*b
    total := c + a*b
		propmed := ide/total
    '

#yes, we see that the gamma -> performance relationship is mediated by elratio
fsem2 <- sem(m2, data=bdfsem, bootstrap=1000, se="bootstrap")
summary(fsem2, fit.measures=TRUE, standardize=TRUE, rsquare=TRUE)




#what about the relationship between trial entropy and overall performance?
#a kind of beta series
bdftest <- bdf %>% dplyr::select(subject, run, trial, entropyH) %>% group_by(subject, run) %>%
    mutate(entropyH=entropyH/mean(entropyH)) %>% ungroup() %>%
    mutate(trial=paste0("e", sprintf("%02d", trial))) %>%
    spread(key=trial, value=entropyH, convert=TRUE) %>% left_join(dplyr::select(bdfruns, subject, run, runreward), by=c("subject", "run")) 

#bdftest <- bdf %>% select(subject, run, trial, entropyFixed) %>% group_by(subject, run) %>%
#    mutate(entropyFixed=entropyFixed/mean(entropyFixed)) %>% ungroup() %>%
#    mutate(trial=paste0("e", sprintf("%02d", trial))) %>%
#    spread(key=trial, value=entropyFixed, convert=TRUE) %>% left_join(select(bdfruns, subject, run, runreward), by=c("subject", "run")) 

cm <- cor(bdftest[,paste0("e", sprintf("%02d", 1:50))])
mcm <- reshape2::melt(cm)
ggplot(mcm, aes(x=Var1, y=Var2, fill=value)) + geom_tile()

#get 50 coefs and ses
cmat <- as.data.frame(matrix(NA_real_, nrow=400, 3, dimnames=list(NULL, c("b", "se", "t")))) #be, se, t
for (r in 1:8) {
  thisrun = filter(bdftest, run==r)
  for (trial in 1:50) {
    #f <- as.formula(paste0("runreward ~ e", sprintf("%02d", trial), " + (1|subject)"))
    f <- as.formula(paste0("runreward ~ e", sprintf("%02d", trial)))
    #m <- lmer(f, bdftest)
    if (sd(thisrun[[paste0("e", sprintf("%02d", trial))]]) < .001) {
      cmat[trial + (r-1)*50,] <- c(0,0,0) #no entropy yet
    } else {
      m <- lm(f, thisrun)
      sm <- summary(m)
      cmat[trial + (r-1)*50,] <- sm$coefficients[2,1:3]
    }
  }
}

cmat$trial <- 1:50
cmat$run <- rep(1:8, each=50)

library(viridis)
pdf("relationship between decay entropy and run rewards by trial.pdf", width=15, height=8)
ggplot(cmat, aes(x=trial, y=b, ymin=b-se, ymax = b+se, color=factor(run))) + geom_pointrange(position=position_dodge(width=1)) + 
    geom_line(position=position_dodge(width=1)) + geom_hline(yintercept=0) + stat_smooth(color="black") +
    scale_color_viridis("Run", discrete=TRUE) + theme_bw(base_size=20)
    #scale_color_brewer("Run", palette="OrRd")
dev.off()

cor.test(~ earlyentropy + lateentropy, bdfagg)

#bdfwide <- bdfruns %>% select(subject, run, elratio, elratioF) %>% spread(key="run", value="elratio")#, elratio, elratioF )

summary(mval <- lmer(runreward ~ earlyentropy + lateentropy + run + (1|subject) + (1|run), filter(bdfruns, run > 1)))
car::Anova(mval)

#why is there no subject variation in MLM? (this looks good...)
ggplot(bdfruns, aes(x=factor(subject), y=earlyentropy)) + geom_boxplot()

#make sure the model is treating subject as factor
bdfruns$subject <- factor(bdfruns$subject)

summary(mval <- lmer(runreward ~ elratio + lateentropy + run + (1|subject), filter(bdfruns, run > 1))) # + (1|run)
car::Anova(mval)


library(ez)
ezANOVA(bdfruns, dv=runreward, wid=subject, within)
bdfruns$runreward <- as.numeric(bdfruns$runreward)
xx <- ezMixed(as.data.frame(bdfruns), dv=.(runreward), random=.(subject), fixed=.(earlyentropy))
print(xx$summary)




summary(nlme::lme(runreward ~ earlyentropy, random=~ 1| subject, filter(bdfruns, run>1)))
summary(nlme::lme(runreward ~ lateentropy, random=~ 1| subject, filter(bdfruns, run>1)))
summary(nlme::lme(runreward ~ lateentropy, random=~ 1| subject, filter(bdfruns, run>1)))
summary(nlme::lme(runreward ~ 1, random=~ 1| subject, bdfruns))
summary(lmer(runreward ~ earlyentropy + (1|subject), bdfruns))
summary(lmer(runreward ~ earlyentropy + (1|subject), filter(bdfruns, run>1)))
summary(lmer(runreward ~ 1 + (1|subject), bdfruns))


summary(mval <- lmer(runreward ~ elratio + (1|subject) + (1|run), filter(bdfruns, run > 1)))
car::Anova(mval)

summary(mval <- lmer(runreward ~ elratio + elratioLag + (1|subject) + (1|run), filter(bdfruns, run > 1)))
car::Anova(mval)




#can we see benefits of early entropy in learning
summary(mval <- lmer(runreward ~ earlyentropy + lateentropy + run + (1|subject) + (1|run), filter(bdfruns, run > 1)))
car::Anova(mval)

summary(mval <- lmer(runreward ~ elratio + run + (1|subject) + (1|run), filter(bdfruns, run > 1)))
car::Anova(mval)

summary(mval <- lmer(runreward ~ elratioF + run + (1|subject) + (1|run), filter(bdfruns, run > 1)))
car::Anova(mval)


summary(mval <- lmer(runreward ~ lateentropy + run + (1|subject) + (1|run), filter(bdfruns, run > 1)))
car::Anova(mval)

summary(mval <- lmer(runreward ~ earlyentropy + run + (1|subject) + (1|run), filter(bdfruns, run > 1)))
car::Anova(mval)


summary(mval <- lmer(runreward ~ earlyentropyF + lateentropyF + run + (1|subject) + (1|run), filter(bdfruns, run > 1)))
car::Anova(mval)



summary(mval <- lmer(runreward ~ lateentropy + allentropy + run + (1|subject) + (1|run), filter(bdfruns, run > 1)))
car::Anova(mval)


summary(mval <- lmer(runreward ~ earlyentropy + lateentropy + (1|subject) + (1|run), bdfruns))
car::Anova(mval)
mval <- lmer(runreward ~ (1|subject), filter(bdfruns, run > 1))

cm <- lmerCellMeans(mval, cont.pts=list(entropyH=c(1,2,3,4), trial=c(1, 10, 25, 40, 50)))
ggplot(cm, aes(x=entropyH, y=ev)) + geom_line() + facet_wrap(~trial)



#compute person means and within-subject centered predictors
bdfcent <- bdf %>% select_(.dots=c("LunaID", "run", "trial", "abstschange", "rtlag", "rt", predictors)) %>%
    #mutate_at(wicenter, funs(bwcent=. - mean(., na.rm=TRUE), bwmean=mean(., na.rm=TRUE))) %>% #not really useful
    group_by(LunaID, run) %>%
    mutate_at(wicenter, funs(wicent=. - mean(., na.rm=TRUE), pmean=mean(., na.rm=TRUE))) %>% #within-person centering and person means
    ungroup() %>% mutate_at(c("trial", "abstschangelag", wicenter, paste0(wicenter, "_pmean")), funs(c=. - mean(., na.rm=TRUE))) #between-person centering of person means

#compute further lags for testing
bdfcent <- bdfcent %>% group_by(LunaID, run) %>% mutate(abstschangelag = lag(abstschange, n=1, order_by=trial),
    abstschangelag2 = lag(abstschange, n=2, order_by=trial),
    abstschangelag3 = lag(abstschange, n=3, order_by=trial),
    abstschangelag4 = lag(abstschange, n=4, order_by=trial),
    abstschangelag5 = lag(abstschange, n=5, order_by=trial),
    abstschangelag6 = lag(abstschange, n=6, order_by=trial),
    abstschangelag7 = lag(abstschange, n=7, order_by=trial),
    abstschangelag8 = lag(abstschange, n=8, order_by=trial),
    abstschangelag9 = lag(abstschange, n=9, order_by=trial),
    abstschangelag10 = lag(abstschange, n=10, order_by=trial),
    entropyHlag_wicentlag2 = lag(entropyHlag_wicent, n=2, order_by=trial),
    entropyHlag_wicentlag3 = lag(entropyHlag_wicent, n=3, order_by=trial),
    entropyHlag_wicentlag4 = lag(entropyHlag_wicent, n=4, order_by=trial),
    entropyHlag_wicentlag5 = lag(entropyHlag_wicent, n=5, order_by=trial),
    entropyHlag_wicentlag6 = lag(entropyHlag_wicent, n=6, order_by=trial),
    entropyHlag_wicentlag7 = lag(entropyHlag_wicent, n=7, order_by=trial),
    entropyHlag_wicentlag8 = lag(entropyHlag_wicent, n=8, order_by=trial),
    entropyHlag_wicentlag9 = lag(entropyHlag_wicent, n=9, order_by=trial),
    entropyHlag_wicentlag10 = lag(entropyHlag_wicent, n=10, order_by=trial)
)


mean(bdfcent$entropyHlag_pmean_c)
mean(bdfcent$entropyHlag_wicent, na.rm=TRUE)
mean(bdfcent$trial_c, na.rm=TRUE)
mean(bdfcent$distfromedgelag_c, na.rm=TRUE)

library(ggplot2)
bdf <- bdf %>% filter(trial_abs > 1) %>% group_by(subject, run) %>% #do({browser()})
    #mutate(entropyH_runnorm = entropyH/max(entropyH), entropyFixed_runnorm=entropyFixed/max(entropyFixed)) %>% ungroup()
    mutate(entropyH_runnorm = scale(entropyH), entropyFixed_runnorm=scale(entropyFixed)) %>% ungroup()
pdf("wild1.pdf", width=24, height=24)
#ggplot(filter(bdf, subject < 15), aes(x=trial, y=abstschange, color=entropyH_runnorm)) + facet_grid(subject ~ run) + stat_smooth() + geom_point() +
#    theme_bw() + scale_color_viridis()


ggplot(filter(bdf, subject < 15), aes(x=trial, color=entropyH_runnorm)) + facet_grid(subject ~ run) + 
    #stat_smooth(aes(y=abstschange), color="blue") + 
    #geom_point() +
    #geom_line(aes(y=abstschange), color="blue") +
    geom_line(aes(y=timestep), color="blue") +
    geom_line(aes(y=entropyFixed_runnorm*5), color="black") +
    geom_line(aes(y=entropyH_runnorm*5), color="red") +
    geom_point(aes(y=timestep, shape=omission, color=vchosen)) +
    scale_y_continuous(sec.axis = sec_axis(~./5, name = "Sel entropy")) +
    theme_bw() + scale_color_viridis()


dev.off()

#effect of multiple lags for RT swings
m1 <- lmer(abstschange ~ distfromedgelag_c*omissionlag*vdevlag_c + trial_c + (1 | LunaID) + (1 | run), bdfcent)
summary(m1)

m2 <- lmer(abstschange ~ distfromedgelag_c*omissionlag*vdevlag_c + trial_c + abstschangelag + (1 | LunaID) + (1 | run), bdfcent)
summary(m2)

m3 <- lmer(abstschange ~ distfromedgelag_c*omissionlag*vdevlag_c + trial_c + abstschangelag + abstschangelag2 + (1 | LunaID) + (1 | run), bdfcent)
summary(m3)

m4 <- lmer(abstschange ~ distfromedgelag_c*omissionlag*vdevlag_c + trial_c + abstschangelag + abstschangelag2 + abstschangelag3 +  (1 | LunaID) + (1 | run), bdfcent)
summary(m4)

m5 <- lmer(abstschange ~ distfromedgelag_c*omissionlag*vdevlag_c + trial_c + abstschangelag + abstschangelag2 + abstschangelag3 + abstschangelag4 + (1 | LunaID) + (1 | run), bdfcent)
summary(m5)

m6 <- lmer(abstschange ~ distfromedgelag_c*omissionlag*vdevlag_c + trial_c + abstschangelag + abstschangelag2 + abstschangelag3 + abstschangelag4 + abstschangelag5 + (1 | LunaID) + (1 | run), bdfcent)
summary(m6)

m7 <- lmer(abstschange ~ distfromedgelag_c*omissionlag*vdevlag_c + trial_c + 
        abstschangelag + abstschangelag2 + abstschangelag3 + abstschangelag4 + abstschangelag5 +
        entropyHlag_pmean_c + (1 | LunaID) + (1 | run), subset(bdfcent, run>1))

m8 <- lmer(abstschange ~ distfromedgelag_c*omissionlag*vdevlag_c + trial_c + 
        abstschangelag + abstschangelag2 + abstschangelag3 + abstschangelag4 + abstschangelag5 +
        entropyHlag_pmean_c + entropyHlag_wicent*trial_c*omissionlag + (1 | LunaID) + (1 | run), subset(bdfcent, run>1))

#summary(m7)
car::Anova(m8)

mcrazy <- lmer(abstschange ~ distfromedgelag_c*omissionlag*vdevlag_c + trial_c + 
        abstschangelag + abstschangelag2 + abstschangelag3 + abstschangelag4 + abstschangelag5 +
        abstschangelag6 + abstschangelag7 + abstschangelag8 + abstschangelag9 + abstschangelag10 +
        entropyHlag_pmean_c + entropyHlag_wicent*trial_c*omissionlag + (1 | LunaID) + (1 | run), subset(bdfcent, run>1))

summary(mcrazy)

car::Anova(mcrazy)

minsane <- lmer(abstschange ~ distfromedgelag_c*omissionlag*vdevlag_c + trial_c + 
        abstschangelag + abstschangelag2 + abstschangelag3 + abstschangelag4 + abstschangelag5 +
        abstschangelag6 + abstschangelag7 + abstschangelag8 + abstschangelag9 + abstschangelag10 +
        entropyHlag_wicentlag2 + entropyHlag_wicentlag3 + entropyHlag_wicentlag4 + entropyHlag_wicentlag5 +
        entropyHlag_wicentlag6 + entropyHlag_wicentlag7 + entropyHlag_wicentlag8 + entropyHlag_wicentlag9 + entropyHlag_wicentlag10 +
        entropyHlag_pmean_c + entropyHlag_wicent*trial_c*omissionlag + (1 | LunaID) + (1 | run), subset(bdfcent, run>1))

summary(minsane)

car::Anova(minsane)


#acf(na.omit(bdfcent$abstschange))
pdf("Abstract art.pdf", width=12, height=10)
ggplot(bdfcent, aes(x=trial, y=abstschange)) + stat_smooth() + geom_jitter(alpha=0.2) + theme_gray(base_size=20)
dev.off()

ggplot(subset(bdf,run>1), aes(x=trial, y=abstschange)) + stat_smooth(method="loess") + theme_gray(base_size=20) + facet_wrap(~msplit) #geom_jitter(alpha=0.2) + 



toresid <- lmer(abstschange ~ distfromedgelag_c*omissionlag*vdevlag_c + trial_c + (1 | LunaID) + (1 | run), bdfcent, na.action=na.exclude)
bdf$leftovers <- resid(toresid)

ggplot(subset(bdf,run>1), aes(x=trial, y=leftovers)) + stat_smooth(method="loess") + theme_gray(base_size=20) + facet_wrap(~msplit) #geom_jitter(alpha=0.2) +






m1 <- lmer(abstschange ~ entropyHlag_pmean_c + entropyHlag_wicent + trial_c + (1 | LunaID), bdfcent)
summary(m1)
car::Anova(m1)

m2 <- lmer(abstschange ~ entropyHlag_pmean_c*distfromedgelag_c + entropyHlag_wicent*distfromedgelag_c + trial_c + (1 | LunaID), bdfcent)
summary(m2)
car::Anova(m2)

m5 <- lmer(abstschange ~ entropyHlag_pmean_c*distfromedgelag_c*omissionlag*vdevlag_c + entropyHlag_wicent*distfromedgelag_c*omissionlag*vdevlag_c + trial_c*entropyHlag_wicent + (1 | LunaID), bdfcent)
summary(m5)
car::Anova(m5)

#current reasonable winner
bdfcent$abstschange_sec <- bdfcent$abstschange*100 
m6 <- lmer(abstschange_sec ~ entropyHlag_pmean_c*distfromedgelag_c*omissionlag*vdevlag_c + entropyHlag_wicent*distfromedgelag_c*omissionlag*vdevlag_c + trial_c*entropyHlag_wicent + entropyHlag_wicent*entropyHlag_pmean_c + abstschangelag_c + (1 | LunaID) + (1 | run), bdfcent)
summary(m6)
car::Anova(m6)

#
bdfcent$abstschange_sec <- bdfcent$abstschange*100 
m6 <- lmer(abstschange_sec ~ entropyHlag_pmean_c*distfromedgelag_c*omissionlag*vdevlag_c + entropyHlag_wicent*distfromedgelag_c*omissionlag*vdevlag_c + trial_c*entropyHlag_wicent + entropyHlag_wicent*entropyHlag_pmean_c + abstschangelag_c + (1 | LunaID) + (1 | run), bdfcent)
summary(m6)
car::Anova(m6)


msimple <- lmer(abstschange_sec ~ entropyHlag_pmean_c*entropyHlag_wicent + trial_c*entropyHlag_wicent + abstschangelag_c + (1 | LunaID) + (1 | run), bdfcent)
summary(msimple)
car::Anova(msimple)

#spot check
bdfcent %>% group_by(LunaID, run) %>% summarize(mean(entropyHlag_wicent, na.rm=TRUE)) %>% print(n=100)

#Here's my current thinking for PLoS paper
#1) Report simple abstschange ~ entropy + trial
#2) Report effect after throwing in a bunch of other stuff
#3) Report effect dissociating within versus between entropy

msimple <- lmer(abstschange_sec ~ entropyHlag_c + (1 | LunaID) + (1 | run), filter(bdfcent, run>1))
summary(msimple)


msimple <- lmer(abstschange_sec ~ entropyHlag_c*trial_c + abstschangelag_c + (1 | LunaID) + (1 | run), filter(bdfcent, run>1))
summary(msimple)
car::Anova(msimple)


pdf("simple entropy for plos.pdf", width=12, height=7)
cm <- lmerCellMeans(msimple, n.cont=10, fixat0=c("trial_c", "abstschangelag_c"))
cm$entropyHlag_c <- cm$entropyHlag_c + mean(bdfcent$entropyHlag, na.rm=TRUE) #uncenter for plotting
ggplot(cm, aes(x=entropyHlag_c, y=abstschange_sec, ymin=abstschange_sec-se, ymax=abstschange_sec+se)) + 
    geom_line(size=2.5) + theme_bw(base_size=24) + geom_pointrange()
dev.off()

msimple2 <- lmer(abstschange_sec ~ entropyHlag_c*trial_c*omissionlag + abstschangelag_c + (1 | LunaID) + (1 | run), filter(bdfcent, run>1))
summary(msimple2)
car::Anova(msimple2)


pdf("simple entropy for plos with omission.pdf", width=12, height=7)
cm <- lmerCellMeans(msimple2, n.cont=10, fixat0=c("trial_c", "abstschangelag_c"))
cm$entropyHlag_c <- cm$entropyHlag_c + mean(bdfcent$entropyHlag, na.rm=TRUE) #uncenter for plotting
ggplot(cm, aes(x=entropyHlag_c, y=abstschange_sec, ymin=abstschange_sec-se, ymax=abstschange_sec+se, color=omissionlag)) + 
    geom_line(size=2.5) + theme_bw(base_size=24) + geom_pointrange()
dev.off()


#controlling for various effects and confounders
mcontrol <- lmer(abstschange_sec ~ entropyHlag_c*trial_c + abstschangelag_c +
        evdevlag_c*omissionlag*vdevlag_c*distfromedgelag_c + (1 | LunaID) + (1 | run), filter(bdfcent, run>1))


summary(mcontrol)
car::Anova(mcontrol)


pdf("wi entropy for plos.pdf", width=5, height=4)
cm <- lmerCellMeans(mcontrol, n.cont=10, fixat0=c("evdevlag_c", "trial_c", "distfromedgelag_c", "vdevlag_c", "abstschangelag_c"))
#average over omissions and rewards
cm$predpoint <- 1:10 #rep 2x
cm <- cm %>% group_by(predpoint) %>% summarize_if(is.numeric, mean)
cm$entropyHlag_c <- cm$entropyHlag_c + mean(bdfcent$entropyHlag, na.rm=TRUE) #uncenter for plotting
ggplot(cm, aes(x=entropyHlag_c, y=abstschange_sec, ymin=abstschange_sec-se, ymax=abstschange_sec+se)) + 
     geom_line(size=1.5) + theme_bw(base_size=20) + geom_pointrange(size=0.8) + ylab("RT swing (ms)") + xlab("Entropy of value distribution") +
     theme(axis.title.y=element_text(margin=margin(r=15)), axis.title.x=element_text(margin=margin(t=10)))
dev.off()

#divide into between and within
mdivide <- lmer(abstschange_sec ~ entropyHlag_pmean_c*trial_c + entropyHlag_wicent*trial_c + abstschangelag_c +
        evdevlag_c*omissionlag*vdevlag_c*distfromedgelag_c + (1 | LunaID) + (1 | run), filter(bdfcent, run>1))

summary(mdivide)
car::Anova(mdivide)

pdf("wi entropy as a function of person avg entropy.pdf", width=12, height=7)
cm <- lmerCellMeans(m6, n.cont=10, divide="entropyHlag_pmean_c", fixat0=c("trial_c", "distfromedgelag_c", "vdevlag_c", "abstschangelag_c"))
ggplot(cm, aes(x=entropyHlag_wicent, y=abstschange_sec, ymin=abstschange_sec-se, ymax=abstschange_sec+se, color=omissionlag)) + 
    facet_wrap(~entropyHlag_pmean_c) + geom_line(size=2.5) + theme_bw(base_size=24) + geom_pointrange()
dev.off()

pdf("wi entropy as a function of trial.pdf", width=12, height=7)
#cm <- lmerCellMeans(m6, n.cont=10, fixat0=c("entropyHlag_pmean_c", "distfromedgelag_c", "vdevlag_c", "abstschangelag_c"), cont.pts=list(trial_c=c(-24, 0, 24)))
cm <- lmerCellMeans(m6, n.cont=10, divide="trial_c", n.divide=3, fixat0=c("entropyHlag_pmean_c", "distfromedgelag_c", "vdevlag_c", "abstschangelag_c"))
levels(cm$trial_c) <- c("Trial 10", "Trial 25", "Trial 40")
#cm$abstschange <- cm$abstschange*100
ggplot(cm, aes(x=entropyHlag_wicent, y=abstschange_sec, ymin=abstschange_sec-se, ymax=abstschange_sec+se, color=omissionlag)) + 
    facet_wrap(~trial_c) + geom_line(size=1.5) + theme_bw(base_size=20) + geom_pointrange() + xlab("Entropy (centered, within run)") + ylab("Change in RT (sec)") +
    scale_color_brewer("Prior outcome", palette="Dark2")
dev.off()

#what about predicting quality of choice as a function of entropy
m6 <- lmer(ev ~ entropyHlag_pmean_c*distfromedgelag_c*omissionlag + entropyHlag_wicent*distfromedgelag_c*omissionlag + trial_c*entropyHlag_wicent + entropyHlag_wicent*entropyHlag_pmean_c + (1 | LunaID) + (1 | run), bdfcent)
summary(m6)
car::Anova(m6)




m8 <- lmer(abstschange ~ entropyHlag*omissionlag + evdevlag*omissionlag*vdevlag + vdevlag*omissionlag + trial_c + abstschangelag_c + (1 | LunaID), bdf)
summary(m8)
car::Anova(m8)

 


anova(m5, m6)


m1 <- lmer(abstschange ~ entropyHlag + trial + (1 | LunaID), bdf)
summary(m1)

m2 <- lmer(abstschange ~ entropyHlag + trial + distfromedgelag + (1 | LunaID), bdf)
summary(m2)

m3 <- lmer(abstschange ~ entropyHlag*distfromedgelag + trial + (1 | LunaID), bdf)
summary(m3)
car::Anova(m3)


m4 <- lmer(abstschange ~ entropyHlag*distfromedgelag*omissionlag + trial + (1 | LunaID), bdf)
summary(m4)
car::Anova(m4)


m5 <- lmer(abstschange ~ entropyHlag*distfromedgelag*omissionlag*vdevlag + trial + (1 | LunaID), bdf)
summary(m5)
car::Anova(m5)


cm <- lmerCellMeans(m4, n.cont=10, divide="distfromedgelag", fixat0="trial")
pdf("me4_entropy_omission_dist.pdf", width=12, height=6)
ggplot(cm, aes(x=entropyHlag, y=abstschange, ymin=abstschange-se, ymax=abstschange+se, color=omissionlag)) + geom_line(size=2.5) + theme_bw(base_size=24) + #geom_linerange(size=1.0, width=250) + 
    ylab("Abs trialwise RT change (cs)") + xlab("Entropy") + scale_color_brewer("Prior Outcome", palette="Set2") + geom_vline(xintercept=0) + geom_hline(yintercept=0) + facet_wrap(~distfromedgelag)
dev.off()



m4 <- lmer(abstschange ~ entropyHlag*omissionlag + trial + (1 | LunaID), bdf)
summary(m4)
car::Anova(m4)

cm <- lmerCellMeans(m4, n.cont=10, fixat0="trial")
pdf("me4_entropy_omission.pdf", width=8, height=6)
ggplot(cm, aes(x=entropyHlag, y=abstschange, ymin=abstschange-se, ymax=abstschange+se, color=omissionlag)) + geom_line(size=2.5) + theme_bw(base_size=24) + #geom_linerange(size=1.0, width=250) + 
    ylab("Abs trialwise RT change (cs)") + xlab("Entropy") + scale_color_brewer("Prior Outcome", palette="Set2") + geom_vline(xintercept=0) + geom_hline(yintercept=0)
dev.off()


anova(m1, m2, m3)

cm <- lmerCellMeans(m1, n.cont=10, fixat0="trial")

pdf("me1_entropyH.pdf", width=8, height=6)
ggplot(cm, aes(x=entropyHlag, y=abstschange, ymin=abstschange-se, ymax=abstschange+se)) + geom_line(size=2.5) + theme_bw(base_size=24) + #geom_linerange(size=1.0, width=250) + 
    ylab("Abs trialwise RT change (cs)") + xlab("Entropy") + scale_color_brewer("Prior Outcome", palette="Set2") + geom_vline(xintercept=0) + geom_hline(yintercept=0)
dev.off()

m4 <- lmer(abstschange ~ entropyHlag + trial + omissionlag + distfromedgelag + (1 | LunaID), bdf)
summary(m4)

cm <- lmerCellMeans(m4, n.cont=10, fixat0="trial")

pdf("me4_entropyH_omission.pdf", width=8, height=6)
ggplot(cm, aes(x=entropyHlag, y=abstschange, ymin=abstschange-se, ymax=abstschange+se, color=omissionlag)) + geom_line(size=2.5) + theme_bw(base_size=24) + #geom_linerange(size=1.0, width=250) + 
    ylab("Abs trialwise RT change (cs)") + xlab("Entropy") + scale_color_brewer("Prior Outcome", palette="Set2") + geom_vline(xintercept=0) + geom_hline(yintercept=0)
dev.off()

m5 <- lmer(abstschange ~ entropyHlag*omissionlag + trial + distfromedgelag + (1 | LunaID), bdf)
summary(m5)

anova(m4, m5)

cm <- lmerCellMeans(m3, n.cont=10, fixat0="trial")

pdf("me3_entropyH_omission.pdf", width=8, height=6)
ggplot(cm, aes(x=entropyHlag, y=abstschange, ymin=abstschange-se, ymax=abstschange+se, color=omissionlag)) + geom_line(size=2.5) + theme_bw(base_size=24) + #geom_linerange(size=1.0, width=250) + 
    ylab("Abs trialwise RT change (cs)") + xlab("Entropy") + scale_color_brewer("Prior Outcome", palette="Set2") + geom_vline(xintercept=0) + geom_hline(yintercept=0)
dev.off()

#add deviation of chosen versus max value on prior trial
m4 <- lmer(abstschange ~ entropyHlag*omissionlag + evdevlag*omissionlag + trial + (1 | LunaID), bdf)
summary(m4)

#anova(m3, m4)

cm <- lmerCellMeans(m4, n.cont=10, divide="evdevlag", fixat0="trial")

pdf("me4_entropyH_omission_evdevlag.pdf", width=14, height=6)
ggplot(cm, aes(x=entropyHlag, y=abstschange, ymin=abstschange-se, ymax=abstschange+se, color=omissionlag)) + geom_line(size=2.5) + theme_bw(base_size=24) + #geom_linerange(size=1.0, width=250) + 
    ylab("Abs trialwise RT change (cs)") + xlab("Entropy") + scale_color_brewer("Prior Outcome", palette="Set2") + geom_vline(xintercept=0) + 
    geom_hline(yintercept=0) + facet_wrap(~evdevlag)
dev.off()


m5 <- lmer(abstschange ~ entropyHlag*omissionlag + evdevlag*omissionlag*vdevlag + vdevlag*omissionlag + trial + (1 | LunaID), bdf)
summary(m5)
car::Anova(m5)

anova(m4, m5)

cm <- lmerCellMeans(m5, n.cont=10, divide=c("evdevlag", "vdevlag"), fixat0="trial")

pdf("me5_entropyH_omission_evdevlag_vdevlag.pdf", width=14, height=14)
ggplot(cm, aes(x=entropyHlag, y=abstschange, ymin=abstschange-se, ymax=abstschange+se, color=omissionlag)) + geom_line(size=2.5) + theme_bw(base_size=24) + #geom_linerange(size=1.0, width=250) + 
    ylab("Abs trialwise RT change (cs)") + xlab("Entropy") + scale_color_brewer("Prior Outcome", palette="Set2") + geom_vline(xintercept=0) + 
    geom_hline(yintercept=0) + facet_grid(vdevlag~evdevlag)
dev.off()

m6 <- lmer(abstschange ~ entropyHlag*omissionlag + evdevlag*omissionlag*vdevlag + vdevlag*omissionlag + trial + (1 | LunaID), subset(bdf, trial > 25))
summary(m6)

m7 <- lmer(abstschange ~ entropyHlag*omissionlag + evdevlag*omissionlag*vdevlag + vdevlag*omissionlag + trial + abstschangelag + (1 | LunaID), subset(bdf, trial > 25))
summary(m7)

m8 <- lmer(abstschange ~ entropyHlag*omissionlag + evdevlag*omissionlag*vdevlag + vdevlag*omissionlag + trial + abstschangelag + (1 | LunaID), bdf)
summary(m8)
car::Anova(m8)

m9 <- lmer(abstschange ~ entropyHlag*omissionlag + evdevlag*omissionlag*vdevlag + vdevlag*omissionlag + trial + abstschangelag + entropyHlag*distfromedgelag + (1 | LunaID), bdf)
summary(m9)
car::Anova(m9)


aggdata = bdf %>% group_by(subject) %>% summarize(totreward=sum(score))
ggplot(aggdata, aes(x=totreward)) + geom_histogram(bins=15)


pdf("rough entropy by trial and performance medsplit.pdf", width=10, height=6)
ggplot(subset(bdf, run>1), aes(x=trial, y=entropyH, color=msplit)) + stat_smooth(alpha=0.1) + theme_bw(base_size=22) + xlab("Trial") + ylab("Entropy of value distribution") +
    scale_color_brewer("Performance", palette="Set1")
dev.off()

pdf("rough entropy by trial and performance medsplit.pdf", width=10, height=6)
ggplot(subset(bdf, run>1), aes(x=trial, y=entropyFixed, color=msplit)) + stat_smooth(alpha=0.1) + theme_bw(base_size=22) + xlab("Trial") + ylab("Entropy of value distribution") +
    scale_color_brewer("Performance", palette="Set1")
dev.off()

entropyagg <- bdf %>% select(LunaID, trial, run, entropyFixed, entropyH, msplit) %>% gather(key="Model", value="entropy", entropyFixed, entropyH) %>%
    mutate(Model=recode(Model, entropyFixed="Fixed LR V", entropyH="Fixed LR V Sel. Maint.")) %>% filter(run>1) %>% group_by(Model, msplit, trial) %>%
    do(data.frame(rbind(Hmisc::smean.cl.boot(.$entropy)))) %>% ungroup()

pdf("rough entropy by trial and performance medsplit.pdf", width=10, height=6)
ggplot(entropyagg, aes(x=trial, y=Mean, ymin=Lower, ymax=Upper, linetype=msplit, color=Model, fill=Model)) + geom_line(size=1.2) + geom_ribbon(alpha=0.1) +
    theme_bw(base_size=22) +
    xlab("Trial") + ylab("Entropy of value distribution") +
    scale_color_brewer("Performance", palette="Dark2") +
    scale_fill_brewer("Performance", palette="Dark2") + 
    scale_linetype("Overall Performance") #+ facet_wrap(~Model)
dev.off()

pdf("rough entropy by trial and performance medsplit decay only.pdf", width=10, height=6)
ggplot(filter(entropyagg, Model=="Fixed LR V Sel. Maint."), aes(x=trial, y=Mean, ymin=Lower, ymax=Upper, linetype=msplit)) + 
    geom_line(size=1.2) + geom_ribbon(alpha=0.1) +
    theme_bw(base_size=22) +
    xlab("Trial") + ylab("Entropy of value distribution") +
    scale_linetype("Overall Performance") #+ facet_wrap(~Model)
dev.off()


entropysplit <- bdf %>% select(LunaID, trial, run, entropyFixed, entropyH, msplit) %>% gather(key="Model", value="entropy", entropyFixed, entropyH) %>%
    mutate(Model=recode(Model, entropyFixed="Fixed LR V", entropyH="Fixed LR V Sel. Maint."))


pdf("rough entropy by trial and performance medsplit.pdf", width=10, height=6)
ggplot(subset(entropysplit, run>1), aes(x=trial, y=entropy, linetype=msplit, color=Model)) + stat_smooth(alpha=0.1, method="loess") + theme_bw(base_size=22) + xlab("Trial") + ylab("Entropy of value distribution") +
    scale_color_brewer("Performance", palette="Dark2") + scale_linetype("Overall Performance") #+ facet_wrap(~Model)
dev.off()

#edesc_aggruns <- edescriptives %>% filter(run > 1) %>% group_by(LunaID, Model, trial) %>% summarize(value=mean(value))


#average magnitude of swings as a function of total reward
m9 <- lmer(abstschange ~ cumreward*trial + abstschangelag + entropyHlag*distfromedgelag + (1 | LunaID), bdf)
summary(m9)
car::Anova(m9)

m9 <- lmer(abstschange ~ totreward*trial + abstschangelag + entropyHlag*distfromedgelag + (1 | LunaID), bdf)
summary(m9)



###
#two ideas: emotion and contingency modulate RT swings and entropy?

bdf %>% filter(trial > 5) %>% group_by(rewFunc, emotion) %>% summarize(mentropyD = mean(entropyH, na.rm=TRUE), mrtswing=mean(abstschange, na.rm=TRUE), 
    sdentropyD = sd(entropyH, na.rm=TRUE), sdrtswing=sd(abstschange, na.rm=TRUE))

summary(lmer(abstschange ~ emotion*rewFunc + (1 | LunaID), bdf))

summary(mblock <- lmer(entropyH ~ emotion + rewFunc + (1 | LunaID), subset(bdf, trial_abs > 3)))
summary(mblock <- lmer(entropyH ~ emotion + rewFunc + (1 | LunaID), bdf))
cm <- lmerCellMeans(mblock)

ggplot(cm, aes(x=rewFunc, y=entropyH, color=emotion, ymin=entropyH-se, ymax=entropyH+se)) + geom_pointrange(size=1.5, position=position_dodge(width=0.2))

summary(mblock <- lmer(entropyH ~ rewFunc + (1 | LunaID), subset(bdf, trial_abs> 3)))


cm <- lmerCellMeans(mblock)

pdf("entropy_by_contingency.pdf", width=10, height=8)
ggplot(cm, aes(x=rewFunc, y=entropyH, ymin=entropyH-se, ymax=entropyH+se)) + geom_pointrange(size=1.5) + theme_bw(base_size=15)
dev.off()

summary(mblock <- lmer(abstschange ~ rewFunc + (1 | LunaID), subset(bdf, trial_abs> 3)))

cm <- lmerCellMeans(mblock)
pdf("rtswing_by_contingency.pdf", width=10, height=8)
ggplot(cm, aes(x=rewFunc, y=abstschange, ymin=abstschange-se, ymax=abstschange+se)) + geom_pointrange(size=1.5) + theme_bw(base_size=15)
dev.off()

summary(mblock <- lmer(abstschange ~ emotion + (1 | LunaID), subset(bdf, trial_abs> 3)))

cm <- lmerCellMeans(mblock)
pdf("rtswing_by_emotion.pdf", width=10, height=8)
ggplot(cm, aes(x=emotion, y=abstschange, ymin=abstschange-se, ymax=abstschange+se)) + geom_pointrange(size=1.5) + theme_bw(base_size=15)
dev.off()


#






#fixed LR model (no decay)
m7 <- lmer(abstschange ~ entropyFlag*omissionlag + evdevlag*omissionlag*vdevlag + vdevlag*omissionlag + trial + (1 | LunaID), subset(bdf, trial > 25))
summary(m7)

m7 <- lmer(abstschange ~ entropyFlag*omissionlag*trial + evdevlag*omissionlag*vdevlag + vdevlag*omissionlag + (1 | LunaID), bdf)
summary(m7)

car::Anova(m7)


m7 <- lmer(abstschange ~ entropyHlag*omissionlag*trial + evdevlag*omissionlag*vdevlag + vdevlag*omissionlag + (1 | LunaID), bdf)
summary(m7)

bdf <- bdf %>% group_by(subject) %>% mutate(medswing = median(abstschange, na.rm=TRUE)) %>% ungroup()
bdf %>% select(subject, run, trial, abstschange, medswing) %>% arrange(subject, run, trial) %>% tail(n=100)

bdf %>% group_by(subject, run) %>% summarize(cent = cor(entropyFixed, entropyH)) %>% summarize()  








pdf("swing_spaghetti.pdf", width=15, height=20)
ggplot(bdf, aes(x=trial, y=abstschange, group=LunaID)) + geom_line(alpha=0.2) + stat_smooth(aes(group=NULL))  + facet_wrap(~run, ncol=1) #subset(bdf, subject < 3)
dev.off()


#model 1: RT change predicted by prior omission and prior RT deviation from V max 
m1 <- lmer(timestepchange ~ omissionlag*vdevlag + trial + (1|LunaID), bdf)
summary(m1)

cm <- lmerCellMeans(m1, n.cont=10, fixat0="trial")
pdf("m1_vdevtest.pdf", width=8, height=6)
ggplot(cm, aes(x=vdevlag, y=timestepchange, color=omissionlag, ymin=timestepchange-se, ymax=timestepchange+se)) + geom_line(size=2.5) + theme_bw(base_size=24) + #geom_linerange(size=1.0, width=250) + 
    ylab("Mean trialwise RT change (ds)") + xlab("RT deviation from maximum value\non prior trial (ms)") + scale_color_brewer("Prior Outcome", palette="Set2") + geom_vline(xintercept=0) + geom_hline(yintercept=0)
dev.off()

#model 2: Does adding magnitude of chosen option versus max help? yes 
m2 <- lmer(timestepchange ~ omissionlag*vdevlag*evdevlag + trial + (1|LunaID), bdf)
summary(m2)

anova(m1, m2)
cm2 <- lmerCellMeans(m2, n.cont=10, divide="evdevlag", fixat0="trial")
pdf("m2_vdevtest_evdiff.pdf", width=15, height=6)
ggplot(cm2, aes(x=vdevlag, y=timestepchange, color=omissionlag, ymin=timestepchange-se, ymax=timestepchange+se)) + geom_line(size=2.5) + theme_bw(base_size=24) + #geom_linerange(size=1.0, width=250) + 
    ylab("Mean trialwise RT change (cs)") + xlab("RT deviation from maximum value\non prior trial (ms)") + scale_color_brewer("Prior Outcome", palette="Set2") + geom_vline(xintercept=0) + geom_hline(yintercept=0) +
    facet_wrap(~evdevlag)
dev.off()

#model 3: adding entropy
#try including entropy
m3 <- lmer(timestepchange ~ omissionlag*vdevlag*entropylag*evdevlag + trial + (1|LunaID), bdf)
summary(m3)
car::Anova(m3)
anova(m2, m3)

cm3 <- lmerCellMeans(m3, divide="entropylag", n.cont=10, fixat0=c("trial", "evdevlag"))

pdf("m3_vdevtest_evdiff_entropy.pdf", width=15, height=6)
ggplot(cm3, aes(x=vdevlag, y=timestepchange, color=omissionlag, ymin=timestepchange-se, ymax=timestepchange+se)) + geom_line(size=2.5) + theme_bw(base_size=24) + #geom_linerange(size=1.0, width=250) + 
    ylab("Mean trialwise RT change (cs)") + xlab("RT deviation from maximum value\non prior trial (cs)") + scale_color_brewer("Prior Outcome", palette="Set2") + geom_vline(xintercept=0) + geom_hline(yintercept=0) +
    facet_wrap(~entropylag, nrow=1)
dev.off()

##current entropy?
test2 <- lmer(timestepchange ~ omissionlag*vdevlag*entropy + (1|LunaID), bdf)
summary(test2)
car::Anova(test2)
anova(test, test2)

cm2 <- lmerCellMeans(test2, divide="entropy", n.cont=10)

pdf("entropytest_cur.pdf", width=15, height=6)
ggplot(cm2, aes(x=vdevlag, y=timestepchange, color=omissionlag, ymin=timestepchange-se, ymax=timestepchange+se)) + geom_line(size=2.5) + theme_bw(base_size=24) + #geom_linerange(size=1.0, width=250) + 
    ylab("Mean trialwise RT change (cs)") + xlab("RT deviation from maximum value\non prior trial (cs)") + scale_color_brewer("Prior Outcome", palette="Set2") + geom_vline(xintercept=0) + geom_hline(yintercept=0) +
    facet_wrap(~entropy, nrow=1)
dev.off()

#incorporate U into models. Does the point of maximal uncertainty predict swings in that direction?
m1 <- lmer(timestepchange ~ omissionlag*vdevlag + (1|LunaID), bdf)
summary(m1)

m2 <- lmer(timestepchange ~ omissionlag*udevlag + (1|LunaID), bdf)
summary(m2)

m2 <- lmer(timestepchange ~ udevlag + (1|LunaID), bdf)
summary(m2)


m3 <- lmer(timestepchange ~ omissionlag*entropylag + (1|LunaID), bdf)
summary(m3)

m4 <- lmer(timestepchange ~ omissionlag*vdevlag + omissionlag*udevlag + (1|LunaID), bdf)
summary(m4)

cm <- lmerCellMeans(m4, divide="vdevlag", n.cont=10)
pdf("udevlag_vdevlag_predicted.pdf", width=10, height=8)
ggplot(cm, aes(x=udevlag, y=timestepchange, ymin=timestepchange-se, ymax=timestepchange+se, color=omissionlag)) + geom_line(size=2) + facet_wrap(~vdevlag) + 
    geom_pointrange(size=1) + 
    theme_bw(base_size=20)
dev.off()
m5 <- lmer(timestepchange ~ omissionlag*vdevlag + omissionlag*udevlag + omissionlag*entropylag + (1|LunaID), bdf)
summary(m5)

#no troubling overall correlations
cor(select(bdf, vdevlag, udevlag, entropylag), use="pairwise.complete.obs")




#try to predict current timestep based on a) prior outcome, b) 
test2 <- lmer(timestep ~ timesteplag + entropy + (1|LunaID), bdf)

library(lme4)
summary(m1 <- lmer(timestep ~ timesteplag + rtvmax + rtumax + (1|LunaID), bdf))
car::Anova(m1)

summary(m2 <- lmer(timestep ~ timesteplag + rtvmax * rtumax + (1|LunaID), bdf))
car::Anova(m2)

summary(m3 <- lmer(timestep ~ timesteplag + rtvmax + rtumax + (1|LunaID) + (1|run), bdf))
car::Anova(m3)

anova(m1, m2)

##Model of reaction times being predicted prior RTs, value, uncertainty
summary(m1 <- lmer(rt ~ rtlag + rtvmaxlag + rtumaxlag + (1|ID) + (1|run), bdf, REML=FALSE))
summary(m2 <- lmer(rt ~ rtlag + rtvmaxlag + rtumaxlag + (1|ID) + (1|ID:run), bdf, REML=FALSE))
summary(m3 <- lmer(rt ~ rtlag + rtvmaxlag + rtumaxlag + (1 |ID) + (1 + rtumaxlag|run) + (1|ID:run), bdf, REML=FALSE))
anova(m1, m2, m3)

summary(m1 <- lmer(rt ~ rtlag + rtvmaxlag + rtumaxlag + (1|LunaID) + (1|run), bdf, REML=FALSE))

car::Anova(m3)
    
cor(coef(m3))

#replicate Alex's analysis of emotion effects -- close, but no cigar... also, was pemax lagged in Alex's analysis?
bdf$emotion <- relevel(bdf$emotion, ref="scram")
bdf$abspe <- abs(bdf$pemax)
summary(memo <- lmer(rt ~ rtlag + rtlag2 + rtlag3 + emotion + rtvmaxlag + abspelag + rewFunc + 
            emotion*omissionlag + emotion*abspelag + omissionlag*abspelag + (1 | ID) + (1 | run) + (1 | run:ID), filter(bdf, rewFunc %in% c("IEV", "DEV") & !ID==11282), REML=FALSE))

car::Anova(memo)

#allow for 4-way interaction
summary(memo <- lmer(rt ~ rtlag + rtlag2 + rtlag3 + emotion*abspelag*rewFunc*omissionlag + rtvmaxlag + (1 | ID) + (1 | run), filter(bdf, rewFunc %in% c("IEV", "DEV") & !ID==11282), REML=FALSE))
car::Anova(memo)

#yes, entropy is in the mix (getting pretty complex!!)
summary(memo <- lmer(rt ~ rtlag + rtlag2 + rtlag3 + emotion*abspelag*rewFunc*omissionlag*entropyHlag + rtvmaxlag + (1 | ID) + (1 | run), filter(bdf, rewFunc %in% c("IEV", "DEV") & !ID==11282), REML=FALSE))
car::Anova(memo)

#dial it back for a second
summary(memo <- lmer(rt ~ rtlag + rtlag2 + rtlag3 + rtvmaxlag + pemaxlag + (1 | ID) + (1 | run), filter(bdf, rewFunc %in% c("IEV", "DEV") & !ID==11282), REML=FALSE))
car::Anova(memo)

summary(memo <- lmer(rt ~ rtlag + rtlag2 + rtlag3 + rtvmaxlag + ppelag + npelag + (1 | ID) + (1 | run), filter(bdf, rewFunc %in% c("IEV", "DEV") & !ID==11282), REML=FALSE))
car::Anova(memo)


summary(memo <- lmer(rt ~ rtlag + rtlag2 + rtlag3 + rtvmaxlag + ppelag*emotion + npelag*emotion + ppelag*rewFunc + npelag*rewFunc + (1 | ID) + (1 | run), filter(bdf, rewFunc %in% c("IEV", "DEV") & !ID==11282), REML=FALSE))
car::Anova(memo)

summary(memo <- lmer(rt ~ rtlag + rtlag2 + rtlag3 + rtvmaxlag + ppelag*emotion*rewFunc + npelag*emotion*rewFunc + (1 | ID) + (1 | run), filter(bdf, rewFunc %in% c("IEV", "DEV") & !ID==11282), REML=FALSE))
car::Anova(memo)

#looks like there's a an entropy x emo effect such that greater entropy in scram and happy associated with longer RTs, but this is reversed in fear
summary(memo <- lmer(rt ~ rtlag + rtlag2 + rtlag3 + rtvmaxlag + abspelag*emotion*rewFunc*entropyHlag + (1 | ID) + (1 | run), filter(bdf, rewFunc %in% c("IEV", "DEV") & !ID==11282), REML=FALSE))
car::Anova(memo)

#equivalent model for RT swing
summary(memo <- lmer(abstschange ~ abstschangelag + abstschangelag2 + abstschangelag3 + rtvmaxlag + abspelag*emotion*rewFunc*entropyHlag + (1 | ID) + (1 | run), filter(bdf, rewFunc %in% c("IEV", "DEV") & !ID==11282), REML=FALSE))
car::Anova(memo)


#######
#w/i versus b/w run effects
summary(m1 <- lmer(rt ~ rtlag + rtvmaxlag_wicent*trial + rtvmaxlag_pmean_c + rtumaxlag_pmean_c + rtumaxlag_wicent*trial + (1|LunaID) + (1|run), bdfcent, REML=FALSE))
car::Anova(m1)
    
summary(m2 <- lmer(timestepchange ~ timesteplag + rtvmax + rtumax + (1|LunaID), bdf))
car::Anova(m2)

#okay, let's bring in 

#cor.test(~rtumax + rtumaxlag, bdf)

#wi-person scale tschange
bdf <- bdf %>% group_by(LunaID, run) %>% mutate(tschange_z = as.vector(scale(timestepchange))) %>% ungroup()

#should we within-person z-score as in frank?
summary(m2 <- lmer(timestepchange ~ timesteplag + rtvmaxlag + rtumaxlag + omissionlag + (1|LunaID), bdf))
car::Anova(m2)

summary(m2 <- lmer(tschange_z ~ timesteplag + rtvmaxlag + rtumaxlag + omissionlag + (1|LunaID), bdf))
car::Anova(m2)

#what about scaling by parameter. grabbed df from sceptic_external_correlates
bdf2 <- select(df, lunaid, fmri_alpha_t, fmri_gamma_t, fmri_beta_t, fixed_uv_tau) %>% inner_join(bdf, c("lunaid" = "LunaID")) 

#summary(m2 <- lmer(timestepchange ~ timesteplag + rtvmaxlag + rtumaxlag*fmri_gamma_t + omissionlag + (1|lunaid), bdf2))
summary(m2 <- lmer(timestep ~ timesteplag + rtvmaxlag + rtumaxlag*fmri_gamma_t + omissionlag + (1|lunaid), bdf2))
car::Anova(m2)

summary(m2 <- lmer(timestep ~ timesteplag + rtvmax + rtumax*fmri_gamma_t*trial + omissionlag + (1|lunaid), bdf2))
car::Anova(m2)

summary(m1 <- lmer(timestep ~ timesteplag + rtvmax + rtumax*trial + omissionlag + (1|lunaid), bdf2))

summary(m2 <- lmer(timestep ~ timesteplag + rtvmax + rtumax*fixed_uv_tau + rtumax*trial + omissionlag + (1|lunaid), bdf2))
car::Anova(m2)

summary(m2 <- lmer(timestep ~ timesteplag + rtvmaxlag + rtumaxlag*fixed_uv_tau + rtumaxlag*trial + omissionlag + (1|lunaid), bdf2))
car::Anova(m2)

summary(m2 <- lmer(timestep ~ timesteplag + rtvmaxlag + rtumaxlag*fixed_uv_tau*trial + omissionlag + (1|lunaid), bdf2))
car::Anova(m2)

summary(m2 <- lmer(timestep ~ timesteplag + rtvmaxlag + rtumaxlag + omissionlag + (1|lunaid), bdf2))
car::Anova(m2)


anova(m1, m2)

#bdf2 <- bdf %>% gather(key=etype, value=entropy_mixlag, entropyHlag, entropyFlag) %>%
#    mutate_at(vars(rtvmaxlag, timesteplag, entropylag, rtumaxlag, entropy_mixlag), funs(cent=. - mean(., na.rm=TRUE))) #%>%

#try standardizing entropy between measures to avoid strange scaling differences
bdf2 <- bdf %>% mutate(entropyHlag=as.vector(scale(entropyHlag)), entropyFlag=as.vector(scale(entropyFlag))) %>% 
    gather(key=etype, value=entropy_mixlag, entropyHlag, entropyFlag) %>%
    mutate(inv_trial = 1/trial) %>%
    mutate_at(vars(trial, inv_trial, timesteplag, rtvmaxlag, entropylag, rtumaxlag, entropy_mixlag), funs(cent=. - mean(., na.rm=TRUE)))

library(lme4)
summary(m3 <- lmer(timestep ~ timesteplag_cent + rtvmaxlag_cent + rtumaxlag_cent*entropy_mixlag_cent*etype*trial_cent + run + (1 + run |LunaID), filter(bdf2, trial_abs > 4))) #omissionlag +
car::Anova(m3)

#break apart for a minute by etype (too many interactions!!)
summary(mfixed <- lmer(timestep ~ timesteplag_cent + rtvmaxlag_cent + rtumaxlag_cent*entropy_mixlag_cent*trial_cent + run + (1 |LunaID), filter(bdf2, etype=="entropyFlag" & trial_abs > 4)))
car::Anova(mfixed)

summary(mdecay <- lmer(timestep ~ timesteplag_cent + rtvmaxlag_cent + rtumaxlag_cent*entropy_mixlag_cent*trial_cent + run + (1 |LunaID), filter(bdf2, etype=="entropyHlag" & trial_abs > 4)))
car::Anova(mdecay)

cm <- lmerCellMeans(mdecay, fixat0=c("rtvmaxlag_cent", "timesteplag_cent", "run"), divide=c("entropy_mixlag_cent", "trial_cent"), n.cont=10)
ggplot(cm, aes(x=rtumaxlag_cent, y=timestep, color=entropy_mixlag_cent, ymin=timestep-se, ymax=timestep+se)) + geom_line() + geom_pointrange() + facet_wrap(~trial_cent)

#better if we use 1/trial? (asymptotic)
summary(mdecay2 <- lmer(timestep ~ timesteplag_cent + rtvmaxlag_cent + rtumaxlag_cent*entropy_mixlag_cent*inv_trial_cent + run + (1 |LunaID), filter(bdf2, etype=="entropyHlag" & trial_abs > 4)))
car::Anova(mdecay2)

anova(mdecay, mdecay2) #inv trial model fits *far* better (~90 AIC points)

#trial on x axis
cm <- lmerCellMeans(mdecay2, fixat0=c("rtvmaxlag_cent", "timesteplag_cent", "run"), divide=c("entropy_mixlag_cent", "rtumaxlag_cent"), n.cont=10)
cm$inv_trial_cent <- 1/(cm$inv_trial_cent + mean(bdf2$inv_trial, na.rm=TRUE))
ggplot(cm, aes(x=inv_trial_cent, y=timestep, color=entropy_mixlag_cent, ymin=timestep-se, ymax=timestep+se)) + 
    geom_line() + geom_pointrange() + facet_wrap(~rtumaxlag_cent)


#more nuanced prediction (without centering, which makes it hard to see)
summary(mdecay2 <- lmer(timestep ~ timesteplag_cent + rtvmaxlag_cent + rtumaxlag*entropy_mixlag_cent*inv_trial + run + (1 |LunaID), filter(bdf2, etype=="entropyHlag" & trial_abs > 4)))
car::Anova(mdecay2)

summary(mdecay2 <- lmer(timestep ~ timesteplag_cent*inv_trial + rtvmaxlag_cent + rtumaxlag*entropy_mixlag_cent*inv_trial + run + (1 |LunaID), filter(bdf2, etype=="entropyHlag" & trial_abs > 4)))
car::Anova(mdecay2)


cm <- lmerCellMeans(mdecay2, fixat0=c("rtvmaxlag_cent", "timesteplag_cent", "run"), divide=c("entropy_mixlag_cent"), n.cont=10,
    cont.pts=list(rtumaxlag=c(5, 15, 25, 35, 45), inv_trial=1/rev(c(3, 5, 15, 25, 35, 45, 47))))

cm$inv_trial <- 1/(cm$inv_trial)
#cm$timestep <- cm$timestep/10
#cm$se <- cm$se/10
pdf("for consideration.pdf", width=20, height=10)
ggplot(cm, aes(x=inv_trial, y=timestep, color=entropy_mixlag_cent, ymin=timestep-se, ymax=timestep+se)) + 
    geom_line(position=position_dodge(width=1.5)) + geom_pointrange(position=position_dodge(width=1.5)) +
    facet_wrap(~rtumaxlag, scales="free_y") + theme_bw(base_size=18) + xlab("Trial") + ylab("Predicted timestep (1-40)") +
    ggtitle("Trial (inverse) on X, Umax in panels")

#get U back on x axis...
ggplot(cm, aes(x=rtumaxlag, y=timestep, color=entropy_mixlag_cent, ymin=timestep-se, ymax=timestep+se)) + 
    geom_line(position=position_dodge(width=1.5)) + geom_pointrange(position=position_dodge(width=1.5)) +
    facet_wrap(~inv_trial, scales="free_y") + theme_bw(base_size=18) + xlab("Timestep of Umax") + ylab("Predicted timestep (1-40)") +
    ggtitle("Umax on X, trial in panels")
dev.off()


#try linear version (although I know it fits worse...)
summary(mdecay2 <- lmer(timestep ~ timesteplag_cent + rtvmaxlag_cent + rtumaxlag*entropy_mixlag_cent*trial + run + (1 |LunaID), filter(bdf2, etype=="entropyHlag" & trial_abs > 4)))
car::Anova(mdecay2)

cm <- lmerCellMeans(mdecay2, fixat0=c("rtvmaxlag_cent", "timesteplag_cent", "run"), divide=c("entropy_mixlag_cent"), n.cont=10,
    cont.pts=list(rtumaxlag=c(5, 25, 45), trial=c(1, 3, 5, 15, 25, 35, 45, 47, 49)))

#cm$inv_trial <- 1/(cm$inv_trial)
pdf("for consideration.pdf", width=9, height=6)
ggplot(cm, aes(x=trial, y=timestep, color=entropy_mixlag_cent, ymin=timestep-se, ymax=timestep+se)) + 
    geom_line() + geom_pointrange() + facet_wrap(~rtumaxlag, scales="free_y") + theme_bw(base_size=18)
dev.off()




#okay, see if we can handle the etype 4-way interaction conceptually...
bdf2$etype <- factor(bdf2$etype) #needed for lmerCellMeans to pick it up properly
summary(mboth <- lmer(timestep ~ timesteplag_cent + rtvmaxlag_cent + rtumaxlag_cent*entropy_mixlag_cent*trial_cent*etype + run + (1 |LunaID), filter(bdf2, trial_abs > 4)))
car::Anova(mboth)
cm <- lmerCellMeans(mboth, fixat0=c("rtvmaxlag_cent", "timesteplag_cent", "run"), divide=c("entropy_mixlag_cent", "trial_cent"), n.cont=10)

ggplot(cm, aes(x=rtumaxlag_cent, y=timestep, color=entropy_mixlag_cent, ymin=timestep-se, ymax=timestep+se)) + geom_line() + geom_pointrange() + facet_grid(etype~trial_cent)



bdf2cent <- bdf %>% gather(key=etype, value=entropy_mixlag, entropyHlag, entropyFlag) %>% group_by(LunaID, run) %>%
    mutate_at(vars(rtvmaxlag, timesteplag, entropylag, rtumaxlag, entropy_mixlag), funs(wicent=. - mean(., na.rm=TRUE), pmean=mean(., na.rm=TRUE))) %>% #within-person centering and person means
    ungroup() %>% mutate_at(vars(timesteplag, trial, rtvmaxlag_pmean, entropylag_pmean, rtumaxlag_pmean, entropy_mixlag_pmean), funs(c=. - mean(., na.rm=TRUE))) #between-person centering of person means
    
    
summary(m3 <- lmer(timestep ~ timesteplag_c + rtvmaxlag_pmean_c + rtvmaxlag_wicent + rtumaxlag_pmean_c*entropy_mixlag_pmean_c + rtumaxlag_wicent*entropy_mixlag_wicent + (1|LunaID), bdf2cent)) #omissionlag +
car::Anova(m3)

summary(m3 <- lmer(timestep ~ timesteplag_c + rtvmaxlag_pmean_c + rtvmaxlag_wicent + rtumaxlag_pmean_c*entropy_mixlag_pmean_c*etype + rtumaxlag_wicent*entropy_mixlag_wicent*etype + run + (1 |LunaID), filter(bdf2cent, trial_abs > 5))) #omissionlag +
car::Anova(m3)

bdf2cent$invtrial <- 1/bdf2cent$trial 
summary(m4 <- lmer(timestep ~ timesteplag_c + rtvmaxlag_pmean_c + rtvmaxlag_wicent + rtumaxlag_pmean_c*entropy_mixlag_pmean_c*etype + rtumaxlag_wicent*entropy_mixlag_wicent*etype*trial + run + (1 + run|LunaID), filter(bdf2cent, trial_abs > 5))) #omissionlag +
car::Anova(m4)

summary(m5 <- lmer(timestep ~ timesteplag_c + rtvmaxlag_pmean_c + rtvmaxlag_wicent + rtumaxlag_pmean_c*entropy_mixlag_pmean_c*etype + rtumaxlag_wicent*entropy_mixlag_wicent*etype*invtrial + run + (1 + run|LunaID), filter(bdf2cent, trial_abs > 5))) #omissionlag +
car::Anova(m5)

anova(m4, m5)


#try to look at effects in plot...
#pull out prior RT and value signals
bdf2cent$timestep_cleanup <- resid(lmer(timestep ~ timesteplag_c + rtvmaxlag_pmean_c + rtvmaxlag_wicent + (1|LunaID), bdf2cent, na.action=na.exclude))
summary(m3 <- lmer(timestep_cleanup ~ rtumaxlag_pmean_c*entropy_mixlag_pmean_c*etype + rtumaxlag_wicent*entropy_mixlag_wicent*etype + (1|LunaID), bdf2cent)) #omissionlag +
car::Anova(m3)



summary(m2 <- lmer(timestep ~ timesteplag_cent + rtvmaxlag_cent + rtumaxlag_cent*entropylag_cent + (1|LunaID), bdf2)) #omissionlag +
car::Anova(m2)

cm <- lmerCellMeans(m2, fixat0=c("rtvmaxlag_cent", "timesteplag_cent"), divide="entropylag_cent")

summary(m3 <- lmer(timestep ~ timesteplag_cent + rtvmaxlag_cent + rtumaxlag_cent*entropyHlag_cent + (1 + run|LunaID), bdf2)) #omissionlag +
car::Anova(m3)
cm3 <- lmerCellMeans(m3, fixat0=c("rtvmaxlag_cent", "timesteplag_cent"), divide="entropyHlag_cent", n.cont=10)

ggplot(cm3, aes(x=rtumaxlag_cent, y=timestep, ymin=timestep-se, ymax=timestep+se, color=entropyHlag_cent)) + geom_line() + geom_pointrange()

summary(m4 <- lmer(timestep ~ timesteplag_cent + rtvmaxlag_cent + rtumaxlag_cent*entropyFlag_cent + (1 + run |LunaID), bdf2)) #omissionlag +
car::Anova(m4)
cm4 <- lmerCellMeans(m4, fixat0=c("rtvmaxlag_cent", "timesteplag_cent"), divide="entropyFlag_cent", n.cont=10)

ggplot(cm4, aes(x=rtumaxlag_cent, y=timestep, ymin=timestep-se, ymax=timestep+se, color=entropyFlag_cent)) + geom_line() + geom_pointrange()

#mm <- merge(cm3, cm4, by="rtumaxlag_cent")
cm3$timestep <- as.vector(cm3$timestep)
cm4$timestep <- as.vector(cm4$timestep)
cm3 <- cm3 %>% rename(entropy = entropyHlag_cent) %>% mutate(model="Decay")
cm4 <- cm4 %>% rename(entropy = entropyFlag_cent) %>% mutate(model="Fixed")

mm <- rbind(cm3, cm4)
mm <- mm %>% mutate(entropy_labeled=sub("(entropyHlag_cent|entropyFlag_cent)", "Entropy", entropy))
pdf("draft fig9.pdf", width=10, height=8)
ggplot(mm, aes(x=rtumaxlag_cent, y=timestep, ymin=timestep-se, ymax=timestep+se, color=entropy_labeled)) + geom_line(position=position_dodge(width=2)) +
    geom_pointrange(position=position_dodge(width=2)) + facet_wrap(~model) + theme_bw(base_size=20)
dev.off()


#cm <- lmerCellMeans(m2, fixat0=c("rtvmaxlag_cent", "timesteplag_cent"), divide="rtumaxlag_cent")
#ggplot(cm, aes(x=entropylag_cent, y=timestep, ymin=timestep-se, ymax=timestep+se, color=omissionlag)) + geom_line() + geom_pointrange() + facet_wrap(~rtumaxlag_cent) 

cm <- lmerCellMeans(m2, fixat0=c("rtvmaxlag_cent", "timesteplag_cent"), divide="entropylag_cent")
ggplot(cm, aes(x=rtumaxlag_cent, y=timestep, ymin=timestep-se, ymax=timestep+se, color=omissionlag)) + geom_line() + geom_pointrange() + facet_wrap(~entropylag_cent) 


summary(m2 <- lmer(timestep ~ timesteplag + entropylag + omissionlag + (1|LunaID), bdf))
car::Anova(m2)



summary(m2 <- lmer(timestep ~ rtumax*trial  + (1|lunaid), bdf2))
car::Anova(m2)

mm <- lmerCellMeans(m2, divide="trial")
pdf("u effects.pdf", width=10, height=8)
ggplot(mm, aes(x=rtumax, y=timestep, color=trial)) + geom_line()
dev.off()


#bdf2 <- bdf2 %>% group_by(lunaid, run) %>% mutate(timestepchangelag=lag(timestepchange, order_by=trial)) %>% ungroup()
#summary(m2 <- lmer(timestepchange ~ timestepchangelag + rtvmaxlag + rtumaxlag*fmri_gamma_t + omissionlag + distfromedgelag + (1|lunaid), bdf2))
#car::Anova(m2)



#what about absolute timestep change?

test2 <- lmer(abstschange ~ omissionlag*absvdevlag*entropylag + (1|LunaID), bdf)

summary(test2)
cm3 <- lmerCellMeans(test2, divide="entropylag", n.cont=10)

pdf("entropytest_abs.pdf", width=15, height=6)
ggplot(cm3, aes(x=absvdevlag, y=abstschange, color=omissionlag, ymin=abstschange-se, ymax=abstschange+se)) + geom_line(size=2.5) + theme_bw(base_size=24) + geom_linerange(size=1.0) + 
    ylab("Mean trialwise RT change (cs)") + xlab("RT deviation from maximum value\non prior trial (cs)") + scale_color_brewer("Prior Outcome", palette="Set2") + geom_vline(xintercept=0) + geom_hline(yintercept=0) +
    facet_wrap(~entropylag, nrow=1)
dev.off()

summary(lmer(abstschange ~ absvdevlag + (1|LunaID), bdf))


#current entropy versus past trial?
test_past <- lmer(abstschange ~ omissionlag*absvdevlag*entropylag + (1|LunaID), bdf) 
#test_cur <- lmer(abstschange ~ omissionlag*absvdevlag*entropy*entropylag + (1|LunaID), bdf)
test_cur <- lmer(abstschange ~ omissionlag*absvdevlag*entropy + (1|LunaID), bdf)

cm4 <- lmerCellMeans(test_cur, divide="entropy", n.cont=10)

pdf("entropytest_cur_abs.pdf", width=15, height=6)
ggplot(cm4, aes(x=absvdevlag, y=abstschange, color=omissionlag, ymin=abstschange-se, ymax=abstschange+se)) + geom_line(size=2.5) + theme_bw(base_size=24) + geom_linerange(size=1.0) + 
    ylab("Mean trialwise RT change (cs)") + xlab("RT deviation from maximum value\non prior trial (cs)") + scale_color_brewer("Prior Outcome", palette="Set2") + geom_vline(xintercept=0) + geom_hline(yintercept=0) +
    facet_wrap(~entropy, nrow=1)
dev.off()

test_cur <- lmer(abstschange ~ omissionlag*absvdevlag*wizentropy*entropy + (1|LunaID), bdf)

summary(lmer(abstschange ~ entropy + (1|LunaID), bdf))

#build up: abs rt change as a function of prior omission and deviation from max value
m1 <- lmer(abstschange ~ omissionlag*absvdevlag + trial + run + (1 |LunaID), bdf)
cm1 <- lmerCellMeans(m1, n.cont=10)

pdf("m1_abschange.pdf", width=15, height=6)
ggplot(cm1, aes(x=absvdevlag, y=abstschange, color=omissionlag, ymin=abstschange-se, ymax=abstschange+se)) + geom_line(size=1.5) + theme_bw(base_size=24) + geom_linerange(size=1.5) + 
    ylab("Mean trialwise RT change (cs)") + xlab("RT deviation from maximum value\non prior trial (cs)") + scale_color_brewer("Prior Outcome", palette="Set2") + geom_vline(xintercept=0) + geom_hline(yintercept=0)
dev.off()

#similar to frank: relative RT swing size as a function of relative entropy (within run)
m2 <- lmer(wizabstschange ~ wizentropy + trial + run + (1 |LunaID), bdf)
cm2 <- lmerCellMeans(m2, n.cont=10, fixat0=c("trial", "run"))

pdf("m2abschange_wizentropy.pdf", width=15, height=6)
ggplot(cm2, aes(x=wizentropy, y=wizabstschange, ymin=wizabstschange-se, ymax=wizabstschange+se)) + geom_line(size=1.5) + theme_bw(base_size=24) + geom_linerange(size=1.5) + 
    ylab("z abs trialwise RT change") + xlab("within-run relative entropy") + geom_vline(xintercept=0) + geom_hline(yintercept=0)
dev.off()


pdf("m2abschange_wizentropy_scatter.pdf", width=15, height=6)
ggplot(bdf, aes(x=wizentropy, y=wizabstschange, group=LunaID)) + geom_point(size=1.0, alpha=0.5) + theme_bw(base_size=24) + stat_smooth(aes(group=NULL), se=TRUE) + 
    ylab("z abs trialwise RT change") + xlab("within-run relative entropy") + geom_vline(xintercept=0) + geom_hline(yintercept=0)
dev.off()

pdf("m2abschange_entropy_scatter.pdf", width=15, height=6)
ggplot(bdf, aes(x=entropy, y=wizabstschange, group=LunaID)) + geom_point(size=1.0, alpha=0.5) + theme_bw(base_size=24) + stat_smooth(aes(group=NULL), se=TRUE) + 
    ylab("z abs trialwise RT change") + xlab("within-run relative entropy") + geom_vline(xintercept=0) + geom_hline(yintercept=0)
dev.off()
