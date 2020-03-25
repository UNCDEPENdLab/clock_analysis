#entropy change betas
library(readr)
library(lme4)
library(tidyverse)
library(psych)
library(emmeans)

base_dir <- "/Users/mnh5174/Box/SCEPTIC_fMRI/specc_betas"
beta_dir <- file.path(base_dir, "sceptic-clock-feedback-v_entropy_change-preconvolve_fse_groupfixed", "v_entropy_change")

setwd(base_dir)

behav <- trial_df <- read.csv("specc_decay_factorize_selective_psequate_fixedparams_ffx_trial_statistics.csv.gz") %>%
  mutate(
    run_trial=case_when(
      trial >= 1 & trial <= 50 ~ trial,
      trial >= 51 & trial <= 100 ~ trial - 50L, #dplyr/rlang has gotten awfully picky about data types!!
      trial >= 101 & trial <= 150 ~ trial - 100L,
      trial >= 151 & trial <= 200 ~ trial - 150L,
      trial >= 201 & trial <= 250 ~ trial - 200L,
      trial >= 251 & trial <= 300 ~ trial - 250L,
      trial >= 301 & trial <= 350 ~ trial - 300L,
      trial >= 351 & trial <= 400 ~ trial - 350L,
      TRUE ~ NA_integer_),
    v_entropy_no5=if_else(run_trial <= 5, NA_real_, v_entropy),
    d_auc_sqrt=if_else(d_auc > 0, NA_real_, sqrt(-1*d_auc)), #only compute the sqrt of d_auc for negative (i.e., reasonable) observations
    v_entropy_sqrt=sqrt(v_entropy),
    rew_om=if_else(score_vba > 0, 1, 0)
  ) %>% #for win/loss maps
  group_by(id, run) %>%
  dplyr::mutate(   #compute rt_swing within run and subject
    rt_csv=rt_csv/1000, #rescale for ease in lmer
    rt_lag = dplyr::lag(rt_csv, 1, order_by=trial),
    rt_vmax_lag = dplyr::lag(rt_vmax, 1, order_by=trial),
    rt_vmax_change = abs(rt_vmax - rt_vmax_lag),
    omission = factor(as.numeric(score_csv > 0), levels=c(0,1), labels=c("Omission", "Reward")),
    omission_lag = lag(omission, order_by=trial),
    v_entropy_lag = dplyr::lag(v_entropy, 1, order_by=trial),
    v_entropy_change = v_entropy - v_entropy_lag, #change in entropy
    v_entropy_change_pos = v_entropy_change*(v_entropy_change > 0),
    v_entropy_change_neg = abs(v_entropy_change*(v_entropy_change < 0)),
    vmax_lag = dplyr::lag(v_max, 1, order_by=trial),
    pe_lag = dplyr::lag(pe_max, 1, order_by=trial),
    ppe_lag = if_else(pe_lag > 0, pe_lag, 0),
    npe_lag = if_else(pe_lag < 0, abs(pe_lag), 0),
    rt_swing = abs( c(NA, diff(rt_csv))),
    rt_swing_sqrt=sqrt(rt_swing)) %>%
  ungroup() %>% select(-dataset, -rt_next, -score_next, -clock_onset, -isi_onset, -feedback_onset, -iti_onset, -image) %>%
  rename(entropy=v_entropy, score=score_csv) %>%
  mutate(emotion=recode(emotion, scram="Scrambled", fear="Fear", happy="Happy"))


design <- read.table(file.path(beta_dir, "v_entropy_change-Intercept-Age_c-BPD_c-BPD_Age_design.txt"), header=TRUE) %>%
  select(feat_input_id, id, SPECC_ID, BPD, Age, Female, Age_c, Age_c, BPD_Age) %>%
  mutate(group=factor(BPD, levels=c(0,1), labels=c("HC", "BPD")),
         Female=factor(Female, levels=c(0,1), labels=c("Male", "Female")))

behav <- behav %>% inner_join(design)

#dust off 2017 analyses
summary(mmm <- lmer(entropy ~ emotion*group + rewFunc*group + (1 | id/run), behav, REML=FALSE))
car::Anova(mmm)


summary(mmm <- lmer(ev ~ emotion*group + rewFunc*group + (1 | id/run), behav %>% filter(rewFunc %in% c("IEV", "DEV")), REML=FALSE))
car::Anova(mmm)

toplot <- summary(lsmeans(mmm, ~emotion | rewFunc))
pdf("figures/ev effects on points.pdf", width=8, height=6)
ggplot(toplot, aes(x=emotion, y=lsmean, ymin=lsmean-SE, ymax=lsmean+SE, color=rewFunc)) + geom_pointrange(size=1.5) +
  ylab("Average EV") + xlab("Emotion") +
  theme_bw(base_size=22) + theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  scale_color_brewer("Contingency", palette="Dark2")
dev.off()

summary(mmm <- lmer(ev ~ emotion*group*Age_c + rewFunc*group*Age_c + (1 | id/run), behav %>% filter(rewFunc %in% c("IEV", "DEV")), REML=FALSE))
car::Anova(mmm, type="3")


#aggregate to run level
bdfruns <- behav %>% group_by(SPECC_ID, run) %>% filter(entropy > 0) %>% arrange(SPECC_ID, run, run_trial) %>% 
  dplyr::summarize(runreward=sum(score), 
                   earlyentropy=mean(entropy[run_trial > 1 & run_trial < 10])/mean(entropy), lateentropy=mean(entropy[run_trial > 40 & run_trial <= 50])/mean(entropy),
                   allentropy=mean(entropy), midentropy=mean(entropy[run_trial > 10 & run_trial <=40]/mean(entropy)),
                   #earlyentropyF=mean(entropyFixed[run_trial > 1 & run_trial < 10])/mean(entropyFixed), lateentropyF=mean(entropyFixed[run_trial > 40 & run_trial <= 50])/mean(entropyFixed),
                   #allentropyF=mean(entropyFixed), midentropyF=mean(entropyFixed[run_trial > 10 & run_trial <=40])/mean(entropyFixed),
                   elratio=earlyentropy/lateentropy, rewFunc=head(rewFunc, n=1), emotion=head(emotion, n=1), group=head(group, n=1), Age=head(Age, n=1)) %>%
  group_by(SPECC_ID) %>% #elratioF=earlyentropyF/lateentropyF, 
  mutate(elratioLag = lag(elratio, order_by=run)) %>% ungroup() %>% #elratioLagF = lag(elratioF, order_by=run
  mutate(emotion=relevel(emotion, ref = "Scrambled"))

bdfruns$Age_c <- bdfruns$Age - mean(bdfruns$Age)

bdfruns$elratio_wins <- winsor(bdfruns$elratio, trim=0.02) #trim a few outliers
bdfruns$runreward_wins <- winsor(bdfruns$runreward, trim=0.02) #trim a few outliers

car::Anova(lmer(runreward_wins ~ Age_c * group * rewFunc * emotion + elratio_wins + (1|SPECC_ID), filter( bdfruns, rewFunc %in% c("IEV", "DEV"))))
summary(lm(runreward_wins ~ elratio_wins + Age_c + group, bdfagg))


summary(lmer(runreward_wins ~ group + elratio_wins + (1|SPECC_ID), filter( bdfruns, rewFunc %in% c("IEV", "DEV"))))
summary(xx <- lmer(runreward_wins ~ group * elratio_wins * emotion + (1|SPECC_ID), filter( bdfruns, rewFunc %in% c("IEV", "DEV"))))

car::Anova(xx, type=3)
summary(lm(runreward_wins ~ elratio_wins + Age_c + group, bdfagg))

summary(lmer(elratio_wins ~ group * emotion + (1|SPECC_ID), filter( bdfruns, rewFunc %in% c("IEV", "DEV"))))
summary(xx <- lmer(runreward_wins ~ group * elratio_wins * emotion + (1|SPECC_ID), filter( bdfruns, rewFunc %in% c("IEV", "DEV"))))


#elratio_wins + 
summary(mm <- lmer(runreward_wins ~ (Age_c * group * emotion * rewFunc)^3 + (1|SPECC_ID), filter(bdfruns, rewFunc %in% c("IEV", "DEV"))))
car::Anova(mm, type="3")

#any evidence of learning failures in the fear condition (Where there are BPD effects)?
#nope, not really
summary(mm <- lmer(runreward_wins ~ (Age_c * group * rewFunc)^3 + (1|SPECC_ID), filter(bdfruns, rewFunc %in% c("IEV", "DEV") & emotion=="Fear")))
car::Anova(mm, type="3")

toplot <- summary(lsmeans(mm, ~emotion | rewFunc))
pdf("figures/emotion effects on points.pdf", width=8, height=6)
ggplot(toplot, aes(x=emotion, y=lsmean, ymin=lsmean-SE, ymax=lsmean+SE, color=rewFunc)) + geom_pointrange(size=1.5) +
  ylab("Average points earned in a single block") + xlab("Emotion") +
  theme_bw(base_size=22) + theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  scale_color_brewer("Contingency", palette="Dark2")
dev.off()

summary(mmtrial2 <- lmer(ev ~ emotion + rewFunc*ppe_lag + rewFunc*npe_lag + run + ppe_lag*emotion + npe_lag*emotion + (1 | SPECC_ID/run), behav))

#this doesn't make sense with so many zeros for ppe or npe
#but funny enough, the effects are essentially identical when we subset down to non-zero trials
summary(mmtrial3 <- lmer(ev ~ emotion*group + rewFunc*group*ppe_lag + rewFunc*group*npe_lag + Age_c*group*emotion + run + 
                           ppe_lag*group*emotion + npe_lag*group*emotion + (1 | SPECC_ID/run), behav))
car::Anova(mmtrial3)


npe_only <- filter(behav, npe_lag > 0 & rewFunc %in% c("DEV", "IEV"))
summary(negmodel <- lmer(ev ~ emotion*group + rewFunc*group*npe_lag + Age_c*group*emotion + run + vmax_lag + 
                           npe_lag*group*emotion + (1 | SPECC_ID/run), npe_only))
car::Anova(negmodel, type="3")

ppe_only <- filter(behav, ppe_lag > 0 & rewFunc %in% c("DEV", "IEV"))
summary(posmodel <- lmer(ev ~ emotion*group + rewFunc*group*ppe_lag + Age_c*group*emotion + run + vmax_lag +
                           ppe_lag*group*emotion + (1 | SPECC_ID/run), ppe_only))

car::Anova(posmodel, type="3")

ppeeff <- summary(lstrends(posmodel, ~ group | emotion, var="ppe_lag", infer=c(TRUE, TRUE), data=ppe_only))
npeeff <- summary(lstrends(negmodel, ~ group | emotion, var="npe_lag", infer=c(TRUE, TRUE), data=npe_only))

both <- rbind(ppeeff %>% rename(trend=ppe_lag.trend) %>% mutate(eff="Positive PE"), npeeff %>% rename(trend=npe_lag.trend) %>% mutate(eff="Negative PE"))

pdf("figures/PE shifts toward high EV BPD diff.pdf", width=11, height=8)
ggplot(both, aes(x=emotion, y=trend, ymin=trend-SE, ymax=trend+SE, color=group)) + facet_wrap(~eff, scales="free_y") + 
  geom_pointrange(position=position_dodge(width=0.5), size=1.6) +
  ylab("Shift toward high-value option after PE") + xlab("Emotion Stimulus") + scale_color_brewer("Group", palette="Set1") +
  theme_bw(base_size=24) + theme(axis.title.x=element_text(margin=margin(t=20, r=0,l=0, b=0)), axis.title.y=element_text(margin=margin(t=0, r=15, l=0, b=0)))
dev.off()


pdf("figures/Negative PE shifts toward high EV BPD diff.pdf", width=9, height=7)
ggplot(filter(both, eff=="Negative PE"), aes(x=emotion, y=trend, ymin=trend-SE, ymax=trend+SE, color=group)) + 
  geom_pointrange(position=position_dodge(width=0.5), size=1.6) +
  ylab("Shift toward high-value option after PE") + xlab("Emotion Stimulus") + scale_color_brewer("Group", palette="Set1") +
  theme_bw(base_size=24) + theme(axis.title.x=element_text(margin=margin(t=20, r=0,l=0, b=0)), 
                                 axis.title.y=element_text(margin=margin(t=0, r=15, l=0, b=0)),
                                 legend.key.size = unit(2.5, 'lines')) +
  geom_hline(yintercept=0)
dev.off()


#similar model, but focusing on RT
summary(negmodel <- lmer(rt_csv ~ rt_lag*emotion*group + rewFunc*group*npe_lag + Age_c*group*emotion + run + vmax_lag + 
                           npe_lag*group*emotion + (1 | SPECC_ID/run), npe_only))

car::Anova(negmodel, type=3)
#BPD shows reduced sensitivity to NPEs: in IEV; also shows reversal of NPE effect in happy
plot_model(negmodel, type="pred", terms=c("rewFunc", "npe_lag", "group"))
plot_model(negmodel, type="pred", terms=c("emotion", "npe_lag", "group"))

