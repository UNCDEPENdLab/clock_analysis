#entropy change betas
library(readr)
library(lme4)
library(dplyr)
library(tidyverse)
library(psych)
library(emmeans)

base_dir <- "/Users/mnh5174/Box/SCEPTIC_fMRI/specc_betas"
beta_dir <- file.path(base_dir, "sceptic-clock-feedback-v_entropy_change-preconvolve_fse_groupfixed", "v_entropy_change")

setwd(base_dir)

behav <- trial_df <- read.csv("specc_decay_factorize_selective_psequate_fixedparams_ffx_trial_statistics.csv.gz") %>%
  mutate(
    trial_rel=case_when(
      trial >= 1 & trial <= 50 ~ trial,
      trial >= 51 & trial <= 100 ~ trial - 50L, #dplyr/rlang has gotten awfully picky about data types!!
      trial >= 101 & trial <= 150 ~ trial - 100L,
      trial >= 151 & trial <= 200 ~ trial - 150L,
      trial >= 201 & trial <= 250 ~ trial - 200L,
      trial >= 251 & trial <= 300 ~ trial - 250L,
      trial >= 301 & trial <= 350 ~ trial - 300L,
      trial >= 351 & trial <= 400 ~ trial - 350L,
      TRUE ~ NA_integer_),
    v_entropy_no5=if_else(trial_rel <= 5, NA_real_, v_entropy),
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
    rt_swing = abs( c(NA, diff(rt_csv))),
    rt_swing_sqrt=sqrt(rt_swing)) %>%
  ungroup() %>% dplyr::select(-dataset, -rt_next, -score_next, -clock_onset, -isi_onset, -feedback_onset, -iti_onset, -image)

design <- read.table(file.path(beta_dir, "v_entropy_change-Intercept-Age_c-BPD_c-BPD_Age_design.txt"), header=TRUE) %>%
    dplyr::select(feat_input_id, id, SPECC_ID, BPD, Age, Female, Age_c, Age_c, BPD_Age)
metadata <- read.csv(file.path(beta_dir, "v_entropy_change_cluster_metadata.csv")) %>%
  filter(model=="Intercept-Age_c-BPD_c-BPD_Age") %>% droplevels()
subj_betas <- read.csv(file.path(beta_dir, "v_entropy_change_subj_betas.csv")) %>%
  dplyr::select(-Intercept, -Age, -BPD, -Age_c, -BPD_c, -BPD_Age, -I_Age_c, -BPD_IAge) %>%
  filter(model=="Intercept-Age_c-BPD_c-BPD_Age") %>% droplevels()
run_betas <- read.csv(file.path(beta_dir, "v_entropy_change_run_betas.csv.gz")) %>% 
  dplyr::select(-Intercept, -Age, -BPD, -Age_c, -BPD_c, -BPD_Age, -I_Age_c, -BPD_IAge) %>%
  filter(model=="Intercept-Age_c-BPD_c-BPD_Age") %>% droplevels()

#overall betas for v_entropy_change
# sanity checks -- all passed
# s2 <- subj_betas %>% inner_join(design)
# s3 <- s2 %>% inner_join(metadata)
# nrow(s2) == nrow(subj_betas)
# nrow(s3) == nrow(subj_betas)

subj_betas <- subj_betas %>% inner_join(design) %>% inner_join(metadata) %>%
  mutate(BPD=factor(BPD, levels=c(0,1), labels=c("HC", "BPD")),
         Female=factor(Female, levels=c(0,1), labels=c("Male", "Female")))

#16 clusters -- but the last 4 are small (43 vox or less, and nothign interesting)
overall_betas <- subj_betas %>% filter(l2_contrast=="overall" & l3_contrast=="BPD_c") %>%
  dplyr::select(-model, -l1_contrast) %>% droplevels() %>% filter(cluster_number <= 12)

for_fa <- overall_betas %>% dplyr::select(id, cluster_number, cope_value) %>%
  spread(key=cluster_number, value = cope_value, sep="_")

for_fa_wins <- apply(for_fa %>% dplyr::select(-id), 2, function(col) { psych::winsor(col, trim=0.02 )})

#clusters 7 and 8 don't want to cooperate
#for_fa_wins <- for_fa_wins[,c(-7,-8)]
for_fa_wins <- for_fa_wins[,c(-7)]

#principal(for_fa, nfactors = 3)
f1 <- fa(for_fa_wins, nfactors = 1, fm="minres")
f2 <- fa(for_fa_wins, nfactors = 2, fm="minres", rotate="oblimin")
f3 <- fa(for_fa_wins, nfactors = 3, fm="minres", rotate = "oblimin")
f4 <- fa(for_fa_wins, nfactors = 4, fm="minres", rotate = "oblimin")
vss(for_fa_wins, rotate="varimax", fm="ml")

#5, 6, and 11 are essentially putamen. And tend to load with each other.
#not sure why 9 likes to go with them.

#3 and 12 are essentially parahippocampal/fusiform. 10 (medial cerebellum) likes to go with them.

#1, 2, 4, 8 go together and are essentially L lateral IPS (2), rACC (1), L adPFC (4), and precuneus/MCC (8)

f3_scores <- factor.scores(for_fa_wins, f3)$scores %>% as.data.frame() %>% setNames(c("c1_putamen", "c2_parahippo", "c3_fp")) %>% bind_cols(for_fa %>% dplyr::select(id))
#f1 is putamen
#f2 is parahippocampal/fusiform
#f3 is something more cog controlly


#join with behavior
behav_w_clus <- behav %>% inner_join(f3_scores) %>% inner_join(design) %>%
  mutate(BPD=factor(BPD, levels=c(0,1), labels=c("HC", "BPD")),
         Female=factor(Female, levels=c(0,1), labels=c("Male", "Female")),
         id=factor(id), run=factor(run)
         ) %>%
  group_by(id, run) %>%
  mutate(v_entropy_change_wi=v_entropy_change - mean(v_entropy_change, na.rm=TRUE), v_entropy_change_pm=mean(v_entropy_change, na.rm=TRUE)) %>%
  ungroup()

  
behav %>% anti_join(f3_scores) %>% pull(id) %>% unique()

summary(m1 <- lmer(rt_csv ~ rt_lag + v_entropy + c1_putamen + c2_parahippo + c3_fp + (1|id/run), behav_w_clus))
summary(m2 <- lmer(rt_csv ~ rt_lag + v_entropy + c1_putamen + (1|id/run), behav_w_clus))
summary(m3 <- lmer(rt_csv ~ rt_lag + v_entropy + c2_parahippo + (1|id/run), behav_w_clus))
summary(m3 <- lmer(rt_csv ~ rt_lag + v_entropy + c3_fp + (1|id/run), behav_w_clus)) #FP cluster appears to speed RTs in general

summary(m3 <- lmer(rt_csv ~ rt_lag + v_entropy * c3_fp + (1|id/run), behav_w_clus))
summary(m3 <- lmer(rt_csv ~ rt_lag + v_entropy * c3_fp * BPD + (1|id/run), behav_w_clus))
summary(m3 <- lmer(rt_csv ~ rt_lag + v_entropy_change * c3_fp * BPD + (1|id/run), behav_w_clus))


summary(m3 <- lmer(rt_csv ~ rt_lag + v_entropy_change_pm + v_entropy_change_wi * c3_fp * BPD + (1|id/run), behav_w_clus))
car::Anova(m3, type="3") #yep

summary(m4 <- lmer(rt_csv ~ rt_lag + v_entropy_change_pm + v_entropy_change_wi * c3_fp * BPD + emotion * v_entropy_change_wi + rewFunc + (1|id/run), behav_w_clus))
car::Anova(m4, type="3") #yep

#add omission lag (doesn't knock out effect)
summary(m4 <- lmer(rt_csv ~ rt_lag + omission_lag*rewFunc + v_entropy_change_pm + 
                     v_entropy_change_wi * c3_fp * BPD + emotion * v_entropy_change_wi + rewFunc + (1|id/run),
                   behav_w_clus, control=lmerControl(optCtrl=list(maxfun=20000), optimizer="Nelder_Mead")))
car::Anova(m4, type="3") #yep

# bobyqa complains about tolerance, though it's super close to the .002 mark
# behav_w_clus_z <- behav_w_clus %>% mutate_if(is.numeric, list(function(x) { as.vector(scale(x)) }) )
# summary(m4 <- lmer(rt_csv ~ rt_lag + omission_lag*rewFunc + v_entropy_change_pm + v_entropy_change_wi * c3_fp * BPD + emotion * v_entropy_change_wi + rewFunc + (1|id/run), behav_w_clus_z))

emtrends(m4, specs="emotion", var="v_entropy_change_wi") #so, in general, in fear condition, entropy increases go with slow down, and decreases with speed-up.
#but in happy and scram, changes don't affect RTs.

plot_model(m4, type="pred", terms=c("v_entropy_change_wi", "c3_fp", "BPD"))

#looks like a cross-over interaction. High values in controls
df <- summary(emtrends(m4, specs= ~ c3_fp*BPD, var="v_entropy_change_wi",
              at=list(c3_fp=quantile(f3_scores$c3_fp, seq(.2, .8, .2)))))

ggplot(df, aes(x=c3_fp, y=v_entropy_change_wi.trend, ymin=v_entropy_change_wi.trend-SE, ymax=v_entropy_change_wi.trend+SE, color=BPD)) + geom_pointrange() + geom_line() +
  geom_rug(data=f3_scores %>% merge(design), aes(x=c3_fp, y=NULL, ymin=NULL, ymax=NULL, color=factor(BPD)))

ggplot(behav_w_clus, aes(y=v_entropy_change_wi, x=c3_fp, color=BPD)) + stat_smooth(method="loess")


hist(behav_w_clus$v_entropy_change_wi)

summary(m5 <- lmer(rt_csv ~ rt_lag + omission_lag*rewFunc + # v_entropy_change_pm + 
                     v_entropy_change_wi * c3_fp * BPD + emotion * v_entropy_change_wi + rewFunc + (1|id/run),
                   behav_w_clus %>% filter(v_entropy_change_wi > 0), control=lmerControl(optCtrl=list(maxfun=20000), optimizer="Nelder_Mead")))
car::Anova(m5, type="3")


summary(m6 <- lmer(rt_csv ~ rt_lag + omission_lag*rewFunc + # v_entropy_change_pm + 
                     v_entropy_change_wi * c3_fp * BPD + emotion * v_entropy_change_wi + rewFunc + (1|id/run),
                   behav_w_clus %>% filter(v_entropy_change_wi < 0), control=lmerControl(optCtrl=list(maxfun=20000), optimizer="Nelder_Mead")))
car::Anova(m6, type="3")



clusters <- f3_scores %>% inner_join(design)
#indeed, all factor scores are higher in BPD
t.test(c1_putamen ~ BPD, clusters)
tapply(clusters$c1_putamen, clusters$BPD, mean)
t.test(c2_parahippo ~ BPD, clusters)
tapply(clusters$c2_parahippo, clusters$BPD, mean)
t.test(c3_fp ~ BPD, clusters)
tapply(clusters$c3_fp, clusters$BPD, mean)

#BPD group shows far higher betas in this region, consistent with selection for cluster
ggplot(clusters, aes(x=factor(BPD), y=c3_fp)) + geom_boxplot()


#####
#does this model apply to other clusters?
#parahippo seems to do nothing here
summary(m5 <- lmer(rt_csv ~ rt_lag + omission_lag*rewFunc + v_entropy_change_pm + 
                     v_entropy_change_wi * c2_parahippo * BPD + emotion * v_entropy_change_wi + rewFunc + (1|id/run),
                   behav_w_clus, control=lmerControl(optCtrl=list(maxfun=20000), optimizer="Nelder_Mead")))
car::Anova(m5, type="3")


#putamen does nothing for RTs per se-- actually, seems to have some control over stickiness (mirror of rt_swing?)
summary(m6 <- lmer(rt_csv ~ rt_lag*c1_putamen + omission_lag*rewFunc + v_entropy_change_pm + rt_vmax_lag + #rt_lag*c1_putamen*BPD
                     v_entropy_change_wi * c1_putamen * BPD + emotion * v_entropy_change_wi + rewFunc + (1|id/run),
                   behav_w_clus, control=lmerControl(optCtrl=list(maxfun=20000), optimizer="Nelder_Mead")))
car::Anova(m6, type="3")

plot_model(m6, type = "pred", terms=c("rt_lag", "c1_putamen")) #, "BPD"))

#deeper dive into putamen
#putamen is controlling magnitude of RT swing -- higher c1 in BPD drives greater RT swings in response to entropy changes

summary(m6 <- lmer(rt_swing ~ omission_lag*rewFunc + v_entropy_change_pm + scale(rt_vmax_lag) +
                     v_entropy_change_wi * c1_putamen * BPD + emotion * v_entropy_change_wi * c1_putamen + rewFunc + (1|id/run),
                   behav_w_clus, control=lmerControl(optCtrl=list(maxfun=20000), optimizer="Nelder_Mead")))
car::Anova(m6, type="3")

# emmip(m6, ~ v_entropy_change_wi * c1_putamen | BPD)
# emtrends(m6, var="v_entropy_change_wi", ~ v_entropy_change_wi * c1_putamen | BPD)

library(sjPlot)

#handy!
plot_model(m6, type = "pred", terms=c("v_entropy_change_wi", "c1_putamen", "BPD"))
plot_model(m6, type = "est")#, terms=c("v_entropy_change_wi", "c1_putamen", "BPD"))

#so, in general, in fear condition, entropy increases go with slow down, and decreases with speed-up.
#but in happy and scram, changes don't affect RTs.


summary(m7 <- lmer(rt_swing ~ omission_lag*rewFunc + v_entropy_change_pm + scale(rt_vmax_lag) + scale(rt_vmax) +
                     v_entropy_change_wi * c1_putamen * BPD + emotion * v_entropy_change_wi * c1_putamen + rewFunc + (1|id/run),
                   behav_w_clus %>% filter(v_entropy_change_wi > 0), control=lmerControl(optCtrl=list(maxfun=20000), optimizer="Nelder_Mead")))
car::Anova(m7, type="3")

plot_model(m7, type = "pred", terms=c("v_entropy_change_wi", "c1_putamen [meansd]", "BPD"), 
           se=FALSE, ci.lvl = .60) + 
  ylab("RT swing (secs)") + scale_color_brewer("Putamen\nbeta", palette="Dark2") +
  xlab("Entropy increase (t vs. t-1)") + ggtitle("RT swings by group and entropy increase") +
  theme_bw(base_size=13)

dat <- ggpredict(m7, type="fe", terms=c("v_entropy_change_wi", "c1_putamen [meansd]", "BPD"))

pdf("RT swing entropy change betas.pdf", width=8, height=4)
ggplot(dat, aes(x = x, y = predicted, colour = group)) +
  stat_smooth(method = "lm", se = FALSE, fullrange = TRUE, size=1.4) +
  facet_wrap(~facet) +
  labs(
    colour = "Putamen\nbeta",
    y = "RT swing (secs)",
    x = "Entropy increase (t vs. t-1)",
    title = "RT swings by group and entropy increase"
  ) + scale_color_brewer(palette="Dark2", breaks=c("-1.01", "0", "1"), labels=c("-1 SD", "M", "+ 1 SD")) +
  theme_bw(base_size=16)
dev.off()
  



summary(m8 <- lmer(rt_swing ~ omission_lag*rewFunc + v_entropy_change_pm + scale(rt_vmax_lag) +
                     v_entropy_change_wi * c1_putamen * BPD + emotion * v_entropy_change_wi * c1_putamen + rewFunc + (1|id/run),
                   behav_w_clus %>% filter(v_entropy_change_wi < 0), control=lmerControl(optCtrl=list(maxfun=20000), optimizer="Nelder_Mead")))
car::Anova(m8, type="3")

plot_model(m8, type = "pred", terms=c("v_entropy_change_wi", "c1_putamen", "BPD"))


###

#replicate above for other clusters?
summary(m9 <- lmer(rt_swing ~ omission_lag*rewFunc + v_entropy_change_pm + scale(rt_vmax_lag) +
                     v_entropy_change_wi * c2_parahippo * BPD + emotion * v_entropy_change_wi * c2_parahippo + rewFunc + (1|id/run),
                   behav_w_clus %>% filter(v_entropy_change_wi > 0), control=lmerControl(optCtrl=list(maxfun=20000), optimizer="Nelder_Mead")))
car::Anova(m9, type="3")

plot_model(m9, type = "pred", terms=c("v_entropy_change_wi", "c2_parahippo", "BPD")) #yes

summary(m10 <- lmer(rt_swing ~ omission_lag*rewFunc + v_entropy_change_pm + scale(rt_vmax_lag) +
                     v_entropy_change_wi * c3_fp * BPD + emotion * v_entropy_change_wi * c3_fp + rewFunc + (1|id/run),
                   behav_w_clus %>% filter(v_entropy_change_wi > 0), control=lmerControl(optCtrl=list(maxfun=20000), optimizer="Nelder_Mead")))
car::Anova(m10, type="3")

plot_model(m10, type = "pred", terms=c("v_entropy_change_wi", "c3_fp", "BPD"))


summary(m11 <- lmer(rt_swing ~ omission_lag*rewFunc + v_entropy_change_pm + scale(rt_vmax_lag) +
                      v_entropy_change_wi * c3_fp * BPD + emotion * v_entropy_change_wi * c3_fp + rewFunc + (1|id/run),
                    behav_w_clus %>% filter(v_entropy_change_wi < 0), control=lmerControl(optCtrl=list(maxfun=20000), optimizer="Nelder_Mead")))
car::Anova(m11, type="3")

plot_model(m11, type = "pred", terms=c("v_entropy_change_wi", "c3_fp", "BPD"))



### entropy as outcome

summary(m10 <- lmer(v_entropy ~ pe_max_pm + c1_vis * BPD * emotion + rt_vmax_lag + emotion * rewFunc + (1|id/run), behav_w_clus))
car::Anova(m10, type="3")

summary(m10 <- lmer(v_entropy ~ emotion * rewFunc * BPD + (1|id/run), behav_w_clus %>% filter(emotion != "scram")))
car::Anova(m10, type="3")

#incl trial?
summary(m10 <- lmer(v_entropy ~ emotion * rewFunc * BPD * trial + (1|id/run), behav_w_clus %>% filter(emotion != "scram")))
car::Anova(m10, type="3")

#fear only
summary(m10 <- lmer(v_entropy ~ rewFunc * BPD * trial + (1|id/run), behav_w_clus %>% filter(emotion == "fear")))
car::Anova(m10, type="3")

#happy only (this is where the action is)
summary(m10 <- lmer(v_entropy ~ rewFunc * BPD * trial + (1|id/run), behav_w_clus %>% filter(emotion == "happy")))
car::Anova(m10, type="3")

library(sjPlot)
plot_model(m10, type="pred", terms=c("trial", "rewFunc", "BPD"))

#look at relative/within run trial
summary(m11 <- lmer(v_entropy ~ emotion * rewFunc * BPD * trial_rel + (1|id/run), behav_w_clus %>% filter(emotion != "scram")))
car::Anova(m11, type="3")

#fear rel
summary(m12 <- lmer(v_entropy ~ rewFunc * BPD * trial_rel + (1|id/run), behav_w_clus %>% filter(emotion == "fear")))
car::Anova(m12, type="3")

summary(m15 <- lmer(v_entropy ~ rewFunc * BPD * trial_rel + (1|id/run), behav_w_clus %>% filter(emotion == "fear" & trial_rel > 5)))
car::Anova(m15, type="3")


#happy rel (this is where action is for rel trial, too.)
summary(m13 <- lmer(v_entropy ~ rewFunc * BPD * trial_rel + (1|id/run), behav_w_clus %>% filter(emotion == "happy")))
car::Anova(m13, type="3")

plot_model(m13, type="pred", terms=c("trial_rel", "rewFunc", "BPD"))

#need to filter early trials
summary(m14 <- lmer(v_entropy ~ rewFunc * BPD * trial_rel + (1|id/run), behav_w_clus %>% filter(emotion == "happy" & trial_rel > 5)))
car::Anova(m14, type="3")

plot_model(m14, type="pred", terms=c("trial_rel", "rewFunc", "BPD"))

summary(m14 <- lmer(v_max ~ BPD * trial_rel + (1|id/run), behav_w_clus %>% filter(emotion == "happy" & trial_rel > 5)))
car::Anova(m14, type="3")

plot_model(m14, type="pred", terms=c("trial_rel", "BPD"))



#mean differences in entropy by group?
behav_bw <- behav_w_clus %>% filter(emotion != "scram") %>% group_by(id, BPD, emotion, rewFunc) %>%
  summarize(emean=mean(v_entropy, na.rm=TRUE)) %>% ungroup() %>% droplevels() %>% mutate(BPD=factor(BPD))
library(afex)
library(emmeans)
mm <- aov_ez(data=behav_bw, dv="emean", id="id", between = c("BPD"), within = c("emotion", "rewFunc"))
summary(mm)

#plot(emmeans(mm, ~rewFunc|BPD))
plot(emmeans(mm, ~rewFunc))


summary(m10 <- lmer(v_entropy ~ emotion * rewFunc * BPD + (1|id/run), behav_w_clus %>% filter(emotion != "scram")))
car::Anova(m10, type="3")



##################

run_betas <- run_betas %>% inner_join(design) %>% inner_join(metadata) %>%
  mutate(BPD=factor(BPD, levels=c(0,1), labels=c("HC", "BPD")),
         Female=factor(Female, levels=c(0,1), labels=c("Male", "Female")))

#run betas
overall_betas <- run_betas %>% filter(l2_contrast=="overall" & l3_contrast=="BPD_c") %>%
  dplyr::select(-model, -l1_contrast) %>% droplevels() %>% filter(cluster_number <= 12)

for_fa <- overall_betas %>% dplyr::select(id, run_num, cluster_number, cope_value) %>%
  spread(key=cluster_number, value = cope_value, sep="_") %>% dplyr::select(-id, -run_num)

lattice::splom(for_fa)

vss(for_fa, rotate="varimax", fm="ml")
f2 <- fa(for_fa, nfactors = 2, fm="ml", rotate="oblimin")
f3 <- fa(for_fa, nfactors = 3, fm="ml", rotate = "oblimin")
f4 <- fa(for_fa, nfactors = 4, fm="ml", rotate = "oblimin")

Sx <- cov(for_fa)
D2 <- mahalanobis(for_fa, colMeans(for_fa), Sx)
D2[which(D2 > quantile(D2, .75) + 1.5*IQR(D2))]

for_fa_wins <- apply(for_fa, 2, function(col) { psych::winsor(col, trim=0.02 )})

lattice::splom(for_fa_wins)

vss(for_fa_wins, rotate="varimax", fm="ml")
f2 <- fa(for_fa_wins, nfactors = 2, fm="ml", rotate="oblimin")
f3 <- fa(for_fa_wins, nfactors = 3, fm="ml", rotate = "oblimin")

