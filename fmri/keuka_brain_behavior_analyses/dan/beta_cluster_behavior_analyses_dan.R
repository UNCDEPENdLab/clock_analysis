# brain-to-behavior analyses with anterior (low entropy) and posterior (PE) hippocampal cluster betas
# first run beta_cluster_import_pca_clean.R if not run once already

library(dplyr)
library(tidyverse)
library(psych)
library(corrplot)
library(lme4)
library(ggpubr)
library(lmerTest)
library(stargazer)
library(car)
library(sjstats)
library(sjPlot)
library(emmeans)
library(cowplot)
# install_github("UNCDEPENdLab/dependlab")
library(dependlab)
source('~/code/Rhelpers/screen.lmerTest.R')
source('~/code/Rhelpers/vif.lme.R')
# library(stringi)

clock_folder <- "~/Data_Analysis/clock_analysis" #michael
# clock_folder <- "~/code/clock_analysis" #alex

# source('~/code/Rhelpers/')
setwd(file.path(clock_folder, 'fmri/keuka_brain_behavior_analyses/dan'))

### load data
# load('trial_df_and_vhdkfpe_clusters.Rdata')
# cleaner version with only H, PE and uncertainty trial vars
load('trial_df_and_vh_pe_hd_clusters_u.Rdata')

############# Main analysis using Schaeffer-based betas

### quick hippocampal sanity check
# check PH PEs extracted at higher threshold
# t <- df %>% select(id, pe_PH_r) %>% unique()
# 
# test <-  lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + rt_vmax_lag_sc + last_outcome + 
#                                 v_max_wi_lag + v_entropy_wi + h_HippAntL + pe_PH_r)^2 + 
#                    rt_lag_sc:last_outcome:h_HippAntL + 
#                    rt_lag_sc:last_outcome:pe_PH_r +
#                    rt_vmax_lag_sc:trial_neg_inv_sc:h_HippAntL + 
#                    rt_vmax_lag_sc:trial_neg_inv_sc:pe_PH_r  +
#                    (1|id/run), df %>% filter(rt_csv<4000))
# screen.lmerTest(test, .05)
# Anova(test, '3')
# 
# mtest <-  lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + rt_vmax_lag_sc + last_outcome + 
#                              v_max_wi_lag + v_entropy_wi + h_HippAntL + pe_PH_r)^2 + 
#                 rt_lag_sc:last_outcome:h_HippAntL + 
#                 rt_lag_sc:last_outcome:pe_PH_r +
#                 rt_vmax_lag_sc:trial_neg_inv_sc:h_HippAntL + 
#                 rt_vmax_lag_sc:trial_neg_inv_sc:pe_PH_r  +
#                 (1|id/run), mdf)# %>% filter(rt_csv<4000))
# screen.lmerTest(mtest, .05)
# Anova(mtest, '3')

# # compare MEG and fMRI
# df$session <- "fMRI"
# mdf$session <- "MEG"
# bdf <- bind_rows(df, mdf)
# ggplot(bdf, aes(run_trial, rt_csv, color = rewFunc, lty = session)) + geom_smooth()
# ggplot(bdf, aes(run_trial, ev, color = rewFunc, lty = session)) + geom_smooth()

############# fMRI
# Effect of DAN on exploration
mb_dan1 <-  lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + rt_vmax_lag_sc + last_outcome + 
                                v_max_wi_lag + v_entropy_wi + dan_parietal + pe_ips)^2 + 
                   rt_lag_sc:last_outcome:dan_parietal + 
                   rt_lag_sc:last_outcome:pe_ips +
                   rt_vmax_lag_sc:trial_neg_inv_sc:dan_parietal + 
                   rt_vmax_lag_sc:trial_neg_inv_sc:pe_ips  +
                   (1|id/run), df %>% filter(rt_csv<4000))
screen.lmerTest(mb_dan1, .05)
summary(mb_dan1)
Anova(mb_dan1, '3')

# basic model with only the general entropy factor:
mb_dan1a <-  lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + rt_vmax_lag_sc + last_outcome + 
                                v_max_wi_lag + v_entropy_wi + general_entropy)^2 + 
                   rt_lag_sc:last_outcome:general_entropy + 
                    rt_lag_sc:trial_neg_inv_sc: v_entropy_wi:general_entropy + 
                   rt_vmax_lag_sc:trial_neg_inv_sc:general_entropy + 
                   (1|id/run), df %>% filter(rt_csv<4000))
screen.lmerTest(mb_dan1a, .05)
summary(mb_dan1a)

# add vlPFC entropy region
mb_dan2 <-  lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + rt_vmax_lag_sc + last_outcome + 
                                v_max_wi_lag + v_entropy_wi + general_entropy + med_par + fef + entropy_vlPFC)^2 + 
                   rt_lag_sc:last_outcome:general_entropy + 
                   rt_lag_sc:last_outcome:med_par +
                   rt_lag_sc:last_outcome:fef +
                   rt_lag_sc:last_outcome:entropy_vlPFC + 
                   rt_vmax_lag_sc:trial_neg_inv_sc:general_entropy + 
                   rt_vmax_lag_sc:trial_neg_inv_sc:med_par  +
                   rt_vmax_lag_sc:trial_neg_inv_sc:fef  +
                   rt_vmax_lag_sc:trial_neg_inv_sc:entropy_vlPFC + 
                   (1|id/run), df %>% filter(rt_csv<4000))
screen.lmerTest(mb_dan2, .05)
summary(mb_dan2)
Anova(mb_dan2, '3')
anova(mb_dan1, mb_dan2)


################
# replication
mmb_dan1 <-  lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + rt_vmax_lag_sc + last_outcome + 
                                v_max_wi_lag + v_entropy_wi + dan_parietal + pe_ips)^2 + 
                   rt_lag_sc:last_outcome:dan_parietal + 
                   rt_lag_sc:last_outcome:pe_ips +
                   rt_vmax_lag_sc:trial_neg_inv_sc:dan_parietal + 
                   rt_vmax_lag_sc:trial_neg_inv_sc:pe_ips  +
                   (1|id/run), mdf %>% filter(rt_csv<4000))
screen.lmerTest(mmb_dan1, .05)
summary(mmb_dan1)
Anova(mmb_dan1, '3')

# add vlPFC entropy region
mmb_dan2 <-  lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + rt_vmax_lag_sc + last_outcome + 
                                 v_max_wi_lag + v_entropy_wi + general_entropy + med_par + fef + entropy_vlPFC)^2 + 
                    rt_lag_sc:last_outcome:general_entropy + 
                    rt_lag_sc:last_outcome:med_par +
                    rt_lag_sc:last_outcome:fef +
                    rt_lag_sc:last_outcome:entropy_vlPFC + 
                    rt_vmax_lag_sc:trial_neg_inv_sc:general_entropy + 
                    rt_vmax_lag_sc:trial_neg_inv_sc:med_par  +
                    rt_vmax_lag_sc:trial_neg_inv_sc:fef  +
                    rt_vmax_lag_sc:trial_neg_inv_sc:entropy_vlPFC + 
                    (1|id/run), mdf %>% filter(rt_csv<4000))
screen.lmerTest(mmb_dan2, .05)
summary(mmb_dan2)
Anova(mmb_dan2, '3')
anova(mmb_dan1, mmb_dan2)

# selection history analyses
mb_dan_lag3 <- lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + rt_lag2_sc + rt_lag3_sc + rt_vmax_lag_sc + last_outcome + 
                                   v_max_wi_lag + v_entropy_wi + general_entropy + med_par + fef)^2 + 
                      rt_lag_sc:last_outcome:general_entropy + 
                      rt_lag_sc:last_outcome:med_par +
                      rt_lag_sc:last_outcome:fef +
                      rt_lag_sc:trial_neg_inv_sc:general_entropy + 
                      rt_lag_sc:trial_neg_inv_sc:med_par +
                      rt_lag_sc:trial_neg_inv_sc:fef +
                      rt_lag2_sc:trial_neg_inv_sc:general_entropy + 
                      rt_lag2_sc:trial_neg_inv_sc:med_par +
                      rt_lag2_sc:trial_neg_inv_sc:fef +
                      rt_lag3_sc:trial_neg_inv_sc:general_entropy + 
                      rt_lag3_sc:trial_neg_inv_sc:med_par +
                      rt_lag3_sc:trial_neg_inv_sc:fef +
                      rt_vmax_lag_sc:trial_neg_inv_sc:general_entropy + 
                      rt_vmax_lag_sc:trial_neg_inv_sc:med_par  +
                      rt_vmax_lag_sc:trial_neg_inv_sc:fef  +
                      (1|id/run), df %>% filter(rt_csv<4000))
screen.lmerTest(mb_dan_lag3, .01)
summary(mb_dan_lag3)
Anova(mb_dan_lag3, '3')


mmb_dan_lag3 <- lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + rt_lag2_sc + rt_lag3_sc + rt_vmax_lag_sc + last_outcome + 
                                    v_max_wi_lag + v_entropy_wi + general_entropy + med_par + fef)^2 + 
                       rt_lag_sc:last_outcome:general_entropy + 
                       rt_lag_sc:last_outcome:med_par +
                       rt_lag_sc:last_outcome:fef +
                       rt_vmax_lag_sc:trial_neg_inv_sc:general_entropy + 
                       rt_vmax_lag_sc:trial_neg_inv_sc:med_par  +
                       rt_vmax_lag_sc:trial_neg_inv_sc:fef  +
                       (1|id/run), mdf %>% filter(rt_csv<4000))
screen.lmerTest(mmb_dan_lag3, .05)
summary(mmb_dan1)
Anova(mmb_dan1, '3')

# see if the last three outcomes are held in the FEF buffer

mb_dan_lag3_rew <- lmer(rt_csv_sc ~  rt_lag_sc*v_entropy_wi*trial_neg_inv_sc +
                           rt_lag_sc*last_outcome*general_entropy + 
                           rt_lag_sc*last_outcome*med_par +
                           rt_lag_sc*last_outcome*fef +
                           rt_lag2_sc*last_outcome_lag2*general_entropy + 
                           rt_lag2_sc*last_outcome_lag2*med_par +
                           rt_lag2_sc*last_outcome_lag2*fef +
                           rt_lag3_sc*last_outcome_lag3*general_entropy + 
                           rt_lag3_sc*last_outcome_lag3*med_par +
                           rt_lag3_sc*last_outcome_lag3*fef +
                           rt_vmax_lag_sc*trial_neg_inv_sc*general_entropy + 
                           rt_vmax_lag_sc*trial_neg_inv_sc*med_par  +
                           rt_vmax_lag_sc*trial_neg_inv_sc*fef  +
                           (1|id/run), df %>% filter(rt_csv<4000))
screen.lmerTest(mb_dan_lag3_rew, .01)
summary(mb_dan_lag3_rew)
Anova(mb_dan_lag3_rew, '3')



mmb_dan_lag3_rew <- lmer(rt_csv_sc ~  rt_lag_sc*v_entropy_wi*trial_neg_inv_sc +
                          rt_lag_sc*last_outcome*general_entropy + 
                          rt_lag_sc*last_outcome*med_par +
                          rt_lag_sc*last_outcome*fef +
                          rt_lag2_sc*last_outcome_lag2*general_entropy + 
                          rt_lag2_sc*last_outcome_lag2*med_par +
                          rt_lag2_sc*last_outcome_lag2*fef +
                          rt_lag3_sc*last_outcome_lag3*general_entropy + 
                          rt_lag3_sc*last_outcome_lag3*med_par +
                          rt_lag3_sc*last_outcome_lag3*fef +
                          rt_vmax_lag_sc*trial_neg_inv_sc*general_entropy + 
                          rt_vmax_lag_sc*trial_neg_inv_sc*med_par  +
                          rt_vmax_lag_sc*trial_neg_inv_sc*fef  +
                          (1|id/run), mdf %>% filter(rt_csv<4000))
screen.lmerTest(mmb_dan_lag3_rew, .01)
summary(mmb_dan_lag3_rew)
Anova(mmb_dan_lag3_rew, '3')

# TEST CONDITION EFFECTS 
mb_dan_lag3_cond <- lmer(rt_csv_sc ~  rt_lag_sc*trial_neg_inv_sc*rewFunc*dan_parietal +
                           rt_lag_sc*trial_neg_inv_sc*rewFunc*pe_ips +
                          rt_lag_sc*last_outcome*dan_parietal + 
                          rt_lag_sc*last_outcome*pe_ips +
                          # rt_lag2_sc*last_outcome_lag2*dan_parietal + 
                          # rt_lag2_sc*last_outcome_lag2*pe_ips +
                          # rt_lag3_sc*last_outcome_lag3*dan_parietal + 
                          # rt_lag3_sc*last_outcome_lag3*pe_ips +
                          rt_vmax_lag_sc*trial_neg_inv_sc*dan_parietal + 
                          rt_vmax_lag_sc*trial_neg_inv_sc*pe_ips  +
                          (1|id/run), df %>% filter(rt_csv<4000))
screen.lmerTest(mb_dan_lag3_cond, .05)
summary(mb_dan_lag3_rew)
Anova(mb_dan_lag3_rew, '3')

mmb_dan_lag3_cond <- lmer(rt_csv_sc ~  rt_lag_sc*trial_neg_inv_sc*rewFunc*dan_parietal +
                           rt_lag_sc*trial_neg_inv_sc*rewFunc*pe_ips +
                           rt_lag_sc*last_outcome*dan_parietal + 
                           rt_lag_sc*last_outcome*pe_ips +
                           # rt_lag2_sc*last_outcome_lag2*dan_parietal + 
                           # rt_lag2_sc*last_outcome_lag2*pe_ips +
                           # rt_lag3_sc*last_outcome_lag3*dan_parietal + 
                           # rt_lag3_sc*last_outcome_lag3*pe_ips +
                           rt_vmax_lag_sc*trial_neg_inv_sc*dan_parietal + 
                           rt_vmax_lag_sc*trial_neg_inv_sc*pe_ips  +
                           (1|id/run), mdf %>% filter(rt_csv<4000))
screen.lmerTest(mmb_dan_lag3_cond, .05)
summary(mmb_dan_lag3_rew)
Anova(mmb_dan_lag3_rew, '3')

###############
# predict score
ms <- lmer(score_csv ~  (rewFunc + trial_neg_inv_sc + pe_ips + dan_parietal)^3 +
             (1|id), df %>% filter(rt_csv<4000))
summary(ms)
car::Anova(ms)

mms <- lmer(score_csv ~  (rewFunc + trial_neg_inv_sc + pe_ips + dan_parietal)^3 +
             (1|id), mdf %>% filter(rt_csv<4000))
summary(mms)
car::Anova(mms)

mev <- lmer(ev ~  (rewFunc + trial_neg_inv_sc + pe_ips + dan_parietal)^3 +
             (1|id), df %>% filter(rt_csv<4000))
summary(mev)
car::Anova(ms)

mmev <- lmer(ev ~  (rewFunc + trial_neg_inv_sc + pe_ips + dan_parietal)^3 +
              (1|id), mdf %>% filter(rt_csv<4000))
summary(mmev)
car::Anova(mms)



ms1 <- lmer(ev ~  (rewFunc + trial_neg_inv_sc + rt_lag_sc + rt_csv_sc + pe_ips + dan_parietal + pe_f3_str)^3 +
             (1|id), df %>% filter(rt_csv<4000))
screen.lmerTest(ms1, .01)

mms1 <- lmer(ev ~  (rewFunc + trial_neg_inv_sc + rt_lag_sc + rt_csv_sc + pe_ips + dan_parietal + pe_f3_str)^3 +
              (1|id), mdf %>% filter(rt_csv<4000))
screen.lmerTest(mms1, .01)


summary(ms1)
car::Anova(ms1)
emtrends(ms, ~rewFunc, var="rt_lag_sc")

################
# Entropy change
################

# Effect of DAN entropy change on exploration
mb_hd1 <-  lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + rt_vmax_lag_sc + last_outcome + 
                                v_max_wi_lag + v_entropy_wi + hd_f2_DAN)^2 + 
                   rt_lag_sc:last_outcome:hd_f2_DAN + 
                   rt_vmax_lag_sc:trial_neg_inv_sc:hd_f2_DAN + 
                   (1|id/run), df %>% filter(rt_csv<4000))
screen.lmerTest(mb_hd1, .05)
summary(mb_hd1)
Anova(mb_hd1, '3')

# replication
mmb_hd1 <-  lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + rt_vmax_lag_sc + last_outcome + 
                               v_max_wi_lag + v_entropy_wi + hd_f2_DAN)^2 + 
                  rt_lag_sc:last_outcome:hd_f2_DAN + 
                  rt_vmax_lag_sc:trial_neg_inv_sc:hd_f2_DAN + 
                  (1|id/run), mdf %>% filter(rt_csv<4000))
screen.lmerTest(mmb_hd1, .05)
summary(mmb_hd1)
Anova(mmb_hd1, '3')

# add other factors:
# Effect of DAN, salience, rlPFC entropy change on exploration
mb_hd2 <-  lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + rt_vmax_lag_sc + last_outcome + 
                               v_max_wi_lag + v_entropy_wi + hd_f2_DAN + hd_f1_salience + hd_f3_rlPFC)^2 + 
                  rt_lag_sc:last_outcome:hd_f2_DAN + 
                  rt_lag_sc:last_outcome:hd_f1_salience + 
                  rt_lag_sc:last_outcome:hd_f3_rlPFC + 
                  rt_vmax_lag_sc:trial_neg_inv_sc:hd_f2_DAN + 
                  rt_vmax_lag_sc:trial_neg_inv_sc:hd_f1_salience + 
                  rt_vmax_lag_sc:trial_neg_inv_sc:hd_f3_rlPFC +
                  (1|id/run), df %>% filter(rt_csv<4000))
screen.lmerTest(mb_hd2, .05)
summary(mb_hd2)
Anova(mb_hd2, '3')

# replication
mmb_hd2 <-  lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + rt_vmax_lag_sc + last_outcome + 
                               v_max_wi_lag + v_entropy_wi + hd_f2_DAN + hd_f1_salience + hd_f3_rlPFC)^2 + 
                  rt_lag_sc:last_outcome:hd_f2_DAN + 
                  rt_lag_sc:last_outcome:hd_f1_salience + 
                  rt_lag_sc:last_outcome:hd_f3_rlPFC + 
                  rt_vmax_lag_sc:trial_neg_inv_sc:hd_f2_DAN + 
                  rt_vmax_lag_sc:trial_neg_inv_sc:hd_f1_salience + 
                  rt_vmax_lag_sc:trial_neg_inv_sc:hd_f3_rlPFC +
                  (1|id/run), mdf %>% filter(rt_csv<4000))
screen.lmerTest(mmb_hd2, .05)
summary(mb_hd2)
Anova(mb_hd2, '3')


###############
## PE clusters
###############

# pe_ips and h_ips
mb_pe_ips_dan <-  lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + rt_vmax_lag_sc + last_outcome + 
                                       v_max_wi_lag + v_entropy_wi + pe_ips + dan_parietal)^3 + 
                          rt_lag_sc:last_outcome:pe_ips:dan_parietal + 
                          rt_vmax_lag_sc:trial_neg_inv_sc:pe_ips:dan_parietal +
                          (1|id/run), df %>% filter(rt_csv<4000))
screen.lmerTest(mb_pe_ips_dan, .01)

cm <- as_tibble(lmer_predict(mb_pe_ips_dan, divide=c("rt_lag_sc", "pe_ips", "dan_parietal"), 
                   fixatmean = c("rt_vmax_lag_sc", "v_entropy_wi", "v_max_wi_lag", "trial_neg_inv_sc")))
ggplot(cm, aes(as.numeric(rt_lag_sc), rt_csv_sc, color = pe_ips, group = pe_ips) ) + geom_line() + facet_wrap(last_outcome~dan_parietal)


mmb_pe_ips_dan <-  lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + rt_vmax_lag_sc + last_outcome + 
                                      v_max_wi_lag + v_entropy_wi + pe_ips + dan_parietal)^3 + 
                         rt_lag_sc:last_outcome:pe_ips:dan_parietal + 
                         rt_vmax_lag_sc:trial_neg_inv_sc:pe_ips:dan_parietal +
                         (1|id/run), mdf %>% filter(rt_csv<4000))
screen.lmerTest(mmb_pe_ips_dan, .01)

cm <- as_tibble(lmer_predict(mmb_pe_ips_dan, divide=c("rt_lag_sc", "pe_ips", "dan_parietal"), 
                             fixat0 = c("rt_vmax_lag_sc", "v_entropy_wi", "v_max_wi_lag", "trial_neg_inv_sc")))
ggplot(cm, aes(as.numeric(rt_lag_sc), rt_csv_sc, color = pe_ips, group = pe_ips) ) + geom_line() + facet_wrap(last_outcome~dan_parietal)


### check if the RTvmax convergence in high-IPS PE beta subjects is early responding in DEV
mb_pe_ips_dan <-  lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + rt_vmax_lag_sc + last_outcome + 
                                      v_max_wi_lag + v_entropy_wi + pe_ips + dan_parietal)^3 + 
                         rt_lag_sc:last_outcome:pe_ips:dan_parietal + 
                         rt_vmax_lag_sc:trial_neg_inv_sc:pe_ips:dan_parietal +
                         (1|id/run), df %>% filter(rt_csv<4000))
screen.lmerTest(mb_pe_ips_dan, .01)

mmb_pe_ips_dan <-  lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + rt_vmax_lag_sc + last_outcome + 
                                       v_max_wi_lag + v_entropy_wi + pe_ips + dan_parietal)^3 + 
                          rt_lag_sc:last_outcome:pe_ips:dan_parietal + 
                          rt_vmax_lag_sc:trial_neg_inv_sc:pe_ips:dan_parietal +
                          (1|id/run), mdf %>% filter(rt_csv<4000))
screen.lmerTest(mmb_pe_ips_dan, .01)





#### model-predicted effects for buffer models

em1 <- as_tibble(emtrends(mb_dan_lag3_rew, var = "rt_lag_sc", specs = c("fef", "last_outcome"), at = list(fef = c(-1.5, 1.17)), options = list()))
em1$study = 'fMRI'
em2 <- as_tibble(emtrends(mmb_dan_lag3_rew, var = "rt_lag_sc", specs = c("fef", "last_outcome"), at = list(fef = c(-1.5, 1.17)), options = list()))
em2$study = 'Replication'
em1 <- rbind(em1, em2)
p1 <- ggplot(em1, aes(x=last_outcome, y=rt_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=as.factor(fef))) + 
  #shape = as.factor(pe_f2_hipp), 
  #p1 <- ggplot(em1, aes(x=last_outcome, y=rt_lag_sc.trend, linetype = as.factor(pe_f2_hipp), ymin=asymp.LCL, ymax=asymp.UCL)) + 
  geom_point(position = position_dodge(width = .6), size=2.5) + 
  #geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), position = position_dodge2(width=0.9), width=0.5) + 
  geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + 
  #geom_crossbar(position=position_dodge(width=0.7), width=0.6) + 
  theme_bw(base_size=12) + facet_wrap(~study)+ ylab("RT swings (AU)\n Small <---------> Large")  + 
  #scale_shape_manual(values=c(15,16), labels = c("10th %ile", "90th %ile")) + 
  #scale_color_brewer("PH RPE\nresponse", palette="Set1", labels = c("10th %ile", "90th %ile")) +
  scale_color_manual("FEF\nresponse", values=c("#1b3840","#4fa3b8"), labels = c("10th %ile", "90th %ile")) +
  labs(shape = "FEF\nresponse") +
  theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) + 
  scale_y_reverse(limits = c(.7, -.1)) 
em12 <- as_tibble(emtrends(mb_dan_lag3_rew, var = "rt_lag2_sc", specs = c("fef", "last_outcome_lag2"), at = list(fef = c(-1.5, 1.17)), options = list()))
em12$study = 'fMRI'
em22 <- as_tibble(emtrends(mmb_dan_lag3_rew, var = "rt_lag2_sc", specs = c("fef", "last_outcome_lag2"), at = list(fef = c(-1.5, 1.17)), options = list()))
em22$study = 'Replication'
em12 <- rbind(em12, em22)
p2 <- ggplot(em12, aes(x=last_outcome_lag2, y=rt_lag2_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=as.factor(fef))) + 
  #shape = as.factor(pe_f2_hipp), 
  #p1 <- ggplot(em1, aes(x=last_outcome, y=rt_lag2_sc, linetype = as.factor(pe_f2_hipp), ymin=asymp.LCL, ymax=asymp.UCL)) + 
  geom_point(position = position_dodge(width = .6), size=2.5) + 
  #geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), position = position_dodge2(width=0.9), width=0.5) + 
  geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + 
  #geom_crossbar(position=position_dodge(width=0.7), width=0.6) + 
  theme_bw(base_size=12) + facet_wrap(~study)+ ylab("RT swings (AU)\n Small <---------> Large")  + 
  #scale_shape_manual(values=c(15,16), labels = c("10th %ile", "90th %ile")) + 
  #scale_color_brewer("PH RPE\nresponse", palette="Set1", labels = c("10th %ile", "90th %ile")) +
  scale_color_manual("FEF\nresponse", values=c("#1b3840","#4fa3b8"), labels = c("10th %ile", "90th %ile")) +
  labs(shape = "FEF\nresponse") +
  theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) + 
  scale_y_reverse(limits = c(.7, -.1)) 

em13 <- as_tibble(emtrends(mb_dan_lag3_rew, var = "rt_lag3_sc", specs = c("fef", "last_outcome_lag3"), at = list(fef = c(-1.5, 1.17)), options = list()))
em13$study = 'fMRI'
em23 <- as_tibble(emtrends(mmb_dan_lag3_rew, var = "rt_lag3_sc", specs = c("fef", "last_outcome_lag3"), at = list(fef = c(-1.5, 1.17)), options = list()))
em23$study = 'Replication'
em13 <- rbind(em13, em23)
p3 <- ggplot(em13, aes(x=last_outcome_lag3, y=rt_lag3_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=as.factor(fef))) + 
  #shape = as.factor(pe_f2_hipp), 
  #p1 <- ggplot(em1, aes(x=last_outcome, y=rt_lag_sc.trend, linetype = as.factor(pe_f2_hipp), ymin=asymp.LCL, ymax=asymp.UCL)) + 
  geom_point(position = position_dodge(width = .6), size=2.5) + 
  #geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), position = position_dodge2(width=0.9), width=0.5) + 
  geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + 
  #geom_crossbar(position=position_dodge(width=0.7), width=0.6) + 
  theme_bw(base_size=12) + facet_wrap(~study)+ ylab("RT swings (AU)\n Small <---------> Large")  + 
  #scale_shape_manual(values=c(15,16), labels = c("10th %ile", "90th %ile")) + 
  #scale_color_brewer("PH RPE\nresponse", palette="Set1", labels = c("10th %ile", "90th %ile")) +
  scale_color_manual("FEF\nresponse", values=c("#1b3840","#4fa3b8"), labels = c("10th %ile", "90th %ile")) +
  labs(shape = "FEF\nresponse") +
  theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) + 
  scale_y_reverse(limits = c(.7, -.1)) 


ggsave("fef1.pdf %>% filter(rt_csv<4000)", p1, width=5, height=4)



# model diagnostics

#residuals plots

plot(mb_dan1)
plot(mmb_dan1)

# try lm

lm_fmri_dan1 <-  lm(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + rt_vmax_lag_sc + last_outcome + 
                                   v_max_wi_lag + v_entropy_wi + h_HippAntL_neg +  pe_f2_hipp + DAN)^2 + 
                      rt_lag_sc:last_outcome:h_HippAntL_neg + 
                      rt_lag_sc:last_outcome:pe_f2_hipp +
                      rt_lag_sc:last_outcome:DAN +
                      rt_vmax_lag_sc:trial_neg_inv_sc:h_HippAntL_neg + 
                      rt_vmax_lag_sc:trial_neg_inv_sc:pe_f2_hipp  +
                      rt_vmax_lag_sc:trial_neg_inv_sc:DAN, df %>% filter(rt_csv<4000), na.action = na.exclude)
plot(lm_fmri_dan1)
# adding contingency does not capture structure

rlm <- resid(lm_fmri_dan1)
rdf <- cbind(df, rlm)
ggplot(rdf, aes(rt_csv, rlm)) + geom_smooth()

# covary out distance from the edge
lm_fmri_dan_sqrt <-  lm(rt_csv^.5 ~ (trial_neg_inv_sc + I(rt_lag^.5) + rt_vmax_lag_sc + last_outcome + 
                                       v_max_wi_lag + v_entropy_wi + h_HippAntL_neg +  pe_f2_hipp + DAN)^2 + 
                          I(rt_lag^.5):last_outcome:h_HippAntL_neg + 
                          I(rt_lag^.5):last_outcome:pe_f2_hipp +
                          I(rt_lag^.5):last_outcome:DAN +
                          rt_vmax_lag_sc:trial_neg_inv_sc:h_HippAntL_neg + 
                          rt_vmax_lag_sc:trial_neg_inv_sc:pe_f2_hipp  +
                          rt_vmax_lag_sc:trial_neg_inv_sc:DAN, df, na.action = na.exclude)
plot(lm_fmri_dan_sqrt)


#### PLOTS #######

ggplot(df, aes(run_trial, rt_csv, color = rewFunc, lty = fef>0)) + geom_smooth(method = 'loess') + facet_wrap(~rewFunc)
ggplot(df, aes(run_trial, rt_csv, color = rewFunc, lty = general_entropy>0)) + geom_smooth(method = 'loess') + facet_wrap(~rewFunc)


ggplot(df, aes(run_trial, rt_swing, color = rewFunc, lty = fef>0)) + geom_smooth(method = 'loess') + facet_wrap(~rewFunc)
ggplot(df, aes(run_trial, rt_swing, color = rewFunc, lty = general_entropy>0)) + geom_smooth(method = 'loess') + facet_wrap(~rewFunc)


## DAN replication forest plot
mterms <- names(fixef(mb_dan1))
setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/dan/plots/')
dan <- plot_models(mb_dan1,mmb_dan1, rm.terms = mterms[c(-25, -31)], m.labels = c("fMRI", "replication"),
                   show.values = T, std.est = "std2", legend.title = "Session", vline.color = "slategray3",
                   wrap.labels = 20, axis.labels = c("RT(t-1) * DAN high entropy response", "RT(Vmax) * DAN high entropy response"), 
                   axis.title = "Less convergence <==> Better convergence")
dan <- dan + ylim(-.01,.6) + geom_hline(yintercept = 0, color = "slategray3")
pdf("dan_beta_models_replication.pdf", height = 3, width = 5)
dan
dev.off()

## PH replication plot
setwd(file.path(clock_folder, 'fmri/keuka_brain_behavior_analyses/plots/'))
ph <- plot_models(mb3hpe_hipp,mmb3hpe_hipp, rm.terms = mterms[c(-22, -39)], m.labels = c("fMRI", "replication"),
                  show.values = T,  std.est = "std2", legend.title = "Session", vline.color = "slategray3",
                  wrap.labels = 15,  axis.labels = c("Previous RT * Omission * Post. hippocampal PE response","Previous RT * Post. hippocampal PE response"),
                  axis.title = "Greater RT swing  <==>  Smaller RT swing")
ph <- ph + ggplot2::ylim(-.1,.1)
pdf("ph_beta_models_replication.pdf", height = 3, width = 5)
ph
dev.off()

############################
## Emtrends plot for betas -> behavior figure
em1 <- as_tibble(emtrends(mb3hpe_hipp, var = "rt_lag_sc", specs = c("pe_f2_hipp", "last_outcome"), at = list(pe_f2_hipp = c(-1.12, 1.07)), options = list()))
em1$study = 'fMRI'
em2 <- as_tibble(emtrends(mmb3hpe_hipp, var = "rt_lag_sc", specs = c("pe_f2_hipp", "last_outcome"), at = list(pe_f2_hipp = c(-1.12, 1.07)), options = list()))
em2$study = 'Replication'
em1 <- rbind(em1, em2)
p1 <- ggplot(em1, aes(x=last_outcome, y=rt_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=as.factor(pe_f2_hipp))) + 
  #shape = as.factor(pe_f2_hipp), 
  #p1 <- ggplot(em1, aes(x=last_outcome, y=rt_lag_sc.trend, linetype = as.factor(pe_f2_hipp), ymin=asymp.LCL, ymax=asymp.UCL)) + 
  geom_point(position = position_dodge(width = .6), size=2.5) + 
  #geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), position = position_dodge2(width=0.9), width=0.5) + 
  geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + 
  #geom_crossbar(position=position_dodge(width=0.7), width=0.6) + 
  theme_bw(base_size=12) + facet_wrap(~study)+ ylab("RT swings (AU)\n Small <---------> Large")  + 
  #scale_shape_manual(values=c(15,16), labels = c("10th %ile", "90th %ile")) + 
  #scale_color_brewer("PH RPE\nresponse", palette="Set1", labels = c("10th %ile", "90th %ile")) +
  scale_color_manual("PH RPE\nresponse", values=c("#1b3840","#4fa3b8"), labels = c("10th %ile", "90th %ile")) +
  labs(shape = "PH RPE\nresponse") +
  theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) + 
  scale_y_reverse(limits = c(.6, 0)) 

p <- ggarrange(p1,p2,p3, ncol = 3, labels = c("Lag1", "Lag2", "Lag3"))
pdf("fef_lags.pdf", width = 12, height = 4)
p
dev.off()
ggsave("p1.pdf", p, width=10, height=4)

em3 <- as_tibble(emtrends(mb3hpe_hipp, var = "rt_lag_sc", specs = c("h_HippAntL_neg", "last_outcome"), at = list(h_HippAntL_neg = c(-.1, .37)), options = list()))
em3$study = 'fMRI'
em4 <- as_tibble(emtrends(mmb3hpe_hipp, var = "rt_lag_sc", specs = c("h_HippAntL_neg", "last_outcome"), at = list(h_HippAntL_neg = c(-.1, .37)), options = list()))
em4$study = 'Replication'
em2 <- rbind(em3, em4)
# p2 <- ggplot(em2, aes(last_outcome, rt_lag_sc.trend, lty = as.factor(h_HippAntL_neg))) + geom_point(position = position_dodge2(width = .9)) + geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), position = position_dodge2()) + 
#   theme_bw() + facet_wrap(~study) +  ylab("RT swings (AU)\n Small <---------> Large") + scale_linetype(labels = c("10th %ile", "90th %ile")) + labs(lty = "AH global max\nresponse") +
#   theme(axis.title.x=element_blank()) + scale_y_reverse(limits = c(.6, 0)) 

p2 <- ggplot(em2, aes(x=last_outcome, y=rt_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=as.factor(h_HippAntL_neg))) + 
  #shape = as.factor(h_HippAntL_neg), 
  #p1 <- ggplot(em1, aes(x=last_outcome, y=rt_lag_sc.trend, linetype = as.factor(pe_f2_hipp), ymin=asymp.LCL, ymax=asymp.UCL)) + 
  geom_point(position = position_dodge(width = .6), size=2.5) + 
  #geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), position = position_dodge2(width=0.9), width=0.5) + 
  geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + 
  #geom_crossbar(position=position_dodge(width=0.7), width=0.6) + 
  theme_bw(base_size=12) + facet_wrap(~study)+ ylab("RT swings (AU)\n Small <---------> Large")  + 
  #scale_shape_manual(values=c(19,23), labels = c("10th %ile", "90th %ile")) + 
  #scale_color_brewer("AH global max\nresponse", palette="Dark2", labels = c("10th %ile", "90th %ile")) +
  scale_color_manual("AH global max\nresponse", values=c("#403202", "#e2b407"), labels = c("10th %ile", "90th %ile")) +
  labs(shape = "AH global max\nresponse") +
  theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) + 
  scale_y_reverse(limits = c(.6, 0)) 

ggarrange(p1,p2)

em5 <- as_tibble(emtrends(mb3hpe_hipp, var = "rt_vmax_lag_sc", specs = c("pe_f2_hipp", "trial_neg_inv_sc"), at = list(pe_f2_hipp = c(-1.12, 1.07), trial_neg_inv_sc = c(-.7, 0.44)), options = list()))
em5$study = 'fMRI'
em6 <- as_tibble(emtrends(mmb3hpe_hipp, var = "rt_vmax_lag_sc", specs = c("pe_f2_hipp", "trial_neg_inv_sc"), at = list(pe_f2_hipp = c(-1.12, 1.07), trial_neg_inv_sc = c(-.7, 0.44)), options = list()))
em6$study = 'Replication'
em3 <- rbind(em5,em6)
# p3 <- ggplot(em3, aes(as.factor(trial_neg_inv_sc), rt_vmax_lag_sc.trend, lty = as.factor(pe_f2_hipp))) + geom_point(position = position_dodge2(width = .9)) + geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), position = position_dodge2()) + 
#   theme_bw() + facet_wrap(~study) +  ylab("Convergence on\nbest RT (AU)") + scale_linetype(labels = c("10th %ile", "90th %ile")) + labs(lty = "PH RPE\nresponse") +
#   scale_x_discrete(name ="Trial", labels=c("-0.7" = "5", "0.44" = "50")) + scale_y_continuous(limits = c(.05, .3))

p3 <- ggplot(em3, aes(x=as.factor(trial_neg_inv_sc), y=rt_vmax_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=as.factor(pe_f2_hipp))) + 
  #shape = as.factor(pe_f2_hipp), 
  #p1 <- ggplot(em1, aes(x=last_outcome, y=rt_lag_sc.trend, linetype = as.factor(pe_f2_hipp), ymin=asymp.LCL, ymax=asymp.UCL)) + 
  geom_point(position = position_dodge(width = .6), size=2.5) + 
  #geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), position = position_dodge2(width=0.9), width=0.5) + 
  geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + 
  #geom_crossbar(position=position_dodge(width=0.7), width=0.6) + 
  theme_bw(base_size=12) + facet_wrap(~study)+ ylab("Convergence on\nbest RT (AU)") +
  #scale_shape_manual(values=c(15,16), labels = c("10th %ile", "90th %ile")) + 
  #scale_color_brewer("PH RPE\nresponse", palette="Set1", labels = c("10th %ile", "90th %ile")) +
  scale_color_manual("PH RPE\nresponse", values=c("#1b3840", "#4fa3b8"), labels = c("10th %ile", "90th %ile")) +
  labs(shape = "PH RPE\nresponse") +
  theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) + 
  scale_x_discrete(name ="Trial", labels=c("-0.7" = "5", "0.44" = "50")) + scale_y_continuous(limits = c(.05, .3))



em7 <- as_tibble(emtrends(mb3hpe_hipp, var = "rt_vmax_lag_sc", specs = c("h_HippAntL_neg", "trial_neg_inv_sc"), at = list(h_HippAntL_neg = c(-.1, .37), trial_neg_inv_sc = c(-.7, 0.44)), options = list()))
em7$study = 'fMRI'
em8 <- as_tibble(emtrends(mmb3hpe_hipp, var = "rt_vmax_lag_sc", specs = c("h_HippAntL_neg", "trial_neg_inv_sc"), at = list(h_HippAntL_neg = c(-.1, .37), trial_neg_inv_sc = c(-.7, 0.44)), options = list()))
em8$study = 'Replication'
em4 <- rbind(em7, em8)
# p4 <- ggplot(em4, aes(as.factor(trial_neg_inv_sc), rt_vmax_lag_sc.trend, lty = as.factor(h_HippAntL_neg))) + geom_point(position = position_dodge2(width = .9)) + geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), position = position_dodge2()) + 
#   theme_bw() + facet_wrap(~study) + ylab("Convergence on\nbest RT (AU)")  + scale_linetype(labels = c("10th %ile", "90th %ile")) + labs(lty = "AH global max\nresponse")  +
#   scale_x_discrete(name ="Trial", labels=c("-0.7" = "5", "0.44" = "50")) + scale_y_continuous(limits = c(.05, .3))

p4 <- ggplot(em4, aes(x=as.factor(trial_neg_inv_sc), y=rt_vmax_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=as.factor(h_HippAntL_neg))) + 
  #shape = as.factor(h_HippAntL_neg), 
  #p1 <- ggplot(em1, aes(x=last_outcome, y=rt_lag_sc.trend, linetype = as.factor(pe_f2_hipp), ymin=asymp.LCL, ymax=asymp.UCL)) + 
  geom_point(position = position_dodge(width = .6), size=2.5) + 
  #geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), position = position_dodge2(width=0.9), width=0.5) + 
  geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + 
  #geom_crossbar(position=position_dodge(width=0.7), width=0.6) + 
  theme_bw(base_size=12) + facet_wrap(~study)+  ylab("Convergence on\nbest RT (AU)") +
  #scale_shape_manual(values=c(19,23), labels = c("10th %ile", "90th %ile")) + 
  #scale_color_brewer("AH global max\nresponse", palette="Dark2", labels = c("10th %ile", "90th %ile")) +
  scale_color_manual("AH global max\nresponse", values=c("#403202", "#e2b407"), labels = c("10th %ile", "90th %ile")) +
  labs(shape = "AH global max\nresponse") +
  theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) + 
  scale_x_discrete(name ="Trial", labels=c("-0.7" = "5", "0.44" = "50")) + scale_y_continuous(limits = c(.05, .3))


setwd("~/OneDrive/collected_letters/papers/sceptic_fmri/hippo/figs/")

#ggarrange(p1,p2,p3,p4)
vv <- cowplot::plot_grid(p1, p2, p3, p4, nrow=2, align="hv")
ggsave("beta_PH_AH_behavior.pdf", vv, height = 4.5, width = 9, useDingbats=FALSE)

#ggsave("beta_PH_AH_behavior_1.pdf", p1, height=3, width=5, useDingbats=FALSE)


# PH PE cluster suppresses the win-stay-lose-switch behaviors
vs2 <- lmer(v_chosen ~ (trial_neg_inv_sc + last_outcome + 
                          h_HippAntL_neg + pe_f2_hipp)^2 + 
              trial_neg_inv_sc*rewFunc + (1|id), df)
screen.lmerTest(vs2, .01)
vemp <- as_tibble(emmeans(vs2, ~last_outcome | pe_f2_hipp, at = list(pe_f2_hipp = c(-2,2)))) %>% mutate(`Chosen value` = emmean)
p1 <- ggplot(vemp, aes(last_outcome, `Chosen value`, color = pe_f2_hipp, group = pe_f2_hipp)) + 
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), position=position_dodge(width=0.5)) + geom_line(position=position_dodge(width=0.5)) + geom_point(position=position_dodge(width=0.5))
vema <- as_tibble(emmeans(vs2, ~last_outcome | h_HippAntL_neg, at = list(h_HippAntL_neg = c(-2,2)))) %>% mutate(`Chosen value` = emmean)
p2 <- ggplot(vema, aes(last_outcome, `Chosen value`, color = h_HippAntL_neg, group = h_HippAntL_neg)) + 
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),position=position_dodge(width=0.5)) + geom_line(position=position_dodge(width=0.5)) + geom_point(position=position_dodge(width=0.5))
setwd('/plots')
pdf('ah_ph_v_chosen_outcome.pdf', height = 3, width = 7)
ggarrange(p1,p2, ncol = 2)
dev.off()


#################
# Sensitivity analysis for the supplement

#####
# Stargazer tables with covariates

# Main analysis
# NB: stargazer only works with lme4, not lmerTest
mb3hpe_hipp_lme4 <-  lme4::lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + rt_vmax_lag_sc + last_outcome + 
                                               v_max_wi_lag + v_entropy_wi + h_HippAntL_neg +  pe_f2_hipp)^2 + 
                                  rt_lag_sc:last_outcome:h_HippAntL_neg + 
                                  rt_lag_sc:last_outcome:pe_f2_hipp +
                                  rt_vmax_lag_sc:trial_neg_inv_sc:h_HippAntL_neg + 
                                  rt_vmax_lag_sc:trial_neg_inv_sc:pe_f2_hipp  +
                                  (1|id/run), df)
mmb3hpe_hipp_lme4 <-  lme4::lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + rt_vmax_lag_sc + last_outcome + 
                                                v_max_wi_lag + v_entropy_wi +h_HippAntL_neg +  pe_f2_hipp)^2 + 
                                   rt_lag_sc:last_outcome:h_HippAntL_neg + 
                                   rt_lag_sc:last_outcome:pe_f2_hipp +
                                   rt_vmax_lag_sc:trial_neg_inv_sc:h_HippAntL_neg + 
                                   rt_vmax_lag_sc:trial_neg_inv_sc:pe_f2_hipp  +
                                   (1|id/run), mdf)
setwd(file.path(clock_folder, 'fmri/keuka_brain_behavior_analyses/tables/'))

stargazer(mb3hpe_hipp_lme4, mmb3hpe_hipp_lme4, type="html", out="hippo_mb_lab.htm", report = "vcs*",
          digits = 1,single.row=TRUE,omit.stat = "bic",
          dep.var.labels = "RT, scaled",
          covariate.labels = c("-1/trial", "RT(t-1)", "RT(Vmax, t-1)", "last outcome: omission vs. reward", "Vmax, within-subject", "entropy, within-subject", "AH low entropy resp.", "PH RPE resp.",  # main effects
                               "-1/trial * RT(t-1)", "-1/trial * last outcome", "-1/trial * RT(Vmax)", "-1/trial * Vmax", "-1/trial * entropy", "-1/trial * AH", "-1/trial * PH", # 2-way
                               "RT(t-1) * RT(Vmax)", "RT(t-1) * last outcome", "RT(t-1) * Vmax", "RT(t-1) * entropy", "RT(t-1) * AH", "RT(t-1) * PH",
                               "RT(Vmax) * last outcome", "RT(Vmax) * Vmax", "RT(Vmax) * entropy", "RT(Vmax) * AH", "RT(Vmax) * PH",
                               "last outcome * Vmax", "last outcome * entropy", "last outcome * AH", "last outcome * PH",
                               "Vmax * entropy", "Vmax * AH", "Vmax * PH", 
                               "Entropy * AH", "Entropy * PH",
                               "AH * PH",
                               "RT(t-1) * last outcome * AH", "RT(t-1) * last outcome * PH",
                               "-1/trial * RT(Vmax) * AH","-1/trial * RT(Vmax) * PH" ), 
          column.labels = c("fMRI session", "Out-of-session replication"),
          star.char = c("*", "**", "***"),
          star.cutoffs = c(0.05, 0.01, 0.001),
          notes = c("* p<0.05; ** p<0.01; *** p<0.001"),
          notes.append = F)

# Sensitivity analyses for fMRI sample
# add trial and contingency
summary(mb4hpe_hipp_lme4 <-  lme4::lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + rt_vmax_lag_sc + last_outcome + 
                                                       v_max_wi_lag + v_entropy_wi + h_HippAntL_neg +  pe_f2_hipp)^2 + 
                                          rt_lag_sc:last_outcome:h_HippAntL_neg + 
                                          rt_lag_sc:last_outcome:pe_f2_hipp +
                                          rt_vmax_lag_sc:trial_neg_inv_sc:h_HippAntL_neg + 
                                          rt_vmax_lag_sc:trial_neg_inv_sc:pe_f2_hipp  + 
                                          trial_neg_inv_sc*rewFunc*h_HippAntL_neg +
                                          trial_neg_inv_sc*rewFunc*pe_f2_hipp +
                                          (1|id/run), df))

# add uncertainty of last choice
df$u_chosen_lag_sc <- scale(df$u_chosen_lag)
summary(mb5hpe_hipp_lme4 <-  lme4::lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + rt_vmax_lag_sc + last_outcome + 
                                                       v_max_wi_lag + v_entropy_wi + h_HippAntL_neg +  pe_f2_hipp)^2 + 
                                          rt_lag_sc:last_outcome:h_HippAntL_neg + 
                                          rt_lag_sc:last_outcome:pe_f2_hipp +
                                          rt_vmax_lag_sc:trial_neg_inv_sc:h_HippAntL_neg + 
                                          rt_vmax_lag_sc:trial_neg_inv_sc:pe_f2_hipp  + 
                                          trial_neg_inv_sc*rewFunc*h_HippAntL_neg +
                                          trial_neg_inv_sc*rewFunc*pe_f2_hipp +
                                          u_chosen_lag_sc*h_HippAntL_neg +
                                          u_chosen_lag_sc*pe_f2_hipp +
                                          (1|id/run), df))

# add subject-level performance
summary(mb6hpe_hipp_lme4 <-  lme4::lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + rt_vmax_lag_sc + last_outcome + 
                                                       v_max_wi_lag + v_entropy_wi + h_HippAntL_neg +  pe_f2_hipp)^2 + 
                                          rt_lag_sc:last_outcome:h_HippAntL_neg + 
                                          rt_lag_sc:last_outcome:pe_f2_hipp +
                                          rt_vmax_lag_sc:trial_neg_inv_sc:h_HippAntL_neg + 
                                          rt_vmax_lag_sc:trial_neg_inv_sc:pe_f2_hipp  + 
                                          trial_neg_inv_sc*rewFunc*h_HippAntL_neg +
                                          trial_neg_inv_sc*rewFunc*pe_f2_hipp +
                                          u_chosen_lag_sc*h_HippAntL_neg +
                                          u_chosen_lag_sc*pe_f2_hipp +
                                          v_entropy_b*h_HippAntL_neg + v_max_b*h_HippAntL_neg +
                                          v_entropy_b*pe_f2_hipp + v_max_b*pe_f2_hipp +
                                          (1|id/run), df))

stargazer(mb3hpe_hipp_lme4, mb4hpe_hipp_lme4, mb5hpe_hipp_lme4,mb6hpe_hipp_lme4, type="html", out="hippo_mb_sens_lab.htm", report = "vcs*",
          digits = 1,single.row=TRUE,omit.stat = c("bic", "LL"),
          dep.var.labels = "RT, scaled",
          column.labels = c("Main analysis", "+ contingency", "+ choice uncertainty", "+ subject-level performance"),
          covariate.labels = c("-1/trial", "RT(t-1)", "RT(Vmax, t-1)", "last outcome: omission vs. reward", "Vmax, within-subject", "entropy, within-subject", "AH low entropy resp.", "PH RPE resp.",  
                               "contingency: CEVR vs. CEV", "contingency: DEV vs. CEV", "contingency: IEV vs. CEV", "uncertainty or last choice", "mean entropy, between-subjects", "mean Vmax, between-subjects",  # main effects
                               "-1/trial * RT(t-1)", "-1/trial * last outcome", "-1/trial * RT(Vmax)", "-1/trial * Vmax", "-1/trial * entropy", "-1/trial * AH", "-1/trial * PH", # 2-way
                               "RT(t-1) * RT(Vmax)", "RT(t-1) * last outcome", "RT(t-1) * Vmax", "RT(t-1) * entropy", "RT(t-1) * AH", "RT(t-1) * PH",
                               "RT(Vmax) * last outcome", "RT(Vmax) * Vmax", "RT(Vmax) * entropy", "RT(Vmax) * AH", "RT(Vmax) * PH",
                               "last outcome * Vmax", "last outcome * entropy", "last outcome * AH", "last outcome * PH",
                               "Vmax * entropy", "Vmax * AH", "Vmax * PH", 
                               "Entropy * AH", "Entropy * PH",
                               "AH * PH",
                               "-1/trial * CEVR","-1/trial * DEV","-1/trial * IEV",
                               "AH * CEVR","AH * DEV","AH * IEV",
                               "PH * CEVR","PH * DEV","PH * IEV",
                               "AH * uncertainty","PH * uncertainty",
                               "AH * mean entropy","AH * mean Vmax",
                               "PH * mean entropy","PH * mean Vmax",
                               "RT(t-1) * last outcome * AH", "RT(t-1) * last outcome * PH",
                               "-1/trial * RT(Vmax) * AH","-1/trial * RT(Vmax) * PH" ,
                               "-1/trial * CEVR * AH","-1/trial * DEV * AH","-1/trial * IEV * AH",
                               "-1/trial * CEVR * PH","-1/trial * DEV * PH","-1/trial * IEV * PH"),
          star.char = c("*", "**", "***"),
          star.cutoffs = c(0.05, 0.01, 0.001),
          notes = c("* p<0.05; ** p<0.01; *** p<0.001"),
          notes.append = F)

# Sensitivity analyses for MEG sample
# add trial and contingency
summary(mmb4hpe_hipp_lme4 <-  lme4::lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + rt_vmax_lag_sc + last_outcome + 
                                                        v_max_wi_lag + v_entropy_wi + h_HippAntL_neg +  pe_f2_hipp)^2 + 
                                           rt_lag_sc:last_outcome:h_HippAntL_neg + 
                                           rt_lag_sc:last_outcome:pe_f2_hipp +
                                           rt_vmax_lag_sc:trial_neg_inv_sc:h_HippAntL_neg + 
                                           rt_vmax_lag_sc:trial_neg_inv_sc:pe_f2_hipp  + 
                                           trial_neg_inv_sc*rewFunc*h_HippAntL_neg +
                                           trial_neg_inv_sc*rewFunc*pe_f2_hipp +
                                           (1|id/run),mdf))

# add uncertainty of last choice
mdf$u_chosen_lag_sc <- scale(mdf$u_chosen_lag)
summary(mmb5hpe_hipp_lme4 <-  lme4::lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + rt_vmax_lag_sc + last_outcome + 
                                                        v_max_wi_lag + v_entropy_wi + h_HippAntL_neg +  pe_f2_hipp)^2 + 
                                           rt_lag_sc:last_outcome:h_HippAntL_neg + 
                                           rt_lag_sc:last_outcome:pe_f2_hipp +
                                           rt_vmax_lag_sc:trial_neg_inv_sc:h_HippAntL_neg + 
                                           rt_vmax_lag_sc:trial_neg_inv_sc:pe_f2_hipp  + 
                                           trial_neg_inv_sc*rewFunc*h_HippAntL_neg +
                                           trial_neg_inv_sc*rewFunc*pe_f2_hipp +
                                           u_chosen_lag_sc*h_HippAntL_neg +
                                           u_chosen_lag_sc*pe_f2_hipp +
                                           (1|id/run),mdf))

# add subject-level performance
summary(mmb6hpe_hipp_lme4 <-  lme4::lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + rt_vmax_lag_sc + last_outcome + 
                                                        v_max_wi_lag + v_entropy_wi + h_HippAntL_neg +  pe_f2_hipp)^2 + 
                                           rt_lag_sc:last_outcome:h_HippAntL_neg + 
                                           rt_lag_sc:last_outcome:pe_f2_hipp +
                                           rt_vmax_lag_sc:trial_neg_inv_sc:h_HippAntL_neg + 
                                           rt_vmax_lag_sc:trial_neg_inv_sc:pe_f2_hipp  + 
                                           trial_neg_inv_sc*rewFunc*h_HippAntL_neg +
                                           trial_neg_inv_sc*rewFunc*pe_f2_hipp +
                                           u_chosen_lag_sc*h_HippAntL_neg +
                                           u_chosen_lag_sc*pe_f2_hipp +
                                           v_entropy_b*h_HippAntL_neg + v_max_b*h_HippAntL_neg +
                                           v_entropy_b*pe_f2_hipp + v_max_b*pe_f2_hipp +
                                           (1|id/run),mdf))

stargazer(mmb3hpe_hipp_lme4, mmb4hpe_hipp_lme4, mmb5hpe_hipp_lme4,mmb6hpe_hipp_lme4, type="html", out="hippo_mmb_sens_lab.htm", report = "vcs*",
          digits = 1,single.row=TRUE,omit.stat = c("bic", "LL"),
          dep.var.labels = "RT, scaled",
          column.labels = c("Main analysis", "+ contingency", "+ choice uncertainty", "+ subject-level performance"),
          covariate.labels = c("-1/trial", "RT(t-1)", "RT(Vmax, t-1)", "last outcome: omission vs. reward", "Vmax, within-subject", "entropy, within-subject", "AH low entropy resp.", "PH RPE resp.",  
                               "contingency: CEVR vs. CEV", "contingency: DEV vs. CEV", "contingency: IEV vs. CEV", "uncertainty or last choice", "mean entropy, between-subjects", "mean Vmax, between-subjects",  # main effects
                               "-1/trial * RT(t-1)", "-1/trial * last outcome", "-1/trial * RT(Vmax)", "-1/trial * Vmax", "-1/trial * entropy", "-1/trial * AH", "-1/trial * PH", # 2-way
                               "RT(t-1) * RT(Vmax)", "RT(t-1) * last outcome", "RT(t-1) * Vmax", "RT(t-1) * entropy", "RT(t-1) * AH", "RT(t-1) * PH",
                               "RT(Vmax) * last outcome", "RT(Vmax) * Vmax", "RT(Vmax) * entropy", "RT(Vmax) * AH", "RT(Vmax) * PH",
                               "last outcome * Vmax", "last outcome * entropy", "last outcome * AH", "last outcome * PH",
                               "Vmax * entropy", "Vmax * AH", "Vmax * PH", 
                               "Entropy * AH", "Entropy * PH",
                               "AH * PH",
                               "-1/trial * CEVR","-1/trial * DEV","-1/trial * IEV",
                               "AH * CEVR","AH * DEV","AH * IEV",
                               "PH * CEVR","PH * DEV","PH * IEV",
                               "AH * uncertainty","PH * uncertainty",
                               "AH * mean entropy","AH * mean Vmax",
                               "PH * mean entropy","PH * mean Vmax",
                               "RT(t-1) * last outcome * AH", "RT(t-1) * last outcome * PH",
                               "-1/trial * RT(Vmax) * AH","-1/trial * RT(Vmax) * PH" ,
                               "-1/trial * CEVR * AH","-1/trial * DEV * AH","-1/trial * IEV * AH",
                               "-1/trial * CEVR * PH","-1/trial * DEV * PH","-1/trial * IEV * PH"),
          star.char = c("*", "**", "***"),
          star.cutoffs = c(0.05, 0.01, 0.001),
          notes = c("* p<0.05; ** p<0.01; *** p<0.001"),
          notes.append = F)


# Model-free (mf) RT analyses -- behavioral variables
mf1 <- lmer(rt_csv ~ (trial_neg_inv_sc + rt_lag_sc + last_outcome )^2  + (1|id/run), df)
summary(mf1)


# ruling out the exploitation account of the PH-guided shifts

# v_chosen change -- the swings in high-PH subjects are toward lower-valued options
mv1 <- lmer(v_chosen_quantile_change ~ (trial_neg_inv_sc + rt_lag_sc + h_HippAntL_neg + rewFunc + last_outcome)^3 +
              (trial_neg_inv_sc + rt_lag_sc + pe_f2_hipp + rewFunc + last_outcome)^3 + v_chosen_quantile_lag + (1|id/run), df)
screen.lmerTest(mv1, .05)
anova(mv1)
em <- as_tibble(emmeans(mv1, ~rt_lag_sc|pe_f2_hipp*rewFunc, at = list(rt_lag_sc = c(-2,0,2), pe_f2_hipp = c(-2,0,2))))
em$`Chosen value quantile change` <- em$emmean
ggplot(em, aes(rt_lag_sc, `Chosen value quantile change`, color = pe_f2_hipp, lty = rewFunc, group = pe_f2_hipp)) + geom_point() + geom_line()
# add reward -- this model is better by a huge margin
df_tmp <- df %>% filter(rewFunc=='IEV' | rewFunc=='DEV') %>% droplevels()
mv2 <- lmer(v_chosen_quantile_change ~ (trial_neg_inv_sc + rt_lag_sc + h_HippAntL_neg + rewFunc + last_outcome)^3 +
              (trial_neg_inv_sc + rt_lag_sc + pe_f2_hipp + rewFunc + last_outcome)^3 + v_chosen_quantile_lag + (1|id/run), df_tmp)
screen.lmerTest(mv2, .05)
Anova(mv2, '3')
# anova(mv1,mv2)
em <- as_tibble(emmeans(mv2, ~rt_lag_sc|pe_f2_hipp*rewFunc, at = list(rt_lag_sc = c(-2,0,2), pe_f2_hipp = c(-2,0,2))))
em$`Chosen value quantile change` <- em$emmean
ggplot(em, aes(rt_lag_sc, `Chosen value quantile change`, color = pe_f2_hipp, group = pe_f2_hipp)) + geom_point() + facet_wrap(~rewFunc)
mrt2 <- lmer(rt_csv ~ (trial_neg_inv_sc + rt_lag_sc + h_HippAntL_neg + rewFunc + last_outcome)^3 +
               (trial_neg_inv_sc + rt_lag_sc + pe_f2_hipp + rewFunc + last_outcome)^3 + v_chosen_quantile_lag + (1|id/run), df)
screen.lmerTest(mrt2, .05)
em <- as_tibble(emtrends(mrt2, var = "rt_lag_sc", specs = c("rewFunc", "pe_f2_hipp"), at = list( pe_f2_hipp = c(-2,0,2))))
ggplot(em, aes(pe_f2_hipp, rt_lag_sc.trend, color = rewFunc)) + geom_point()

ggplot(df, aes(run_trial, rt_csv, color = rewFunc, lty =  pe_f2_hipp_resp)) + geom_smooth()
# monotonicity of reward function -- no rew*rewFunc*hipp
mf2 <-  lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + last_outcome +
                            h_HippAntL_neg + pe_f2_hipp)^2 +
               rt_lag_sc:last_outcome:h_HippAntL_neg +
               rt_lag_sc:last_outcome:pe_f2_hipp + trial_neg_inv_sc*rewFunc +
               rewFunc:last_outcome:h_HippAntL_neg + rewFunc:last_outcome:pe_f2_hipp + (1|id/run), df)
screen.lmerTest(mf2)
ggplot(df, aes(rt_lag, rt_csv, color = pe_f2_hipp_resp)) + geom_smooth(method = 'gam', formula = y~splines::ns(x,3)) + facet_wrap(~last_outcome)
ggplot(mdf, aes(rt_lag, rt_csv, color = pe_f2_hipp_resp)) + geom_smooth(method = 'gam', formula = y~splines::ns(x,4))

# v_chosen
mf2a <-  lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + last_outcome +
                             h_HippAntL_neg + pe_f2_hipp)^2 +
                rt_lag_sc:last_outcome:h_HippAntL_neg +
                rt_lag_sc:last_outcome:pe_f2_hipp + trial_neg_inv_sc*rewFunc +
                rewFunc:last_outcome:h_HippAntL_neg + rewFunc:last_outcome:pe_f2_hipp + (1|id/run), df)
screen.lmerTest(mf2)


##
## MEG data for out-of-session replication
mmf3hpe <- lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + last_outcome + 
                               h_HippAntL_neg + pe_f2_hipp)^2 + 
                  rt_lag_sc:last_outcome:h_HippAntL_neg + 
                  rt_lag_sc:last_outcome:pe_f2_hipp + trial_neg_inv_sc*rewFunc + (1|id/run), mdf)
screen.lmerTest(mmf3hpe)


##################
# model-free sensitivity analyses
mf3hpe <-  lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + last_outcome + 
                               h_HippAntL_neg + pe_f2_hipp)^2 + 
                  rt_lag_sc:last_outcome:h_HippAntL_neg + 
                  rt_lag_sc:last_outcome:pe_f2_hipp + trial_neg_inv_sc*rewFunc + (1|id/run), df)
# summary(mf3hpe)
screen.lmerTest(mf3hpe)
summary(mf3hpe)
##
## MEG data for out-of-session replication
mmf3hpe <- lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + last_outcome + 
                               h_HippAntL_neg + pe_f2_hipp)^2 + 
                  rt_lag_sc:last_outcome:h_HippAntL_neg + 
                  rt_lag_sc:last_outcome:pe_f2_hipp + trial_neg_inv_sc*rewFunc + (1|id/run), mdf)
screen.lmerTest(mmf3hpe)
summary(mmf3hpe)
Anova(mmf3hpe,'3')
summary(mmf3hpe)

# reduced model without last outcome

mf2hpe <-  lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc +  
                               h_HippAntL_neg + pe_f2_hipp)^2 + 
                  rt_lag_sc:h_HippAntL_neg + 
                  rt_lag_sc:pe_f2_hipp + trial_neg_inv_sc*rewFunc + (1|id/run), df)
# summary(mf3hpe)
screen.lmerTest(mf2hpe)
summary(mf2hpe)
##
## MEG data for out-of-session replication
mmf2hpe <- lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc +  
                               h_HippAntL_neg + pe_f2_hipp)^2 + 
                  rt_lag_sc:h_HippAntL_neg + 
                  rt_lag_sc:pe_f2_hipp + trial_neg_inv_sc*rewFunc + (1|id/run), mdf)
screen.lmerTest(mmf2hpe)
summary(mmf2hpe)
Anova(mmf2hpe,'3')
summary(mmf2hpe)


# include RTvmax

mf4hpe <-  lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + last_outcome + 
                               h_HippAntL_neg + pe_f2_hipp)^2 + 
                  rt_lag_sc:last_outcome:h_HippAntL_neg + 
                  rt_lag_sc:last_outcome:pe_f2_hipp + 
                  rt_vmax_lag_sc*trial_neg_inv_sc*h_HippAntL_neg + 
                  rt_vmax_lag_sc*trial_neg_inv_sc*pe_f2_hipp  + 
                  trial_neg_inv_sc*rewFunc + (1|id/run), df)
screen.lmerTest(mf4hpe)
anova(mf4hpe)
##
## MEG data for out-of-session replication
mmf4hpe <-  lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + last_outcome + 
                                h_HippAntL_neg + pe_f2_hipp)^2 + 
                   rt_lag_sc:last_outcome:h_HippAntL_neg + 
                   rt_lag_sc:last_outcome:pe_f2_hipp + 
                   rt_vmax_lag_sc*trial_neg_inv_sc*h_HippAntL_neg + 
                   rt_vmax_lag_sc*trial_neg_inv_sc*pe_f2_hipp  + 
                   trial_neg_inv_sc*rewFunc + (1|id/run), mdf)
summary(mmf4hpe)
anova(mf4hpe)
Anova(mmf4hpe,'3')

# include RTvmax without outcome

mf2ahpe <-  lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + 
                                h_HippAntL_neg + pe_f2_hipp)^2 + 
                   rt_vmax_lag_sc*trial_neg_inv_sc*h_HippAntL_neg + 
                   rt_vmax_lag_sc*trial_neg_inv_sc*pe_f2_hipp  + 
                   (1|id/run), df)
summary(mf2ahpe)
anova(mf2ahpe)
mmf2ahpe <-  lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + 
                                 h_HippAntL_neg + pe_f2_hipp)^2 + 
                    rt_vmax_lag_sc*trial_neg_inv_sc*h_HippAntL_neg + 
                    rt_vmax_lag_sc*trial_neg_inv_sc*pe_f2_hipp  + 
                    (1|id/run), mdf)
summary(mmf2ahpe)
anova(mmf2ahpe)


## ascertain replication -- visual check
# plot_models(mmf3hpe,mf3hpe)
# ascertain replication overall 
p1 <- plot_model(mf3hpe, show.values = T)
p2 <- plot_model(mmf3hpe, show.values = T)
pdf("model_free_beta_replication.pdf", height = 6, width = 12)
ggarrange(p1,p2,ncol = 2, labels  = c("fMRI", "MEG"))
dev.off()

# save model statistics for supplement
mf3hpe_lme4 <-  lme4::lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + last_outcome + 
                                          h_HippAntL_neg + pe_f2_hipp)^2 + 
                             rt_lag_sc:last_outcome:h_HippAntL_neg + 
                             rt_lag_sc:last_outcome:pe_f2_hipp + trial_neg_inv_sc*rewFunc + (1|id/run), df)
# summary(mf3hpe)
screen.lmerTest(mf3hpe)

###############
# compare performance in fMRI and MEG sessions: better in MEG, particularly in IEV
p1 <- ggplot(df, aes(run_trial, rt_csv, color = rewFunc)) + geom_smooth(method = "gam", formula = y~splines::ns(x,4)) +  coord_cartesian(ylim=c(1300,2300))
p2 <- ggplot(mdf, aes(run_trial, rt_csv, color = rewFunc)) + geom_smooth(method = "gam", formula = y~splines::ns(x,4)) +  coord_cartesian(ylim=c(1300,2300))
ggarrange(p1,p2, labels = c("fMRI", "MEG"), align = "hv")
## MEG data for out-of-session replication
#
mmf3hpe_lme4 <-  lme4::lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + last_outcome + 
                                           h_HippAntL_neg + pe_f2_hipp)^2 + 
                              rt_lag_sc:last_outcome:h_HippAntL_neg + 
                              rt_lag_sc:last_outcome:pe_f2_hipp + trial_neg_inv_sc*rewFunc + (1|id/run), mdf)

setwd(file.path(clock_folder, 'fmri/keuka_brain_behavior_analyses/tables/'))
# NB: stargazer only runs with lmer, not lmerTest objects
stargazer(mf3hpe_lme4, mmf3hpe_lme4, type="html", out="hippo_mf.htm", report = "vcs*",
          digits = 1,single.row=TRUE,omit.stat = "bic",
          column.labels = c("fMRI session", "Out-of-session replication"),
          star.char = c("*", "**", "***"),
          star.cutoffs = c(0.05, 0.01, 0.001),
          notes = c("* p<0.05; ** p<0.01; *** p<0.001"),
          notes.append = F)

##############
# Sensitivity analyses (cont.):
# without covariates: effects are unchanged
summary(m0 <-  lme4::lmer(rt_csv_sc ~  rt_lag_sc*last_outcome*pe_f2_hipp +
                            rt_vmax_lag_sc*trial_neg_inv_sc*h_HippAntL_neg + 
                            (1|id/run), df))
summary(mm0 <-  lme4::lmer(rt_csv_sc ~  rt_lag_sc*last_outcome*pe_f2_hipp +
                             rt_vmax_lag_sc*trial_neg_inv_sc*h_HippAntL_neg + 
                             (1|id/run), mdf))
setwd(file.path(clock_folder, 'fmri/keuka_brain_behavior_analyses/tables/'))
stargazer(m0, mm0, type="html", out="hippo_mb_no_covariates.htm", report = "vcs*",
          digits = 1,single.row=TRUE,omit.stat = "bic",
          column.labels = c("fMRI session", "Out-of-session replication"),
          star.char = c("*", "**", "***"),
          star.cutoffs = c(0.05, 0.01, 0.001),
          notes = c("* p<0.05; ** p<0.01; *** p<0.001"),
          notes.append = F)


##############
# sensitivity analyses (cont.): controlling for uncertainty
mf3hpe_u <-  lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + last_outcome + scale(u_chosen_lag) +
                                 h_HippAntL_neg + pe_f2_hipp)^2 + 
                    rt_lag_sc:last_outcome:h_HippAntL_neg + 
                    rt_lag_sc:last_outcome:pe_f2_hipp + trial_neg_inv_sc*rewFunc + (1|id/run), df)
while (any(grepl("failed to converge", mf3hpe_u@optinfo$conv$lme4$messages) )) {
  print(mf3hpe_u@optinfo$conv$lme4$conv)
  ss <- getME(mf3hpe_u,c("theta","fixef"))
  mf3hpe_u <- update(mf3hpe_u, start=ss)}

# summary(mf3hpe)
screen.lmerTest(mf3hpe_u)
##
## MEG data for out-of-session replication
mmf3hpe_u <-  lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + last_outcome + scale(u_chosen_lag) +
                                  h_HippAntL_neg + pe_f2_hipp)^2 + 
                     rt_lag_sc:last_outcome:h_HippAntL_neg + 
                     rt_lag_sc:last_outcome:pe_f2_hipp + trial_neg_inv_sc*rewFunc + (1|id/run), mdf)
while (any(grepl("failed to converge", mmf3hpe_u@optinfo$conv$lme4$messages) )) {
  print(mmf3hpe_u@optinfo$conv$lme4$conv)
  ss <- getME(mmf3hpe_u,c("theta","fixef"))
  mmf3hpe_u <- update(mmf3hpe_u, start=ss)}

# summary(mmf3hpe)
screen.lmerTest(mmf3hpe_u)

##############
# Sensitivity analyses (cont.):
# R vs. L PH 
#
mb3hpe_hipp_rl <-  lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + rt_vmax_lag_sc + last_outcome + 
                                       v_max_wi_lag + v_entropy_wi + h_HippAntL_neg +  pe_PH)^2 + 
                          rt_lag_sc:last_outcome:h_HippAntL_neg + 
                          rt_lag_sc:last_outcome:pe_PH +
                          rt_vmax_lag_sc:trial_neg_inv_sc:h_HippAntL_neg + 
                          rt_vmax_lag_sc:trial_neg_inv_sc:pe_PH  +
                          (1|id/run), df)
summary(mb3hpe_hipp_rl)
screen.lmerTest(mb3hpe_hipp_rl, .05)
# Anova(mmb3hpe_hipp, '3')


mb3hpe_hipp_rl <-  lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + rt_vmax_lag_sc + last_outcome + 
                                       v_max_wi_lag + v_entropy_wi + h_HippAntL_neg +  pe_PH)^2 + 
                          rt_lag_sc:last_outcome:h_HippAntL_neg + 
                          rt_lag_sc:last_outcome:pe_PH +
                          rt_vmax_lag_sc:trial_neg_inv_sc:h_HippAntL_neg + 
                          rt_vmax_lag_sc:trial_neg_inv_sc:pe_PH  +
                          (1|id/run), df)
# summary(mb3hpe_hipp)
screen.lmerTest(mb3hpe_hipp_rl, .05)
# Anova(mmb3hpe_hipp, '3')

mmb3hpe_hipp_rl <-  lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + rt_vmax_lag_sc + last_outcome + 
                                        v_max_wi_lag + v_entropy_wi +h_HippAntL_neg +  pe_PH)^2 + 
                           rt_lag_sc:last_outcome:h_HippAntL_neg + 
                           rt_lag_sc:last_outcome:pe_PH +
                           rt_vmax_lag_sc:trial_neg_inv_sc:h_HippAntL_neg + 
                           rt_vmax_lag_sc:trial_neg_inv_sc:pe_PH  +
                           (1|id/run), mdf)
while (any(grepl("failed to converge", mmb3hpe_hipp_rl@optinfo$conv$lme4$messages) )) {
  print(mmb3hpe_hipp_rl@optinfo$conv$lme4$conv)
  ss <- getME(mmb3hpe_hipp_rl,c("theta","fixef"))
  mmb3hpe_hipp_rl <- update(mmb3hpe_hipp_rl, start=ss)}
screen.lmerTest(mmb3hpe_hipp_rl, .05)
summary(mmb3hpe_hipp_rl)
Anova(mmb3hpe_hipp_rl, '3')



screen.lmerTest(mmb3hpe_hipp_rl, .05)
summary(mmb3hpe_hipp_rl)
Anova(mmb3hpe_hipp_rl, '3')

mb3hpe_hipp_r <-  lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + rt_vmax_lag_sc + last_outcome + 
                                      v_max_wi_lag + v_entropy_wi + h_HippAntL_neg +  pe_PH_r)^2 + 
                         rt_lag_sc:last_outcome:h_HippAntL_neg + 
                         rt_lag_sc:last_outcome:pe_PH_r +
                         rt_vmax_lag_sc:trial_neg_inv_sc:h_HippAntL_neg + 
                         rt_vmax_lag_sc:trial_neg_inv_sc:pe_PH_r  +
                         (1|id/run), df)
# summary(mb3hpe_hipp)
screen.lmerTest(mb3hpe_hipp_r, .05)
# Anova(mmb3hpe_hipp, '3')

mmb3hpe_hipp_r <-  lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + rt_vmax_lag_sc + last_outcome + 
                                       v_max_wi_lag + v_entropy_wi +h_HippAntL_neg +  pe_PH_r)^2 + 
                          rt_lag_sc:last_outcome:h_HippAntL_neg + 
                          rt_lag_sc:last_outcome:pe_PH_r +
                          rt_vmax_lag_sc:trial_neg_inv_sc:h_HippAntL_neg + 
                          rt_vmax_lag_sc:trial_neg_inv_sc:pe_PH_r  +
                          (1|id/run), mdf)
screen.lmerTest(mmb3hpe_hipp_r, .05)
summary(mmb3hpe_hipp_r)
Anova(mmb3hpe_hipp_r, '3')

# bump refinement vs. more global search -- add rt_lag*rt_vmax interactions
mb4hpe_hipp_rl <-  lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + rt_vmax_lag_sc + last_outcome + 
                                       v_max_wi_lag + v_entropy_wi + h_HippAntL_neg +  pe_PH)^2 + 
                          rt_lag_sc:last_outcome:h_HippAntL_neg + 
                          rt_lag_sc:last_outcome:pe_PH +
                          rt_vmax_lag_sc:trial_neg_inv_sc:h_HippAntL_neg + 
                          rt_vmax_lag_sc:trial_neg_inv_sc:pe_PH  +
                          (1|id/run), df)
# summary(mb3hpe_hipp)
screen.lmerTest(mb4hpe_hipp_rl, .05)
# Anova(mmb3hpe_hipp, '3')

mmb4hpe_hipp_rl <-  lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + rt_vmax_lag_sc + last_outcome + 
                                        v_max_wi_lag + v_entropy_wi + h_HippAntL_neg +  pe_PH)^2 + 
                           rt_lag_sc:last_outcome:h_HippAntL_neg + 
                           rt_lag_sc:last_outcome:pe_PH +
                           rt_vmax_lag_sc:trial_neg_inv_sc:h_HippAntL_neg + 
                           rt_vmax_lag_sc:trial_neg_inv_sc:pe_PH  +
                           (1|id/run), mdf)
screen.lmerTest(mmb4hpe_hipp_rl, .05)

## PH replication plot for mean of R and left
mterms <- names(fixef(mb3hpe_hipp_rl))
setwd(file.path(clock_folder, 'fmri/keuka_brain_behavior_analyses/plots/'))
ph <- plot_models(mb3hpe_hipp_rl,mmb3hpe_hipp_rl, rm.terms = mterms[c(-22, -39)], m.labels = c("fMRI", "replication"),
                  show.values = T,  std.est = "std2", legend.title = "Session", vline.color = "slategray3",
                  wrap.labels = 15,  axis.labels = c("Previous RT * Omission * Post. hippocampal PE response","Previous RT * Post. hippocampal PE response"),
                  axis.title = "Greater RT swing  <==>  Smaller RT swing")
ph <- ph + ggplot2::ylim(-.1,.1)
pdf("ph_beta_models_replication_mean_LR.pdf", height = 3, width = 5)
ph
dev.off()

# visual sanity checks
# p1 <- ggplot(df, aes(rt_lag, rt_csv, color = pe_f2_hipp_resp, lty = last_outcome)) + geom_smooth(method = 'glm') #+ facet_wrap(~rewFunc)
# p2 <- ggplot(mdf, aes(rt_lag, rt_csv, color = pe_f2_hipp_resp, lty = last_outcome)) + geom_smooth(method = 'glm') #+ facet_wrap(~rewFunc)
# ggarrange(p1,p2, ncol = 1, nrow = 2, labels = c("fMRI", "MEG"))
# 
# p1 <- ggplot(df, aes(rt_vmax_lag, rt_csv, color = pe_f2_hipp_resp)) + geom_smooth(method = 'glm') + facet_wrap(~learning_epoch)
# p2 <- ggplot(mdf, aes(rt_vmax_lag, rt_csv, color = pe_f2_hipp_resp)) + geom_smooth(method = 'glm') + facet_wrap(~learning_epoch)
# ggarrange(p1,p2, ncol = 1, nrow = 2, labels = c("fMRI", "MEG"))

# p1 <- ggplot(df, aes(rt_lag, rt_csv, color = h_HippAntL_resp, lty = last_outcome)) + geom_smooth(method = 'glm') #+ facet_wrap(~learning_epoch)
# p2 <- ggplot(mdf, aes(rt_lag, rt_csv, color = h_HippAntL_resp, lty = last_outcome)) + geom_smooth(method = 'glm') #+ facet_wrap(~learning_epoch)
# ggarrange(p1,p2, ncol = 1, nrow = 2, labels = c("fMRI", "MEG"))
# 
# p1 <- ggplot(df, aes(rt_lag, rt_csv, color = h_HippAntL_resp, lty = last_outcome)) + geom_smooth(method = 'glm') #+ facet_wrap(~learning_epoch)
# p2 <- ggplot(mdf, aes(rt_lag, rt_csv, color = h_HippAntL_resp, lty = last_outcome)) + geom_smooth(method = 'glm') #+ facet_wrap(~learning_epoch)
# ggarrange(p1,p2, ncol = 1, nrow = 2, labels = c("fMRI", "MEG"))
# 
# 
# p1 <- ggplot(df, aes(rt_vmax_lag, rt_csv, color = h_HippAntL_resp)) + geom_smooth(method = 'glm') + facet_wrap(~learning_epoch)
# p2 <- ggplot(mdf, aes(rt_vmax_lag, rt_csv, color = h_HippAntL_resp)) + geom_smooth(method = 'glm') + facet_wrap(~learning_epoch)
# ggarrange(p1,p2, ncol = 1, nrow = 2, labels = c("fMRI", "MEG"))
# 

# # understand rt_vmax_change effect
# ggplot(df, aes(rt_vmax_change, rt_csv, color = pe_f2_hipp_resp)) + geom_smooth(method = "glm")

######################
# Neural specificity
######################

# Add: v_f1_neg_cog"             "v_f2_paralimb"            "h_f1_fp"                  "h_f2_neg_paralimb"       
# "pe_f1_cort_str"           "pe_f2_hipp"               "pe_PH"         
mn1 <-  lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + rt_vmax_lag_sc + last_outcome + 
                            v_max_wi_lag + v_entropy_wi + h_HippAntL_neg + pe_f1_cort_str + pe_f2_hipp + v_f2_paralimb + h_f1_fp)^2 + 
               rt_lag_sc:last_outcome:h_HippAntL_neg + 
               rt_lag_sc:last_outcome:pe_f2_hipp +
               rt_lag_sc:last_outcome:pe_f1_cort_str + 
               rt_lag_sc:last_outcome:v_f2_paralimb + 
               rt_lag_sc:last_outcome:h_f1_fp + 
               rt_vmax_lag_sc:trial_neg_inv_sc:h_HippAntL_neg + 
               rt_vmax_lag_sc:trial_neg_inv_sc:pe_f2_hipp  +
               rt_vmax_lag_sc:trial_neg_inv_sc:pe_f1_cort_str  +
               rt_vmax_lag_sc:trial_neg_inv_sc:v_f2_paralimb  +
               rt_vmax_lag_sc:trial_neg_inv_sc:h_f1_fp  +
               (1|id/run), df)
summary(mn1)
Anova(mn1, '3')
mmn1 <-  lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + rt_vmax_lag_sc + last_outcome + 
                             v_max_wi_lag + v_entropy_wi + h_HippAntL_neg + pe_f1_cort_str + pe_f2_hipp + v_f2_paralimb + h_f1_fp)^2 + 
                rt_lag_sc:last_outcome:h_HippAntL_neg + 
                rt_lag_sc:last_outcome:pe_f2_hipp +
                rt_lag_sc:last_outcome:pe_f1_cort_str + 
                rt_lag_sc:last_outcome:v_f2_paralimb + 
                rt_lag_sc:last_outcome:h_f1_fp + 
                rt_vmax_lag_sc:trial_neg_inv_sc:h_HippAntL_neg + 
                rt_vmax_lag_sc:trial_neg_inv_sc:pe_f2_hipp  +
                rt_vmax_lag_sc:trial_neg_inv_sc:pe_f1_cort_str  +
                rt_vmax_lag_sc:trial_neg_inv_sc:v_f2_paralimb  +
                rt_vmax_lag_sc:trial_neg_inv_sc:h_f1_fp  +
                (1|id/run), mdf)
summary(mmn1)
Anova(mmn1, '3')

# supplementary figure illustrating v_chosen by PH response to prove that PH explorers were rewarded.

ggplot(df, aes(run_trial, score_csv, color = pe_f2_hipp_resp)) + geom_smooth(method = 'gam', formula = y ~ splines::ns(x,4))
setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/hippo/figs/final/supp')
pdf('score_by_PH.pdf', height = 3, width = 3)
ggplot(df , aes(run_trial, score_csv, color = pe_f2_hipp_resp)) + geom_smooth(method = 'loess')
dev.off()
######################
# Uncertainty models #
######################
# archival: these questions are now better addressed in coxme analyses
# Sanity checks
# ggplot(df, aes(run_trial, u_chosen, group = interaction(rewFunc, rt_lag>2000), color = rewFunc, lty = rt_lag>2000)) + geom_smooth()
# ggplot(df, aes(run_trial, u_chosen, group = rewFunc, color = rewFunc)) + geom_smooth()
# ggplot(df, aes(run_trial, u_chosen_change, group = rewFunc, color = rewFunc)) + geom_smooth(method = 'gam', formula = y ~ splines::ns(x, 3))

# Inspect correlations to estimate the uncertainty/value confound -- not huge
vars <- df %>% select(u_chosen, u_chosen_change, u_chosen_quantile, u_chosen_quantile_change, 
                      v_chosen, v_chosen_quantile, v_chosen_quantile_change, 
                      run_trial, magnitude, probability,rt_csv, rt_lag, rt_vmax_lag)
u_cor <- corr.test(vars,method = 'pearson', adjust = 'none')

setwd(file.path(clock_folder, 'fmri/keuka_brain_behavior_analyses/plots'))
pdf("u_corr_reset.pdf", width=12, height=12)
corrplot(u_cor$r, cl.lim=c(-1,1),
         method = "circle", tl.cex = 1.5, type = "upper", tl.col = 'black',
         order = "hclust", diag = FALSE,
         addCoef.col="black", addCoefasPercent = FALSE,
         p.mat = u_cor$p, sig.level=0.05, insig = "blank")
dev.off()

# timecourses of choice uncertainty by HIPP response
# ggplot(df, aes(run_trial, rt_csv, color = rewFunc, lty = h_HippAntL_resp, group = interaction(rewFunc, h_HippAntL_resp))) + geom_smooth(method = 'gam', formula = y ~ splines::ns(x, 3)) #+ facet_wrap(~run)
p1 <- ggplot(df, aes(run_trial, u_chosen, color = rt_lag>2000, lty = pe_f2_hipp_resp, group = interaction(rt_lag>2000, pe_f2_hipp_resp))) + geom_smooth(method = 'gam', formula = y ~ splines::ns(x, 3))  + facet_wrap(last_outcome~rewFunc)
p2 <- ggplot(mdf, aes(run_trial, u_chosen, color = rt_lag>2000, lty = pe_f2_hipp_resp, group = interaction(rt_lag>2000, pe_f2_hipp_resp))) + geom_smooth(method = 'gam', formula = y ~ splines::ns(x, 3))  + facet_wrap(last_outcome~rewFunc)
ggarrange(p1,p2)

# RT timecourses by hipp response
p1 <- ggplot(df, aes(run_trial, rt_csv, color = rewFunc, lty = pe_f2_hipp_resp, group = interaction(rewFunc, pe_f2_hipp_resp))) + geom_smooth(method = 'gam', formula = y ~ splines::ns(x, 4))
p2 <- ggplot(mdf, aes(run_trial, rt_csv, color = rewFunc, lty = pe_f2_hipp_resp, group = interaction(rewFunc, pe_f2_hipp_resp))) + geom_smooth(method = 'gam', formula = y ~ splines::ns(x, 4))
ggarrange(p1,p2)

p1 <- ggplot(df, aes(run_trial, rt_csv, color = rewFunc, lty = h_HippAntL_resp, group = interaction(rewFunc, h_HippAntL_resp))) + geom_smooth(method = 'gam', formula = y ~ splines::ns(x, 4))
p2 <- ggplot(mdf, aes(run_trial, rt_csv, color = rewFunc, lty = h_HippAntL_resp, group = interaction(rewFunc, h_HippAntL_resp))) + geom_smooth(method = 'gam', formula = y ~ splines::ns(x, 4))
ggarrange(p1,p2)

# the differences in timecourse of RT swings between high and low PH are visible in fMRI, subtle in MEG
# p1 <- ggplot(df %>% filter(!is.na(last_outcome)), aes(run_trial, rt_swing, color = rewFunc, lty = pe_f2_hipp_resp, group = interaction(rewFunc, pe_f2_hipp_resp))) + geom_smooth(method = 'loess') + facet_wrap(~last_outcome)
# p2 <- ggplot(mdf %>% filter(!is.na(last_outcome)), aes(run_trial, rt_swing, color = rewFunc, lty = pe_f2_hipp_resp, group = interaction(rewFunc, pe_f2_hipp_resp))) + geom_smooth(method = 'loess')+ facet_wrap(~last_outcome)
p1 <- ggplot(df %>% filter(!is.na(last_outcome) & rewFunc!="CEV" & rewFunc!="CEVR"), aes(run_trial, rt_swing, color = pe_f2_hipp_resp)) + geom_smooth(method = 'loess') 
p2 <- ggplot(mdf %>% filter(!is.na(last_outcome) & rewFunc!="CEV" & rewFunc!="CEVR"), aes(run_trial, rt_swing, color = pe_f2_hipp_resp)) + geom_smooth(method = 'loess') 
pdf('rt_swings_by_PH_resp.pdf', height = 5, width = 10)
ggarrange(p1,p2, labels = c("fMRI", "MEG"))
dev.off()

p1 <- ggplot(df %>% filter(!is.na(last_outcome) & rewFunc!="CEV" & rewFunc!="CEVR"), aes(run_trial, rt_swing, color = h_HippAntL_resp)) + geom_smooth(method = 'loess') 
p2 <- ggplot(mdf %>% filter(!is.na(last_outcome) & rewFunc!="CEV" & rewFunc!="CEVR"), aes(run_trial, rt_swing, color = h_HippAntL_resp)) + geom_smooth(method = 'loess') 
pdf('rt_swings_by_AH_resp.pdf', height = 5, width = 10)
ggarrange(p1,p2, labels = c("fMRI", "MEG"))
dev.off()

##
# Builidng the uncertainty model for interactions with betas
# more sanity checks on quantiles (relative uncertainty) -- looks right
# ggplot(df, aes(run_trial, u_chosen)) + geom_smooth()+ facet_wrap(~rewFunc)
ggplot(df, aes(run_trial, u_chosen_quantile)) + geom_smooth()+ facet_wrap(~rewFunc)
ggplot(df, aes(run_trial, u_chosen_quantile_change)) + geom_smooth()+ facet_wrap(~rewFunc)

# Regress value out of uncertainty for plotting
m <- lmer(u_chosen_quantile ~ v_chosen + (1|ID), df)
df$u_chosen_v_partialed_out <- resid(m)

m <- lmer(u_chosen_quantile ~ v_chosen + (1|ID), mdf)
mdf$u_chosen_v_partialed_out <- resid(m)

#
mdf <- mdf %>% group_by(ID, run) %>%  arrange(run_trial) %>% mutate(v_chosen_lag = lag(v_chosen)) %>% ungroup()
ldf <- df %>% filter(rewFunc=='IEV' | rewFunc=='DEV')
lmdf <- mdf %>% filter(rewFunc=='IEV' | rewFunc=='DEV')

p1 <- ggplot(ldf, aes(run_trial, u_chosen_quantile, color = pe_f2_hipp_resp)) + geom_smooth(method = 'gam', formula = y ~ splines::ns(x,2))+ facet_wrap(~rewFunc)
p2 <- ggplot(ldf, aes(run_trial, u_chosen_v_partialed_out, color = pe_f2_hipp_resp)) + geom_smooth(method = 'gam', formula = y ~ splines::ns(x,2))+ facet_wrap(~rewFunc)
p3 <- ggplot(lmdf, aes(run_trial, u_chosen_quantile, color = pe_f2_hipp_resp)) + geom_smooth(method = 'gam', formula = y ~ splines::ns(x,2))+ facet_wrap(~rewFunc)
p4 <- ggplot(lmdf, aes(run_trial, u_chosen_v_partialed_out, color = pe_f2_hipp_resp)) + geom_smooth(method = 'gam', formula = y ~ splines::ns(x,2))+ facet_wrap(~rewFunc)
setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/hippo/figs/beta_uncertainty/')
pdf('u_chosen_by_PH.pdf', width = 9, height = 6)
ggarrange(p1,p2, p3, p4, ncol = 2, nrow = 2, labels = c("fMRI", "fMRI, value partialed out", "replication", "replication, value partialed out"))
dev.off()
## same for value
p1 <- ggplot(ldf, aes(run_trial, v_chosen, color = pe_f2_hipp_resp)) + geom_smooth(method = 'loess')+ facet_wrap(~rewFunc)
p2 <- ggplot(lmdf, aes(run_trial, v_chosen, color = pe_f2_hipp_resp)) + geom_smooth(method = 'loess')+ facet_wrap(~rewFunc)
p3 <- ggplot(ldf, aes(run_trial, v_chosen, color = h_HippAntL_resp)) + geom_smooth(method = 'loess')+ facet_wrap(~rewFunc)
p4 <- ggplot(lmdf, aes(run_trial, v_chosen, color = h_HippAntL_resp)) + geom_smooth(method = 'loess')+ facet_wrap(~rewFunc)
setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/hippo/figs/beta_uncertainty/')
pdf('v_chosen_by_PH_AH.pdf', width = 9, height = 8)
ggarrange(p1,p2, p3, p4, ncol = 2, nrow = 2, labels = c("fMRI", "replication","fMRI", "replication"))
dev.off()

# most interpretable set of uncertainty models
umb1 <- lmer(u_chosen_quantile ~ (trial_neg_inv_sc + rt_lag_sc + last_outcome + v_entropy_wi + h_HippAntL_neg)^2 +
               (trial_neg_inv_sc + rt_lag_sc + last_outcome + v_entropy_wi + pe_f2_hipp)^2 +
               scale(u_chosen_quantile_lag) + rt_lag_sc*rewFunc + (1|id/run), df)
screen.lmerTest(umb1, .05)

umb1v <- lmer(u_chosen_quantile ~ (trial_neg_inv_sc + rt_lag_sc + last_outcome + v_entropy_wi + h_HippAntL_neg)^2 +
                (trial_neg_inv_sc + rt_lag_sc + last_outcome + v_entropy_wi + pe_f2_hipp)^2 + scale(u_chosen_quantile_lag) + v_chosen_quantile_change + rt_lag_sc*rewFunc + (1|id/run), df)
screen.lmerTest(umb1v, .05)
summary(umb1v)

# attempt MEG replication
mumb1 <- lmer(u_chosen ~ (trial_neg_inv_sc + rt_lag_sc + last_outcome + v_entropy_wi + h_HippAntL_neg)^2 +
                (trial_neg_inv_sc + rt_lag_sc + last_outcome + v_entropy_wi + pe_f2_hipp)^2 +
                scale(u_chosen_lag) + rt_lag_sc*rewFunc + (1|id/run), mdf)
screen.lmerTest(mumb1, .05)

mumb1v <- lmer(u_chosen_change ~ (trial_neg_inv_sc + rt_lag_sc + last_outcome + v_entropy_wi + h_HippAntL_neg)^2 +
                 (trial_neg_inv_sc + rt_lag_sc + last_outcome + v_entropy_wi + pe_f2_hipp)^2 +
                 scale(u_chosen_lag) + v_chosen_change + rt_lag_sc*rewFunc + (1|id/run), mdf)
screen.lmerTest(mumb1v, .05)

# remove entropy to minimize circularity
umb2 <- lmer(u_chosen_quantile ~ trial_neg_inv_sc + rt_lag_sc * rt_vmax_lag_sc * h_HippAntL_neg + 
               rt_lag_sc * rt_vmax_lag_sc * pe_f2_hipp + last_outcome + 
               scale(u_chosen_quantile_lag) +  rt_lag_sc*rewFunc + (1|id/run), df)
screen.lmerTest(umb2, .05)
# summary(umb2)
umb2v <- lmer(u_chosen_quantile ~ trial_neg_inv_sc + rt_lag_sc * rt_vmax_lag_sc * h_HippAntL_neg + 
                rt_lag_sc * rt_vmax_lag_sc * pe_f2_hipp + last_outcome + 
                scale(u_chosen_quantile_lag) +  v_chosen_quantile_change + rt_lag_sc*rewFunc + (1|id/run), df)
screen.lmerTest(umb2v, .05)
# summary(umb2v)

# use U instead of quantile(U) to enable MEG replication
umb3 <- lmer(u_chosen_change ~ trial_neg_inv_sc + rt_lag_sc * rt_vmax_lag_sc * h_HippAntL_neg + 
               rt_lag_sc * rt_vmax_lag_sc * pe_f2_hipp + last_outcome + 
               scale(u_chosen_lag) +  rt_lag_sc*rewFunc + (1|id/run), df)
screen.lmerTest(umb3, .05)
# summary(umb3)
umb3v <- lmer(u_chosen_change ~ trial_neg_inv_sc + rt_lag_sc * rt_vmax_lag_sc * h_HippAntL_neg + 
                rt_lag_sc * rt_vmax_lag_sc * pe_f2_hipp + last_outcome + 
                scale(u_chosen_lag) +  v_chosen_change + rt_lag_sc*rewFunc + (1|id/run), df)
screen.lmerTest(umb3v, .05)
# summary(umb3v)

# MEG replication
mumb3 <- lmer(u_chosen_change ~ trial_neg_inv_sc + rt_lag_sc * rt_vmax_lag_sc * h_HippAntL_neg + 
                rt_lag_sc * rt_vmax_lag_sc * pe_f2_hipp + last_outcome + 
                scale(u_chosen_lag) +  rt_lag_sc*rewFunc + (1|id/run), mdf)
screen.lmerTest(mumb3, .05)
# summary(mumb3)
mumb3v <- lmer(u_chosen_change ~ trial_neg_inv_sc + rt_lag_sc * rt_vmax_lag_sc * h_HippAntL_neg + 
                 rt_lag_sc * rt_vmax_lag_sc * pe_f2_hipp + last_outcome + 
                 scale(u_chosen_lag) +  v_chosen_change + rt_lag_sc*rewFunc + (1|id/run), mdf)
screen.lmerTest(mumb3v, .05)
# summary(mumb3v)

# plot uncertainty change following long vs short RTs
setwd(file.path(clock_folder, 'fmri/keuka_brain_behavior_analyses/plots'))
p1 <- ggplot(df, aes(rt_lag, u_chosen_quantile, color = pe_f2_hipp_resp)) + geom_smooth(method = 'loess') #+ facet_wrap(~rewFunc)
p2 <- ggplot(mdf, aes(rt_lag, u_chosen_quantile, color = pe_f2_hipp_resp)) + geom_smooth(method = 'loess')

# pdf('uncertainty_change_prev_RT_by_PH.pdf', height = 6, width = 12)
ggarrange(p1,p2, labels = c("fMRI", "MEG"))
# dev.off()
# summary(ub3v)
p1 <- ggplot(df %>% filter(rewFunc=='CEVR'), aes(rt_csv, u_chosen_quantile, color = pe_f2_hipp_resp)) + geom_smooth() #+ facet_wrap(~rewFunc)
p2 <- ggplot(mdf %>% filter(rewFunc=='CEVR'), aes(rt_csv, u_chosen_quantile, color = pe_f2_hipp_resp)) + geom_smooth()
# pdf('uncertainty_change_prev_RT_by_PH.pdf', height = 6, width = 12)
ggarrange(p1,p2, labels = c("fMRI", "MEG"))


# sanity check in simple models: they both predict uncertainty-aversion
us1 <- lmer(u_chosen_quantile ~ (trial_neg_inv_sc + h_HippAntL_neg + rewFunc + last_outcome)^3 +
              (trial_neg_inv_sc + pe_f2_hipp + rewFunc + last_outcome )^3 + (1|id/run), df)
anova(us1)
screen.lmerTest(us1, .05)
emu <- as_tibble(emmeans(us1, ~trial_neg_inv_sc | last_outcome * pe_f2_hipp, at = list(pe_f2_hipp = c(-2,2), trial_neg_inv_sc = c(-2,2)))) %>% mutate(`Chosen uncertainty quantile` = emmean)
ggplot(emu, aes(trial_neg_inv_sc, `Chosen uncertainty quantile`, color = as.factor(pe_f2_hipp), lty = last_outcome)) + geom_line()
us1v <- lmer(u_chosen_quantile ~ (trial_neg_inv_sc + h_HippAntL_neg + rewFunc+ last_outcome)^3 +
               (trial_neg_inv_sc + pe_f2_hipp + rewFunc+ last_outcome)^3 + v_chosen_quantile + (1|id/run), df)
screen.lmerTest(us1v, .05)


mus1 <- lmer(u_chosen ~ (trial_neg_inv_sc + h_HippAntL_neg + rewFunc+ last_outcome + rt_swing)^3 +
               (trial_neg_inv_sc + pe_f2_hipp + rewFunc+ last_outcome + rt_swing)^3 + (1|id/run), mdf)
screen.lmerTest(mus1, .05)
emu <- as_tibble(emmeans(mus1, ~trial_neg_inv_sc | last_outcome * pe_f2_hipp, at = list(pe_f2_hipp = c(-2,2), trial_neg_inv_sc = c(-2,2)))) %>% mutate(`Chosen uncertainty quantile` = emmean)
ggplot(emu, aes(trial_neg_inv_sc, `Chosen uncertainty quantile`, color = as.factor(pe_f2_hipp), lty = last_outcome)) + geom_line()

mus1v <- lmer(u_chosen ~ (trial_neg_inv_sc + h_HippAntL_neg + rewFunc+ last_outcome+ rt_swing)^3 +
                (trial_neg_inv_sc + pe_f2_hipp + rewFunc+ last_outcome+ rt_swing)^3 + v_chosen + (1|id/run), mdf)
screen.lmerTest(mus1v, .05)



# demonstrate that this holds across contingencies and early/late learning
pdf('AH_entropy_uncertainty_aversion_by_cond.pdf', height = 6, width = 8)
ggplot(df %>% filter(!is.na(v_entropy_wi)), aes(run_trial, u_chosen_quantile, lty = v_entropy_wi>0, color = h_HippAntL_resp)) + 
  geom_smooth(method = 'gam', formula = y ~ splines::ns(x,3)) + facet_wrap(~rewFunc)
dev.off()


# plot the effect of HIPP on uncertainty sensitivity
p1 <- ggplot(df, aes(run_trial, v_chosen, lty = pe_f2_hipp_resp, group = pe_f2_hipp_resp)) + geom_smooth() + facet_wrap(~rewFunc)
p2 <- ggplot(df, aes(run_trial, u_chosen_quantile, lty = pe_f2_hipp_resp, group = (pe_f2_hipp_resp))) + geom_smooth() + facet_wrap(~rewFunc)
p3 <- ggplot(df, aes(run_trial, v_chosen, lty = h_HippAntL_resp, group = h_HippAntL_resp)) + geom_smooth()+ facet_wrap(~rewFunc)
p4 <- ggplot(df, aes(run_trial, u_chosen_quantile, lty = h_HippAntL_resp, group = h_HippAntL_resp)) + geom_smooth() + facet_wrap(~rewFunc)
pdf("PH_AH_on_u_sensitivity_rewFunc.pdf", height = 8, width = 8)
ggarrange(p1,p2,p3,p4,ncol = 2, nrow = 2)
dev.off()

p1 <- ggplot(df, aes(run_trial, u_chosen_quantile, lty = pe_f1_cort_str_resp, group = interaction(pe_f1_cort_str_resp, last_outcome), color = last_outcome)) + geom_smooth() #+ facet_wrap(~run)
p2 <- ggplot(df, aes(run_trial, u_chosen_quantile, lty = pe_f2_hipp_resp, group = interaction(pe_f2_hipp_resp, last_outcome), color = last_outcome)) + geom_smooth() #+ facet_wrap(~run)
p3 <- ggplot(df, aes(run_trial, u_chosen_quantile, lty = h_f1_fp_resp, group = interaction(h_f1_fp_resp, last_outcome), color = last_outcome)) + geom_smooth() #+ facet_wrap(~run)
p4 <- ggplot(df, aes(run_trial, u_chosen_quantile, lty = h_HippAntL_resp, group = interaction(h_HippAntL_resp, last_outcome), color = last_outcome)) + geom_smooth() #+ facet_wrap(~run)

pdf("clusters_u_sensitivity_reward.pdf", height = 4, width = 8)
ggarrange(p1,p2,p3, p4, ncol = 2, nrow = 2)
dev.off()

####### Nat Comm R1: timecourses of PE and H
# p1 <- ggplot(df, aes(run_trial, pe_max, color = rewFunc)) + geom_smooth(method = 'gam', formula = y ~ splines::ns(x, 3))
# p2 <- ggplot(df, aes(run_trial, v_entropy_wi, color = rewFunc)) + geom_smooth(method = 'gam', formula = y ~ splines::ns(x, 3))
p1 <- ggplot(df, aes(run_trial, pe_max, color = rewFunc)) + geom_smooth(method = 'loess')
p2 <- ggplot(df, aes(run_trial, v_entropy_wi, color = rewFunc)) + geom_smooth(method = 'loess')

pdf("timecourses_of_pe_h_by_cond.pdf", height = 6, width = 4)
ggarrange(p1, p2, ncol = 1, nrow = 2)
dev.off()
#
save(file = 'vhd_u_meg_models.Rdata', list = ls(all.names = TRUE))
# load(file = 'vhd_u_meg_models.Rdata')