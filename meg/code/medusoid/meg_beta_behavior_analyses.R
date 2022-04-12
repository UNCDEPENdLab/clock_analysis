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
# library(dependlab)
source('~/code/Rhelpers/screen.lmerTest.R')
source('~/code/Rhelpers/vif.lme.R')
# library(stringi)

# clock_folder <- "~/Data_Analysis/clock_analysis" #michael
clock_folder <- "~/code/clock_analysis" #alex
# source('~/code/Rhelpers/')
setwd(file.path(clock_folder, 'fmri/keuka_brain_behavior_analyses/dan'))

### load data
source("~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/get_trial_data.R")
#get_trial_data <- function(repo_directory=NULL, dataset="mmclock_fmri", groupfixed=TRUE) 
df <- get_trial_data(repo_directory = clock_folder, dataset = "mmclock_meg", groupfixed = T)

# add meg data
# wbetas <- readRDS("~/OneDrive/collected_letters/papers/meg/plots/wholebrain/betas/MEG_betas_wide_echange_vmax_reward_Nov30_2021.RDS") %>% 
wbetas <- readRDS("~/code/clock_analysis/meg/data/MEG_betas_entropy_change_v_max_reward_signed_pe_rs_abs_pe_Mar_14_2022.RDS") %>% 
  mutate(omission_early_theta = - reward_early_theta,
         omission_late_delta = - reward_late_delta) %>% 
  mutate(entropy_change_early_beta_supp = -  entropy_change_early_beta_entropy_change,
         entropy_change_late_beta_supp = - entropy_change_late_beta_entropy_change,
         # entropy_change_early_beta_supp_ec = -  entropy_change_early_beta_entropy_change_ec_sensors,
         # entropy_change_late_beta_supp_ec = - entropy_change_late_beta_entropy_change_ec_sensors,
         # abs_pe_late_beta_supp_ec = - abspe_ec_late_beta,
         pe_late_beta_supp = - pe_late_beta,
         neg_pe_early_theta = - pe_early_theta,
         vmax_late_alpha = vmax_late_beta) %>%
  select(c(id, omission_early_theta, omission_late_delta, 
           entropy_change_early_beta_supp, entropy_change_late_beta_supp,
           # entropy_change_early_beta_supp_ec, entropy_change_late_beta_supp_ec, 
           vmax_late_alpha, vmax_late_beta,
           neg_pe_early_theta, pe_late_beta_supp,
           abspe_early_theta))
           # abs_pe_late_beta_supp_ec, abs_pe_late_beta_supp))
# merge
df <- df %>% inner_join(wbetas, by = "id")


############# MEG
# Effect of MEG betas on exploration
# only the effects of interest
meg_late_beta_only <-  
  lme4::lmer(rt_csv_sc ~ 
               rt_lag_sc * last_outcome * entropy_change_late_beta_supp + 
               rt_vmax_lag_sc * trial_neg_inv_sc * entropy_change_late_beta_supp + 
               (1|id/run), df %>% filter(rt_csv<4000))
# screen.lmerTest(meg1, .01)
summary(meg_late_beta_only)
Anova(meg_late_beta_only, '3')


# because post-omission theta and delta are correlated at 0.75, will not enter delta here
meg1_om_theta <-  
  lme4::lmer(rt_csv_sc ~ 
               rt_lag_sc * last_outcome * entropy_change_early_beta_supp + 
               rt_lag_sc * last_outcome * entropy_change_late_beta_supp + 
               rt_lag_sc * last_outcome * abspe_early_theta +
               rt_lag_sc * last_outcome * omission_early_theta +
               rt_vmax_lag_sc * trial_neg_inv_sc * entropy_change_early_beta_supp + 
               rt_vmax_lag_sc * trial_neg_inv_sc * entropy_change_late_beta_supp + 
               rt_vmax_lag_sc * trial_neg_inv_sc * abspe_early_theta  +
               rt_vmax_lag_sc * trial_neg_inv_sc * omission_early_theta  +
               (1|id/run), df %>% filter(rt_csv<4000))
# screen.lmerTest(meg1, .01)
summary(meg1_om_theta)
Anova(meg1_om_theta, '3')

meg1_pe_theta <-  
  lme4::lmer(rt_csv_sc ~ 
               rt_lag_sc * last_outcome * entropy_change_early_beta_supp + 
               rt_lag_sc * last_outcome * entropy_change_late_beta_supp + 
               rt_lag_sc * last_outcome * vmax_late_alpha +
               rt_lag_sc * last_outcome * neg_pe_early_theta +
               rt_vmax_lag_sc * trial_neg_inv_sc * entropy_change_early_beta_supp + 
               rt_vmax_lag_sc * trial_neg_inv_sc * entropy_change_late_beta_supp + 
               rt_vmax_lag_sc * trial_neg_inv_sc * vmax_late_alpha  +
               rt_vmax_lag_sc * trial_neg_inv_sc * neg_pe_early_theta  +
               (1|id/run), df %>% filter(rt_csv<4000))
# screen.lmerTest(meg1, .01)
summary(meg1_pe_theta)
Anova(meg1_pe_theta, '3')
anova(meg1_om_theta, meg1_pe_theta)
# note: omission predicts much better than PE early theta

meg_ec <-  
  lme4::lmer(rt_csv_sc ~ 
               rt_lag_sc * last_outcome * entropy_change_early_beta_supp + 
               rt_lag_sc * last_outcome * entropy_change_late_beta_supp + 
               rt_vmax_lag_sc * trial_neg_inv_sc * entropy_change_early_beta_supp + 
               rt_vmax_lag_sc * trial_neg_inv_sc * entropy_change_late_beta_supp + 
               (1|id/run), df %>% filter(rt_csv<4000))
# screen.lmerTest(meg1, .01)
summary(meg_ec)
Anova(meg_ec, '3')

meg_ec_pe <-  
  lme4::lmer(rt_csv_sc ~ 
               rt_lag_sc * last_outcome * entropy_change_early_beta_supp + 
               rt_lag_sc * last_outcome * entropy_change_late_beta_supp + 
               rt_lag_sc * last_outcome * pe_late_beta_supp +
               rt_lag_sc * last_outcome * neg_pe_early_theta + 
               rt_vmax_lag_sc * trial_neg_inv_sc * entropy_change_early_beta_supp + 
               rt_vmax_lag_sc * trial_neg_inv_sc * entropy_change_late_beta_supp + 
               rt_vmax_lag_sc * trial_neg_inv_sc * pe_late_beta_supp +
               rt_vmax_lag_sc * trial_neg_inv_sc * neg_pe_early_theta +
               (1|id/run), df %>% filter(rt_csv<4000))
# screen.lmerTest(meg1, .01)
summary(meg_ec_pe)
Anova(meg_ec_pe, '3')

# hypothesis: best model combines omission early theta and PE late beta with EC responses
# correct!
meg_hybrid <-  
  lme4::lmer(rt_csv_sc ~ 
               # rt_lag_sc * last_outcome * entropy_change_early_beta_supp + 
               rt_lag_sc * last_outcome * entropy_change_late_beta_supp + 
               rt_lag_sc * last_outcome * omission_early_theta + 
               rt_vmax_lag_sc * trial_neg_inv_sc * entropy_change_late_beta_supp + 
               rt_vmax_lag_sc * trial_neg_inv_sc * omission_early_theta +
               (1|id/run), df %>% filter(rt_csv<4000))
# screen.lmerTest(meg1, .01)
summary(meg_hybrid)
Anova(meg_hybrid, '3')

# compare models
anova(meg_hybrid, meg_ec_pe, meg_ec, meg1_om_theta)

# add random slopes of RT and RT_Vmax
meg_hybrid_rs <-  
  lme4::lmer(rt_csv_sc ~ 
               # rt_lag_sc * last_outcome * entropy_change_early_beta_supp + 
               rt_lag_sc * last_outcome * entropy_change_late_beta_supp + 
               rt_lag_sc * last_outcome * omission_early_theta + 
               rt_vmax_lag_sc * trial_neg_inv_sc * entropy_change_late_beta_supp + 
               rt_vmax_lag_sc * trial_neg_inv_sc * omission_early_theta +
               (rt_lag_sc + rt_vmax_lag_sc|id/run), df %>% filter(rt_csv<4000))
# screen.lmerTest(meg1, .01)
summary(meg_hybrid_rs)
Anova(meg_hybrid_rs, '3')

# examine interaction with condition: everything holds AND late beta enhances condition effects over trials
meg_cond <-  
  lme4::lmer(rt_csv_sc ~ 
               # rt_lag_sc * last_outcome * entropy_change_early_beta_supp + 
               (rt_lag_sc + last_outcome +  entropy_change_late_beta_supp + rewFunc)^3 + 
               (rt_lag_sc + last_outcome + omission_early_theta + rewFunc)^3 + 
               # rt_vmax_lag_sc * trial_neg_inv_sc * entropy_change_early_beta_supp + 
               (rt_vmax_lag_sc + trial_neg_inv_sc + entropy_change_late_beta_supp + rewFunc)^3 + 
               (rt_vmax_lag_sc + trial_neg_inv_sc + omission_early_theta + rewFunc)^3 + 
               (1|id/run), df %>% filter(rt_csv<4000))
# screen.lmerTest(meg1, .01)
summary(meg_cond)
Anova(meg_cond, '3')
vif(meg_cond)

# load fmri data
fdf <- get_trial_data(repo_directory = clock_folder, dataset = "mmclock_fmri", groupfixed = T)
fdf <- fdf %>% mutate(id = as.character(id)) %>% inner_join(wbetas, by = "id")

# now that we have done model selection in the MEG study, replicate only the best model
fmri_hybrid <-  
  lme4::lmer(rt_csv_sc ~ 
               # rt_lag_sc * last_outcome * entropy_change_early_beta_supp + 
               rt_lag_sc * last_outcome * entropy_change_late_beta_supp + 
               rt_lag_sc * last_outcome * omission_early_theta + 
               # rt_vmax_lag_sc * trial_neg_inv_sc * entropy_change_early_beta_supp + 
               rt_vmax_lag_sc * trial_neg_inv_sc * entropy_change_late_beta_supp + 
               rt_vmax_lag_sc * trial_neg_inv_sc * omission_early_theta +
               (1|id/run), fdf %>% filter(rt_csv<4000))
# screen.lmerTest(meg1, .01)
summary(fmri_hybrid)
Anova(fmri_hybrid, '3')

# add random slopes of RT and RT_Vmax
fmri_hybrid_rs <-  
  lme4::lmer(rt_csv_sc ~ 
               # rt_lag_sc * last_outcome * entropy_change_early_beta_supp + 
               rt_lag_sc * last_outcome * entropy_change_late_beta_supp + 
               rt_lag_sc * last_outcome * pe_late_beta_supp +
               rt_lag_sc * last_outcome * omission_early_theta + 
               # rt_vmax_lag_sc * trial_neg_inv_sc * entropy_change_early_beta_supp + 
               rt_vmax_lag_sc * trial_neg_inv_sc * entropy_change_late_beta_supp + 
               rt_vmax_lag_sc * trial_neg_inv_sc * pe_late_beta_supp +
               rt_vmax_lag_sc * trial_neg_inv_sc * omission_early_theta +
               (rt_lag_sc + rt_vmax_lag_sc|id/run), fdf %>% filter(rt_csv<4000))
# screen.lmerTest(meg1, .01)
summary(fmri_hybrid_rs)
Anova(fmri_hybrid_rs, '3')

# only the effects of interest
fmri_late_beta_only <-  
  lme4::lmer(rt_csv_sc ~ 
               rt_lag_sc * last_outcome * entropy_change_late_beta_supp + 
               rt_vmax_lag_sc * trial_neg_inv_sc * entropy_change_late_beta_supp + 
               (1|id/run), fdf %>% filter(rt_csv<4000))
# screen.lmerTest(meg1, .01)
summary(fmri_late_beta_only)
Anova(fmri_late_beta_only, '3')


fmri1theta <- lme4::lmer(rt_csv_sc ~ 
                         rt_lag_sc * last_outcome * entropy_change_early_beta_supp + 
                         rt_lag_sc * last_outcome * entropy_change_late_beta_supp + 
                         rt_lag_sc * last_outcome * abspe_early_theta +
                         rt_lag_sc * last_outcome * omission_early_theta +
                         rt_vmax_lag_sc * trial_neg_inv_sc * entropy_change_early_beta_supp + 
                         rt_vmax_lag_sc * trial_neg_inv_sc * entropy_change_late_beta_supp + 
                         rt_vmax_lag_sc * trial_neg_inv_sc * abspe_early_theta  +
                         rt_vmax_lag_sc * trial_neg_inv_sc * omission_early_theta  +
                         (1|id/run), fdf %>% filter(rt_csv<4000))
# screen.lmerTest(fmri1, .01)
summary(fmri1theta)
Anova(fmri1theta, '3')


fmri_cond <-  
  lme4::lmer(rt_csv_sc ~ 
               # rt_lag_sc * last_outcome * entropy_change_early_beta_supp + 
               (rt_lag_sc + last_outcome +  entropy_change_late_beta_supp + rewFunc)^3 + 
               (rt_lag_sc + last_outcome + omission_early_theta + rewFunc)^3 + 
               # rt_vmax_lag_sc * trial_neg_inv_sc * entropy_change_early_beta_supp + 
               (rt_vmax_lag_sc + trial_neg_inv_sc + entropy_change_late_beta_supp + rewFunc)^3 + 
               (rt_vmax_lag_sc + trial_neg_inv_sc + omission_early_theta + rewFunc)^3 + 
               (1|id/run), fdf %>% filter(rt_csv<4000))
# screen.lmerTest(meg1, .01)
summary(fmri_cond)
Anova(fmri_cond, '3')
vif(fmri_cond)


# fmri1delta <- lme4::lmer(rt_csv_sc ~ 
#                               rt_lag_sc * last_outcome * entropy_change_early_beta_supp + 
#                               rt_lag_sc * last_outcome * entropy_change_late_beta_supp + 
#                               rt_lag_sc * last_outcome * vmax_late_alpha +
#                               rt_lag_sc * last_outcome * omission_late_delta +
#                               rt_vmax_lag_sc * trial_neg_inv_sc * entropy_change_early_beta_supp + 
#                               rt_vmax_lag_sc * trial_neg_inv_sc * entropy_change_late_beta_supp + 
#                               rt_vmax_lag_sc * trial_neg_inv_sc * vmax_late_alpha  +
#                               rt_vmax_lag_sc * trial_neg_inv_sc * omission_late_delta  +
#                               (1|id/run), fdf %>% filter(rt_csv<4000))
# # screen.lmerTest(fmri1, .01)
# summary(fmri1delta)
# Anova(fmri1delta, '3')

# # late betas suppression to EC from only the posterior "EC" sensors predicts as well as from all sensors
# fmri_ec <-  
#   lme4::lmer(rt_csv_sc ~ 
#                rt_lag_sc * last_outcome * entropy_change_early_beta_supp + 
#                rt_lag_sc * last_outcome * entropy_change_late_beta_supp + 
#                rt_vmax_lag_sc * trial_neg_inv_sc * entropy_change_early_beta_supp + 
#                rt_vmax_lag_sc * trial_neg_inv_sc * entropy_change_late_beta_supp + 
#                (1|id/run), fdf %>% filter(rt_csv<4000))
# # screen.lmerTest(fmri1, .01)
# summary(fmri_ec)
# Anova(fmri_ec, '3')
# 
# # add abs_pe_late_beta_supp_ec
# fmri_ec_pe <-  
#   lme4::lmer(rt_csv_sc ~ 
#                rt_lag_sc * last_outcome * entropy_change_early_beta_supp + 
#                rt_lag_sc * last_outcome * entropy_change_late_beta_supp + 
#                rt_lag_sc * last_outcome * pe_late_beta_supp +
#                rt_lag_sc * last_outcome * neg_pe_early_theta + 
#                rt_vmax_lag_sc * trial_neg_inv_sc * entropy_change_early_beta_supp + 
#                rt_vmax_lag_sc * trial_neg_inv_sc * entropy_change_late_beta_supp + 
#                rt_vmax_lag_sc * trial_neg_inv_sc * pe_late_beta_supp +
#                rt_vmax_lag_sc * trial_neg_inv_sc * neg_pe_early_theta +
#                (1|id/run), fdf %>% filter(rt_csv<4000))
# # screen.lmerTest(meg1, .01)
# summary(fmri_ec_pe)
# Anova(fmri_ec_pe, '3')

###########################
## RT lag, reward effects
###########################

# late beta to entropy change
em1 <- as_tibble(emtrends(meg_hybrid, data = df,  var = "rt_lag_sc", specs = c("entropy_change_late_beta_supp", "last_outcome"), at = list(entropy_change_late_beta_supp = c(-.071, .221)), options = list()))
em1$study = "1. MEG"
em2 <- as_tibble(emtrends(fmri_hybrid, data = fdf, var = "rt_lag_sc", specs = c("entropy_change_late_beta_supp", "last_outcome"), at = list(entropy_change_late_beta_supp = c(-.071, .221)), options = list()))
em2$study = '2. fMRI replication'
em1 <- rbind(em2, em1)
ec_late_beta <- ggplot(em1, aes(x=last_outcome, y=rt_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=as.factor(entropy_change_late_beta_supp))) + 
  #shape = as.factor(pe_f2_hipp), 
  #p1 <- ggplot(em1, aes(x=last_outcome, y=rt_lag_sc.trend, linetype = as.factor(pe_f2_hipp), ymin=asymp.LCL, ymax=asymp.UCL)) + 
  geom_point(position = position_dodge(width = .6), size=2.5) + 
  #geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), position = position_dodge2(width=0.9), width=0.5) + 
  geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + 
  #geom_crossbar(position=position_dodge(width=0.7), width=0.6) + 
  theme_bw(base_size=12) + facet_wrap(~study)+ ylab("RT swings (AU)\n Small <---------> Large")  + 
  #scale_shape_manual(values=c(15,16), labels = c("10th %ile", "90th %ile")) + 
  #scale_color_brewer("PH RPE\nresponse", palette="Set1", labels = c("10th %ile", "90th %ile")) +
  scale_color_manual("Entropy change\nlate beta\nsuppression", values=c("#1b3840","#4fa3b8"), labels = c("10th %ile", "90th %ile")) +
  labs(shape = "Entropy change\nlate beta\nsuppression") +
  theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) + 
  scale_y_reverse(limits = c(.65, 0)) 

# early theta to reward omission
rem1 <- as_tibble(emtrends(meg_hybrid, data = df,  var = "rt_lag_sc", specs = c("omission_early_theta", "last_outcome"), at = list(omission_early_theta = c(-.040, .610)), options = list()))
rem1$study = "1. MEG"
rem2 <- as_tibble(emtrends(fmri_hybrid, data = fdf, var = "rt_lag_sc", specs = c("omission_early_theta", "last_outcome"), at = list(omission_early_theta = c(-.040, .610)), options = list()))
rem2$study = '2. fMRI replication'
rem1 <- rbind(rem2, rem1)
reward_early_theta <- ggplot(rem1, aes(x=last_outcome, y=rt_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=as.factor(omission_early_theta))) + 
  #shape = as.factor(pe_f2_hipp), 
  #p1 <- ggplot(em1, aes(x=last_outcome, y=rt_lag_sc.trend, linetype = as.factor(pe_f2_hipp), ymin=asymp.LCL, ymax=asymp.UCL)) + 
  geom_point(position = position_dodge(width = .6), size=2.5) + 
  #geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), position = position_dodge2(width=0.9), width=0.5) + 
  geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + 
  #geom_crossbar(position=position_dodge(width=0.7), width=0.6) + 
  theme_bw(base_size=12) + facet_wrap(~study)+ ylab("RT swings (AU)\n Small <---------> Large")  + 
  #scale_shape_manual(values=c(15,16), labels = c("10th %ile", "90th %ile")) + 
  #scale_color_brewer("PH RPE\nresponse", palette="Set1", labels = c("10th %ile", "90th %ile")) +
  scale_color_manual("Reward omission\nearly theta\nsynchronization", values=c("#1b3840","#4fa3b8"), labels = c("10th %ile", "90th %ile")) +
  labs(shape = "Reward omission\nearly theta\nsynchronization") +
  theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) + 
  scale_y_reverse(limits = c(.6, 0)) 

setwd('~/OneDrive/collected_letters/papers/meg/plots/meg_to_behavior/')
pdf("meg_session_level_to_beh.pdf", height = 4, width = 11)
ggarrange(reward_early_theta, ec_late_beta)
dev.off()

###### same with condition
emc1 <- as_tibble(emtrends(meg_cond, data = df,  var = "rt_lag_sc", specs = c("entropy_change_late_beta_supp", "last_outcome", "rewFunc"), at = list(entropy_change_late_beta_supp = c(-.071, .221)), options = list()))
emc1$study = "1. MEG"
emc2 <- as_tibble(emtrends(fmri_cond, data = fdf, var = "rt_lag_sc", specs = c("entropy_change_late_beta_supp", "last_outcome", "rewFunc"), at = list(entropy_change_late_beta_supp = c(-.071, .221)), options = list()))
emc2$study = '2. fMRI replication'
emc1 <- rbind(emc2, emc1)
ec_late_beta_cond <- ggplot(emc1, aes(x=last_outcome, y=rt_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=as.factor(entropy_change_late_beta_supp))) + 
  #shape = as.factor(pe_f2_hipp), 
  #p1 <- ggplot(emc1, aes(x=last_outcome, y=rt_lag_sc.trend, linetype = as.factor(pe_f2_hipp), ymin=asymp.LCL, ymax=asymp.UCL)) + 
  geom_point(position = position_dodge(width = .6), size=2.5) + 
  #geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), position = position_dodge2(width=0.9), width=0.5) + 
  geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + 
  #geom_crossbar(position=position_dodge(width=0.7), width=0.6) + 
  theme_bw(base_size=12) + facet_grid(rewFunc~study)+ ylab("RT swings (AU)\n Small <---------> Large")  + 
  #scale_shape_manual(values=c(15,16), labels = c("10th %ile", "90th %ile")) + 
  #scale_color_brewer("PH RPE\nresponse", palette="Set1", labels = c("10th %ile", "90th %ile")) +
  scale_color_manual("Entropy change\nlate beta\nsuppression", values=c("#1b3840","#4fa3b8"), labels = c("10th %ile", "90th %ile")) +
  labs(shape = "Entropy change\nlate beta\nsuppression") +
  theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) + 
  scale_y_reverse(limits = c(.7, 0)) 

# early theta to reward omission
remc1 <- as_tibble(emtrends(meg_cond, data = df,  var = "rt_lag_sc", specs = c("omission_early_theta", "last_outcome", "rewFunc"), at = list(omission_early_theta = c(-.040, .610)), options = list()))
remc1$study = "1. MEG"
remc2 <- as_tibble(emtrends(fmri_cond, data = fdf, var = "rt_lag_sc", specs = c("omission_early_theta", "last_outcome", "rewFunc"), at = list(omission_early_theta = c(-.040, .610)), options = list()))
remc2$study = '2. fMRI replication'
remc1 <- rbind(remc2, remc1)
reward_early_theta_cond <- ggplot(remc1, aes(x=last_outcome, y=rt_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=as.factor(omission_early_theta))) + 
  #shape = as.factor(pe_f2_hipp), 
  #p1 <- ggplot(emc1, aes(x=last_outcome, y=rt_lag_sc.trend, linetype = as.factor(pe_f2_hipp), ymin=asymp.LCL, ymax=asymp.UCL)) + 
  geom_point(position = position_dodge(width = .6), size=2.5) + 
  #geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), position = position_dodge2(width=0.9), width=0.5) + 
  geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + 
  #geom_crossbar(position=position_dodge(width=0.7), width=0.6) + 
  theme_bw(base_size=12) + facet_grid(rewFunc~study) + ylab("RT swings (AU)\n Small <---------> Large")  + 
  #scale_shape_manual(values=c(15,16), labels = c("10th %ile", "90th %ile")) + 
  #scale_color_brewer("PH RPE\nresponse", palette="Set1", labels = c("10th %ile", "90th %ile")) +
  scale_color_manual("Reward omission\nearly theta\nsynchronization", values=c("#1b3840","#4fa3b8"), labels = c("10th %ile", "90th %ile")) +
  labs(shape = "Reward omission\nearly theta\nsynchronization") +
  theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) + 
  scale_y_reverse(limits = c(.7, 0)) 

setwd('~/OneDrive/collected_letters/papers/meg/plots/meg_to_behavior/')
pdf("meg_session_level_to_beh_rewFunc.pdf", height = 8, width = 11)
ggarrange(reward_early_theta_cond, ec_late_beta_cond)
dev.off()
# rlem1 <- as_tibble(emtrends(meg1delta, data = df,  var = "rt_lag_sc", specs = c("omission_late_delta", "last_outcome"), at = list(omission_late_delta = c(-.3, .31)), options = list()))
# rlem1$study = "1. MEG"
# rlem2 <- as_tibble(emtrends(fmri1delta, data = fdf, var = "rt_lag_sc", specs = c("omission_late_delta", "last_outcome"), at = list(omission_late_delta = c(-.3, .31)), options = list()))
# rlem2$study = '2. fMRI replication'
# rlem1 <- rbind(rlem2, rlem1)


# reward_late_delta <- ggplot(rlem1, aes(x=last_outcome, y=rt_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=as.factor(omission_late_delta))) + 
#   #shape = as.factor(pe_f2_hipp), 
#   #p1 <- ggplot(em1, aes(x=last_outcome, y=rt_lag_sc.trend, linetype = as.factor(pe_f2_hipp), ymin=asymp.LCL, ymax=asymp.UCL)) + 
#   geom_point(position = position_dodge(width = .6), size=2.5) + 
#   #geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), position = position_dodge2(width=0.9), width=0.5) + 
#   geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + 
#   #geom_crossbar(position=position_dodge(width=0.7), width=0.6) + 
#   theme_bw(base_size=12) + facet_wrap(~study)+ ylab("RT swings (AU)\n Small <---------> Large")  + 
#   #scale_shape_manual(values=c(15,16), labels = c("10th %ile", "90th %ile")) + 
#   #scale_color_brewer("PH RPE\nresponse", palette="Set1", labels = c("10th %ile", "90th %ile")) +
#   scale_color_manual("Reward omission\nlate delta\nsynchronization", values=c("#1b3840","#4fa3b8"), labels = c("10th %ile", "90th %ile")) +
#   labs(shape = "Reward omission\nlate delta\nsynchronization") +
#   theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
#         axis.text=element_text(size=8.5, color="grey10")) + 
#   scale_y_reverse(limits = c(.7, -.1)) 
# vem1 <- as_tibble(emtrends(meg1theta, data = df,  var = "rt_lag_sc", specs = c("vmax_late_alpha", "last_outcome"), at = list(vmax_late_alpha = c(-.08, .08)), options = list()))
# vem1$study = "1. MEG"
# vem2 <- as_tibble(emtrends(fmri1theta, data = fdf, var = "rt_lag_sc", specs = c("vmax_late_alpha", "last_outcome"), at = list(vmax_late_alpha = c(-.08, .08)), options = list()))
# vem2$study = '2. fMRI replication'
# vem1 <- rbind(vem2, vem1)
# vmax_late_alpha <- ggplot(vem1, aes(x=last_outcome, y=rt_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=as.factor(vmax_late_alpha))) + 
#   #shape = as.factor(pe_f2_hipp), 
#   #p1 <- ggplot(em1, aes(x=last_outcome, y=rt_lag_sc.trend, linetype = as.factor(pe_f2_hipp), ymin=asymp.LCL, ymax=asymp.UCL)) + 
#   geom_point(position = position_dodge(width = .6), size=2.5) + 
#   #geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), position = position_dodge2(width=0.9), width=0.5) + 
#   geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + 
#   #geom_crossbar(position=position_dodge(width=0.7), width=0.6) + 
#   theme_bw(base_size=12) + facet_wrap(~study)+ ylab("RT swings (AU)\n Small <---------> Large")  + 
#   #scale_shape_manual(values=c(15,16), labels = c("10th %ile", "90th %ile")) + 
#   #scale_color_brewer("PH RPE\nresponse", palette="Set1", labels = c("10th %ile", "90th %ile")) +
#   scale_color_manual("Vmax\nlate alpha\nresponse", values=c("#1b3840","#4fa3b8"), labels = c("10th %ile", "90th %ile")) +
#   labs(shape = "Vmax\nlate alpha\nresponse") +
#   theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
#         axis.text=element_text(size=8.5, color="grey10")) + 
#   scale_y_reverse(limits = c(.7, -.1)) 


ec1 <- as_tibble(emtrends(meg_hybrid, data = df,  var = "rt_vmax_lag_sc", specs = c("entropy_change_late_beta_supp", "trial_neg_inv_sc"), at = list(entropy_change_late_beta_supp =  c(-.071, .221), trial_neg_inv_sc = c(-.88, 0.39)), options = list()))
ec1$study = "1. MEG"
ec2 <- as_tibble(emtrends(fmri_hybrid, data = fdf, var = "rt_vmax_lag_sc", specs = c("entropy_change_late_beta_supp", "trial_neg_inv_sc"), at = list(entropy_change_late_beta_supp =  c(-.071, .221), trial_neg_inv_sc = c(-.88, 0.39)), options = list()))
ec2$study = '2. fMRI replication'
ec <- rbind(ec2, ec1)
ec_late_beta_rtvmax <- ggplot(ec, aes(x=as.factor(trial_neg_inv_sc), y=rt_vmax_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=as.factor(entropy_change_late_beta_supp))) + 
  geom_point(position = position_dodge(width = .6), size=2.5) + 
  geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + 
  theme_bw(base_size=12) + facet_wrap(~study)+  ylab("Convergence on\nbest RT (AU)") +
  scale_color_manual("Entropy change\nlate beta\nsuppression", values=c("#403202", "#e2b407"), labels = c("10th %ile", "90th %ile")) +
  labs(shape = "Entropy change\nlate beta\nsuppression") +
  theme(panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) + 
  scale_x_discrete(name ="Trial", labels=c("-0.88" = "5", "0.39" = "50")) +  scale_y_continuous(limits = c(0, .25)) 




# PE late beta suppression: behavioral effects
pe1 <- as_tibble(emtrends(meg_hybrid, data = df,  var = "rt_vmax_lag_sc", specs = c("pe_late_beta_supp", "trial_neg_inv_sc"), at = list(pe_late_beta_supp = c(-.00431, .00597), trial_neg_inv_sc = c(-.88, 0.39)), options = list()))
pe1$study = "1. MEG"
pe2 <- as_tibble(emtrends(fmri_hybrid, data = fdf, var = "rt_vmax_lag_sc", specs = c("pe_late_beta_supp", "trial_neg_inv_sc"), at = list(pe_late_beta_supp = c(-.00431, .00597), trial_neg_inv_sc = c(-.88, 0.39)), options = list()))
pe2$study = '2. fMRI replication'
pe <- rbind(pe2, pe1)
pe_late_beta <- ggplot(pe, aes(x=as.factor(trial_neg_inv_sc), y=rt_vmax_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=as.factor(pe_late_beta_supp))) + 
  geom_point(position = position_dodge(width = .6), size=2.5) + 
  geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + 
  theme_bw(base_size=12) + facet_wrap(~study)+  ylab("Convergence on\nbest RT (AU)") +
  scale_color_manual("Prediction error\nlate beta\nsuppression", values=c("#403202", "#e2b407"), labels = c("10th %ile", "90th %ile")) +
  labs(shape = "Prediction error\nlate beta\nsuppression") +
  theme(panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) + 
  scale_x_discrete(name ="Trial", labels=c("-0.88" = "5", "0.39" = "50")) + scale_y_continuous(limits = c(0, .25))

pe3 <- as_tibble(emtrends(meg_hybrid, data = df,  var = "rt_lag_sc", specs = c("pe_late_beta_supp", "last_outcome"), at = list(pe_late_beta_supp = c(-.00431, .00597)), options = list()))
pe3$study = "1. MEG"
pe4 <- as_tibble(emtrends(fmri_hybrid, data = fdf, var = "rt_lag_sc", specs = c("pe_late_beta_supp", "last_outcome"), at = list(pe_late_beta_supp = c(-.00431, .00597)), options = list()))
pe4$study = '2. fMRI replication'
pe5 <- rbind(pe4, pe3)
pe_late_beta_wsls <- ggplot(pe5, aes(x=last_outcome, y=rt_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=as.factor(pe_late_beta_supp))) + 
  #shape = as.factor(pe_f2_hipp), 
  #p1 <- ggplot(em1, aes(x=last_outcome, y=rt_lag_sc.trend, linetype = as.factor(pe_f2_hipp), ymin=asymp.LCL, ymax=asymp.UCL)) + 
  geom_point(position = position_dodge(width = .6), size=2.5) + 
  #geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), position = position_dodge2(width=0.9), width=0.5) + 
  geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + 
  #geom_crossbar(position=position_dodge(width=0.7), width=0.6) + 
  theme_bw(base_size=12) + facet_wrap(~study)+ ylab("RT swings (AU)\n Small <---------> Large")  + 
  #scale_shape_manual(values=c(15,16), labels = c("10th %ile", "90th %ile")) + 
  #scale_color_brewer("PH RPE\nresponse", palette="Set1", labels = c("10th %ile", "90th %ile")) +
  scale_color_manual("Prediction error\nlate beta\nsuppression", values=c("#1b3840","#4fa3b8"), labels = c("10th %ile", "90th %ile")) +
  labs(shape = "Prediction error\nlate beta\nsuppression") +
  theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) + 
  scale_y_reverse(limits = c(.6, 0)) 



ggarrange(ec_late_beta, ec_late_beta_rtvmax, reward_early_theta, pe_late_beta, pe_late_beta_wsls, ncol = 2, nrow = 3)

setwd("~/OneDrive/collected_letters/papers/meg/plots/wholebrain")
pdf(file = "MEG_to_behavior_ec_pe_reward.pdf", height = 9, width = 12)
ggarrange(ec_late_beta, ec_late_beta_rtvmax, reward_early_theta, pe_late_beta, pe_late_beta_wsls, ncol = 2, nrow = 3)
dev.off()

##################
# RT Vmax effects 
#################

# not replicable or interpretable for late beta suppression
emodel_condition <- as_tibble(emtrends(fmri1delta, var = "rt_vmax_lag_sc", specs = c("omission_late_delta", "trial_neg_inv_sc"), at = list(omission_late_delta = c(-.3, .31), trial_neg_inv_sc = c(-.88, 0.39)), options = list()))
emodel_condition$study = '2. fMRI replication'
emodel_RT_Vmax <- as_tibble(emtrends(meg1delta, var = "rt_vmax_lag_sc", specs = c("omission_late_delta", "trial_neg_inv_sc"), at = list(omission_late_delta = c(-.3, .31), trial_neg_inv_sc = c(-.88, 0.39)), options = list()))
emodel_RT_Vmax$study = "1. MEG"
em4 <- rbind(emodel_condition, emodel_RT_Vmax)
omission_late_delta_rtvmax <- ggplot(em4, aes(x=as.factor(trial_neg_inv_sc), y=rt_vmax_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=as.factor(omission_late_delta))) + 
  geom_point(position = position_dodge(width = .6), size=2.5) + 
  geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + 
  theme_bw(base_size=12) + facet_wrap(~study)+  ylab("Convergence on\nbest RT (AU)") +
  scale_color_manual("Reward omission\nlate delta\nsynchronization", values=c("#403202", "#e2b407"), labels = c("10th %ile", "90th %ile")) +
  labs(shape = "Reward omission\nlate delta\nsynchronization") +
  theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) + 
  scale_x_discrete(name ="Trial", labels=c("-0.88" = "5", "0.39" = "50")) + scale_y_continuous(limits = c(.0, .25))

setwd("~/OneDrive/collected_letters/papers/meg/plots/wholebrain")
pdf(file = "MEG_to_behavior_late_delta_RT_Vmax.pdf", height = 3, width = 6)
print(omission_late_delta_rtvmax)
dev.off()

################
# Add entropy
################
meg2 <-  lme4::lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + rt_vmax_lag_sc + last_outcome + 
                                      v_entropy_wi + entropy_change_late_beta_supp + vmax_late_beta + 
                                      omission_early_theta)^2 + 
                         rt_lag_sc:last_outcome:entropy_change_late_beta_supp + 
                         rt_lag_sc:last_outcome:vmax_late_beta +
                         rt_lag_sc:last_outcome:omission_early_theta +
                         rt_lag_sc:v_entropy_wi:entropy_change_late_beta_supp + 
                         rt_lag_sc:v_entropy_wi:vmax_late_beta +
                         rt_lag_sc:v_entropy_wi:omission_early_theta +
                         rt_vmax_lag_sc:trial_neg_inv_sc:entropy_change_late_beta_supp + 
                         rt_vmax_lag_sc:trial_neg_inv_sc:vmax_late_beta  +
                         rt_vmax_lag_sc:trial_neg_inv_sc:omission_early_theta  +
                         (1|id/run), df %>% filter(rt_csv<4000))
# screen.lmerTest(meg2, .01)
summary(meg2)
Anova(meg2, '3')

# load fmri data
fmri2 <-  lme4::lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + rt_vmax_lag_sc + last_outcome + 
                                       v_entropy_wi + entropy_change_late_beta_supp + vmax_late_beta + 
                                       omission_early_theta)^2 + 
                          rt_lag_sc:last_outcome:entropy_change_late_beta_supp + 
                          rt_lag_sc:last_outcome:vmax_late_beta +
                          rt_lag_sc:last_outcome:omission_early_theta +
                          rt_lag_sc:v_entropy_wi:entropy_change_late_beta_supp + 
                          rt_lag_sc:v_entropy_wi:vmax_late_beta +
                          rt_lag_sc:v_entropy_wi:omission_early_theta +
                          rt_vmax_lag_sc:trial_neg_inv_sc:entropy_change_late_beta_supp + 
                          rt_vmax_lag_sc:trial_neg_inv_sc:vmax_late_beta  +
                          rt_vmax_lag_sc:trial_neg_inv_sc:omission_early_theta  +
                          (1|id/run), fdf %>% filter(rt_csv<4000))
# screen.lmerTest(fmri2, .01)
summary(fmri2)
Anova(fmri2, '3')

# linear trial
meg3 <-  lme4::lmer(rt_csv_sc ~ (scale(run_trial) + rt_lag_sc + rt_vmax_lag_sc + last_outcome + 
                                      v_entropy_wi + entropy_change_late_beta_supp + vmax_late_beta + 
                                      omission_early_theta)^2 + 
                         rt_lag_sc:last_outcome:entropy_change_late_beta_supp + 
                         rt_lag_sc:last_outcome:vmax_late_beta +
                         rt_lag_sc:last_outcome:omission_early_theta +
                         rt_vmax_lag_sc:scale(run_trial):entropy_change_late_beta_supp + 
                         rt_vmax_lag_sc:scale(run_trial):vmax_late_beta  +
                         rt_vmax_lag_sc:scale(run_trial):omission_early_theta  +
                         (1|id/run), df %>% filter(rt_csv<4000))
# screen.lmerTest(meg3, .01)
summary(meg3)
Anova(meg3, '3')

# fMRI
fmri3 <-  lme4::lmer(rt_csv_sc ~ (scale(run_trial) + rt_lag_sc + rt_vmax_lag_sc + last_outcome + 
                                       v_entropy_wi + entropy_change_late_beta_supp + vmax_late_beta + 
                                       omission_early_theta)^2 + 
                          rt_lag_sc:last_outcome:entropy_change_late_beta_supp + 
                          rt_lag_sc:last_outcome:vmax_late_beta +
                          rt_lag_sc:last_outcome:omission_early_theta +
                          rt_vmax_lag_sc:scale(run_trial):entropy_change_late_beta_supp + 
                          rt_vmax_lag_sc:scale(run_trial):vmax_late_beta  +
                          rt_vmax_lag_sc:scale(run_trial):omission_early_theta  +
                          (1|id/run), fdf %>% filter(rt_csv<4000))
# screen.lmerTest(fmri3, .01)
summary(fmri3)
Anova(fmri3, '3')


# add random slopes
meg1_rs <-  lme4::lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + rt_vmax_lag_sc + last_outcome + 
                                         v_entropy_wi + entropy_change_late_beta_supp + vmax_late_beta + 
                                         omission_early_theta)^2 + 
                            rt_lag_sc:last_outcome:entropy_change_late_beta_supp + 
                            rt_lag_sc:last_outcome:vmax_late_beta +
                            rt_lag_sc:last_outcome:omission_early_theta +
                            rt_vmax_lag_sc:trial_neg_inv_sc:entropy_change_late_beta_supp + 
                            rt_vmax_lag_sc:trial_neg_inv_sc:vmax_late_beta  +
                            rt_vmax_lag_sc:trial_neg_inv_sc:omission_early_theta  +
                            (rt_lag_sc + rt_vmax_lag_sc|id/run), df %>% filter(rt_csv<4000))
# screen.lmerTest(meg1_rs, .01)
summary(meg1_rs)
Anova(meg1_rs, '3')

# fMRI

fmri1_rs <-  lme4::lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + rt_vmax_lag_sc + last_outcome + 
                                          v_entropy_wi + entropy_change_late_beta_supp + vmax_late_beta + 
                                          omission_early_theta)^2 + 
                             rt_lag_sc:last_outcome:entropy_change_late_beta_supp + 
                             rt_lag_sc:last_outcome:vmax_late_beta +
                             rt_lag_sc:last_outcome:omission_early_theta +
                             rt_vmax_lag_sc:trial_neg_inv_sc:entropy_change_late_beta_supp + 
                             rt_vmax_lag_sc:trial_neg_inv_sc:vmax_late_beta  +
                             rt_vmax_lag_sc:trial_neg_inv_sc:omission_early_theta  +
                             (rt_lag_sc + rt_vmax_lag_sc|id/run), fdf %>% filter(rt_csv<4000))
# screen.lmerTest(fmri1_rs, .01)
summary(fmri1_rs)
Anova(fmri1_rs, '3')

# late beta to entropy change
em1 <- as_tibble(emtrends(meg1_rs, data = df,  var = "rt_lag_sc", specs = c("entropy_change_late_beta_supp", "last_outcome"), at = list(entropy_change_late_beta_supp = c(-.14, .12)), options = list()))
em1$study = "1. MEG"
em2 <- as_tibble(emtrends(fmri1_rs, data = fdf, var = "rt_lag_sc", specs = c("entropy_change_late_beta_supp", "last_outcome"), at = list(entropy_change_late_beta_supp = c(-.14, .12)), options = list()))
em2$study = '2. fMRI replication'
em1 <- rbind(em1, em2)
p1 <- ggplot(em1, aes(x=last_outcome, y=rt_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=as.factor(entropy_change_late_beta_supp))) + 
  #shape = as.factor(pe_f2_hipp), 
  #p1 <- ggplot(em1, aes(x=last_outcome, y=rt_lag_sc.trend, linetype = as.factor(pe_f2_hipp), ymin=asymp.LCL, ymax=asymp.UCL)) + 
  geom_point(position = position_dodge(width = .6), size=2.5) + 
  #geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), position = position_dodge2(width=0.9), width=0.5) + 
  geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + 
  #geom_crossbar(position=position_dodge(width=0.7), width=0.6) + 
  theme_bw(base_size=12) + facet_wrap(~study)+ ylab("RT swings (AU)\n Small <---------> Large")  + 
  #scale_shape_manual(values=c(15,16), labels = c("10th %ile", "90th %ile")) + 
  #scale_color_brewer("PH RPE\nresponse", palette="Set1", labels = c("10th %ile", "90th %ile")) +
  scale_color_manual("Entropy change\nlate beta\nresponse", values=c("#1b3840","#4fa3b8"), labels = c("10th %ile", "90th %ile")) +
  labs(shape = "Entropy change\nlate beta\nresponse") +
  theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) + 
  scale_y_reverse(limits = c(.7, -.1)) 

# early theta to reward omission
rem1 <- as_tibble(emtrends(meg1_rs, data = df,  var = "rt_lag_sc", specs = c("omission_early_theta", "last_outcome"), at = list(omission_early_theta = c(-.33, .3)), options = list()))
rem1$study = "1. MEG"
rem2 <- as_tibble(emtrends(fmri1_rs, data = fdf, var = "rt_lag_sc", specs = c("omission_early_theta", "last_outcome"), at = list(omission_early_theta = c(-.33, .3)), options = list()))
rem2$study = '2. fMRI replication'
rem1 <- rbind(rem1, rem2)
rp1 <- ggplot(rem1, aes(x=last_outcome, y=rt_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=as.factor(omission_early_theta))) + 
  #shape = as.factor(pe_f2_hipp), 
  #p1 <- ggplot(em1, aes(x=last_outcome, y=rt_lag_sc.trend, linetype = as.factor(pe_f2_hipp), ymin=asymp.LCL, ymax=asymp.UCL)) + 
  geom_point(position = position_dodge(width = .6), size=2.5) + 
  #geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), position = position_dodge2(width=0.9), width=0.5) + 
  geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + 
  #geom_crossbar(position=position_dodge(width=0.7), width=0.6) + 
  theme_bw(base_size=12) + facet_wrap(~study)+ ylab("RT swings (AU)\n Small <---------> Large")  + 
  #scale_shape_manual(values=c(15,16), labels = c("10th %ile", "90th %ile")) + 
  #scale_color_brewer("PH RPE\nresponse", palette="Set1", labels = c("10th %ile", "90th %ile")) +
  scale_color_manual("Reward omission\nearly theta\nresponse", values=c("#1b3840","#4fa3b8"), labels = c("10th %ile", "90th %ile")) +
  labs(shape = "Reward omission\nearly theta\nresponse") +
  theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) + 
  scale_y_reverse(limits = c(.7, -.1)) 
ggarrange(p1, rp1)
setwd("~/OneDrive/collected_letters/papers/meg/plots/wholebrain")
pdf(file = "MEG_to_behavior_rs.pdf", height = 4, width = 12)
print(ggarrange(p1, rp1))
dev.off()


