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
wbetas <- readRDS("~/OneDrive/collected_letters/papers/meg/plots/wholebrain/betas/MEG_betas_wide_echange_vmax_reward_Nov30_2021.RDS") %>% 
  mutate(omission_early_theta = - avg_reward_early_theta,
         omission_late_delta = - avg_reward_late_delta) %>% 
  mutate(entropy_change_early_beta_supp = -  avg_entropy_change_early_beta,
         entropy_change_late_beta_supp = - avg_entropy_change_late_beta,
         vmax_late_alpha = avg_vmax_late_beta) %>%
  select(c(id, omission_early_theta, omission_late_delta, entropy_change_early_beta_supp, entropy_change_late_beta_supp, vmax_late_alpha))
# merge
df <- df %>% inner_join(wbetas, by = "id")
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
mb_meg1 <-  
  lme4::lmer(rt_csv_sc ~ 
               (trial_neg_inv_sc + rt_lag_sc + rt_vmax_lag_sc + last_outcome + 
                  v_entropy_wi + entropy_change_early_beta_supp + entropy_change_late_beta_supp + vmax_late_alpha + 
                  omission_early_theta + omission_late_delta)^2 + 
               rt_lag_sc:last_outcome:entropy_change_early_beta_supp + 
               rt_lag_sc:last_outcome:entropy_change_late_beta_supp + 
               rt_lag_sc:last_outcome:vmax_late_alpha +
               rt_lag_sc:last_outcome:omission_early_theta +
               rt_lag_sc:last_outcome:omission_late_delta +
               rt_vmax_lag_sc:trial_neg_inv_sc:entropy_change_early_beta_supp + 
               rt_vmax_lag_sc:trial_neg_inv_sc:entropy_change_late_beta_supp + 
               rt_vmax_lag_sc:trial_neg_inv_sc:vmax_late_alpha  +
               rt_vmax_lag_sc:trial_neg_inv_sc:omission_early_theta  +
               rt_vmax_lag_sc:trial_neg_inv_sc:omission_late_delta  +
               (1|id/run), df %>% filter(rt_csv<4000))
# screen.lmerTest(mb_meg1, .01)
summary(mb_meg1)
Anova(mb_meg1, '3')

# load fmri data
fdf <- get_trial_data(repo_directory = clock_folder, dataset = "mmclock_fmri", groupfixed = T)
fdf <- fdf %>% mutate(id = as.character(id)) %>% inner_join(wbetas, by = "id")

fmb_meg1 <- lme4::lmer(rt_csv_sc ~ 
                         (trial_neg_inv_sc + rt_lag_sc + rt_vmax_lag_sc + last_outcome + 
                            v_entropy_wi + entropy_change_early_beta_supp + entropy_change_late_beta_supp + vmax_late_alpha + 
                            omission_early_theta + omission_late_delta)^2 + 
                         rt_lag_sc:last_outcome:entropy_change_early_beta_supp + 
                         rt_lag_sc:last_outcome:entropy_change_late_beta_supp + 
                         rt_lag_sc:last_outcome:vmax_late_alpha +
                         rt_lag_sc:last_outcome:omission_early_theta +
                         rt_lag_sc:last_outcome:omission_late_delta +
                         rt_vmax_lag_sc:trial_neg_inv_sc:entropy_change_early_beta_supp + 
                         rt_vmax_lag_sc:trial_neg_inv_sc:entropy_change_late_beta_supp + 
                         rt_vmax_lag_sc:trial_neg_inv_sc:vmax_late_alpha  +
                         rt_vmax_lag_sc:trial_neg_inv_sc:omission_early_theta  +
                         rt_vmax_lag_sc:trial_neg_inv_sc:omission_late_delta  +
                         (1|id/run), fdf %>% filter(rt_csv<4000))
# screen.lmerTest(fmb_meg1, .01)
summary(fmb_meg1)
Anova(fmb_meg1, '3')

# late beta to entropy change
em1 <- as_tibble(emtrends(mb_meg1, data = df,  var = "rt_lag_sc", specs = c("entropy_change_late_beta_supp", "last_outcome"), at = list(entropy_change_late_beta_supp = c(-.13, .12)), options = list()))
em1$study = 'MEG'
em2 <- as_tibble(emtrends(fmb_meg1, data = fdf, var = "rt_lag_sc", specs = c("entropy_change_late_beta_supp", "last_outcome"), at = list(entropy_change_late_beta_supp = c(-.13, .12)), options = list()))
em2$study = 'fMRI replication'
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
  scale_y_reverse(limits = c(.7, -.1)) 

# early theta to reward omission
rem1 <- as_tibble(emtrends(mb_meg1, data = df,  var = "rt_lag_sc", specs = c("omission_early_theta", "last_outcome"), at = list(omission_early_theta = c(-.33, .3)), options = list()))
rem1$study = 'MEG'
rem2 <- as_tibble(emtrends(fmb_meg1, data = fdf, var = "rt_lag_sc", specs = c("omission_early_theta", "last_outcome"), at = list(omission_early_theta = c(-.33, .3)), options = list()))
rem2$study = 'fMRI replication'
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
  scale_y_reverse(limits = c(.7, -.1)) 
rlem1 <- as_tibble(emtrends(mb_meg1, data = df,  var = "rt_lag_sc", specs = c("omission_late_delta", "last_outcome"), at = list(omission_late_delta = c(-.33, .3)), options = list()))
rlem1$study = 'MEG'
rlem2 <- as_tibble(emtrends(fmb_meg1, data = fdf, var = "rt_lag_sc", specs = c("omission_late_delta", "last_outcome"), at = list(omission_late_delta = c(-.33, .3)), options = list()))
rlem2$study = 'fMRI replication'
rlem1 <- rbind(rlem2, rlem1)
reward_late_delta <- ggplot(rlem1, aes(x=last_outcome, y=rt_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=as.factor(omission_late_delta))) + 
  #shape = as.factor(pe_f2_hipp), 
  #p1 <- ggplot(em1, aes(x=last_outcome, y=rt_lag_sc.trend, linetype = as.factor(pe_f2_hipp), ymin=asymp.LCL, ymax=asymp.UCL)) + 
  geom_point(position = position_dodge(width = .6), size=2.5) + 
  #geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), position = position_dodge2(width=0.9), width=0.5) + 
  geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + 
  #geom_crossbar(position=position_dodge(width=0.7), width=0.6) + 
  theme_bw(base_size=12) + facet_wrap(~study)+ ylab("RT swings (AU)\n Small <---------> Large")  + 
  #scale_shape_manual(values=c(15,16), labels = c("10th %ile", "90th %ile")) + 
  #scale_color_brewer("PH RPE\nresponse", palette="Set1", labels = c("10th %ile", "90th %ile")) +
  scale_color_manual("Reward omission\nlate delta\nsynchronization", values=c("#1b3840","#4fa3b8"), labels = c("10th %ile", "90th %ile")) +
  labs(shape = "Reward omission\nlate delta\nsynchronization") +
  theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) + 
  scale_y_reverse(limits = c(.7, -.1)) 
vem1 <- as_tibble(emtrends(mb_meg1, data = df,  var = "rt_lag_sc", specs = c("vmax_late_alpha", "last_outcome"), at = list(vmax_late_alpha = c(-.08, .08)), options = list()))
vem1$study = 'MEG'
vem2 <- as_tibble(emtrends(fmb_meg1, data = fdf, var = "rt_lag_sc", specs = c("vmax_late_alpha", "last_outcome"), at = list(vmax_late_alpha = c(-.08, .08)), options = list()))
vem2$study = 'fMRI replication'
vem1 <- rbind(vem2, vem1)
vmax_late_alpha <- ggplot(vem1, aes(x=last_outcome, y=rt_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=as.factor(vmax_late_alpha))) + 
  #shape = as.factor(pe_f2_hipp), 
  #p1 <- ggplot(em1, aes(x=last_outcome, y=rt_lag_sc.trend, linetype = as.factor(pe_f2_hipp), ymin=asymp.LCL, ymax=asymp.UCL)) + 
  geom_point(position = position_dodge(width = .6), size=2.5) + 
  #geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), position = position_dodge2(width=0.9), width=0.5) + 
  geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + 
  #geom_crossbar(position=position_dodge(width=0.7), width=0.6) + 
  theme_bw(base_size=12) + facet_wrap(~study)+ ylab("RT swings (AU)\n Small <---------> Large")  + 
  #scale_shape_manual(values=c(15,16), labels = c("10th %ile", "90th %ile")) + 
  #scale_color_brewer("PH RPE\nresponse", palette="Set1", labels = c("10th %ile", "90th %ile")) +
  scale_color_manual("Vmax\nlate alpha\nresponse", values=c("#1b3840","#4fa3b8"), labels = c("10th %ile", "90th %ile")) +
  labs(shape = "Vmax\nlate alpha\nresponse") +
  theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) + 
  scale_y_reverse(limits = c(.7, -.1)) 

ggarrange(ec_late_beta, reward_early_theta, reward_late_delta, vmax_late_alpha)

setwd("~/OneDrive/collected_letters/papers/meg/plots/wholebrain")
pdf(file = "MEG_to_behavior.pdf", height = 7, width = 12)
print(ggarrange(ec_late_beta, reward_early_theta, reward_late_delta, vmax_late_alpha))
dev.off()

################
# Add entropy
################
mb_meg2 <-  lme4::lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + rt_vmax_lag_sc + last_outcome + 
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
# screen.lmerTest(mb_meg2, .01)
summary(mb_meg2)
Anova(mb_meg2, '3')

# load fmri data
fmb_meg2 <-  lme4::lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + rt_vmax_lag_sc + last_outcome + 
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
# screen.lmerTest(fmb_meg2, .01)
summary(fmb_meg2)
Anova(fmb_meg2, '3')

# linear trial
mb_meg3 <-  lme4::lmer(rt_csv_sc ~ (scale(run_trial) + rt_lag_sc + rt_vmax_lag_sc + last_outcome + 
                                      v_entropy_wi + entropy_change_late_beta_supp + vmax_late_beta + 
                                      omission_early_theta)^2 + 
                         rt_lag_sc:last_outcome:entropy_change_late_beta_supp + 
                         rt_lag_sc:last_outcome:vmax_late_beta +
                         rt_lag_sc:last_outcome:omission_early_theta +
                         rt_vmax_lag_sc:scale(run_trial):entropy_change_late_beta_supp + 
                         rt_vmax_lag_sc:scale(run_trial):vmax_late_beta  +
                         rt_vmax_lag_sc:scale(run_trial):omission_early_theta  +
                         (1|id/run), df %>% filter(rt_csv<4000))
# screen.lmerTest(mb_meg3, .01)
summary(mb_meg3)
Anova(mb_meg3, '3')

# fMRI
fmb_meg3 <-  lme4::lmer(rt_csv_sc ~ (scale(run_trial) + rt_lag_sc + rt_vmax_lag_sc + last_outcome + 
                                       v_entropy_wi + entropy_change_late_beta_supp + vmax_late_beta + 
                                       omission_early_theta)^2 + 
                          rt_lag_sc:last_outcome:entropy_change_late_beta_supp + 
                          rt_lag_sc:last_outcome:vmax_late_beta +
                          rt_lag_sc:last_outcome:omission_early_theta +
                          rt_vmax_lag_sc:scale(run_trial):entropy_change_late_beta_supp + 
                          rt_vmax_lag_sc:scale(run_trial):vmax_late_beta  +
                          rt_vmax_lag_sc:scale(run_trial):omission_early_theta  +
                          (1|id/run), fdf %>% filter(rt_csv<4000))
# screen.lmerTest(fmb_meg3, .01)
summary(fmb_meg3)
Anova(fmb_meg3, '3')


# add random slopes
mb_meg1_rs <-  lme4::lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + rt_vmax_lag_sc + last_outcome + 
                                         v_entropy_wi + entropy_change_late_beta_supp + vmax_late_beta + 
                                         omission_early_theta)^2 + 
                            rt_lag_sc:last_outcome:entropy_change_late_beta_supp + 
                            rt_lag_sc:last_outcome:vmax_late_beta +
                            rt_lag_sc:last_outcome:omission_early_theta +
                            rt_vmax_lag_sc:trial_neg_inv_sc:entropy_change_late_beta_supp + 
                            rt_vmax_lag_sc:trial_neg_inv_sc:vmax_late_beta  +
                            rt_vmax_lag_sc:trial_neg_inv_sc:omission_early_theta  +
                            (rt_lag_sc + rt_vmax_lag_sc|id/run), df %>% filter(rt_csv<4000))
# screen.lmerTest(mb_meg1_rs, .01)
summary(mb_meg1_rs)
Anova(mb_meg1_rs, '3')

# fMRI

fmb_meg1_rs <-  lme4::lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + rt_vmax_lag_sc + last_outcome + 
                                          v_entropy_wi + entropy_change_late_beta_supp + vmax_late_beta + 
                                          omission_early_theta)^2 + 
                             rt_lag_sc:last_outcome:entropy_change_late_beta_supp + 
                             rt_lag_sc:last_outcome:vmax_late_beta +
                             rt_lag_sc:last_outcome:omission_early_theta +
                             rt_vmax_lag_sc:trial_neg_inv_sc:entropy_change_late_beta_supp + 
                             rt_vmax_lag_sc:trial_neg_inv_sc:vmax_late_beta  +
                             rt_vmax_lag_sc:trial_neg_inv_sc:omission_early_theta  +
                             (rt_lag_sc + rt_vmax_lag_sc|id/run), fdf %>% filter(rt_csv<4000))
# screen.lmerTest(fmb_meg1_rs, .01)
summary(fmb_meg1_rs)
Anova(fmb_meg1_rs, '3')

# late beta to entropy change
em1 <- as_tibble(emtrends(mb_meg1_rs, data = df,  var = "rt_lag_sc", specs = c("entropy_change_late_beta_supp", "last_outcome"), at = list(entropy_change_late_beta_supp = c(-.14, .12)), options = list()))
em1$study = 'MEG'
em2 <- as_tibble(emtrends(fmb_meg1_rs, data = fdf, var = "rt_lag_sc", specs = c("entropy_change_late_beta_supp", "last_outcome"), at = list(entropy_change_late_beta_supp = c(-.14, .12)), options = list()))
em2$study = 'fMRI replication'
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
rem1 <- as_tibble(emtrends(mb_meg1_rs, data = df,  var = "rt_lag_sc", specs = c("omission_early_theta", "last_outcome"), at = list(omission_early_theta = c(-.33, .3)), options = list()))
rem1$study = 'MEG'
rem2 <- as_tibble(emtrends(fmb_meg1_rs, data = fdf, var = "rt_lag_sc", specs = c("omission_early_theta", "last_outcome"), at = list(omission_early_theta = c(-.33, .3)), options = list()))
rem2$study = 'fMRI replication'
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


