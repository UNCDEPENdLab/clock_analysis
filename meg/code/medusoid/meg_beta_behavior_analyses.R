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
plots <- F
# source('~/code/Rhelpers/')
setwd(file.path(clock_folder, 'fmri/keuka_brain_behavior_analyses/dan'))


perform_checks = F # initial checks ensuring that the extracted betas are valid

### load data
source("~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/get_trial_data.R")
#get_trial_data <- function(repo_directory=NULL, dataset="mmclock_fmri", groupfixed=TRUE) 
df <- get_trial_data(repo_directory = clock_folder, dataset = "mmclock_meg", groupfixed = T)

# add meg data
# wbetas <- readRDS("~/OneDrive/collected_letters/papers/meg/plots/wholebrain/betas/MEG_betas_wide_echange_vmax_reward_Nov30_2021.RDS") %>% 
wbetas <- readRDS("~/code/clock_analysis/meg/data/MEG_betas_entropy_change_v_max_reward_abs_pe_Mar_14_2022.RDS") %>% 
  mutate(omission_early_theta = - reward_early_theta,
         omission_late_delta = - reward_late_delta) %>% 
  mutate(entropy_change_early_beta_supp = -  entropy_change_early_beta_entropy_change,
         entropy_change_late_beta_supp = - entropy_change_late_beta_entropy_change,
         abspe_late_beta_supp = - abspe_late_beta
         # ,
         # entropy_change_early_beta_supp_ec = -  entropy_change_early_beta_entropy_change_ec_sensors,
         # entropy_change_late_beta_supp_ec = - entropy_change_late_beta_entropy_change_ec_sensors,
         # abs_pe_late_beta_supp_ec = - abspe_ec_late_beta,
         # abs_pe_late_beta_supp = - abspe_late_beta,
         # vmax_late_alpha = vmax_late_beta
  ) %>%
  select(c(id, omission_early_theta, omission_late_delta, 
           entropy_change_early_beta_supp, entropy_change_late_beta_supp, abspe_late_beta_supp
           # ,
           # entropy_change_early_beta_supp_ec, entropy_change_late_beta_supp_ec, 
           # vmax_late_alpha, 
           # abs_pe_late_beta_supp_ec, abs_pe_late_beta_supp
  ))
# merge
df <- df %>% inner_join(wbetas, by = "id")


if (perform_checks) {

#############
# Check correlations with earnings in learnable blocks: not the slightest
sdf <- df %>% group_by(id, run) %>% filter(rewFunc=="IEV" | rewFunc=="DEV") %>% ungroup() %>%
  select(id,  score_csv, rewFunc) %>% group_by(id, rewFunc) %>% summarize(run_earnings = mean(score_csv), .groups = "keep") %>% unique() %>% inner_join(wbetas)
car::Anova(lm(run_earnings ~ rewFunc * omission_early_theta + rewFunc * omission_late_delta + rewFunc * abspe_late_beta_supp + rewFunc * entropy_change_early_beta_supp + 
                rewFunc * entropy_change_late_beta_supp, sdf))
               # rewFunc * entropy_change_late_beta_supp + rewFunc * vmax_late_alpha, sdf))
# Performance: EV as outcome
ldf <- df %>% filter(rewFunc=="IEV" | rewFunc=="DEV")
ev_meg2 <-  
  lmerTest::lmer(ev ~ 
               (trial_neg_inv_sc + entropy_change_late_beta_supp + abspe_late_beta_supp +
                  omission_early_theta + rewFunc)^2 +
               (1|id), ldf %>% filter(rt_csv<4000))
screen.lmerTest(ev_meg2, .01)
summary(ev_meg2)
Anova(ev_meg2, '3')
}

if (perform_checks) {
  #############
  # Check correlations with earnings in learnable blocks: not the slightest
  sdf <- df %>% group_by(id, run) %>% filter(rewFunc=="IEV" | rewFunc=="DEV") %>% ungroup() %>%
    select(id,  score_csv, rewFunc) %>% group_by(id, rewFunc) %>% summarize(run_earnings = mean(score_csv), .groups = "keep") %>% unique() %>% inner_join(wbetas)
  car::Anova(lm(run_earnings ~ rewFunc * omission_early_theta + rewFunc * omission_late_delta + rewFunc * abspe_late_beta_supp + rewFunc * entropy_change_early_beta_supp + 
                  rewFunc * entropy_change_late_beta_supp, sdf))
  # rewFunc * entropy_change_late_beta_supp + rewFunc * vmax_late_alpha, sdf))
  # Performance: EV as outcome
  ldf <- df %>% filter(rewFunc=="IEV" | rewFunc=="DEV")
  ev_meg2 <-  
    lmerTest::lmer(ev ~ 
                     (trial_neg_inv_sc + entropy_change_late_beta_supp + abspe_late_beta_supp +
                        omission_early_theta + rewFunc)^2 +
                     (1|id), ldf %>% filter(rt_csv<4000))
  screen.lmerTest(ev_meg2, .01)
  summary(ev_meg2)
  Anova(ev_meg2, '3')
  
  ev_meg3 <-  
    lmerTest::lmer(ev ~ 
                     (trial_neg_inv_sc + entropy_change_late_beta_supp + 
                        omission_early_theta + rewFunc)^3 +
                     (1|id), ldf %>% filter(rt_csv<4000))
  screen.lmerTest(ev_meg3, .01)
  summary(ev_meg3)
  Anova(ev_meg3, '3')
  anova(ev_meg2, ev_meg3)
  
  # Condition effects:
  cond_meg2 <-  
    lmerTest::lmer(rt_csv ~ 
                     (trial_neg_inv_sc + entropy_change_late_beta_supp + 
                        omission_early_theta + rewFunc)^2 +
                     (1|id), df %>% filter(rt_csv<4000))
  screen.lmerTest(cond_meg2)
  summary(cond_meg2)
  Anova(cond_meg2, '3')
  
  cond_meg3 <-  
    lmerTest::lmer(ev ~ 
                     (trial_neg_inv_sc + entropy_change_late_beta_supp + 
                        omission_early_theta + rewFunc)^3 +
                     (1|id), ldf %>% filter(rt_csv<4000))
  screen.lmerTest(cond_meg3, .01)
  summary(cond_meg3)
  Anova(cond_meg2, '3')
  anova(cond_meg2, cond_meg3)
  
  ggplot(ldf, aes(trial_neg_inv_sc, ev, lty = omission_early_theta>0, color = entropy_change_late_beta_supp>0)) + geom_smooth(method = "gam") + facet_wrap(~rewFunc)
  ggplot(ldf, aes(trial_neg_inv_sc, rt_csv, lty = omission_early_theta>0, color = entropy_change_late_beta_supp>0)) + geom_smooth(method = "gam") + facet_wrap(~rewFunc)
  ggplot(ldf, aes(trial_neg_inv_sc, rt_swing, lty = omission_early_theta>0, color = entropy_change_late_beta_supp>0)) + geom_smooth(method = "gam") + facet_wrap(~rewFunc)
  ggplot(ldf, aes(trial_neg_inv_sc, abs(rt_csv - rt_vmax), lty = omission_early_theta>0, color = entropy_change_late_beta_supp>0)) + geom_smooth(method = "gam") + facet_wrap(~rewFunc)
  ggplot(df, aes(trial_neg_inv_sc, abs(rt_csv - rt_vmax), lty = omission_early_theta>0, color = entropy_change_late_beta_supp>0)) + geom_smooth(method = "gam") + facet_wrap(~rewFunc)
  
  
  ggplot(df %>% filter(!is.na(reward_lag)), aes(trial_neg_inv_sc, rt_swing, lty = omission_early_theta>0, color = rewFunc)) + geom_smooth(method = "gam") + facet_wrap(~reward_lag)
  ggplot(df %>% filter(!is.na(reward_lag) & abs(omission_early_theta)>.3), aes(rt_swing, color = omission_early_theta>0)) + geom_boxplot(notch = T, varwidth = T) + facet_grid(rewFunc~reward_lag)
  
  ggplot(df %>% filter (!is.na(reward_lag)), aes(omission_early_theta, rt_swing, lty = reward_lag)) + geom_smooth(method = "gam") + facet_grid(~rewFunc)
  ggplot(ldf %>% filter (!is.na(reward_lag)), aes(omission_early_theta, ev, lty = reward_lag)) + geom_smooth(method = "gam") + facet_grid(rt_lag<2~rewFunc)
}



ggplot(df %>% filter(!is.na(reward_lag)), aes(trial_neg_inv_sc, rt_swing, lty = omission_early_theta>0, color = rewFunc)) + geom_smooth(method = "gam") + facet_wrap(~reward_lag)
ggplot(df %>% filter(!is.na(reward_lag) & abs(omission_early_theta)>.3), aes(rt_swing, color = omission_early_theta>0)) + geom_boxplot(notch = T, varwidth = T) + facet_grid(rewFunc~reward_lag)

ggplot(df %>% filter (!is.na(reward_lag)), aes(omission_early_theta, rt_swing, lty = reward_lag)) + geom_smooth(method = "gam") + facet_grid(~rewFunc)
ggplot(ldf %>% filter (!is.na(reward_lag)), aes(omission_early_theta, ev, lty = reward_lag)) + geom_smooth(method = "gam") + facet_grid(rt_lag<2~rewFunc)


# ggplot(fdf %>% filter (!is.na(reward_lag)), aes(omission_early_theta, rt_swing, lty = reward_lag)) + geom_smooth(method = "gam") + facet_grid(rt_lag<2~rewFunc)


# ev_meg3 <-  
#   lmerTest::lmer(ev ~ 
#                (trial_neg_inv_sc + last_outcome + rewFunc + entropy_change_early_beta_supp)^3 +
#                (trial_neg_inv_sc + last_outcome + rewFunc + entropy_change_late_beta_supp)^3 +
#                (trial_neg_inv_sc + last_outcome + rewFunc + vmax_late_alpha)^3 +
#                (trial_neg_inv_sc + last_outcome + rewFunc + omission_early_theta)^3 +
#                (trial_neg_inv_sc + last_outcome + rewFunc + omission_late_delta)^3 +
#                (1|id/run), ldf %>% filter(rt_csv<4000))
# screen.lmerTest(ev_meg3, .01)
# summary(ev_meg3)
# Anova(ev_meg3, '3')
# 
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
# mdf$session <- "1. MEG"
# bdf <- bind_rows(df, mdf)
# ggplot(bdf, aes(run_trial, rt_csv, color = rewFunc, lty = session)) + geom_smooth()
# ggplot(bdf, aes(run_trial, ev, color = rewFunc, lty = session)) + geom_smooth()

############# MEG

# Effect of MEG betas on exploration
# because post-omission theta and delta are correlated at 0.75, will test in different models
mb_meg1theta <-  
  lme4::lmer(rt_csv_sc ~ 
               # rt_lag_sc * last_outcome * entropy_change_early_beta_supp + 
               rt_lag_sc * last_outcome * entropy_change_late_beta_supp + 
               # rt_lag_sc * last_outcome * vmax_late_alpha +
               rt_lag_sc * last_outcome * omission_early_theta +
               # rt_vmax_lag_sc * trial_neg_inv_sc * entropy_change_early_beta_supp + 
               rt_vmax_lag_sc * trial_neg_inv_sc * entropy_change_late_beta_supp + 
               # rt_vmax_lag_sc * trial_neg_inv_sc * vmax_late_alpha  +
               rt_vmax_lag_sc * trial_neg_inv_sc * omission_early_theta  +
               (1|id/run), df %>% filter(rt_csv<4000))
# screen.lmerTest(mb_meg1, .01)
summary(mb_meg1theta)
Anova(mb_meg1theta, '3')

# abspe_late_beta has an independent pro-RT_lag effect, but weaker than entropy_change_late_beta in both MEG and fMRI
mb_meg1theta1 <-  
  lme4::lmer(rt_csv_sc ~ 
               # rt_lag_sc * last_outcome * entropy_change_early_beta_supp + 
               rt_lag_sc * last_outcome * entropy_change_late_beta_supp + 
               # rt_lag_sc * last_outcome * vmax_late_alpha +
               rt_lag_sc * last_outcome * omission_early_theta +
               rt_lag_sc * last_outcome * abspe_late_beta_supp +
               # rt_vmax_lag_sc * trial_neg_inv_sc * entropy_change_early_beta_supp + 
               rt_vmax_lag_sc * trial_neg_inv_sc * entropy_change_late_beta_supp + 
               rt_vmax_lag_sc * trial_neg_inv_sc * abspe_late_beta_supp + 
               # rt_vmax_lag_sc * trial_neg_inv_sc * vmax_late_alpha  +
               rt_vmax_lag_sc * trial_neg_inv_sc * omission_early_theta  +
               (1|id/run), df %>% filter(rt_csv<4000))
# screen.lmerTest(mb_meg1, .01)
summary(mb_meg1theta1)
Anova(mb_meg1theta1, '3')


mb_meg1theta1_rewFunc <-  
  lme4::lmer(rt_csv_sc ~ 
               # rt_lag_sc * last_outcome * entropy_change_early_beta_supp + 
               rt_lag_sc * last_outcome * entropy_change_late_beta_supp * rewFunc + 
               # rt_lag_sc * last_outcome * vmax_late_alpha +
               rt_lag_sc * last_outcome * omission_early_theta * rewFunc +
               rt_lag_sc * last_outcome * abspe_late_beta_supp * rewFunc +
               # rt_vmax_lag_sc * trial_neg_inv_sc * entropy_change_early_beta_supp + 
               rt_vmax_lag_sc * trial_neg_inv_sc * entropy_change_late_beta_supp * rewFunc + 
               rt_vmax_lag_sc * trial_neg_inv_sc * abspe_late_beta_supp * rewFunc + 
               # rt_vmax_lag_sc * trial_neg_inv_sc * vmax_late_alpha  +
               rt_vmax_lag_sc * trial_neg_inv_sc * omission_early_theta * rewFunc  +
               (1|id/run), df %>% filter(rt_csv<4000))
# screen.lmerTest(mb_meg1, .01)
summary(mb_meg1theta1_rewFunc)
Anova(mb_meg1theta1_rewFunc, '3')


# mb_meg1delta <-  
#   lme4::lmer(rt_csv_sc ~ 
#                rt_lag_sc * last_outcome * entropy_change_early_beta_supp + 
#                rt_lag_sc * last_outcome * entropy_change_late_beta_supp + 
#                rt_lag_sc * last_outcome * vmax_late_alpha +
#                rt_lag_sc * last_outcome * omission_late_delta +
#                rt_vmax_lag_sc * trial_neg_inv_sc * entropy_change_early_beta_supp + 
#                rt_vmax_lag_sc * trial_neg_inv_sc * entropy_change_late_beta_supp + 
#                rt_vmax_lag_sc * trial_neg_inv_sc * vmax_late_alpha  +
#                rt_vmax_lag_sc * trial_neg_inv_sc * omission_late_delta  +
#                (1|id/run), df %>% filter(rt_csv<4000))
# # screen.lmerTest(mb_meg1, .01)
# summary(mb_meg1delta)
# Anova(mb_meg1delta, '3')

mb_meg_ec_allsensors <-  
  lme4::lmer(rt_csv_sc ~ 
               rt_lag_sc * last_outcome * entropy_change_early_beta_supp + 
               rt_lag_sc * last_outcome * entropy_change_late_beta_supp + 
               rt_vmax_lag_sc * trial_neg_inv_sc * entropy_change_early_beta_supp + 
               rt_vmax_lag_sc * trial_neg_inv_sc * entropy_change_late_beta_supp + 
               (1|id/run), df %>% filter(rt_csv<4000))
# screen.lmerTest(mb_meg1, .01)
summary(mb_meg_ec_allsensors)
Anova(mb_meg_ec_allsensors, '3')

mb_meg_ec_late_only <-  
  lme4::lmer(rt_csv_sc ~ 
               rt_lag_sc * last_outcome * entropy_change_late_beta_supp + 
               rt_vmax_lag_sc * trial_neg_inv_sc * entropy_change_late_beta_supp + 
               (1|id/run), df %>% filter(rt_csv<4000))
# screen.lmerTest(mb_meg1, .01)
summary(mb_meg_ec_late_only)
Anova(mb_meg_ec_late_only, '3')

mb_meg_ec_late_only_rs <-  
  lme4::lmer(rt_csv_sc ~ 
               rt_lag_sc * last_outcome * entropy_change_late_beta_supp + 
               rt_vmax_lag_sc * trial_neg_inv_sc * entropy_change_late_beta_supp + 
               (1 + rt_lag_sc + rt_vmax_lag_sc + last_outcome |id) + (1|id:run), df %>% filter(rt_csv<4000))
# screen.lmerTest(mb_meg1, .01)
while (!is.null(mb_meg_ec_late_only_rs@optinfo$conv$lme4$messages)) {
  ss <- getME(mb_meg_ec_late_only_rs,c("theta","fixef"))
  mb_meg_ec_late_only_rs <- update(mb_meg_ec_late_only_rs, start=ss, control=lmerControl(optimizer = "bobyqa",optCtr=list(maxfun=2e5)))}

summary(mb_meg_ec_late_only_rs)
Anova(mb_meg_ec_late_only_rs, '3')


##############
# fMRI out-of-session replication
# load fmri data frame (fdf)
fdf <- get_trial_data(repo_directory = clock_folder, dataset = "mmclock_fmri", groupfixed = T)
fdf <- fdf %>% mutate(id = as.character(id)) %>% inner_join(wbetas, by = "id")

if (perform_checks) {
  # select learnable conditions for EV analyses
  lfdf <- fdf %>% filter(rewFunc=="IEV" | rewFunc=="DEV")
  ev_fmri2 <-  
    lmerTest::lmer(ev ~ 
                     (trial_neg_inv_sc + entropy_change_late_beta_supp + 
                        omission_early_theta + rewFunc)^2 +
                     (1|id), lfdf %>% filter(rt_csv<4000))
  screen.lmerTest(ev_fmri2, .01)
  summary(ev_fmri2)
  Anova(ev_fmri2, '3')
  
  ev_fmri3 <-  
    lmerTest::lmer(ev ~ 
                     (trial_neg_inv_sc + entropy_change_late_beta_supp + 
                        omission_early_theta + rewFunc)^3 +
                     (1|id), lfdf %>% filter(rt_csv<4000))
  screen.lmerTest(ev_fmri3, .01)
  summary(ev_fmri3)
  Anova(ev_fmri3, '3')
  anova(ev_fmri2, ev_fmri3)
  
  # Condition effects:
  cond_fmri2 <-  
    lmerTest::lmer(rt_csv ~ 
                     (trial_neg_inv_sc + entropy_change_late_beta_supp + 
                        omission_early_theta + rewFunc)^2 +
                     (1|id), fdf %>% filter(rt_csv<4000))
  screen.lmerTest(cond_fmri2)
  summary(cond_fmri2)
  Anova(cond_fmri2, '3')
  
  cond_fmri3 <-  
    lmerTest::lmer(ev ~ 
                     (trial_neg_inv_sc + entropy_change_late_beta_supp + 
                        omission_early_theta + rewFunc)^3 +
                     (1|id), fdf %>% filter(rt_csv<4000))
  screen.lmerTest(cond_fmri3, .01)
  summary(cond_fmri3)
  Anova(ev_fmri2, '3')
  anova(ev_fmri2, ev_fmri3)
}

fmb_meg1theta <- lme4::lmer(rt_csv_sc ~ 
                              rt_lag_sc * last_outcome * entropy_change_early_beta_supp + 
                              rt_lag_sc * last_outcome * entropy_change_late_beta_supp + 
                              # rt_lag_sc * last_outcome * vmax_late_alpha +
                              rt_lag_sc * last_outcome * omission_early_theta +
                              rt_vmax_lag_sc * trial_neg_inv_sc * entropy_change_early_beta_supp + 
                              rt_vmax_lag_sc * trial_neg_inv_sc * entropy_change_late_beta_supp + 
                              # rt_vmax_lag_sc * trial_neg_inv_sc * vmax_late_alpha  +
                              rt_vmax_lag_sc * trial_neg_inv_sc * omission_early_theta  +
                              (1|id/run), fdf %>% filter(rt_csv<4000))
# screen.lmerTest(fmb_meg1, .01)
summary(fmb_meg1theta)
Anova(fmb_meg1theta, '3')

fmb_meg1theta1 <-  
  lme4::lmer(rt_csv_sc ~ 
               # rt_lag_sc * last_outcome * entropy_change_early_beta_supp + 
               rt_lag_sc * last_outcome * entropy_change_late_beta_supp + 
               # rt_lag_sc * last_outcome * vmax_late_alpha +
               rt_lag_sc * last_outcome * omission_early_theta +
               rt_lag_sc * last_outcome * abspe_late_beta_supp +
               # rt_vmax_lag_sc * trial_neg_inv_sc * entropy_change_early_beta_supp + 
               rt_vmax_lag_sc * trial_neg_inv_sc * entropy_change_late_beta_supp + 
               rt_vmax_lag_sc * trial_neg_inv_sc * abspe_late_beta_supp + 
               # rt_vmax_lag_sc * trial_neg_inv_sc * vmax_late_alpha  +
               rt_vmax_lag_sc * trial_neg_inv_sc * omission_early_theta  +
               (1|id/run), fdf %>% filter(rt_csv<4000))
# screen.lmerTest(mb_meg1, .01)
summary(fmb_meg1theta1)
Anova(fmb_meg1theta1, '3')


fmb_fmri_ec_late_only_rs <-  
  lme4::lmer(rt_csv_sc ~ 
               rt_lag_sc * last_outcome * entropy_change_late_beta_supp + 
               rt_vmax_lag_sc * trial_neg_inv_sc * entropy_change_late_beta_supp + 
               (1 + rt_lag_sc + rt_vmax_lag_sc + last_outcome|id) + (1|id:run), fdf %>% filter(rt_csv<4000))
while (!is.null(fmb_fmri_ec_late_only_rs@optinfo$conv$lme4$messages)) {
  ss <- getME(fmb_fmri_ec_late_only_rs,c("theta","fixef"))
  fmb_fmri_ec_late_only_rs <- update(fmb_fmri_ec_late_only_rs, start=ss, control=lmerControl(optimizer = "bobyqa",optCtr=list(maxfun=2e5)))}
summary(fmb_fmri_ec_late_only_rs)
Anova(fmb_fmri_ec_late_only_rs, '3')
Anova(fmb_fmri_ec_late_only_rs, '2')

fmb_fmri_ec_early_late_rs <-  
  lme4::lmer(rt_csv_sc ~ 
               rt_lag_sc * last_outcome * entropy_change_early_beta_supp + 
               rt_lag_sc * last_outcome * entropy_change_late_beta_supp + 
               rt_vmax_lag_sc * trial_neg_inv_sc * entropy_change_early_beta_supp + 
               rt_vmax_lag_sc * trial_neg_inv_sc * entropy_change_late_beta_supp + 
               (1 + rt_lag_sc + rt_vmax_lag_sc + last_outcome|id) + (1|id:run), fdf %>% filter(rt_csv<4000))
while (!is.null(fmb_fmri_ec_early_late_rs@optinfo$conv$lme4$messages)) {
  ss <- getME(fmb_fmri_ec_early_late_rs,c("theta","fixef"))
  fmb_fmri_ec_early_late_rs <- update(fmb_fmri_ec_early_late_rs, start=ss, control=lmerControl(optimizer = "bobyqa",optCtr=list(maxfun=2e5)))}
summary(fmb_fmri_ec_early_late_rs)
Anova(fmb_fmri_ec_early_late_rs, '3')
Anova(fmb_fmri_ec_early_late_rs, '2')



if (plots) {
  # late beta to entropy change
  em1 <- as_tibble(emtrends(mb_meg1theta, data = df,  var = "rt_lag_sc", specs = c("entropy_change_late_beta_supp", "last_outcome"), at = list(entropy_change_late_beta_supp = c(-.07, .22)), options = list()))
  # em1 <- as_tibble(emtrends(mb_meg_ec_late_only_rs, data = df,  var = "rt_lag_sc", specs = c("entropy_change_late_beta_supp", "last_outcome"), at = list(entropy_change_late_beta_supp = c(-.07, .22)), options = list()))
  em1$study = "1. MEG"
  em2 <- as_tibble(emtrends(fmb_meg1theta, data = fdf, var = "rt_lag_sc", specs = c("entropy_change_late_beta_supp", "last_outcome"), at = list(entropy_change_late_beta_supp = c(-.07, .22)), options = list()))
  # em2 <- as_tibble(emtrends(fmb_fmri_ec_late_only_rs, data = df,  var = "rt_lag_sc", specs = c("entropy_change_late_beta_supp", "last_outcome"), at = list(entropy_change_late_beta_supp = c(-.07, .22)), options = list()))
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
  rem1 <- as_tibble(emtrends(mb_meg1theta, data = df,  var = "rt_lag_sc", specs = c("omission_early_theta", "last_outcome"), at = list(omission_early_theta = c(-.31, .27)), options = list()))
  rem1$study = "1. MEG"
  rem2 <- as_tibble(emtrends(fmb_meg1theta, data = fdf, var = "rt_lag_sc", specs = c("omission_early_theta", "last_outcome"), at = list(omission_early_theta = c(-.31, .27)), options = list()))
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
    scale_y_reverse(limits = c(.7, -.1)) 
  rlem1 <- as_tibble(emtrends(mb_meg1delta, data = df,  var = "rt_lag_sc", specs = c("omission_late_delta", "last_outcome"), at = list(omission_late_delta = c(-.3, .31)), options = list()))
  rlem1$study = "1. MEG"
  rlem2 <- as_tibble(emtrends(fmb_meg1delta, data = fdf, var = "rt_lag_sc", specs = c("omission_late_delta", "last_outcome"), at = list(omission_late_delta = c(-.3, .31)), options = list()))
  rlem2$study = '2. fMRI replication'
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
  vem1 <- as_tibble(emtrends(mb_meg1theta, data = df,  var = "rt_lag_sc", specs = c("vmax_late_alpha", "last_outcome"), at = list(vmax_late_alpha = c(-.08, .08)), options = list()))
  vem1$study = "1. MEG"
  vem2 <- as_tibble(emtrends(fmb_meg1theta, data = fdf, var = "rt_lag_sc", specs = c("vmax_late_alpha", "last_outcome"), at = list(vmax_late_alpha = c(-.08, .08)), options = list()))
  vem2$study = '2. fMRI replication'
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
  
  # absolute PE late beta suppression effects
  # late beta to entropy change
  am1 <- as_tibble(emtrends(mb_meg_ec_abspe_ecsensors, data = df,  var = "rt_lag_sc", specs = c("abs_pe_late_beta_supp_ec", "last_outcome"), at = list(abs_pe_late_beta_supp_ec = c(-.01, .01)), options = list()))
  am1$study = "1. MEG"
  am2 <- as_tibble(emtrends(mb_fmri_ec_abspe_ecsensors, data = fdf, var = "rt_lag_sc", specs = c("abs_pe_late_beta_supp_ec", "last_outcome"), at = list(abs_pe_late_beta_supp_ec = c(-.01, .01)), options = list()))
  am2$study = '2. fMRI replication'
  am1 <- rbind(am2, am1)
  abspe_late_beta_ec <- ggplot(am1, aes(x=last_outcome, y=rt_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=as.factor(abs_pe_late_beta_supp_ec))) + 
    #shape = as.factor(pe_f2_hipp), 
    #p1 <- ggplot(am1, aes(x=last_outcome, y=rt_lag_sc.trend, linetype = as.factor(pe_f2_hipp), ymin=asymp.LCL, ymax=asymp.UCL)) + 
    geom_point(position = position_dodge(width = .6), size=2.5) + 
    #geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), position = position_dodge2(width=0.9), width=0.5) + 
    geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + 
    #geom_crossbar(position=position_dodge(width=0.7), width=0.6) + 
    theme_bw(base_size=12) + facet_wrap(~study)+ ylab("RT swings (AU)\n Small <---------> Large")  + 
    #scale_shape_manual(values=c(15,16), labels = c("10th %ile", "90th %ile")) + 
    #scale_color_brewer("PH RPE\nresponse", palette="Set1", labels = c("10th %ile", "90th %ile")) +
    scale_color_manual("Absolute PE\nlate beta\nsuppression", values=c("#1b3840","#4fa3b8"), labels = c("10th %ile", "90th %ile")) +
    labs(shape = "Absolute PE\nlate beta\nsuppression") +
    theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
          axis.text=element_text(size=8.5, color="grey10")) + 
    scale_y_reverse(limits = c(.7, -.1)) 
  
  
  ggarrange(ec_late_beta, abspe_late_beta_ec, reward_early_theta, reward_late_delta, vmax_late_alpha)
  
  setwd("~/OneDrive/collected_letters/papers/meg/plots/wholebrain")
  pdf(file = "MEG_to_behavior.pdf", height = 7, width = 12)
  print(ggarrange(ec_late_beta, abspe_late_beta_ec, reward_early_theta))
  dev.off()
}
##################
# RT Vmax effects 
#################
# 
# # not replicable or interpretable for late beta suppression
# emodel_condition <- as_tibble(emtrends(fmb_meg1delta, var = "rt_vmax_lag_sc", specs = c("omission_late_delta", "trial_neg_inv_sc"), at = list(omission_late_delta = c(-.3, .31), trial_neg_inv_sc = c(-.88, 0.39)), options = list()))
# emodel_condition$study = '2. fMRI replication'
# emodel_RT_Vmax <- as_tibble(emtrends(mb_meg1delta, var = "rt_vmax_lag_sc", specs = c("omission_late_delta", "trial_neg_inv_sc"), at = list(omission_late_delta = c(-.3, .31), trial_neg_inv_sc = c(-.88, 0.39)), options = list()))
# emodel_RT_Vmax$study = "1. MEG"
# em4 <- rbind(emodel_condition, emodel_RT_Vmax)
# omission_late_delta_rtvmax <- ggplot(em4, aes(x=as.factor(trial_neg_inv_sc), y=rt_vmax_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=as.factor(omission_late_delta))) + 
#   geom_point(position = position_dodge(width = .6), size=2.5) + 
#   geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + 
#   theme_bw(base_size=12) + facet_wrap(~study)+  ylab("Convergence on\nbest RT (AU)") +
#   scale_color_manual("Reward omission\nlate delta\nsynchronization", values=c("#403202", "#e2b407"), labels = c("10th %ile", "90th %ile")) +
#   labs(shape = "Reward omission\nlate delta\nsynchronization") +
#   theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
#         axis.text=element_text(size=8.5, color="grey10")) + 
#   scale_x_discrete(name ="Trial", labels=c("-0.88" = "5", "0.39" = "50")) + scale_y_continuous(limits = c(.0, .25))
# 
# setwd("~/OneDrive/collected_letters/papers/meg/plots/wholebrain")
# pdf(file = "MEG_to_behavior_late_delta_RT_Vmax.pdf", height = 3, width = 6)
# print(omission_late_delta_rtvmax)
# dev.off()
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
                                         v_entropy_wi + entropy_change_late_beta_supp + #vmax_late_beta + 
                                         omission_early_theta)^2 + 
                            rt_lag_sc:last_outcome:entropy_change_late_beta_supp + 
                            # rt_lag_sc:last_outcome:vmax_late_beta +
                            rt_lag_sc:last_outcome:omission_early_theta +
                            rt_vmax_lag_sc:trial_neg_inv_sc:entropy_change_late_beta_supp + 
                            # rt_vmax_lag_sc:trial_neg_inv_sc:vmax_late_beta  +
                            rt_vmax_lag_sc:trial_neg_inv_sc:omission_early_theta  +
                            (rt_lag_sc + rt_vmax_lag_sc + last_outcome|id) + (1|id:run), df %>% filter(rt_csv<4000))
# screen.lmerTest(mb_meg1_rs, .01)
summary(mb_meg1_rs)
Anova(mb_meg1_rs, '3')

# fMRI

fmb_meg1_rs <-  lme4::lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + rt_vmax_lag_sc + last_outcome + 
                                          v_entropy_wi + entropy_change_late_beta_supp + #vmax_late_beta + 
                                          omission_early_theta)^2 + 
                             rt_lag_sc:last_outcome:entropy_change_late_beta_supp + 
                             # rt_lag_sc:last_outcome:vmax_late_beta +
                             rt_lag_sc:last_outcome:omission_early_theta +
                             rt_vmax_lag_sc:trial_neg_inv_sc:entropy_change_late_beta_supp + 
                             # rt_vmax_lag_sc:trial_neg_inv_sc:vmax_late_beta  +
                             rt_vmax_lag_sc:trial_neg_inv_sc:omission_early_theta  +
                             (rt_lag_sc + rt_vmax_lag_sc|id/run), fdf %>% filter(rt_csv<4000))
# screen.lmerTest(fmb_meg1_rs, .01)
summary(fmb_meg1_rs)
Anova(fmb_meg1_rs, '3')

# late beta to entropy change
em1 <- as_tibble(emtrends(mb_meg1_rs, data = df,  var = "rt_lag_sc", specs = c("entropy_change_late_beta_supp", "last_outcome"), at = list(entropy_change_late_beta_supp = c(-.14, .12)), options = list()))
em1$study = "1. MEG"
em2 <- as_tibble(emtrends(fmb_meg1_rs, data = fdf, var = "rt_lag_sc", specs = c("entropy_change_late_beta_supp", "last_outcome"), at = list(entropy_change_late_beta_supp = c(-.14, .12)), options = list()))
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
rem1 <- as_tibble(emtrends(mb_meg1_rs, data = df,  var = "rt_lag_sc", specs = c("omission_early_theta", "last_outcome"), at = list(omission_early_theta = c(-.33, .3)), options = list()))
rem1$study = "1. MEG"
rem2 <- as_tibble(emtrends(fmb_meg1_rs, data = fdf, var = "rt_lag_sc", specs = c("omission_early_theta", "last_outcome"), at = list(omission_early_theta = c(-.33, .3)), options = list()))
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

# coxme models
# MEG
surv_mdf <- inner_join(mbb %>% mutate(id = as.character(id)), wbetas, by = "id")
mb_cox  <- coxme(Surv(t1,t2,response) ~ value_wi + value_wi:entropy_change_late_beta_supp + uncertainty_wi  + 
                   (1|ID), surv_mdf)
summary(mb_cox)
mb_cox1  <- coxme(Surv(t1,t2,response) ~ value_wi + value_wi:entropy_change_late_beta_supp + uncertainty_wi  + uncertainty_wi:entropy_change_late_beta_supp +
                    (1|ID), surv_mdf)
summary(mb_cox1)
# uncertainty random slope:
mb_cox1_rsu  <- coxme(Surv(t1,t2,response) ~ value_wi + value_wi:entropy_change_late_beta_supp + uncertainty_wi  + uncertainty_wi:entropy_change_late_beta_supp +
                        (1 + uncertainty_wi|ID), surv_mdf)
summary(mb_cox1_rsu)

# value random slope:
mb_cox1_rsv  <- coxme(Surv(t1,t2,response) ~ value_wi + value_wi:entropy_change_late_beta_supp + uncertainty_wi  + uncertainty_wi:entropy_change_late_beta_supp +
                        (1 + value_wi|ID), surv_mdf)
summary(mb_cox1_rsv)



# fMRI
surv_fdf <- inner_join(bb %>% mutate(id = as.character(ID)), wbetas, by = "id")
fb_cox  <- coxme(Surv(t1,t2,response) ~ value_wi + value_wi:entropy_change_late_beta_supp + uncertainty_wi  + 
                   (1|ID), surv_fdf)
summary(fb_cox)
fb_cox1  <- coxme(Surv(t1,t2,response) ~ value_wi + value_wi:entropy_change_late_beta_supp + uncertainty_wi  + uncertainty_wi:entropy_change_late_beta_supp +
                    (1|ID), surv_fdf)
summary(fb_cox1)


## old models

# 
# fmb_meg1delta <- lme4::lmer(rt_csv_sc ~ 
#                               rt_lag_sc * last_outcome * entropy_change_early_beta_supp + 
#                               rt_lag_sc * last_outcome * entropy_change_late_beta_supp + 
#                               rt_lag_sc * last_outcome * vmax_late_alpha +
#                               rt_lag_sc * last_outcome * omission_late_delta +
#                               rt_vmax_lag_sc * trial_neg_inv_sc * entropy_change_early_beta_supp + 
#                               rt_vmax_lag_sc * trial_neg_inv_sc * entropy_change_late_beta_supp + 
#                               rt_vmax_lag_sc * trial_neg_inv_sc * vmax_late_alpha  +
#                               rt_vmax_lag_sc * trial_neg_inv_sc * omission_late_delta  +
#                               (1|id/run), fdf %>% filter(rt_csv<4000))
# # screen.lmerTest(fmb_meg1, .01)
# summary(fmb_meg1delta)
# Anova(fmb_meg1delta, '3')
# 
# # late betas suppression to EC from only the posterior "EC" sensors predicts as well as from all sensors
# mb_fmri_ec_allsensors <-  
#   lme4::lmer(rt_csv_sc ~ 
#                rt_lag_sc * last_outcome * entropy_change_early_beta_supp + 
#                rt_lag_sc * last_outcome * entropy_change_late_beta_supp + 
#                rt_vmax_lag_sc * trial_neg_inv_sc * entropy_change_early_beta_supp + 
#                rt_vmax_lag_sc * trial_neg_inv_sc * entropy_change_late_beta_supp + 
#                (1|id/run), fdf %>% filter(rt_csv<4000))
# # screen.lmerTest(mb_fmri1, .01)
# summary(mb_fmri_ec_allsensors)
# Anova(mb_fmri_ec_allsensors, '3')
# 
# # add abs_pe_late_beta_supp_ec
# mb_fmri_ec_abspe_allsensors <-  
#   lme4::lmer(rt_csv_sc ~ 
#                rt_lag_sc * last_outcome * entropy_change_early_beta_supp + 
#                rt_lag_sc * last_outcome * entropy_change_late_beta_supp + 
#                rt_lag_sc * last_outcome * abs_pe_late_beta_supp + 
#                rt_vmax_lag_sc * trial_neg_inv_sc * entropy_change_early_beta_supp + 
#                rt_vmax_lag_sc * trial_neg_inv_sc * entropy_change_late_beta_supp + 
#                rt_vmax_lag_sc * trial_neg_inv_sc * abs_pe_late_beta_supp + 
#                (1|id/run), fdf %>% filter(rt_csv<4000))
# # screen.lmerTest(mb_fmri1, .01)
# summary(mb_fmri_ec_abspe_allsensors)
# Anova(mb_fmri_ec_abspe_allsensors, '3')
# 
# mb_fmri_ec_ecsensors <-  
#   lme4::lmer(rt_csv_sc ~ 
#                rt_lag_sc * last_outcome * entropy_change_early_beta_supp_ec + 
#                rt_lag_sc * last_outcome * entropy_change_late_beta_supp_ec + 
#                rt_vmax_lag_sc * trial_neg_inv_sc * entropy_change_early_beta_supp_ec + 
#                rt_vmax_lag_sc * trial_neg_inv_sc * entropy_change_late_beta_supp_ec + 
#                (1|id/run), fdf %>% filter(rt_csv<4000))
# # screen.lmerTest(mb_fmri1, .01)
# summary(mb_fmri_ec_ecsensors)
# Anova(mb_fmri_ec_ecsensors, '3')
# 
# # add abs_pe_late_beta_supp_ec
# mb_fmri_ec_abspe_ecsensors <-  
#   lme4::lmer(rt_csv_sc ~ 
#                rt_lag_sc * last_outcome * entropy_change_early_beta_supp_ec + 
#                rt_lag_sc * last_outcome * entropy_change_late_beta_supp_ec + 
#                rt_lag_sc * last_outcome * abs_pe_late_beta_supp_ec + 
#                rt_vmax_lag_sc * trial_neg_inv_sc * entropy_change_early_beta_supp_ec + 
#                rt_vmax_lag_sc * trial_neg_inv_sc * entropy_change_late_beta_supp_ec + 
#                rt_vmax_lag_sc * trial_neg_inv_sc * abs_pe_late_beta_supp_ec + 
#                (1|id/run), fdf %>% filter(rt_csv<4000))
# # screen.lmerTest(mb_fmri1, .01)
# summary(mb_fmri_ec_abspe_ecsensors)
# Anova(mb_fmri_ec_abspe_ecsensors, '3')
