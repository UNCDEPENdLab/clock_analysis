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
plot_dir <- "~/OneDrive/collected_letters/papers/meg/plots/meg_to_behavior/"

### load data
source("~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/get_trial_data.R")
# get_trial_data <- function(repo_directory=NULL, dataset="mmclock_fmri", groupfixed=TRUE) 
df <- get_trial_data(repo_directory = clock_folder, dataset = "mmclock_meg", groupfixed = T)

# add wide-format, condition-wise MEG betas, keep all the processing transparent here
cond_wbetas <- readRDS("~/code/clock_analysis/meg/data/MEG_betas_ec_rewfunc_rt_next_reward_rewfunc_April_5_2022.RDS") %>% 
  mutate(entropy_change_early_beta_supp = -  entropy_change_early_beta_ec_rewfunc,
         entropy_change_late_beta_supp = - entropy_change_late_beta_ec_rewfunc,
         rt_shorten_late_beta_supp = - rt_next_late_beta_rt_next,
         omission_early_theta = - reward_early_theta_reward_rewfunc) %>%
  select(id, rewFunc, entropy_change_early_beta_supp, entropy_change_late_beta_supp, rt_shorten_late_beta_supp, 
         omission_early_theta) 
avg <- cond_wbetas %>% group_by(id) %>% summarise(omission_early_theta_avg = mean(omission_early_theta),
                                                  entropy_change_late_beta_avg = mean(entropy_change_late_beta_supp)) %>%
  ungroup()

cond_wbetas <- cond_wbetas %>%  merge(avg, by = "id") %>% mutate(
  ec_lbeta_wi = entropy_change_late_beta_supp - entropy_change_late_beta_avg,
  om_theta_wi = omission_early_theta - omission_early_theta_avg)

# inspect

qmeg <- cond_wbetas %>% group_by(rewFunc) %>% summarize(omission_early_theta = quantile(omission_early_theta, c(.1, .9), names = F),
                                                        q = c("10th %ile", "90th %ile"),
                                                        qcolor = c("#1b3840","#4fa3b8"),
                                                        theta_color = c("orange4", "orange"),
                                                        entropy_change_late_beta_supp = quantile(entropy_change_late_beta_supp, c(.1, .9), names = F),
                                                        entropy_change_late_beta_avg = quantile(entropy_change_late_beta_avg, c(.1, .9), names = F),
                                                        omission_early_theta_avg = quantile(omission_early_theta_avg, c(.1, .9), names = F),
                                                        ec_lbeta_wi = quantile(ec_lbeta_wi, c(.1, .9), names = F),
                                                        om_theta_wi = quantile(om_theta_wi, c(.1, .9), names = F),
                                                        ) %>% ungroup()

# inspect
setwd("~/OneDrive/collected_letters/papers/meg/plots/meg_to_fmri/")
betas_to_inspect <- cond_wbetas %>% select(is.numeric)
cormat <- corr.test(betas_to_inspect)
corrplot(cormat$r, p.mat = cormat$p, order = "hclust", tl.cex = .8, insig = 'blank', method = 'number')

# compare with session-level betas
wbetas <- readRDS("~/code/clock_analysis/meg/data/MEG_betas_entropy_change_entropy_change_ec_sensors_v_max_reward_abspe_ec_sensors_abs_pe_Dec15_2021.RDS") %>% 
  mutate(omission_early_theta_session = - reward_early_theta) %>% 
  mutate(entropy_change_early_beta_supp = -  entropy_change_early_beta_entropy_change,
         entropy_change_late_beta_supp = - entropy_change_late_beta_entropy_change
  ) %>%
  select(c(id, omission_early_theta_session, 
           entropy_change_early_beta_supp, entropy_change_late_beta_supp
  ))

# merge
df <- df %>% inner_join(cond_wbetas, by = c("id", "rewFunc"))


pdf("omission_early_theta_by_cond.pdf", height = 3, width = 5)
ggplot(cond_wbetas, (aes(rewFunc, omission_early_theta, color = rewFunc))) + 
  geom_violin(draw_quantiles = .5) + geom_jitter(alpha = .3)
dev.off()
summary(lm(omission_early_theta ~ rewFunc, cond_wbetas))

ggplot(cond_wbetas, aes(rewFunc, rt_shorten_late_beta_supp, color = rewFunc)) + 
  geom_violin(draw_quantiles = .5) + geom_jitter(alpha = .3)
summary(lm(rt_shorten_late_beta_supp ~ rewFunc, cond_wbetas))
# intercept NS, IEV<CEV,CEVR,DEV NS

pdf("ec_late_beta_by_cond.pdf", height = 3, width = 5)
ggplot(cond_wbetas, aes(rewFunc, entropy_change_late_beta_supp, color = rewFunc)) + 
  geom_violin(draw_quantiles = .5) + geom_boxplot() + geom_jitter(alpha = .3)
dev.off()
Anova(lmer(entropy_change_late_beta_supp ~ rewFunc + (1|id), cond_wbetas))


# only learnable
ldf <- df %>% filter(rewFunc=="IEV" | rewFunc=="DEV")
ev_meg2_rewFunc <-  
  lmerTest::lmer(ev ~ trial_neg_inv_sc * entropy_change_late_beta_supp * rewFunc + 
                   trial_neg_inv_sc * omission_early_theta * rewFunc +
                   (1|id), ldf %>% filter(rt_csv<4000))
screen.lmerTest(ev_meg2_rewFunc, .01)
summary(ev_meg2_rewFunc)
Anova(ev_meg2_rewFunc, '3')

ev_meg3_rewFunc <-  
  lmerTest::lmer(ev ~ 
                   (trial_neg_inv_sc + entropy_change_late_beta_supp + 
                      omission_early_theta + rewFunc)^3 +
                   (1|id), ldf %>% filter(rt_csv<4000))
screen.lmerTest(ev_meg3_rewFunc, .01)
summary(ev_meg3_rewFunc)
Anova(ev_meg3_rewFunc, '3')
anova(ev_meg2_rewFunc, ev_meg3_rewFunc)
ec1 <- ggplot(ldf, aes(trial_neg_inv_sc, ev, color = entropy_change_late_beta_supp > 0.19136, lty = reward_lag)) + geom_smooth(method = "gam") +
  facet_wrap(~rewFunc) + scale_color_manual("Entropy change late beta suppression", values = c("1b3840", "4fa3b8"), labels = c("low", "high"))  
o1 <- ggplot(ldf, aes(trial_neg_inv_sc, ev, color = omission_early_theta > 0.46125, lty = reward_lag)) + geom_smooth(method = "gam") +
  facet_wrap(~rewFunc) + scale_color_manual("Omission early theta", values = c("orange4", "orange"), labels = c("low", "high"))  
ec2 <- ggplot(df, aes(trial_neg_inv_sc, rt_csv, color = entropy_change_late_beta_supp > 0.19136, lty = reward_lag)) + geom_smooth(method = "gam") +
  facet_wrap(~rewFunc) + scale_color_manual("Entropy change late beta suppression", values = c("1b3840", "4fa3b8"), labels = c("low", "high"))  
o2 <- ggplot(df, aes(trial_neg_inv_sc, rt_csv, color = omission_early_theta > 0.46125, lty = reward_lag)) + geom_smooth(method = "gam") +
  facet_wrap(~rewFunc) + scale_color_manual("Omission early theta", values = c("orange4", "orange"), labels = c("low", "high"))  
ec3 <- ggplot(df, aes(trial_neg_inv_sc, rt_swing, color = entropy_change_late_beta_supp > 0.19136, lty = reward_lag)) + geom_smooth(method = "gam") +
  facet_wrap(~rewFunc) + scale_color_manual("Entropy change late beta suppression", values = c("1b3840", "4fa3b8"), labels = c("low", "high"))  
o3 <- ggplot(df, aes(trial_neg_inv_sc, rt_swing, color = omission_early_theta > 0.46125, lty = reward_lag)) + geom_smooth(method = "gam") +
  facet_wrap(~rewFunc) + scale_color_manual("Omission early theta", values = c("orange4", "orange"), labels = c("low", "high"))  
# setwd(plot_dir)
# pdf("performance_by_meg_beta_rewfunc.pdf", height = 8, width = 18)
# ggarrange(o1, o2, o3, ec1, ec2, ec3, nrow = 2, ncol = 3)
# dev.off()

# look at convergence on RT_Vmax
df <- df %>% mutate(rt_vmax_delta  = abs(rt_csv_sc - rt_vmax_lag_sc))
ggplot(df, aes(run_trial, rt_vmax_delta, color = omission_early_theta > 0.46125, lty = reward_lag)) + geom_smooth(method = "gam") +
  facet_wrap(~rewFunc) + scale_color_manual("Omission early theta", values = c("orange4", "orange"), labels = c("low", "high"))  

ggplot(df, aes(run_trial, rt_vmax_delta, color = entropy_change_late_beta_supp > 0.19136, lty = reward_lag)) + geom_smooth(method = "gam") +
  facet_wrap(~rewFunc) + scale_color_manual("Entropy change late beta suppression", values = c("1b3840", "4fa3b8"), labels = c("low", "high"))  

# example subjects: not revealing
# theta
# library(ggh4x)
# setwd(plot_dir)
# pdf("theta_rtswing_subject_megaplot.pdf", height = 30, width = 60)
# ggplot(ldf %>% filter(run_trial>1)) + geom_line( aes(run_trial, log(rt_swing), color = id, group = interaction(id, run))) + geom_smooth(aes(run_trial, log(rt_swing)), method = "gam") +
#   facet_nested((omission_early_theta > 0.46125)  ~ rewFunc + reward_lag) 
# dev.off()
# 
# pdf("beta_rtswing_subject_megaplot.pdf", height = 30, width = 60)
# ggplot(ldf %>% filter(run_trial>1)) + geom_line( aes(run_trial, log(rt_swing), color = id, group = interaction(id, run))) + geom_smooth(aes(run_trial, log(rt_swing)), method = "gam") +
#   facet_nested((entropy_change_late_beta_supp > 0.19136)  ~ rewFunc + reward_lag) 
# dev.off()
# 
# pdf("theta_rt_csv_subject_megaplot.pdf", height = 30, width = 60)
# ggplot(ldf %>% filter(run_trial>1)) + geom_line( aes(run_trial, log(rt_csv), color = id, group = interaction(id, run))) + geom_smooth(aes(run_trial, log(rt_csv)), method = "gam") +
#   facet_nested((omission_early_theta > 0.46125)  ~ rewFunc + reward_lag) 
# dev.off()
# 
# pdf("beta_rt_csv_subject_megaplot.pdf", height = 30, width = 60)
# ggplot(ldf %>% filter(run_trial>1)) + geom_line( aes(run_trial, log(rt_csv), color = id, group = interaction(id, run))) + geom_smooth(aes(run_trial, log(rt_csv)), method = "gam") +
#   facet_nested((entropy_change_late_beta_supp > 0.19136)  ~ rewFunc + reward_lag) 
# dev.off()


# Condition effects on RT:
cond_meg2_rewFunc <-  
  lmerTest::lmer(rt_csv ~ 
                   (trial_neg_inv_sc + omission_early_theta + rewFunc)^2 +
                   (trial_neg_inv_sc + entropy_change_late_beta_supp + rewFunc)^2 +
                 (1|id), df %>% filter(rt_csv<4000))
screen.lmerTest(cond_meg2_rewFunc)
summary(cond_meg2_rewFunc)
Anova(cond_meg2_rewFunc, '3')

# based on deviations
cond_meg2_rewFunc_wi <-  
  lmerTest::lmer(rt_csv ~ 
                   (trial_neg_inv_sc + omission_early_theta + rewFunc)^2 +
                   (trial_neg_inv_sc + om_theta_wi + rewFunc)^2 +
                   (trial_neg_inv_sc + entropy_change_late_beta_supp + rewFunc)^2 +
                   (trial_neg_inv_sc + ec_lbeta_wi + rewFunc)^2 +
                   (1|id), df %>% filter(rt_csv<4000))
screen.lmerTest(cond_meg2_rewFunc_wi)
anova(cond_meg2_rewFunc, cond_meg2_rewFunc_wi)
summary(cond_meg2_rewFunc)
Anova(cond_meg2_rewFunc, '3')


cond_meg3_rewFunc <-  
  lmerTest::lmer(ev ~ 
                   (trial_neg_inv_sc + entropy_change_late_beta_supp + 
                      omission_early_theta + rewFunc)^3 +
                   (1|id), ldf %>% filter(rt_csv<4000))
screen.lmerTest(cond_meg3_rewFunc, .01)
summary(cond_meg3_rewFunc)
Anova(ev_meg2_rewFunc, '3')
anova(ev_meg2_rewFunc, ev_meg3_rewFunc)


############# MEG
# Effect of MEG betas on exploration
# only the effects of interest
beh <-  
  lme4::lmer(rt_csv_sc ~ 
               rt_lag_sc * last_outcome * rewFunc + 
               rt_vmax_lag_sc * trial_neg_inv_sc * rewFunc + 
               (1|id/run), df %>% filter(rt_csv<4000))
# screen.lmerTest(meg1, .01)
summary(beh)
Anova(beh, '3')

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


meg_late_beta_rewFunc <-  
  lme4::lmer(rt_csv_sc ~ 
               rt_lag_sc * last_outcome * entropy_change_late_beta_supp * rewFunc + 
               rt_vmax_lag_sc * trial_neg_inv_sc * entropy_change_late_beta_supp * rewFunc + 
               (1|id/run), df %>% filter(rt_csv<4000))
# screen.lmerTest(meg1, .01)
summary(meg_late_beta_rewFunc)
Anova(meg_late_beta_rewFunc, '3')

# do RT swings decrease with trial?
meg_late_beta_rewFunc_trial <-  
  lme4::lmer(rt_csv_sc ~ 
               (rt_lag_sc + last_outcome + entropy_change_late_beta_supp + rewFunc + trial_neg_inv_sc)^3 +
               rt_vmax_lag_sc * trial_neg_inv_sc * entropy_change_late_beta_supp * rewFunc + 
               (1|id/run), df %>% filter(rt_csv<4000))
# screen.lmerTest(meg1, .01)
summary(meg_late_beta_rewFunc)
Anova(meg_late_beta_rewFunc, '3')

meg_rts_beta_rewFunc <-  
  lme4::lmer(rt_csv_sc ~ 
               rt_lag_sc * last_outcome * rt_shorten_late_beta_supp * rewFunc + 
               rt_vmax_lag_sc * trial_neg_inv_sc * rt_shorten_late_beta_supp * rewFunc + 
               (1|id/run), df %>% filter(rt_csv<4000))
# screen.lmerTest(meg1, .01)
summary(meg_rts_beta_rewFunc)
Anova(meg_rts_beta_rewFunc, '3')

meg_all_late_beta_rewFunc <-  
  lme4::lmer(rt_csv_sc ~ 
               rt_lag_sc * last_outcome * entropy_change_late_beta_supp * rewFunc + 
               rt_vmax_lag_sc * trial_neg_inv_sc * entropy_change_late_beta_supp * rewFunc + 
               rt_lag_sc * last_outcome * rt_shorten_late_beta_supp * rewFunc + 
               rt_vmax_lag_sc * trial_neg_inv_sc * rt_shorten_late_beta_supp * rewFunc + 
               (1|id/run), df %>% filter(rt_csv<4000))
# screen.lmerTest(meg1, .01)
summary(meg_all_late_beta_rewFunc)
Anova(meg_all_late_beta_rewFunc, '3')

meg_om_theta_lbeta_rewFunc <-
  lmerTest::lmer(rt_csv_sc ~
               rt_lag_sc * last_outcome * entropy_change_late_beta_supp * rewFunc +
               rt_lag_sc * last_outcome * omission_early_theta * rewFunc +
               rt_vmax_lag_sc * trial_neg_inv_sc * entropy_change_late_beta_supp * rewFunc +
               rt_vmax_lag_sc * trial_neg_inv_sc * omission_early_theta * rewFunc  +
               (1|id/run), df %>% filter(rt_csv<4000))
screen.lmerTest(meg_om_theta_lbeta_rewFunc, .01)
summary(meg_om_theta_lbeta_rewFunc)
Anova(meg_om_theta_lbeta_rewFunc, '3')

# decomposing into subject and condition
meg_om_theta_lbeta_decomposed <-
  lme4::lmer(rt_csv_sc ~
               rt_lag_sc * last_outcome * entropy_change_late_beta_avg * rewFunc +
               rt_lag_sc * last_outcome * ec_lbeta_wi * rewFunc +
               rt_lag_sc * last_outcome * omission_early_theta_avg * rewFunc +
               rt_lag_sc * last_outcome * om_theta_wi * rewFunc +
               # rt_vmax_lag_sc * trial_neg_inv_sc * entropy_change_late_beta_avg * rewFunc +
               # rt_vmax_lag_sc * trial_neg_inv_sc * ec_lbeta_wi * rewFunc +
               # rt_vmax_lag_sc * trial_neg_inv_sc * omission_early_theta_avg * rewFunc  +
               # rt_vmax_lag_sc * trial_neg_inv_sc * om_theta_wi * rewFunc  +
               (1|id/run), df %>% filter(rt_csv<4000))
# screen.lmerTest(meg_om_theta_lbeta_decomposed, .01)
# anova(meg_om_theta_lbeta_rewFunc, meg_om_theta_lbeta_decomposed)

summary(meg_om_theta_lbeta_decomposed)
Anova(meg_om_theta_lbeta_decomposed, '3')


# understand divergence from previous results: average omission theta across condition has NS pro-RT swing effect
meg_om_theta_avg <-
  lme4::lmer(rt_csv_sc ~
               rt_lag_sc * last_outcome * omission_early_theta_avg  +
               rt_vmax_lag_sc * trial_neg_inv_sc * omission_early_theta_avg  +
               (1|id/run), df %>% filter(rt_csv<4000))
# screen.lmerTest(meg1, .01)
summary(meg_om_theta_avg)
Anova(meg_om_theta_avg, '3')

# with condition
meg_om_theta_avg_rewFunc <-
  lme4::lmer(rt_csv_sc ~
               rt_lag_sc * last_outcome * omission_early_theta_avg * rewFunc  +
               rt_vmax_lag_sc * trial_neg_inv_sc * omission_early_theta_avg * rewFunc  +
               (1|id/run), df %>% filter(rt_csv<4000))
# screen.lmerTest(meg1, .01)
summary(meg_om_theta_avg_rewFunc)
Anova(meg_om_theta_avg_rewFunc, '3')

# 
# meg1_pe_theta <-  
#   lme4::lmer(rt_csv_sc ~ 
#                rt_lag_sc * last_outcome * entropy_change_early_beta_supp + 
#                rt_lag_sc * last_outcome * entropy_change_late_beta_supp + 
#                rt_lag_sc * last_outcome * vmax_late_alpha +
#                rt_lag_sc * last_outcome * neg_pe_early_theta +
#                rt_vmax_lag_sc * trial_neg_inv_sc * entropy_change_early_beta_supp + 
#                rt_vmax_lag_sc * trial_neg_inv_sc * entropy_change_late_beta_supp + 
#                rt_vmax_lag_sc * trial_neg_inv_sc * vmax_late_alpha  +
#                rt_vmax_lag_sc * trial_neg_inv_sc * neg_pe_early_theta  +
#                (1|id/run), df %>% filter(rt_csv<4000))
# # screen.lmerTest(meg1, .01)
# summary(meg1_pe_theta)
# Anova(meg1_pe_theta, '3')
# anova(meg1_om_theta, meg1_pe_theta)
# # note: omission predicts much better than PE early theta
# 
# meg_ec <-  
#   lme4::lmer(rt_csv_sc ~ 
#                rt_lag_sc * last_outcome * entropy_change_early_beta_supp + 
#                rt_lag_sc * last_outcome * entropy_change_late_beta_supp + 
#                rt_vmax_lag_sc * trial_neg_inv_sc * entropy_change_early_beta_supp + 
#                rt_vmax_lag_sc * trial_neg_inv_sc * entropy_change_late_beta_supp + 
#                (1|id/run), df %>% filter(rt_csv<4000))
# # screen.lmerTest(meg1, .01)
# summary(meg_ec)
# Anova(meg_ec, '3')
# 
# meg_ec_pe <-  
#   lme4::lmer(rt_csv_sc ~ 
#                rt_lag_sc * last_outcome * entropy_change_early_beta_supp + 
#                rt_lag_sc * last_outcome * entropy_change_late_beta_supp + 
#                rt_lag_sc * last_outcome * pe_late_beta_supp +
#                rt_lag_sc * last_outcome * neg_pe_early_theta + 
#                rt_vmax_lag_sc * trial_neg_inv_sc * entropy_change_early_beta_supp + 
#                rt_vmax_lag_sc * trial_neg_inv_sc * entropy_change_late_beta_supp + 
#                rt_vmax_lag_sc * trial_neg_inv_sc * pe_late_beta_supp +
#                rt_vmax_lag_sc * trial_neg_inv_sc * neg_pe_early_theta +
#                (1|id/run), df %>% filter(rt_csv<4000))
# # screen.lmerTest(meg1, .01)
# summary(meg_ec_pe)
# Anova(meg_ec_pe, '3')
# 
# # hypothesis: best model combines omission early theta and PE late beta with EC responses
# # correct!
# meg_hybrid <-  
#   lme4::lmer(rt_csv_sc ~ 
#                # rt_lag_sc * last_outcome * entropy_change_early_beta_supp + 
#                rt_lag_sc * last_outcome * entropy_change_late_beta_supp + 
#                rt_lag_sc * last_outcome * pe_late_beta_supp +
#                rt_lag_sc * last_outcome * omission_early_theta + 
#                # rt_vmax_lag_sc * trial_neg_inv_sc * entropy_change_early_beta_supp + 
#                rt_vmax_lag_sc * trial_neg_inv_sc * entropy_change_late_beta_supp + 
#                rt_vmax_lag_sc * trial_neg_inv_sc * pe_late_beta_supp +
#                rt_vmax_lag_sc * trial_neg_inv_sc * omission_early_theta +
#                (1|id/run), df %>% filter(rt_csv<4000))
# # screen.lmerTest(meg1, .01)
# summary(meg_hybrid)
# Anova(meg_hybrid, '3')
# 
# # compare models
# anova(meg_hybrid, meg_ec_pe, meg_ec, meg1_om_theta)
# 
# # add random slopes of RT and RT_Vmax
# meg_hybrid_rs <-  
#   lme4::lmer(rt_csv_sc ~ 
#                # rt_lag_sc * last_outcome * entropy_change_early_beta_supp + 
#                rt_lag_sc * last_outcome * entropy_change_late_beta_supp + 
#                rt_lag_sc * last_outcome * pe_late_beta_supp +
#                rt_lag_sc * last_outcome * omission_early_theta + 
#                # rt_vmax_lag_sc * trial_neg_inv_sc * entropy_change_early_beta_supp + 
#                rt_vmax_lag_sc * trial_neg_inv_sc * entropy_change_late_beta_supp + 
#                rt_vmax_lag_sc * trial_neg_inv_sc * pe_late_beta_supp +
#                rt_vmax_lag_sc * trial_neg_inv_sc * omission_early_theta +
#                (rt_lag_sc + rt_vmax_lag_sc|id/run), df %>% filter(rt_csv<4000))
# # screen.lmerTest(meg1, .01)
# summary(meg_hybrid_rs)
# Anova(meg_hybrid_rs, '3')
# 

# load fmri data
fdf <- get_trial_data(repo_directory = clock_folder, dataset = "mmclock_fmri", groupfixed = T)
fdf <- fdf %>% mutate(id = as.character(id)) %>% inner_join(cond_wbetas, by = c("id", "rewFunc"))

lfdf <-fdf %>% filter(rewFunc=="IEV" | rewFunc=="DEV")

# descriptives
fec1 <- ggplot(lfdf, aes(trial_neg_inv_sc, ev, color = entropy_change_late_beta_supp > 0.19136, lty = reward_lag)) + geom_smooth(method = "gam") +
  facet_wrap(~rewFunc) + scale_color_manual("Entropy change late beta suppression", values = c("1b3840", "4fa3b8"), labels = c("low", "high"))  
fo1 <- ggplot(lfdf, aes(trial_neg_inv_sc, ev, color = omission_early_theta > 0.46125, lty = reward_lag)) + geom_smooth(method = "gam") +
  facet_wrap(~rewFunc) + scale_color_manual("Omission early theta", values = c("orange4", "orange"), labels = c("low", "high"))  
fec2 <- ggplot(fdf, aes(trial_neg_inv_sc, rt_csv, color = entropy_change_late_beta_supp > 0.19136, lty = reward_lag)) + geom_smooth(method = "gam") +
  facet_wrap(~rewFunc) + scale_color_manual("Entropy change late beta suppression", values = c("1b3840", "4fa3b8"), labels = c("low", "high"))  
fo2 <- ggplot(fdf, aes(trial_neg_inv_sc, rt_csv, color = omission_early_theta > 0.46125, lty = reward_lag)) + geom_smooth(method = "gam") +
  facet_wrap(~rewFunc) + scale_color_manual("Omission early theta", values = c("orange4", "orange"), labels = c("low", "high"))  
fec3 <- ggplot(fdf, aes(trial_neg_inv_sc, rt_swing, color = entropy_change_late_beta_supp > 0.19136, lty = reward_lag)) + geom_smooth(method = "gam") +
  facet_wrap(~rewFunc) + scale_color_manual("Entropy change late beta suppression", values = c("1b3840", "4fa3b8"), labels = c("low", "high"))  
fo3 <- ggplot(fdf, aes(trial_neg_inv_sc, rt_swing, color = omission_early_theta > 0.46125, lty = reward_lag)) + geom_smooth(method = "gam") +
  facet_wrap(~rewFunc) + scale_color_manual("Omission early theta", values = c("orange4", "orange"), labels = c("low", "high"))  


thetaplot <- ggarrange(o1, o2, o3, 
                       fo1, fo2, fo3,  nrow = 2, ncol = 3, common.legend = T, labels = c("MEG","", "", "fMRI", "", ""), label.y = 1.1)
betaplot <- ggarrange(ec1, ec2, ec3, 
                      fec1, fec2, fec3,  nrow = 2, ncol = 3, common.legend = T)

# convergence on RT_Vmax
fdf <- fdf %>% mutate(rt_vmax_delta  = abs(rt_csv_sc - rt_vmax_lag_sc))
ggplot(fdf, aes(run_trial, rt_vmax_delta, color = omission_early_theta > 0.46125, lty = reward_lag)) + geom_smooth(method = "gam") +
  facet_wrap(~rewFunc) + scale_color_manual("Omission early theta", values = c("orange4", "orange"), labels = c("low", "high"))  
ggplot(fdf, aes(run_trial, rt_vmax_delta, color = entropy_change_late_beta_supp > 0.19136, lty = reward_lag)) + geom_smooth(method = "gam") +
  facet_wrap(~rewFunc) + scale_color_manual("Entropy change late beta suppression", values = c("1b3840", "4fa3b8"), labels = c("low", "high"))  



setwd(plot_dir)
pdf("performance_by_meg_beta_rewfunc_replication.pdf", height = 6, width = 22)
ggarrange(thetaplot, betaplot, nrow = 1, ncol = 2)
dev.off()



ev_fmri2_rewFunc <-  
  lmerTest::lmer(ev ~ 
                   (trial_neg_inv_sc + entropy_change_late_beta_supp + 
                      omission_early_theta + rewFunc)^2 +
                   (1|id), lfdf %>% filter(rt_csv<4000))
screen.lmerTest(ev_fmri2_rewFunc, .01)
summary(ev_fmri2_rewFunc)
Anova(ev_fmri2_rewFunc, '3')

ev_fmri3_rewFunc <-  
  lmerTest::lmer(ev ~ 
                   (trial_neg_inv_sc + entropy_change_late_beta_supp + 
                      omission_early_theta + rewFunc)^3 +
                   (1|id), lfdf %>% filter(rt_csv<4000))
screen.lmerTest(ev_fmri3_rewFunc, .01)
summary(ev_fmri3_rewFunc)
Anova(ev_fmri3_rewFunc, '3')
anova(ev_fmri2_rewFunc, ev_fmri3_rewFunc)

# Condition effects:
cond_fmri2_rewFunc <-  
  lmerTest::lmer(rt_csv ~ 
                   (trial_neg_inv_sc + entropy_change_late_beta_supp + 
                      omission_early_theta + rewFunc)^2 +
                   (1|id),fdf %>% filter(rt_csv<4000))
screen.lmerTest(cond_fmri2_rewFunc)
summary(cond_fmri2_rewFunc)
Anova(cond_fmri2_rewFunc, '3')

cond_fmri3_rewFunc <-  
  lmerTest::lmer(ev ~ 
                   (trial_neg_inv_sc + entropy_change_late_beta_supp + 
                      omission_early_theta + rewFunc)^3 +
                   (1|id),fdf %>% filter(rt_csv<4000))
screen.lmerTest(cond_fmri3_rewFunc, .01)
summary(cond_fmri3_rewFunc)
Anova(ev_fmri2_rewFunc, '3')
anova(ev_fmri2_rewFunc, ev_fmri3_rewFunc)

# decompose betas into subject- and condition-level
fmri_om_theta_lbeta_decomposed <-
  lme4::lmer(rt_csv_sc ~
                   rt_lag_sc * last_outcome * entropy_change_late_beta_avg * rewFunc +
                   rt_lag_sc * last_outcome * ec_lbeta_wi * rewFunc +
                   rt_lag_sc * last_outcome * omission_early_theta_avg * rewFunc +
                   rt_lag_sc * last_outcome * om_theta_wi * rewFunc +
                   rt_vmax_lag_sc * trial_neg_inv_sc * entropy_change_late_beta_avg * rewFunc +
                   rt_vmax_lag_sc * trial_neg_inv_sc * ec_lbeta_wi * rewFunc +
                   rt_vmax_lag_sc * trial_neg_inv_sc * omission_early_theta_avg * rewFunc  +
                   rt_vmax_lag_sc * trial_neg_inv_sc * om_theta_wi * rewFunc  +
                   (1|id/run), fdf %>% filter(rt_csv<4000))
# screen.lmerTest(fmri_om_theta_lbeta_decomposed, .01)
anova(fmri_om_theta_lbeta_rewFunc, fmri_om_theta_lbeta_decomposed)

summary(fmri_om_theta_lbeta_decomposed)
Anova(fmri_om_theta_lbeta_decomposed, '3')

# # now that we have done model selection in the MEG study, replicate only the best model
# fmri_hybrid <-  
#   lme4::lmer(rt_csv_sc ~ 
#                # rt_lag_sc * last_outcome * entropy_change_early_beta_supp + 
#                rt_lag_sc * last_outcome * entropy_change_late_beta_supp + 
#                rt_lag_sc * last_outcome * pe_late_beta_supp +
#                rt_lag_sc * last_outcome * omission_early_theta + 
#                # rt_vmax_lag_sc * trial_neg_inv_sc * entropy_change_early_beta_supp + 
#                rt_vmax_lag_sc * trial_neg_inv_sc * entropy_change_late_beta_supp + 
#                rt_vmax_lag_sc * trial_neg_inv_sc * pe_late_beta_supp +
#                rt_vmax_lag_sc * trial_neg_inv_sc * omission_early_theta +
#                (1|id/run),fdf %>% filter(rt_csv<4000))
# # screen.lmerTest(meg1, .01)
# summary(fmri_hybrid)
# Anova(fmri_hybrid, '3')
# 
# # add random slopes of RT and RT_Vmax
# fmri_hybrid_rs <-  
#   lme4::lmer(rt_csv_sc ~ 
#                # rt_lag_sc * last_outcome * entropy_change_early_beta_supp + 
#                rt_lag_sc * last_outcome * entropy_change_late_beta_supp + 
#                rt_lag_sc * last_outcome * pe_late_beta_supp +
#                rt_lag_sc * last_outcome * omission_early_theta + 
#                # rt_vmax_lag_sc * trial_neg_inv_sc * entropy_change_early_beta_supp + 
#                rt_vmax_lag_sc * trial_neg_inv_sc * entropy_change_late_beta_supp + 
#                rt_vmax_lag_sc * trial_neg_inv_sc * pe_late_beta_supp +
#                rt_vmax_lag_sc * trial_neg_inv_sc * omission_early_theta +
#                (rt_lag_sc + rt_vmax_lag_sc|id/run),fdf %>% filter(rt_csv<4000))
# # screen.lmerTest(meg1, .01)
# summary(fmri_hybrid_rs)
# Anova(fmri_hybrid_rs, '3')

# only the effects of interest
fmri_late_beta_only <-  
  lme4::lmer(rt_csv_sc ~ 
               rt_lag_sc * last_outcome * entropy_change_late_beta_supp + 
               rt_vmax_lag_sc * trial_neg_inv_sc * entropy_change_late_beta_supp + 
               (1|id/run),fdf %>% filter(rt_csv<4000))
# screen.lmerTest(meg1, .01)
summary(fmri_late_beta_only)
Anova(fmri_late_beta_only, '3')

# effect of EC late beta only in fMRI sample, by condition
fmri_late_beta_rewFunc <-  
  lme4::lmer(rt_csv_sc ~ 
               rt_lag_sc * last_outcome * entropy_change_late_beta_supp * rewFunc + 
               rt_vmax_lag_sc * trial_neg_inv_sc * entropy_change_late_beta_supp * rewFunc + 
               (1|id/run),fdf %>% filter(rt_csv<4000))
# screen.lmerTest(meg1, .01)
summary(fmri_late_beta_rewFunc)
Anova(fmri_late_beta_rewFunc, '3')

# effect of RT shortening late beta only in fMRI sample, by condition
fmri_rts_beta_rewFunc <-  
  lme4::lmer(rt_csv_sc ~ 
               rt_lag_sc * last_outcome * rt_shorten_late_beta_supp * rewFunc + 
               rt_vmax_lag_sc * trial_neg_inv_sc * rt_shorten_late_beta_supp * rewFunc + 
               (1|id/run),fdf %>% filter(rt_csv<4000))
# screen.lmerTest(meg1, .01)
summary(fmri_rts_beta_rewFunc)
Anova(fmri_rts_beta_rewFunc, '3')

# effect of both
fmri_all_late_beta_rewFunc <-  
  lme4::lmer(rt_csv_sc ~ 
               rt_lag_sc * last_outcome * entropy_change_late_beta_supp * rewFunc + 
               rt_vmax_lag_sc * trial_neg_inv_sc * entropy_change_late_beta_supp * rewFunc + 
               rt_lag_sc * last_outcome * rt_shorten_late_beta_supp * rewFunc + 
               rt_vmax_lag_sc * trial_neg_inv_sc * rt_shorten_late_beta_supp * rewFunc + 
               (1|id/run),fdf %>% filter(rt_csv<4000))
# screen.lmerTest(meg1, .01)
summary(fmri_all_late_beta_rewFunc)
Anova(fmri_all_late_beta_rewFunc, '3')


fmri_om_theta_lbeta_rewFunc <- lme4::lmer(rt_csv_sc ~
                                            rt_lag_sc * last_outcome * entropy_change_late_beta_supp * rewFunc +
                                            rt_lag_sc * last_outcome * omission_early_theta * rewFunc +
                                            rt_vmax_lag_sc * trial_neg_inv_sc * entropy_change_late_beta_supp * rewFunc +
                                            rt_vmax_lag_sc * trial_neg_inv_sc * omission_early_theta * rewFunc  +
                                            (1|id/run),fdf %>% filter(rt_csv<4000))
# screen.lmerTest(fmri1, .01)
summary(fmri_om_theta_lbeta_rewFunc)
Anova(fmri_om_theta_lbeta_rewFunc, '3')

# fmri1delta <- lme4::lmer(rt_csv_sc ~ 
#                               rt_lag_sc * last_outcome * entropy_change_early_beta_supp + 
#                               rt_lag_sc * last_outcome * entropy_change_late_beta_supp + 
#                               rt_lag_sc * last_outcome * vmax_late_alpha +
#                               rt_lag_sc * last_outcome * omission_late_delta +
#                               rt_vmax_lag_sc * trial_neg_inv_sc * entropy_change_early_beta_supp + 
#                               rt_vmax_lag_sc * trial_neg_inv_sc * entropy_change_late_beta_supp + 
#                               rt_vmax_lag_sc * trial_neg_inv_sc * vmax_late_alpha  +
#                               rt_vmax_lag_sc * trial_neg_inv_sc * omission_late_delta  +
#                               (1|id/run),fdf %>% filter(rt_csv<4000))
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
#                (1|id/run),fdf %>% filter(rt_csv<4000))
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
#                (1|id/run),fdf %>% filter(rt_csv<4000))
# # screen.lmerTest(meg1, .01)
# summary(fmri_ec_pe)
# Anova(fmri_ec_pe, '3')



# late beta to entropy change
em1 <- as_tibble(emtrends(meg_late_beta_rewFunc, data = df,  var = "rt_lag_sc", specs = c("entropy_change_late_beta_supp", "last_outcome", "rewFunc"), at = list(entropy_change_late_beta_supp = qmeg$entropy_change_late_beta_supp))) %>% 
  inner_join(qmeg, by = c("entropy_change_late_beta_supp", "rewFunc"))
em1$study = "1. MEG"
em2 <- as_tibble(emtrends(fmri_late_beta_rewFunc, data =fdf, var = "rt_lag_sc", specs = c("entropy_change_late_beta_supp", "last_outcome", "rewFunc"), at = list(entropy_change_late_beta_supp = qmeg$entropy_change_late_beta_supp))) %>%  
  inner_join(qmeg, by = c("entropy_change_late_beta_supp", "rewFunc"))
em2$study = '2. fMRI replication'
em1 <- rbind(em2, em1)
ec_late_beta <- 
  ggplot(em1, aes(x=last_outcome, y=rt_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=q)) + 
  geom_point(position = position_dodge(width = .6), size=2.5) + 
  geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + 
  theme_bw(base_size=12) + facet_grid(rewFunc~study)+ ylab("RT swings (AU)\n Small <---------> Large")  + 
  scale_color_manual("Entropy change\nlate beta\nsuppression", values = em1$qcolor)  +
  labs(shape = "Entropy change\nlate beta\nsuppression") +
  theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) + 
  scale_y_reverse(limits = c(.75, -.05)) 

setwd("~/OneDrive/collected_letters/papers/meg/plots/meg_to_fmri/")
pdf("meg_ec_late_beta_to_behavior_rewFunc.pdf", height = 8, width = 8)
print(ec_late_beta)
dev.off()

# no rewFunc
em3 <- as_tibble(emtrends(meg_late_beta_rewFunc, data = df,  var = "rt_lag_sc", specs = c("entropy_change_late_beta_supp", "last_outcome"), at = list(entropy_change_late_beta_supp = c(-0.183, 0.317)), options = list()))
em3$study = "1. MEG"
em4 <- as_tibble(emtrends(fmri_late_beta_rewFunc, data =fdf, var = "rt_lag_sc", specs = c("entropy_change_late_beta_supp", "last_outcome"), at = list(entropy_change_late_beta_supp = c(-0.183, 0.317)), options = list()))
em4$study = '2. fMRI replication'
em5 <- rbind(em3, em4)
ec_late_beta_no_rewFunc <- 
  ggplot(em5, aes(x=last_outcome, y=rt_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=as.factor(entropy_change_late_beta_supp))) + 
  geom_point(position = position_dodge(width = .6), size=2.5) + 
  geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + 
  theme_bw(base_size=12) + facet_grid(~study)+ ylab("RT swings (AU)\n Small <---------> Large")  + 
  scale_color_manual("Entropy change\nlate beta\nsuppression", values = unique(qmeg$qcolor), labels = unique(qmeg$q)) +
  labs(shape = "Entropy change\nlate beta\nsuppression") +
  theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) + 
  scale_y_reverse(limits = c(.7, -.01)) 

setwd("~/OneDrive/collected_letters/papers/meg/plots/meg_to_fmri/")
pdf("meg_ec_late_beta_to_behavior_all_conditions.pdf", height = 3, width = 5)
print(ec_late_beta_no_rewFunc)
dev.off()

# late beta to RT shortening
# df$rt_shorten_late_beta_supp 
# Mean      Gmd      .05      .10      .25      .50      .75      .90 
# -0.03118   0.2789 -0.40689 -0.29657 -0.14373 -0.01340  0.08042  0.22498 
rtm1 <- as_tibble(emtrends(meg_all_late_beta_rewFunc, data = df,  var = "rt_lag_sc", specs = c("rt_shorten_late_beta_supp", "last_outcome", "rewFunc"), at = list(rt_shorten_late_beta_supp  = c(-0.297,  0.225)), options = list()))
rtm1$study = "1. MEG"
rtm2 <- as_tibble(emtrends(fmri_all_late_beta_rewFunc, data =fdf, var = "rt_lag_sc", specs = c("rt_shorten_late_beta_supp", "last_outcome", "rewFunc"), at = list(rt_shorten_late_beta_supp  = c(-0.297,  0.225)), options = list()))
rtm2$study = '2. fMRI replication'
rtm1 <- rbind(rtm2, rtm1)
rt_late_beta <- 
  ggplot(rtm1, aes(x=last_outcome, y=rt_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=as.factor(rt_shorten_late_beta_supp ))) + 
  geom_point(position = position_dodge(width = .6), size=2.5) + 
  geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + 
  theme_bw(base_size=12) + facet_grid(rewFunc~study) + ylab("RT swings (AU)\n Small <---------> Large")  + 
  scale_color_manual("RT shortening\nlate beta\nsuppression", values=c("#1b3840","#4fa3b8"), labels = c("10th %ile", "90th %ile")) +
  labs(shape = "RT shortening\nlate beta\nsuppression") +
  theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) + 
  scale_y_reverse(limits = c(.7, -.05)) 

setwd("~/OneDrive/collected_letters/papers/meg/plots/meg_to_fmri/")
pdf("meg_rt_shortening_late_beta_to_behavior_rewFunc.pdf", height = 8, width = 8)
print(rt_late_beta)
dev.off()

# EC late beta controlling for RT-shortening late beta
emr1 <- as_tibble(emtrends(meg_all_late_beta_rewFunc, data = df,  var = "rt_lag_sc", specs = c("entropy_change_late_beta_supp", "last_outcome", "rewFunc"), at = list(omission_early_theta = qmeg$omission_early_theta))) %>%   inner_join(qmeg, by = c("omission_early_theta", "rewFunc"))
emr1$study = "1. MEG"
emr2 <- as_tibble(emtrends(fmri_all_late_beta_rewFunc, data =fdf, var = "rt_lag_sc", specs = c("entropy_change_late_beta_supp", "last_outcome", "rewFunc"), at = list(omission_early_theta = qmeg$omission_early_theta))) %>%   inner_join(qmeg, by = c("omission_early_theta", "rewFunc"))
emr2$study = '2. fMRI replication'
emr1 <- rbind(emr2, emr1)
ec_late_beta_rt <- 
  ggplot(emr1, aes(x=last_outcome, y=rt_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=as.factor(entropy_change_late_beta_supp))) + 
  geom_point(position = position_dodge(width = .6), size=2.5) + 
  geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + 
  theme_bw(base_size=12) + facet_grid(rewFunc~study) + ylab("RT swings (AU)\n Small <---------> Large")  + 
  scale_color_manual("Entropy change\nlate beta\nsuppression", values=c("#1b3840","#4fa3b8"), labels = c("10th %ile", "90th %ile")) +
  labs(shape = "Entropy change\nlate beta\nsuppression") +
  theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) + 
  scale_y_reverse(limits = c(.7, -.01)) 

setwd("~/OneDrive/collected_letters/papers/meg/plots/meg_to_fmri/")
pdf("meg_ec_late_beta_to_behavior_rtlbeta_adjusted_rewFunc.pdf", height = 8, width = 8)
print(ec_late_beta_rt)
dev.off()


# early theta to reward omission, with condition
oem1 <- as_tibble(emtrends(meg_om_theta_lbeta_rewFunc, data = df,  var = "rt_lag_sc", specs = c("omission_early_theta", "last_outcome", "rewFunc"), at = list(omission_early_theta = qmeg$omission_early_theta))) %>% 
  inner_join(qmeg, by = c("omission_early_theta", "rewFunc"))
oem1$study = "1. MEG"
oem2 <- as_tibble(emtrends(fmri_om_theta_lbeta_rewFunc, data =fdf, var = "rt_lag_sc", specs = c("omission_early_theta", "last_outcome", "rewFunc"), at = list(omission_early_theta = qmeg$omission_early_theta))) %>%  
  inner_join(qmeg, by = c("omission_early_theta", "rewFunc"))
oem2$study = '2. fMRI replication'
oem1 <- rbind(oem2, oem1)
omission_early_theta_rewFunc <- ggplot(oem1, aes(x=last_outcome, y=rt_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=q)) + 
  #shape = as.factor(pe_f2_hipp), 
  #p1 <- ggplot(em1, aes(x=last_outcome, y=rt_lag_sc.trend, linetype = as.factor(pe_f2_hipp), ymin=asymp.LCL, ymax=asymp.UCL)) + 
  geom_point(position = position_dodge(width = .6), size=2.5) + 
  #geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), position = position_dodge2(width=0.9), width=0.5) + 
  geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + 
  #geom_crossbar(position=position_dodge(width=0.7), width=0.6) + 
  theme_bw(base_size=12) + facet_grid(rewFunc~study)+ ylab("RT swings (AU)\n Small <---------> Large")  + 
  #scale_shape_manual(values=c(15,16), labels = c("10th %ile", "90th %ile")) + 
  #scale_color_brewer("PH RPE\nresponse", palette="Set1", labels = c("10th %ile", "90th %ile")) +
  scale_color_manual("Reward omission\nearly theta\nsynchronization", values = oem1$theta_color, labels = oem1$q) +
  labs(shape = "Reward omission\nearly theta\nsynchronization") +
  theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) + 
  scale_y_reverse(limits = c(.65, -.02)) 
pdf("meg_om_early_theta_to_behavior_lbeta_adjusted_rewFunc.pdf", height = 6, width = 5)
print(omission_early_theta_rewFunc)
dev.off()

# compile new results accounting for condition
setwd("~/OneDrive/collected_letters/papers/meg/plots/meg_to_behavior/")
pdf("meg_betas_to_behavior_cond_level_rewFunc.pdf", height = 6, width = 11)
print(ggarrange(omission_early_theta_rewFunc, ec_late_beta, nrow = 1, ncol = 2))
dev.off()


# early theta to reward omission, without condition - not much of a 2-way interaction indicating RT swings
rem1 <- as_tibble(emtrends(meg_om_theta_lbeta_rewFunc, data = df,  var = "rt_lag_sc", specs = c("omission_early_theta", "last_outcome"), at = list(omission_early_theta = c(-0.120, 0.661)), options = list()))
rem1$study = "1. MEG"
rem2 <- as_tibble(emtrends(fmri_om_theta_lbeta_rewFunc, data =fdf, var = "rt_lag_sc", specs = c("omission_early_theta", "last_outcome"), at = list(omission_early_theta = c(-0.120, 0.661)), options = list()))
rem2$study = '2. fMRI replication'
rem1 <- rbind(rem2, rem1)
omission_early_theta <- ggplot(rem1, aes(x=last_outcome, y=rt_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=as.factor(omission_early_theta))) + 
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
  scale_y_reverse(limits = c(.54, 0.075)) 
# rlem1 <- as_tibble(emtrends(meg1delta, data = df,  var = "rt_lag_sc", specs = c("omission_late_delta", "last_outcome"), at = list(omission_late_delta = c(-.3, .31)), options = list()))
# rlem1$study = "1. MEG"
# rlem2 <- as_tibble(emtrends(fmri1delta, data =fdf, var = "rt_lag_sc", specs = c("omission_late_delta", "last_outcome"), at = list(omission_late_delta = c(-.3, .31)), options = list()))
# rlem2$study = '2. fMRI replication'
# rlem1 <- rbind(rlem2, rlem1)


###############
# decomposition analyses

# lbeta
dem1 <- as_tibble(emtrends(meg_om_theta_lbeta_decomposed, data = df,  var = "rt_lag_sc", specs = c("entropy_change_late_beta_avg", "last_outcome", "rewFunc"), 
                           at = list(entropy_change_late_beta_avg = unique(qmeg$entropy_change_late_beta_avg)))) %>%  
  inner_join(qmeg, by = c("entropy_change_late_beta_avg", "rewFunc"))
dem1$study = "1. MEG"
dem2 <- as_tibble(emtrends(fmri_om_theta_lbeta_decomposed, data =fdf, var = "rt_lag_sc", specs = c("entropy_change_late_beta_avg", "last_outcome", "rewFunc"), 
                           at = list(entropy_change_late_beta_avg = unique(qmeg$entropy_change_late_beta_avg)))) %>%  
  inner_join(qmeg, by = c("entropy_change_late_beta_avg", "rewFunc"))
dem2$study = '2. fMRI replication'
dem1 <- rbind(dem2, dem1)
lbeta_decomposed <- 
  ggplot(dem1, aes(x=last_outcome, y=rt_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=q)) + 
  geom_point(position = position_dodge(width = .6), size=2.5) + 
  geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + 
  theme_bw(base_size=12) + facet_grid(rewFunc~study)+ ylab("RT swings (AU)\n Small <---------> Large")  + 
  scale_color_manual("Entropy change\nlate beta\nsuppression\n(across conditions)", values = dem1$qcolor)  +
  labs(shape = "Entropy change\nlate beta\nsuppression\n(across conditions)") +
  theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) + 
  scale_y_reverse(limits = c(.75, -.09)) 

wdem1 <- as_tibble(emtrends(meg_om_theta_lbeta_decomposed, data = df,  var = "rt_lag_sc", specs = c("ec_lbeta_wi", "last_outcome", "rewFunc"), 
                           at = list(ec_lbeta_wi = unique(qmeg$ec_lbeta_wi)))) %>%  
  inner_join(qmeg, by = c("ec_lbeta_wi", "rewFunc"))
wdem1$study = "1. MEG"
wdem2 <- as_tibble(emtrends(fmri_om_theta_lbeta_decomposed, data =fdf, var = "rt_lag_sc", specs = c("ec_lbeta_wi", "last_outcome", "rewFunc"), 
                           at = list(ec_lbeta_wi = unique(qmeg$ec_lbeta_wi)))) %>%  
  inner_join(qmeg, by = c("ec_lbeta_wi", "rewFunc"))
wdem2$study = '2. fMRI replication'
wdem1 <- rbind(wdem2, wdem1)
w_lbeta_decomposed <- 
  ggplot(wdem1, aes(x=last_outcome, y=rt_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=q)) + 
  geom_point(position = position_dodge(width = .6), size=2.5) + 
  geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + 
  theme_bw(base_size=12) + facet_grid(rewFunc~study)+ ylab("RT swings (AU)\n Small <---------> Large")  + 
  scale_color_manual("Entropy change\nlate beta\nsuppression\n(within condition)", values = dem1$qcolor)  +
  labs(shape = "Entropy change\nlate beta\nsuppression\n(within condition)") +
  theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) + 
  scale_y_reverse(limits = c(.9, -.09)) 


setwd("~/OneDrive/collected_letters/papers/meg/plots/meg_to_fmri/")
pdf("decomposed_meg_ec_lbeta_behavior.pdf", height = 8, width = 14)
ggarrange(lbeta_decomposed, w_lbeta_decomposed)
dev.off()


# early theta to reward omission, with condition
dom1 <- as_tibble(emtrends(meg_om_theta_lbeta_decomposed, data = df,  var = "rt_lag_sc", specs = c("omission_early_theta_avg", "last_outcome", "rewFunc"), 
                           at = list(omission_early_theta_avg = unique(qmeg$omission_early_theta_avg)))) %>%  
  inner_join(qmeg, by = c("omission_early_theta_avg", "rewFunc"))
dom1$study = "1. MEG"
dom2 <- as_tibble(emtrends(fmri_om_theta_lbeta_decomposed, data =fdf, var = "rt_lag_sc", specs = c("omission_early_theta_avg", "last_outcome", "rewFunc"), 
                           at = list(omission_early_theta_avg = unique(qmeg$omission_early_theta_avg)))) %>%  
  inner_join(qmeg, by = c("omission_early_theta_avg", "rewFunc"))
dom2$study = '2. fMRI replication'
dom1 <- rbind(dom2, dom1)
etheta_decomposed <- 
  ggplot(dom1, aes(x=last_outcome, y=rt_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=q)) + 
  geom_point(position = position_dodge(width = .6), size=2.5) + 
  geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + 
  theme_bw(base_size=12) + facet_grid(rewFunc~study)+ ylab("RT swings (AU)\n Small <---------> Large")  + 
  scale_color_manual("Omission\nearly theta\nsynchronization\n(across conditions)", values = dom1$theta_color)  +
  labs(shape = "Omission\nearly theta\nsynchronization\n(across conditions)") +
  theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) + 
  scale_y_reverse(limits = c(.75, -.09)) 

wdom1 <- as_tibble(emtrends(meg_om_theta_lbeta_decomposed, data = df,  var = "rt_lag_sc", specs = c("om_theta_wi", "last_outcome", "rewFunc"), 
                            at = list(om_theta_wi = unique(qmeg$om_theta_wi)))) %>%  
  inner_join(qmeg, by = c("om_theta_wi", "rewFunc"))
wdom1$study = "1. MEG"
wdom2 <- as_tibble(emtrends(fmri_om_theta_lbeta_decomposed, data =fdf, var = "rt_lag_sc", specs = c("om_theta_wi", "last_outcome", "rewFunc"), 
                            at = list(om_theta_wi = unique(qmeg$om_theta_wi)))) %>%  
  inner_join(qmeg, by = c("om_theta_wi", "rewFunc"))
wdom2$study = '2. fMRI replication'
wdom1 <- rbind(wdom2, wdom1)
w_etheta_decomposed <- 
  ggplot(wdom1, aes(x=last_outcome, y=rt_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=q)) + 
  geom_point(position = position_dodge(width = .6), size=2.5) + 
  geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + 
  theme_bw(base_size=12) + facet_grid(rewFunc~study)+ ylab("RT swings (AU)\n Small <---------> Large")  + 
  scale_color_manual("Omission\nearly theta\nsynchronization\n(within condition)", values = dom1$theta_color)  +
  labs(shape = "Omission\nearly theta\nsynchronization\n(within condition)") +
  theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) + 
  scale_y_reverse(limits = c(.9, -.09)) 

setwd("~/OneDrive/collected_letters/papers/meg/plots/meg_to_fmri/")
pdf("decomposed_meg_omission_theta_behavior.pdf", height = 8, width = 14)
ggarrange(etheta_decomposed, w_etheta_decomposed)
dev.off()

# early theta across conditions
across_conditions_1 <- as_tibble(emtrends(meg_om_theta_lbeta_decomposed, data = df,  var = "rt_lag_sc", specs = c("omission_early_theta_avg", "last_outcome"), 
                           at = list(omission_early_theta_avg = unique(qmeg$omission_early_theta_avg)))) %>%  
  inner_join(qmeg, by = c("omission_early_theta_avg"))
across_conditions_1$study = "1. MEG"
across_conditions_2 <- as_tibble(emtrends(fmri_om_theta_lbeta_decomposed, data =fdf, var = "rt_lag_sc", specs = c("omission_early_theta_avg", "last_outcome"), 
                           at = list(omission_early_theta_avg = unique(qmeg$omission_early_theta_avg)))) %>%  
  inner_join(qmeg, by = c("omission_early_theta_avg"))
across_conditions_2$study = '2. fMRI replication'
across_conditions_1 <- rbind(across_conditions_2, across_conditions_1)
etheta_decomposed <- 
  ggplot(across_conditions_1, aes(x=last_outcome, y=rt_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=q)) + 
  geom_point(position = position_dodge(width = .6), size=2.5) + 
  geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + 
  theme_bw(base_size=12) + facet_grid(rewFunc~study)+ ylab("RT swings (AU)\n Small <---------> Large")  + 
  scale_color_manual("Omission\nearly theta\nsynchronization\n(across conditions)", values = across_conditions_1$theta_color)  +
  labs(shape = "Omission\nearly theta\nsynchronization\n(across conditions)") +
  theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) + 
  scale_y_reverse(limits = c(.75, -.09)) 

wacross_conditions_1 <- as_tibble(emtrends(meg_om_theta_lbeta_decomposed, data = df,  var = "rt_lag_sc", specs = c("om_theta_wi", "last_outcome"), 
                            at = list(om_theta_wi = unique(qmeg$om_theta_wi)))) %>%  
  inner_join(qmeg, by = c("om_theta_wi"))
wacross_conditions_1$study = "1. MEG"
wacross_conditions_2 <- as_tibble(emtrends(fmri_om_theta_lbeta_decomposed, data =fdf, var = "rt_lag_sc", specs = c("om_theta_wi", "last_outcome"), 
                            at = list(om_theta_wi = unique(qmeg$om_theta_wi)))) %>%  
  inner_join(qmeg, by = c("om_theta_wi"))
wacross_conditions_2$study = '2. fMRI replication'
wacross_conditions_1 <- rbind(wacross_conditions_2, wacross_conditions_1)
w_etheta_decomposed <- 
  ggplot(wacross_conditions_1, aes(x=last_outcome, y=rt_lag_sc.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=q)) + 
  geom_point(position = position_dodge(width = .6), size=2.5) + 
  geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + 
  theme_bw(base_size=12) + facet_grid(rewFunc~study)+ ylab("RT swings (AU)\n Small <---------> Large")  + 
  scale_color_manual("Omission\nearly theta\nsynchronization\n(within condition)", values = across_conditions_1$theta_color)  +
  labs(shape = "Omission\nearly theta\nsynchronization\n(within condition)") +
  theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) + 
  scale_y_reverse(limits = c(.9, -.09)) 

setwd("~/OneDrive/collected_letters/papers/meg/plots/meg_to_fmri/")
pdf("decomposed_meg_omission_theta_behavior.pdf", height = 8, width = 14)
ggarrange(etheta_decomposed, w_etheta_decomposed)
dev.off()


setwd("~/OneDrive/collected_letters/papers/meg/plots/meg_to_fmri/")
pdf("decomposed_meg_betas_behavior.pdf", height = 8, width = 22)
ggarrange(etheta_decomposed, w_etheta_decomposed, lbeta_decomposed, w_lbeta_decomposed, ncol = 4, nrow = 1)
dev.off()


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
# vem2 <- as_tibble(emtrends(fmri1theta, data =fdf, var = "rt_lag_sc", specs = c("vmax_late_alpha", "last_outcome"), at = list(vmax_late_alpha = c(-.08, .08)), options = list()))
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


ec1 <- as_tibble(emtrends(meg_hybrid, data = df,  var = "rt_vmax_lag_sc", specs = c("entropy_change_late_beta_supp", "trial_neg_inv_sc"), at = list(entropy_change_late_beta_supp =  c(-0.183, 0.317), trial_neg_inv_sc = c(-.88, 0.39)), options = list()))
ec1$study = "1. MEG"
ec2 <- as_tibble(emtrends(fmri_hybrid, data =fdf, var = "rt_vmax_lag_sc", specs = c("entropy_change_late_beta_supp", "trial_neg_inv_sc"), at = list(entropy_change_late_beta_supp =  c(-0.183, 0.317), trial_neg_inv_sc = c(-.88, 0.39)), options = list()))
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
pe2 <- as_tibble(emtrends(fmri_hybrid, data =fdf, var = "rt_vmax_lag_sc", specs = c("pe_late_beta_supp", "trial_neg_inv_sc"), at = list(pe_late_beta_supp = c(-.00431, .00597), trial_neg_inv_sc = c(-.88, 0.39)), options = list()))
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
pe4 <- as_tibble(emtrends(fmri_hybrid, data =fdf, var = "rt_lag_sc", specs = c("pe_late_beta_supp", "last_outcome"), at = list(pe_late_beta_supp = c(-.00431, .00597)), options = list()))
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
                                   v_entropy_wi + entropy_change_late_beta_supp + 
                                   omission_early_theta)^2 + 
                      rt_lag_sc:last_outcome:entropy_change_late_beta_supp + 
                      rt_lag_sc:last_outcome:omission_early_theta +
                      rt_lag_sc:v_entropy_wi:entropy_change_late_beta_supp + 
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
                       (1|id/run),fdf %>% filter(rt_csv<4000))
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
                       (1|id/run),fdf %>% filter(rt_csv<4000))
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
                          (rt_lag_sc + rt_vmax_lag_sc|id/run),fdf %>% filter(rt_csv<4000))
# screen.lmerTest(fmri1_rs, .01)
summary(fmri1_rs)
Anova(fmri1_rs, '3')

# late beta to entropy change
em1 <- as_tibble(emtrends(meg1_rs, data = df,  var = "rt_lag_sc", specs = c("entropy_change_late_beta_supp", "last_outcome"), at = list(entropy_change_late_beta_supp = c(-.14, .12)), options = list()))
em1$study = "1. MEG"
em2 <- as_tibble(emtrends(fmri1_rs, data =fdf, var = "rt_lag_sc", specs = c("entropy_change_late_beta_supp", "last_outcome"), at = list(entropy_change_late_beta_supp = c(-.14, .12)), options = list()))
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
rem2 <- as_tibble(emtrends(fmri1_rs, data =fdf, var = "rt_lag_sc", specs = c("omission_early_theta", "last_outcome"), at = list(omission_early_theta = c(-.33, .3)), options = list()))
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



ggplot(df, aes(outcome, pe_max, color = rewFunc)) + geom_boxplot()
