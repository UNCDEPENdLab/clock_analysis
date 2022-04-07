# relating MEG to run-wise fMRI betas at network level by condition

library(tidyverse)
library(lme4)
library(ggpubr)
library(car)
library(viridis)
library(ggnewscale)
library(RColorBrewer)
library(emmeans)
source("~/code/Rhelpers/theme_black.R")
# install_github("UNCDEPENdLab/dependlab")
# library(dependlab)
source('~/code/Rhelpers/screen.lmerTest.R')
source('~/code/Rhelpers/vif.lme.R')
# library(stringi)

# clock_folder <- "~/Data_Analysis/clock_analysis" #michael
clock_folder <- "~/code/clock_analysis" #alex
# source('~/code/Rhelpers/')
fmri_dir <- '/Volumes/GoogleDrive/.shortcut-targets-by-id/1ukjK6kTlaR-LXIqX6nylYOPWu1j3XGyF/SCEPTIC_fMRI/wholebrain_betas'
source("~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/get_trial_data.R")

# get design
design <- get_trial_data(repo_directory = clock_folder, dataset = "mmclock_fmri", groupfixed = T) %>% select(id, run, rewFunc) %>% unique() %>%
  rename(run_number = run) %>% mutate(id = as.character(id))

# load meg data
# wbetas <- readRDS("~/OneDrive/collected_letters/papers/meg/plots/wholebrain/betas/MEG_betas_wide_echange_vmax_reward_Nov30_2021.RDS") %>% 
# add meg data, keep all the processing transparent here
# wbetas <- readRDS("~/OneDrive/collected_letters/papers/meg/plots/wholebrain/betas/MEG_betas_wide_echange_vmax_reward_Nov30_2021.RDS") %>% 
wbetas <- readRDS("~/code/clock_analysis/meg/data/MEG_betas_ec_rewfunc_rt_next_reward_rewfunc_April_5_2022.RDS") %>% 
  mutate(entropy_change_early_beta_supp = -  entropy_change_early_beta_ec_rewfunc,
         entropy_change_late_beta_supp = - entropy_change_late_beta_ec_rewfunc,
         rt_shorten_late_beta_supp = - rt_next_late_beta_rt_next,
         omission_early_theta = - reward_early_theta_reward_rewfunc) %>%
  select(id, rewFunc, entropy_change_early_beta_supp, entropy_change_late_beta_supp, rt_shorten_late_beta_supp, 
         omission_early_theta) 
omission_ave <- wbetas %>% group_by(id) %>% summarise(omission_early_theta_avg = mean(omission_early_theta)) %>%
  ungroup()
wbetas <- wbetas %>%  merge(omission_ave, by = "id")

# load fMRI betas: reward_omission, entropy, signed_pe
# reward/omission
setwd(file.path(fmri_dir, 'L1m-rew_om'))
rew_betas <- read_csv("Schaefer2018_200Parcels_7Networks_order_fonov_2.3mm_ants_cope_l1.csv.gz") %>% filter(l1_cope_name=="EV_rew_om") %>% 
  dplyr::select(id, run_number, l1_model, mask_value, value) %>% mutate(id = as.character(id),
                                                                        mask_value = as.factor(mask_value),
                                                                        run_mc  = scale(run_number, center = T, scale = F))

# merge
df <- rew_betas   %>% inner_join(design, by = c("id", "run_number")) %>% inner_join(wbetas, by = c("id", "rewFunc"))

labels <- read_delim("~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/Schaefer2018_200Parcels_7Networks_order_manual.txt", 
                     delim = "\t", escape_double = FALSE, 
                     col_names = FALSE, trim_ws = TRUE) %>% dplyr::select(1:2) %>% rename(mask_value = X1, label = X2) %>%
  mutate(hemi = stringr::str_extract(label, "_([^_]+)_"), 
         hemi = stringr::str_extract(hemi, "[^_]"),
         network = substr(label, 14, 17),
         mask_value = as.factor(mask_value))
df <- df %>% inner_join(labels, by = "mask_value")

# inspect
fbetas <- df %>% select(id, rewFunc, network, value) %>% group_by(id, rewFunc, network) %>% summarize(beta = mean(value)) %>%
  pivot_wider(names_from = c(rewFunc, network), values_from = beta) %>% ungroup() %>% select(where(is.numeric))
cormat <- corr.test(fbetas, method = "pearson")
# note negative correlations between CEVR and other betas
setwd("~/OneDrive/collected_letters/papers/meg/plots/meg_to_fmri/")
pdf("rew_om_betas_by_rewFunc.pdf", height = 10, width = 10)
corrplot(cormat$r, p.mat = cormat$p, order = "hclust", tl.cex = .8, insig = 'blank', method = 'number')
dev.off()

# basic model with ec late beta
m_rewom_lbeta <- lmer(value ~ rewFunc * network * entropy_change_late_beta_supp + (1|id), df)
while (any(grepl("failed to converge", m_rewom_lbeta@optinfo$conv$lme4$messages) )) {
ss <- getME(m_rewom_lbeta,c("theta","fixef"))
m_rewom_lbeta <- update(m_rewom_lbeta,start=ss,control=lmerControl(optCtrl=list(maxfun=2e4)))}
summary(m_rewom_lbeta)
Anova(m_rewom_lbeta, '3')

# plot
em_rewom_lbeta <- as_tibble(emmeans(m_rewom_lbeta, data = df, ~network|entropy_change_late_beta_supp*rewFunc , at = list(entropy_change_late_beta_supp = c(-0.1842, 0.3160))))
rewom_ec_lbeta <- ggplot(em_rewom_lbeta, aes(x=network, y=emmean, ymin=asymp.LCL, ymax=asymp.UCL, color=as.factor(entropy_change_late_beta_supp))) +
  geom_point(position = position_dodge(width = .6), size=2.5) +
  geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + facet_wrap(~rewFunc) +
  theme_bw(base_size=12) +  ylab("BOLD response to reward > omission")  +
  scale_color_manual("Entropy change\nlate beta\nsuppression", values=c("#1b3840","#4fa3b8"), labels = c("10th %ile", "90th %ile")) +
  labs(shape = "Entropy change\nlate beta\nsuppression") +
  theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) # +

# basic model for early theta
m_rewom_etheta <- lmer(value ~ rewFunc * network * omission_early_theta + (1|id), df)
while (any(grepl("failed to converge", m_rewom_etheta@optinfo$conv$lme4$messages) )) {
  ss <- getME(m_rewom_etheta,c("theta","fixef"))
  m_rewom_etheta <- update(m_rewom_etheta,start=ss,control=lmerControl(optCtrl=list(maxfun=2e4)))}
summary(m_rewom_etheta)
Anova(m_rewom_etheta, '3')

em_rewom_etheta <- as_tibble(emmeans(m_rewom_etheta, data = df, ~network|omission_early_theta*rewFunc , at = list(omission_early_theta = c(-0.13878, 0.66119)))) 
rewom_rewom_theta <- ggplot(em_rewom_etheta, aes(x=network, y=emmean, ymin=asymp.LCL, ymax=asymp.UCL, color=as.factor(omission_early_theta))) +
  geom_point(position = position_dodge(width = .6), size=2.5) +
  geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + facet_wrap(~rewFunc) +
  theme_bw(base_size=12) +  ylab("BOLD response to reward > omission")  +
  scale_color_manual("Reward omission\nearly theta\nsynchronization", values=c("#1b3840","#4fa3b8"), labels = c("10th %ile", "90th %ile")) +
  labs(shape = "Reward omission\nearly theta\nsynchronization") +
  theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) # +


# Prediction errors
setwd(file.path(fmri_dir, 'L1m-pe'))
pe_betas <- read_csv("Schaefer2018_200Parcels_7Networks_order_fonov_2.3mm_ants_cope_l1.csv.gz") %>% filter(l1_cope_name=="EV_pe") %>% 
  dplyr::select(id, run_number, l1_model, mask_value, value) %>% mutate(id = as.character(id),
                                                                        mask_value = as.factor(mask_value),
                                                                        run_mc  = scale(run_number, center = T, scale = F))
# merge
df <- pe_betas   %>% inner_join(design, by = c("id", "run_number")) %>% 
  inner_join(wbetas, by = c("id", "rewFunc")) %>% inner_join(labels, by = "mask_value")

# inspect
fbetas <- df %>% select(id, rewFunc, network, value) %>% group_by(id, rewFunc, network) %>% summarize(beta = mean(value)) %>%
  pivot_wider(names_from = c(rewFunc, network), values_from = beta) %>% ungroup() %>% select(where(is.numeric))
cormat <- corr.test(fbetas, method = "pearson")
# note negative correlations between CEVR and other betas
setwd("~/OneDrive/collected_letters/papers/meg/plots/meg_to_fmri/")
pdf("pe_betas_by_rewFunc.pdf", height = 10, width = 10)
corrplot(cormat$r, p.mat = cormat$p, order = "hclust", tl.cex = .8, insig = 'blank', method = 'number')
dev.off()
# plotbetas <- df %>% select(id, rewFunc, network, value) %>% group_by(id, rewFunc, network) %>% summarize(beta = mean(value))
# ggplot(plotbetas, aes(network, beta, color = rewFunc)) + geom_violin() + geom_jitter(alpha = .1)
# PE beta prediction
m_pe_lbeta <- lmer(value ~ rewFunc * network * entropy_change_late_beta_supp + (1|id), df)
while (any(grepl("failed to converge", m_pe_lbeta@optinfo$conv$lme4$messages) )) {
  ss <- getME(m_pe_lbeta,c("theta","fixef"))
  m_pe_lbeta <- update(m_pe_lbeta,start=ss,control=lmerControl(optCtrl=list(maxfun=2e4)))}
summary(m_pe_lbeta)
Anova(m_pe_lbeta, '3')
em_pe_lbeta <- as_tibble(emmeans(m_pe_lbeta, data = df, ~network|entropy_change_late_beta_supp*rewFunc , at = list(entropy_change_late_beta_supp = c(-0.1842, 0.3160))))
pe_ec_lbeta <- ggplot(em_pe_lbeta, aes(x=network, y=emmean, ymin=asymp.LCL, ymax=asymp.UCL, color=as.factor(entropy_change_late_beta_supp))) +
  geom_point(position = position_dodge(width = .6), size=2.5) +
  geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + facet_wrap(~rewFunc) +
  theme_bw(base_size=12) +  ylab("BOLD response to reward prediction error")  +
  scale_color_manual("Entropy change\nlate beta\nsuppression", values=c("#1b3840","#4fa3b8"), labels = c("10th %ile", "90th %ile")) +
  labs(shape = "Entropy change\nlate beta\nsuppression") +
  theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) # +

# predicting PE BOLD with early theta
# basic model for early theta
m_pe_etheta <- lmer(value ~ rewFunc * network * omission_early_theta + (1|id), df)
while (any(grepl("failed to converge", m_pe_etheta@optinfo$conv$lme4$messages) )) {
  ss <- getME(m_pe_etheta,c("theta","fixef"))
  m_pe_etheta <- update(m_pe_etheta,start=ss,control=lmerControl(optCtrl=list(maxfun=2e4)))}
summary(m_pe_etheta)
Anova(m_pe_etheta, '3')

em_pe_etheta <- as_tibble(emmeans(m_pe_etheta, data = df, ~network|omission_early_theta*rewFunc , at = list(omission_early_theta = c(-0.13878, 0.66119)))) 
pe_rewom_theta <- ggplot(em_pe_etheta, aes(x=network, y=emmean, ymin=asymp.LCL, ymax=asymp.UCL, color=as.factor(omission_early_theta))) +
  geom_point(position = position_dodge(width = .6), size=2.5) +
  geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + facet_wrap(~rewFunc) +
  theme_bw(base_size=12) +  ylab("BOLD response to reward prediction error")  +
  scale_color_manual("Reward omission\nearly theta\nsynchronization", values=c("#1b3840","#4fa3b8"), labels = c("10th %ile", "90th %ile")) +
  labs(shape = "Reward omission\nearly theta\nsynchronization") +
  theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) # +



# read in entropy change betas
setwd(file.path(fmri_dir, 'L1m-entropy_echange'))
ec_betas <- read_csv("Schaefer2018_200Parcels_7Networks_order_fonov_2.3mm_ants_cope_l1.csv.gz") %>% filter(l1_cope_name=="EV_entropy_change_feedback") %>% 
  dplyr::select(id, run_number, l1_model, mask_value, value) %>% mutate(id = as.character(id),
                                                                        mask_value = as.factor(mask_value),
                                                                        run_mc  = scale(run_number, center = T, scale = F))
# merge
df <- ec_betas   %>% inner_join(design, by = c("id", "run_number"))  %>% 
  inner_join(wbetas, by = c("id", "rewFunc")) %>% inner_join(labels, by = "mask_value")
setwd("~/OneDrive/collected_letters/papers/meg/plots/meg_to_fmri/")

# inspect
fbetas <- df %>% select(id, rewFunc, network, value) %>% group_by(id, rewFunc, network) %>% summarize(beta = mean(value)) %>%
  pivot_wider(names_from = c(rewFunc, network), values_from = beta) %>% ungroup() %>% select(where(is.numeric))
cormat <- corr.test(fbetas, method = "pearson")

pdf("ec_betas_by_rewFunc.pdf", height = 10, width = 10)
corrplot(cormat$r, p.mat = cormat$p, order = "hclust", tl.cex = .8, insig = 'blank', method = 'number')
dev.off()

# not much for late beta outside of DMN
m_ec_lbeta <- lmer(value ~ rewFunc *  network * entropy_change_late_beta_supp + (1|id), df)
while (any(grepl("failed to converge", m_ec_lbeta@optinfo$conv$lme4$messages) )) {
  ss <- getME(m_ec_lbeta,c("theta","fixef"))
  m_ec_lbeta <- update(m_ec_lbeta,start=ss,control=lmerControl(optCtrl=list(maxfun=2e4)))}

summary(m_ec_lbeta)
Anova(m_ec_lbeta, '3')


# more for early beta in DMN and somatomotor, visual
m_ec_ebeta <- lmer(value ~ rewFunc * network * entropy_change_early_beta_supp + (1|id), df)
while (any(grepl("failed to converge", m_ec_ebeta@optinfo$conv$lme4$messages) )) {
  ss <- getME(m_ec_ebeta,c("theta","fixef"))
  m_ec_ebeta <- update(m_ec_ebeta,start=ss,control=lmerControl(optCtrl=list(maxfun=2e4)))}
summary(m_ec_ebeta)
Anova(m_ec_ebeta, '3')
# df$entropy_change_late_beta_supp 
# n  missing distinct     Info     Mean      Gmd      .05      .10      .25      .50      .75      .90      .95 
# 92800        0       58        1 -0.00563   0.1265 -0.18217 -0.14210 -0.07790 -0.01109  0.08370  0.15427  0.21376 

em_ec_lbeta <- as_tibble(emmeans(m_ec_lbeta, data = df, ~network|entropy_change_late_beta_supp*rewFunc , at = list(entropy_change_late_beta_supp = c(-0.1842, 0.3160))))
ec_ec_lbeta <- ggplot(em_ec_lbeta, aes(x=network, y=emmean, ymin=asymp.LCL, ymax=asymp.UCL, color=as.factor(entropy_change_late_beta_supp))) +
  geom_point(position = position_dodge(width = .6), size=2.5) +
  geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + facet_wrap(~rewFunc) +
  theme_bw(base_size=12) +  ylab("BOLD response to entropy change")  +
  scale_color_manual("Entropy change\nlate beta\nsuppression", values=c("#1b3840","#4fa3b8"), labels = c("10th %ile", "90th %ile")) +
  labs(shape = "Entropy change\nlate beta\nsuppression") +
  theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) # +

# basic model for early theta
m_ec_etheta <- lmer(value ~ rewFunc * network * omission_early_theta + (1|id), df)
while (any(grepl("failed to converge", m_ec_etheta@optinfo$conv$lme4$messages) )) {
  ss <- getME(m_ec_etheta,c("theta","fixef"))
  m_ec_etheta <- update(m_ec_etheta,start=ss,control=lmerControl(optCtrl=list(maxfun=2e4)))}
summary(m_ec_etheta)
Anova(m_ec_etheta, '3')

em_ec_etheta <- as_tibble(emmeans(m_ec_etheta, data = df, ~network|omission_early_theta*rewFunc , at = list(omission_early_theta = c(-0.13878, 0.66119)))) 
ec_rewom_theta <- ggplot(em_ec_etheta, aes(x=network, y=emmean, ymin=asymp.LCL, ymax=asymp.UCL, color=as.factor(omission_early_theta))) +
  geom_point(position = position_dodge(width = .6), size=2.5) +
  geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + facet_wrap(~rewFunc) +
  theme_bw(base_size=12) +  ylab("BOLD response to entropy change")  +
  scale_color_manual("Reward omission\nearly theta\nsynchronization", values=c("#1b3840","#4fa3b8"), labels = c("10th %ile", "90th %ile")) +
  labs(shape = "Reward omission\nearly theta\nsynchronization") +
  theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) # +



# read in entropy EV_entropy_wiz_clock
setwd(file.path(fmri_dir, 'L1m-entropy_wiz'))
e_betas <- read_csv("Schaefer2018_200Parcels_7Networks_order_fonov_2.3mm_ants_cope_l1.csv.gz") %>% filter(l1_cope_name=="EV_entropy_wiz_clock") %>% 
  dplyr::select(id, run_number, l1_model, mask_value, value) %>% mutate(id = as.character(id),
                                                                        mask_value = as.factor(mask_value),
                                                                        run_mc  = scale(run_number, center = T, scale = F))
# merge
df <- e_betas   %>% inner_join(design, by = c("id", "run_number"))  %>% 
  inner_join(wbetas, by = c("id", "rewFunc")) %>% inner_join(labels, by = "mask_value")

# late beta: again, mostly DMN
m_e_lbeta <- lmer(value ~ rewFunc * network * entropy_change_late_beta_supp + (1|id), df)
while (any(grepl("failed to converge", m_e_lbeta@optinfo$conv$lme4$messages) )) {
  ss <- getME(m_e_lbeta,c("theta","fixef")) 
  m_e_lbeta <- update(m_e_lbeta,start=ss,control=lmerControl(optCtrl=list(maxfun=2e4)))}

summary(m_e_lbeta)
Anova(m_e_lbeta, '3')
em_e_lbeta <- as_tibble(emmeans(m_e_lbeta, data = df, ~network|entropy_change_late_beta_supp*rewFunc , at = list(entropy_change_late_beta_supp = c(-0.1842, 0.3160))))
e_ec_lbeta <- ggplot(em_e_lbeta, aes(x=network, y=emmean, ymin=asymp.LCL, ymax=asymp.UCL, color=as.factor(entropy_change_late_beta_supp))) +
  geom_point(position = position_dodge(width = .6), size=2.5) +
  geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + facet_wrap(~rewFunc) +
  theme_bw(base_size=12) +  ylab("BOLD response to entropy, clock-aligned")  +
  scale_color_manual("Entropy change\nlate beta\nsuppression", values=c("#1b3840","#4fa3b8"), labels = c("10th %ile", "90th %ile")) +
  labs(shape = "Entropy change\nlate beta\nsuppression") +
  theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) # +

# basic model for early theta
m_e_etheta <- lmer(value ~ rewFunc * network * omission_early_theta + (1|id), df)
while (any(grepl("failed to converge", m_e_etheta@optinfo$conv$lme4$messages) )) {
  ss <- getME(m_e_etheta,c("theta","fixef"))
  m_e_etheta <- update(m_e_etheta,start=ss,control=lmerControl(optCtrl=list(maxfun=2e4)))}
summary(m_e_etheta)
Anova(m_e_etheta, '3')

em_e_etheta <- as_tibble(emmeans(m_e_etheta, data = df, ~network|omission_early_theta*rewFunc , at = list(omission_early_theta = c(-0.13878, 0.66119)))) 
e_rewom_theta <- ggplot(em_e_etheta, aes(x=network, y=emmean, ymin=asymp.LCL, ymax=asymp.UCL, color=as.factor(omission_early_theta))) +
  geom_point(position = position_dodge(width = .6), size=2.5) +
  geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + facet_wrap(~rewFunc) +
  theme_bw(base_size=12) +  ylab("BOLD response to entropy")  +
  scale_color_manual("Reward omission\nearly theta\nsynchronization", values=c("#1b3840","#4fa3b8"), labels = c("10th %ile", "90th %ile")) +
  labs(shape = "Reward omission\nearly theta\nsynchronization") +
  theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) # +

setwd("~/OneDrive/collected_letters/papers/meg/plots/meg_to_fmri/")


pdf("MEG_to_fMRI_runwise_betas_rewFunc.pdf", height = 10, width = 32)
ggarrange(rewom_rewom_theta, pe_rewom_theta, ec_rewom_theta, e_rewom_theta,
  rewom_ec_lbeta, pe_ec_lbeta, ec_ec_lbeta, e_ec_lbeta,
          ncol = 4, nrow = 2)
dev.off()


# prototype the emtrends version, first for EC/EC_lbeta

# emt6 <- as_tibble(emtrends(m6, data = df, var = "entropy_change_late_beta_supp",   ~network|run_number, at = list(run_number = c(1,8)))) %>% 
#   mutate(run_number = as.character(run_number))
# ec_ec_lbeta <- ggplot(emt6, aes(x=network, y=entropy_change_late_beta_supp.trend, group = run_number, ymin=asymp.LCL, ymax=asymp.UCL, color = run_number)) +
#   geom_point(position = position_dodge(width = .6), size=2.5) + geom_line(position = position_dodge(width = .6)) +
#   geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + 
#   theme_bw(base_size=12) +  ylab("effect of oscilation response on BOLD response to entropy change")  +
#   scale_color_manual("Run", values=c("#1b3840","#4fa3b8"), labels = c("1", "8")) +
#   labs(shape = "Run") +
#   theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
#         axis.text=element_text(size=8.5, color="grey10")) # +


