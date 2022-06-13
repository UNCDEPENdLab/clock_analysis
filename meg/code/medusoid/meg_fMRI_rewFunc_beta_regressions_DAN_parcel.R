# relating MEG to run-wise DAN fMRI betas at parcel level, by condition

library(tidyverse)
library(lme4)
library(ggpubr)
library(car)
library(viridis)
library(ggnewscale)
library(RColorBrewer)
library(emmeans)
library(readxl)
library(readr)
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

decompose = T # whether to decompose MEG betas into session and condition levels

# get design
design <- get_trial_data(repo_directory = clock_folder, dataset = "mmclock_fmri", groupfixed = T) %>% select(id, run, rewFunc) %>% unique() %>%
  rename(run_number = run) %>% mutate(id = as.character(id),
                                      learnable = case_when(
                                        rewFunc == "IEV" | rewFunc == "DEV" ~ "learnable",
                                        rewFunc == "CEV" | rewFunc == "CEVR" ~ "unlearnable"
                                      ))

design_face <- get_trial_data(repo_directory = clock_folder, dataset = "mmclock_fmri", groupfixed = T) %>% select(id, run, rewFunc, emotion) %>% unique() %>%
  rename(run_number = run) %>% mutate(id = as.character(id),
                                      learnable = case_when(
                                        rewFunc == "IEV" | rewFunc == "DEV" ~ "learnable",
                                        rewFunc == "CEV" | rewFunc == "CEVR" ~ "unlearnable"
                                      ))

# get fMRI parcel labels for 400
setwd("~/code/schaefer_wb_parcellation")
schaefer_7 <- read.csv("labels/Schaefer2018_400Parcels_7Networks_order.csv") %>%
  mutate(network=factor(network), net_num = as.numeric(network)) %>%
  rename(network7=network, net_num7=net_num)

# this has the spatial coordinate, spatial_roi_num
schaefer_7_lookup <- read.csv("labels/Schaefer_400_7networks_labels.csv")

schaefer_7 <- schaefer_7 %>% inner_join(schaefer_7_lookup, by="roi_num") %>%
  rename(roi_num7=roi_num, subregion7=subregion)

schaefer_17 <- read.csv("labels/Schaefer2018_400Parcels_17Networks_order.csv") %>%
  mutate(network=factor(network), net_num = as.numeric(network)) %>%
  rename(network17=network, net_num17=net_num) %>%
  select(-hemi) # mirrored in 7

# this has the spatial coordinate, spatial_roi_num
schaefer_17_lookup <- read.csv("labels/Schaefer_400_17networks_labels.csv") %>%
  select(roi_num, spatial_roi_num) # x,y,z and labels already duplicated in 7-network lookup

schaefer_17 <- schaefer_17 %>% inner_join(schaefer_17_lookup, by="roi_num") %>%
  rename(roi_num17=roi_num, subregion17=subregion)

both <- inner_join(schaefer_7, schaefer_17, by="spatial_roi_num") %>%
  select(spatial_roi_num, roi_num7, roi_num17, network7, network17, net_num7, net_num17, subregion7, subregion17, everything())
setDT(both)
labels <- both %>% filter(network7=="DorsAttn" & (network17=="DorsAttnA" | network17=="DorsAttnB")) %>% 
  mutate(roi_num7 = as.factor(roi_num7)) %>% 
  # label lobes
  mutate(lobe = case_when(
    str_detect(subregion17, "Temp") ~ "temporal",
    str_detect(subregion17, "Par") | str_detect(subregion17, "SPL") | str_detect(subregion17, "PostC") |
      str_detect(subregion17, "IPS") | str_detect(subregion17, "IPL") | str_detect(subregion17, "pCun") ~ "parietal",
    str_detect(subregion17, "PFC") | str_detect(subregion17, "FEF") | str_detect(subregion17, "PrCv") ~ "frontal"),
    vm_gradient17 = case_when(
      lobe == "temporal" ~ "MT+",
      lobe == "parietal" & network17 == "DorsAttnA" ~ "PPCcaudal",
      lobe == "parietal" & network17 == "DorsAttnB" ~ "PPCrostral",
      lobe == "frontal" ~ "premotor",
      TRUE ~ as.character(network17))
  )


# load meg data
# wbetas <- readRDS("~/OneDrive/collected_letters/papers/meg/plots/wholebrain/betas/MEG_betas_wide_echange_vmax_reward_Nov30_2021.RDS") %>% 

if (decompose) {
  wbetas <- readRDS("~/code/clock_analysis/meg/data/MEG_betas_ec_rewfunc_rt_next_reward_rewfunc_April_5_2022.RDS") %>% 
    mutate(entropy_change_early_beta_supp = -  entropy_change_early_beta_ec_rewfunc,
           entropy_change_late_beta_supp = - entropy_change_late_beta_ec_rewfunc,
           rt_shorten_late_beta_supp = - rt_next_late_beta_rt_next,
           omission_early_theta = - reward_early_theta_reward_rewfunc) %>%
    select(id, rewFunc, entropy_change_early_beta_supp, entropy_change_late_beta_supp, rt_shorten_late_beta_supp, 
           omission_early_theta) 
  avg <- wbetas %>% group_by(id) %>% summarise(omission_early_theta_avg = mean(omission_early_theta),
                                                    entropy_change_late_beta_avg = mean(entropy_change_late_beta_supp)) %>%
    ungroup()
  
  wbetas <- wbetas %>%  merge(avg, by = "id") %>% mutate(
    ec_lbeta_wi = entropy_change_late_beta_supp - entropy_change_late_beta_avg,
    om_theta_wi = omission_early_theta - omission_early_theta_avg)
  
  # inspect
  
  qmeg <- wbetas %>% group_by(rewFunc) %>% summarize(omission_early_theta = quantile(omission_early_theta, c(.1, .9), names = F),
                                                          q = c("10th %ile", "90th %ile"),
                                                          qcolor = c("#1b3840","#4fa3b8"),
                                                          theta_color = c("orange4", "orange"),
                                                          entropy_change_late_beta_supp = quantile(entropy_change_late_beta_supp, c(.1, .9), names = F),
                                                          entropy_change_late_beta_avg = quantile(entropy_change_late_beta_avg, c(.1, .9), names = F),
                                                          omission_early_theta_avg = quantile(omission_early_theta_avg, c(.1, .9), names = F),
                                                          ec_lbeta_wi = quantile(ec_lbeta_wi, c(.1, .9), names = F),
                                                          om_theta_wi = quantile(om_theta_wi, c(.1, .9), names = F),
  ) %>% ungroup()
} else {
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
  
  # get MEG quantiles by condition
  qmeg <- wbetas %>% group_by(rewFunc) %>% summarize(omission_early_theta = quantile(omission_early_theta, c(.1, .9)),
                                                     q = c("10th %ile", "90th %ile"),
                                                     qcolor = c("#1b3840","#4fa3b8"),
                                                     theta_color = c("orange4", "orange"),
                                                     entropy_change_late_beta_supp = quantile(entropy_change_late_beta_supp, c(.1, .9))) %>% ungroup()
}
# load fMRI betas: reward_omission, entropy, signed_pe
# reward/omission

setwd(file.path(fmri_dir, 'L1m-rew_om'))
rew_betas <- read_csv("Schaefer2018_200Parcels_7Networks_order_fonov_2.3mm_ants_cope_l1.csv.gz") %>% filter(l1_cope_name=="EV_rew_om") %>% 
  dplyr::select(id, run_number, l1_model, mask_value, value) %>% mutate(id = as.character(id),
                                                                        mask_value = as.character(mask_value),
                                                                        run_mc  = scale(run_number, center = T, scale = F))
# merge
df <- rew_betas %>% inner_join(wbetas, by = "id", "rewFunc") %>% inner_join(design_face, by = c("id", "run_number", "rewFunc")) 

labels <- read_delim("~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/Schaefer2018_200Parcels_7Networks_order_manual.txt", 
                     delim = "\t", escape_double = FALSE, 
                     col_names = FALSE, trim_ws = TRUE) %>% dplyr::select(1:2) %>% rename(mask_value = X1, label = X2) %>%
  mutate(hemi = stringr::str_extract(label, "_([^_]+)_"), 
         hemi = stringr::str_extract(hemi, "[^_]"),
         network = substr(label, 14, 17),
         mask_value = as.factor(mask_value))
dan_labels <- read_excel("~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/MNH Dan Labels.xlsx") %>% 
  mutate(mask_value = as.character(roinum),
         parcel = word(MNHLabel, 2, sep = "_")) %>% select(mask_value, MNHLabel, parcel, Stream, Visuomotor_Gradient, lobe)
df <- df %>% inner_join(labels, by = "mask_value") %>% filter(network=="Dors") %>% inner_join(dan_labels, by = "mask_value")

# inspect
rew_dist <-  ggplot(df, aes(rewFunc, value, color = rewFunc)) + geom_violin(draw_quantiles = .5) + stat_summary(fun=mean, geom="point") + ylab("Reward>omission fMRI beta")
rew_means <- ggplot(df, aes(rewFunc, value, color = rewFunc)) + stat_summary(fun=mean, geom="point") + ylab("Reward>omission fMRI beta, mean")

summary(lm(value ~ rewFunc * Stream, df %>% filter(emotion == "scram")))
summary(lm(value ~ rewFunc * Stream * run_mc, df))

# models with early theta
m_rewom_etheta <- lmer(value ~ rewFunc *  Stream * omission_early_theta  +  (1|id), df)
summary(m_rewom_etheta)
Anova(m_rewom_etheta, '3')

m_rewom_etheta_l <- lmer(value ~ learnable *  Stream * omission_early_theta  +  (1|id), df)
summary(m_rewom_etheta)
Anova(m_rewom_etheta_l, '3')
anova(m_rewom_etheta, m_rewom_etheta_l) # individual conditions are much stronger

# em_rewom_etheta <- as_tibble(emmeans(m_rewom_etheta, data = df, ~Stream|omission_early_theta*rewFunc, at = list(omission_early_theta = c(-0.13878, 0.66119))))
# sample MEG at condition-wise quantiles
em_rewom_etheta <- as_tibble(emmeans(m_rewom_etheta, data = df, ~Stream|omission_early_theta*rewFunc, at = list(omission_early_theta = qmeg$omission_early_theta))) %>%
  inner_join(qmeg, by = c("omission_early_theta", "rewFunc"))

rewom_rewom_etheta <- ggplot(em_rewom_etheta, aes(x = Stream, y=emmean, ymin=asymp.LCL, ymax=asymp.UCL, lty = rewFunc, color=as.factor(omission_early_theta))) +
  geom_point(position = position_dodge(width = .6), size=2.5) +
  geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + facet_wrap(~(rewFunc=="DEV" | rewFunc=="IEV")) +
  theme_bw(base_size=12) +  ylab("BOLD response to reward > omission")  +
  scale_color_manual("Reward omission\nearly theta\nsynchronization", values = em_rewom_etheta$theta_color, labels = em_rewom_etheta$q) +
  labs(shape = "Reward omission\nearly theta\nsynchronization") +
  theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) # +

ggplot(em_rewom_etheta, aes(x = Stream, y=emmean, ymin=asymp.LCL, ymax=asymp.UCL, color = rewFunc, lty=as.factor(q))) +
  geom_point(position = position_dodge(width = .6), size=2.5) +
  geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + facet_wrap(~(rewFunc=="DEV" | rewFunc=="IEV")) +
  theme_bw(base_size=12) +  ylab("BOLD response to reward > omission")  +
  # scale_color_manual("Reward omission\nearly theta\nsynchronization", values = em_rewom_etheta$theta_color, labels = em_rewom_etheta$q) +
  labs(lty = "Reward omission\nearly theta\nsynchronization") +
  theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) # +



# models with late beta
# streams
m_rewom_lbeta <- lmer(value ~ rewFunc *  Stream * entropy_change_late_beta_supp + (1|id), df)
summary(m_rewom_lbeta)
Anova(m_rewom_lbeta, '3')
# # thankfully, no major face effects
# m_rewom_lbetas <- lmer(value ~ rewFunc *  Stream * entropy_change_late_beta_supp + (1|id), df %>% filter(emotion == "scram"))
# summary(m_rewom_lbetas)
# Anova(m_rewom_lbetas, '3')
# # control for emotion
# m_rewom_lbetae <- lmer(value ~ rewFunc *  Stream * entropy_change_late_beta_supp + emotion * Stream * entropy_change_late_beta_supp + (1|id), df)
# summary(m_rewom_lbetae)
# Anova(m_rewom_lbetae, '3')

# interaction with run?
m_rewom_lbeta_run <- lmer(value ~ rewFunc *  Stream * entropy_change_late_beta_supp * run_mc + (1|id), df)
summary(m_rewom_lbeta_run)
Anova(m_rewom_lbeta_run, '3')


# em_rewom_lbeta <- as_tibble(emmeans(m_rewom_lbeta, data = df, ~Stream|entropy_change_late_beta_supp*rewFunc, at = list(entropy_change_late_beta_supp = c(-0.07553, 0.22084))))
em_rewom_lbeta <- as_tibble(emmeans(m_rewom_lbeta, data = df, ~Stream|entropy_change_late_beta_supp*rewFunc, 
                                    at = list(entropy_change_late_beta_supp = qmeg$entropy_change_late_beta_supp))) %>%
  inner_join(qmeg, by = c("entropy_change_late_beta_supp", "rewFunc"))


rewom_ec_lbeta <- ggplot(em_rewom_lbeta, aes(x = Stream, y=emmean, ymin=asymp.LCL, ymax=asymp.UCL, color=as.factor(entropy_change_late_beta_supp))) +
  geom_point(position = position_dodge(width = .6), size=2.5) +
  geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + facet_wrap(~rewFunc) +
  theme_bw(base_size=12) +  ylab("BOLD response to reward > omission")  +
  scale_color_manual("Entropy change\nlate beta\nsuppression", values = em_rewom_lbeta$qcolor, labels = em_rewom_lbeta$q)  +
  labs(shape = "Entropy change\nlate beta\nsuppression") +
  theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) # +

# # scrambled only
# em_rewom_lbetas <- as_tibble(emmeans(m_rewom_lbetas, data = df, ~Stream|entropy_change_late_beta_supp*rewFunc, at = list(entropy_change_late_beta_supp = c(-0.07553, 0.22084))))
# rewom_ec_lbeta_s <- ggplot(em_rewom_lbetas, aes(x = Stream, y=emmean, ymin=asymp.LCL, ymax=asymp.UCL, color=as.factor(entropy_change_late_beta_supp))) +
#   geom_point(position = position_dodge(width = .6), size=2.5) +
#   geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + facet_wrap(~rewFunc) +
#   theme_bw(base_size=12) +  ylab("BOLD response to reward > omission")  +
#   scale_color_manual("Entropy change\nlate beta\nsuppression\nexcluding faces", values=c("#1b3840","#4fa3b8"), labels = c("10th %ile", "90th %ile")) +
#   labs(shape = "Entropy change\nlate beta\nsuppression\nexcluding faces") +
#   theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
#         axis.text=element_text(size=8.5, color="grey10")) # +
# # co-vary for emotion
# em_rewom_lbetae <- as_tibble(emmeans(m_rewom_lbetae, data = df, ~Stream|entropy_change_late_beta_supp*rewFunc, at = list(entropy_change_late_beta_supp = c(-0.07553, 0.22084))))
# rewom_ec_lbeta_e <- ggplot(em_rewom_lbetae, aes(x = Stream, y=emmean, ymin=asymp.LCL, ymax=asymp.UCL, color=as.factor(entropy_change_late_beta_supp))) +
#   geom_point(position = position_dodge(width = .6), size=2.5) +
#   geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + facet_wrap(~rewFunc) +
#   theme_bw(base_size=12) +  ylab("BOLD response to reward > omission")  +
#   scale_color_manual("Entropy change\nlate beta\nsuppression", values=c("#1b3840","#4fa3b8"), labels = c("10th %ile", "90th %ile")) +
#   labs(shape = "Entropy change\nlate beta\nsuppression") +
#   theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
#         axis.text=element_text(size=8.5, color="grey10")) # +
# 
# pdf("effects_of_emotion_rewom_lbeta.pdf", height = 6, width = 18)
# ggarrange(rewom_ec_lbeta, rewom_ec_lbeta_s, rewom_ec_lbeta_e, 
#           labels = c("original", "scrambled only", "covary for emotion"), nrow = 1, ncol = 3)
# dev.off()
setwd("~/OneDrive/collected_letters/papers/meg/plots/meg_to_fmri/")
pdf("MEG_on_rewom_fMRI_beta_only.pdf", height = 12, width = 10)
ggarrange(rewom_rewom_etheta, rewom_ec_lbeta, nrow = 2, ncol = 1)
dev.off()


# Prediction errors
setwd(file.path(fmri_dir, 'L1m-pe'))
pe_betas <- read_csv("Schaefer2018_200Parcels_7Networks_order_fonov_2.3mm_ants_cope_l1.csv.gz") %>% filter(l1_cope_name=="EV_pe") %>% 
  dplyr::select(id, run_number, l1_model, mask_value, value) %>% mutate(id = as.character(id),
                                                                        mask_value = as.factor(mask_value),
                                                                        run_mc  = scale(run_number, center = T, scale = F),
                                                                        beta_winsor = psych::winsor(value, trim = .05))
# merge
df <- pe_betas %>% inner_join(design, by = c("id", "run_number")) %>%
  inner_join(wbetas, by = c("id", "rewFunc"))
# hist(df$beta_winsor)

labels <- read_delim("~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/Schaefer2018_200Parcels_7Networks_order_manual.txt", 
                     delim = "\t", escape_double = FALSE, 
                     col_names = FALSE, trim_ws = TRUE) %>% dplyr::select(1:2) %>% rename(mask_value = X1, label = X2) %>%
  mutate(hemi = stringr::str_extract(label, "_([^_]+)_"), 
         hemi = stringr::str_extract(hemi, "[^_]"),
         network = substr(label, 14, 17),
         mask_value = as.factor(mask_value))
df <- df %>% inner_join(labels, by = "mask_value") %>% filter(network=="Dors") %>% inner_join(dan_labels, by = "mask_value")

# inspect
pe_dist <-  ggplot(df, aes(rewFunc, value, color = rewFunc)) + geom_violin(draw_quantiles = .5) + ylab("PE fMRI beta")
pe_means <- ggplot(df, aes(rewFunc, value, color = rewFunc)) + stat_summary(fun=mean, geom="point") + ylab("PE fMRI beta, mean")

# revisit MEG betas
ggplot(df, aes(rewFunc, omission_early_theta, color = rewFunc)) + geom_violin(draw_quantiles = .5) + ylab("MEG omission early theta")


# PE beta prediction
# early theta
m_pe_etheta <- lmer(value ~ rewFunc *  Stream * omission_early_theta  +  (1|id), df)
summary(m_pe_etheta)
Anova(m_pe_etheta, '3')

# em_pe_etheta <- as_tibble(emmeans(m_pe_etheta, data = df, ~Stream|omission_early_theta*rewFunc, at = list(omission_early_theta = c(-0.13878, 0.66119))))
# sample MEG at condition-wise quantiles
em_pe_etheta <- as_tibble(emmeans(m_pe_etheta, data = df, ~Stream|omission_early_theta*rewFunc, at = list(omission_early_theta = qmeg$omission_early_theta))) %>%
  inner_join(qmeg, by = c("omission_early_theta", "rewFunc"))

pe_pe_etheta <- ggplot(em_pe_etheta, aes(x = Stream, y=emmean, ymin=asymp.LCL, ymax=asymp.UCL, color=as.factor(omission_early_theta))) +
  geom_point(position = position_dodge(width = .6), size=2.5) +
  geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + facet_wrap(~rewFunc) +
  theme_bw(base_size=12) +  ylab("BOLD response to reward prediction error")  +
  scale_color_manual("Reward omission\nearly theta\nsynchronization", values = em_pe_etheta$theta_color, labels = em_pe_etheta$q) +
  labs(shape = "Reward omission\nearly theta\nsynchronization") +
  theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) # +

ggplot(em_pe_etheta, aes(x = Stream, y=emmean, ymin=asymp.LCL, ymax=asymp.UCL, color = rewFunc, lty=as.factor(q))) +
  geom_point(position = position_dodge(width = .6), size=2.5) +
  geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + facet_wrap(~(rewFunc=="DEV" | rewFunc=="IEV")) +
  theme_bw(base_size=12) +  ylab("BOLD response to reward prediction error")  +
  # scale_color_manual("Reward omission\nearly theta\nsynchronization", values = em_rewom_etheta$theta_color, labels = em_rewom_etheta$q) +
  labs(lty = "Reward omission\nearly theta\nsynchronization", color = "Condition") +
  theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) # +

# effects of early theta decomposed into condition and session levels
if (decompose) {
  m_pe_etheta_dec <- lmer(value ~ rewFunc *  Stream * omission_early_theta_avg  +
                            rewFunc *  Stream * om_theta_wi  +(1|id), df)
  summary(m_pe_etheta_dec)
  anova(m_pe_etheta, m_pe_etheta_dec)
  Anova(m_pe_etheta_dec, '3')
}



# entropy change late beta
m_pe_lbeta <- lmer(value ~ rewFunc * Stream * entropy_change_late_beta_supp + (1|id), df)
# summary(m_pe_lbeta)
Anova(m_pe_lbeta, '3')
em_pe_lbeta <- as_tibble(emmeans(m_pe_lbeta, data = df, ~Stream|entropy_change_late_beta_supp*rewFunc, 
                                 at = list(entropy_change_late_beta_supp = qmeg$entropy_change_late_beta_supp))) %>%
  inner_join(qmeg, by = c("entropy_change_late_beta_supp", "rewFunc"))


pe_ec_lbeta <- ggplot(em_pe_lbeta, aes(x = Stream, y=emmean, ymin=asymp.LCL, ymax=asymp.UCL, color=as.factor(entropy_change_late_beta_supp))) +
  geom_point(position = position_dodge(width = .6), size=2.5) +
  geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + facet_wrap(~rewFunc) +
  theme_bw(base_size=12) +  ylab("BOLD response toreward prediction error")  +
  scale_color_manual("Entropy change\nlate beta\nsuppression", values = em_pe_lbeta$qcolor, labels = em_pe_lbeta$q)  +
  labs(shape = "Entropy change\nlate beta\nsuppression") +
  theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) # +

ggplot(em_pe_lbeta, aes(x = Stream, y=emmean, ymin=asymp.LCL, ymax=asymp.UCL, color = rewFunc, lty=as.factor(q))) +
  geom_point(position = position_dodge(width = .6), size=2.5) +
  geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + facet_wrap(~(rewFunc=="DEV" | rewFunc=="IEV")) +
  theme_bw(base_size=12) +  ylab("BOLD response to reward > omission")  +
  # scale_color_manual("Reward omission\nearly theta\nsynchronization", values = em_rewom_etheta$theta_color, labels = em_rewom_etheta$q) +
  labs(lty = "Entropy change\nlate beta\nsuppression", color = "Condition") +
  theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) # +


setwd("~/OneDrive/collected_letters/papers/meg/plots/meg_to_fmri/")
pdf("MEG_on_pe_fMRI_beta_only.pdf", height = 12, width = 10)
ggarrange(pe_pe_etheta, pe_ec_lbeta, nrow = 2, ncol = 1)
dev.off()



# read in entropy change betas
# setwd(file.path(fmri_dir, 'L1m-entropy_echange'))
setwd(file.path(fmri_dir, 'L1m-echange'))
ec_betas <- read_csv("Schaefer_444_final_2009c_2.3mm_cope_l1.csv.gz") %>% filter(l1_cope_name=="EV_entropy_change_feedback") %>% 
  dplyr::select(id, run_number, l1_model, mask_value, value) %>% rename(roi_num7 = "mask_value")  %>% 
  mutate(id = as.character(id),
  roi_num7 = as.factor(roi_num7),
  run_mc  = scale(run_number, center = T, scale = F),
  beta_winsor = psych::winsor(value, trim = .05)) 
df <- ec_betas %>% inner_join(design, by = c("id", "run_number")) %>% inner_join(wbetas, by = c("id", "rewFunc")) %>% inner_join(labels, by = "roi_num7") 
  
# ec_betas <- read_csv("Schaefer2018_200Parcels_7Networks_order_fonov_2.3mm_ants_cope_l1.csv.gz") %>% filter(l1_cope_name=="EV_entropy_change_feedback") %>% 
#   dplyr::select(id, run_number, l1_model, mask_value, value) %>% mutate(id = as.character(id),
#                                                                         mask_value = as.factor(mask_value),
#                                                                         run_mc  = scale(run_number, center = T, scale = F),
#                                                                         beta_winsor = psych::winsor(value, trim = .05)) 

# df <- ec_betas %>% inner_join(design, by = c("id", "run_number")) %>% inner_join(wbetas, by = c("id", "rewFunc")) %>% inner_join(labels, by = "mask_value") %>% filter(network=="Dors") %>% 
#   inner_join(dan_labels, by = "mask_value") #%>%
#  filter(Stream!='visual-motion')

# summary(lm(value ~ rewFunc * Stream * run_mc, df))
summary(lm(value ~ rewFunc * vm_gradient17 * run_mc, df))
anova(lm(value ~ rewFunc * vm_gradient17 * run_mc, df))


# inspect
ec_dist <-  ggplot(df, aes(rewFunc, value, color = rewFunc)) + geom_violin(draw_quantiles = .5) + ylab("Entropy change change fMRI beta")
ec_means <- ggplot(df, aes(rewFunc, value, color = rewFunc)) + stat_summary(fun=mean, geom="point") + ylab("Entropy change fMRI beta,\nmean")# not much for late beta outside of DMN


# winsorizing does not help
# m6 <- lmer(beta_winsor ~ rewFunc * Stream * entropy_change_late_beta_supp + (1|id), df)

m_ec_etheta <- lmer(value ~ rewFunc *  vm_gradient17 * omission_early_theta  +  (1|id), df)
summary(m_ec_etheta)
Anova(m_ec_etheta, '3')


# em_ec_etheta <- as_tibble(emmeans(m_ec_etheta, data = df, ~vm_gradient17|omission_early_theta*rewFunc, at = list(omission_early_theta = c(-0.13878, 0.66119))))
# sample MEG at condition-wise quantiles
em_ec_etheta <- as_tibble(emmeans(m_ec_etheta, data = df, ~vm_gradient17|omission_early_theta*rewFunc, at = list(omission_early_theta = qmeg$omission_early_theta))) %>%
  inner_join(qmeg, by = c("omission_early_theta", "rewFunc"))

ec_etheta <- ggplot(em_ec_etheta, aes(x = vm_gradient17, y=emmean, ymin=asymp.LCL, ymax=asymp.UCL, color=as.factor(omission_early_theta))) +
  geom_point(position = position_dodge(width = .6), size=2.5) +
  geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + facet_wrap(~rewFunc) +
  theme_bw(base_size=12) +  ylab("BOLD response to entropy change")  +
  scale_color_manual("Reward omission\nearly theta\nsynchronization", values = em_ec_etheta$theta_color, labels = em_ec_etheta$q) +
  labs(shape = "Reward omission\nearly theta\nsynchronization") +
  theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) # +


# entropy change late beta
m_ec_lbeta <- lmer(value ~ rewFunc * vm_gradient17 * entropy_change_late_beta_supp + (1|id), df)
summary(m_ec_lbeta)
Anova(m_ec_lbeta, '3')

em_ec_lbeta <- as_tibble(emmeans(m_ec_lbeta, data = df, ~vm_gradient17|entropy_change_late_beta_supp*rewFunc, 
                                 at = list(entropy_change_late_beta_supp = qmeg$entropy_change_late_beta_supp))) %>%
  inner_join(qmeg, by = c("entropy_change_late_beta_supp", "rewFunc")) 

ec_lbeta <- ggplot(em_ec_lbeta, aes(x = vm_gradient17, y=emmean, ymin=asymp.LCL, ymax=asymp.UCL, color=as.factor(entropy_change_late_beta_supp))) +
  geom_point(position = position_dodge(width = .6), size=2.5) +
  geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + facet_wrap(~rewFunc) +
  theme_bw(base_size=12) +  ylab("BOLD response to entropy change")  +
  scale_color_manual("Entropy change\nlate beta\nsuppression", values = em_ec_lbeta$qcolor, labels = em_ec_lbeta$q)  +
  labs(shape = "Entropy change\nlate beta\nsuppression") +
  theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) # +

em_ec_lbeta_learn <- em_ec_lbeta %>% filter(rewFunc=="DEV" | rewFunc=="IEV")

ec_lbeta_learn <- ggplot(em_ec_lbeta_learn, 
                         aes(x = vm_gradient17, y=emmean, ymin=asymp.LCL, ymax=asymp.UCL, color=qcolor)) +
  geom_point(position = position_dodge(width = .6), size=2.5) +
  geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + facet_grid(~ rewFunc) +
  theme_bw(base_size=12) +  ylab("BOLD response to entropy change")  +
  scale_color_identity("Entropy change\nlate beta\nsuppression", guide = "legend", labels = unique(em_ec_lbeta_learn$q))  +
  labs(shape = "Entropy change\nlate beta\nsuppression") +
  theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) # +


em_ec_lbeta_unlearn <- em_ec_lbeta %>% filter(rewFunc=="CEV" | rewFunc=="CEVR")

ec_lbeta_unlearn <- ggplot(em_ec_lbeta_unlearn, 
                           aes(x = vm_gradient17, y=emmean, ymin=asymp.LCL, ymax=asymp.UCL, color=qcolor)) +
  geom_point(position = position_dodge(width = .6), size=2.5) +
  geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + facet_grid(~ rewFunc) +
  theme_bw(base_size=12) +  ylab("BOLD response to entropy change")  +
  scale_color_identity("Entropy change\nlate beta\nsuppression", guide = "legend", labels = unique(em_ec_lbeta_learn$q))  +
  labs(shape = "Entropy change\nlate beta\nsuppression") + 
  theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) # +

setwd("~/OneDrive/collected_letters/papers/meg/plots/meg_to_fmri/")
pdf("lbeta_learnable_unlearnable.pdf", height = 6, width = 9, onefile = F)
ggarrange(NULL, ec_lbeta_learn, NULL, ec_lbeta_unlearn, 
          ncol = 1, nrow = 4, heights = c(.1, 1, .1, 1), 
          common.legend = T, legend = "right", align = "h", 
          labels = c("", "Learnable conditions"," ", "Unlearnable conditions"), hjust = -1.6, vjust = -.2)
dev.off()
ggplot(em_ec_lbeta, aes(x = vm_gradient17, y=emmean, ymin=asymp.LCL, ymax=asymp.UCL, color = rewFunc, lty=as.factor(q))) +
  geom_point(position = position_dodge(width = .6), size=2.5) +
  geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + facet_wrap(~(rewFunc=="DEV" | rewFunc=="IEV")) +
  theme_bw(base_size=12) +  ylab("BOLD response to entropy change")  +
  # scale_color_manual("Reward omission\nearly theta\nsynchronization", values = em_rewom_etheta$theta_color, labels = em_rewom_etheta$q) +
  labs(lty = "Entropy change\nlate beta\nsuppression", color = "Condition") +
  theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) # +


setwd("~/OneDrive/collected_letters/papers/meg/plots/meg_to_fmri/")
pdf("MEG_on_ec_fMRI_beta_only.pdf", height = 12, width = 10)
ggarrange(ec_etheta, ec_lbeta, nrow = 2, ncol = 1)
dev.off()

pdf("MEG_lbeta_on_ec_fMRI_beta_learnable.pdf", height = 6, width = 10)
ggplot(em_ec_lbeta, aes(x = vm_gradient17, y=emmean, ymin=asymp.LCL, ymax=asymp.UCL, color = rewFunc, lty=as.factor(q))) +
  geom_point(position = position_dodge(width = .6), size=2.5) +
  geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + facet_wrap(~(rewFunc=="DEV" | rewFunc=="IEV")) +
  theme_bw(base_size=12) +  ylab("BOLD response to entropy change")  +
  # scale_color_manual("Reward omission\nearly theta\nsynchronization", values = em_rewom_etheta$theta_color, labels = em_rewom_etheta$q) +
  labs(lty = "Entropy change\nlate beta\nsuppression", color = "Condition") +
  theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) # +
dev.off()
setwd("~/OneDrive/collected_letters/papers/meg/plots/meg_to_fmri/")
pdf("MEG_on_ec_pe_fMRI_betas.pdf", height = 10, width = 18)
ggarrange(pe_pe_etheta, ec_etheta, pe_ec_lbeta,ec_lbeta, nrow = 2, ncol = 2)
dev.off()

ggarrange(pe_pe_etheta, pe_ec_lbeta, nrow = 2, ncol = 1)

# parcel-level late beta effects

# entropy change late beta
df$roi_num17 <- factor(df$roi_num17)
m_ec_lbeta_parcel <- lmer(value ~ rewFunc * roi_num17 * entropy_change_late_beta_supp + (1|id), df)
summary(m_ec_lbeta_parcel)
Anova(m_ec_lbeta_parcel, '3')


# interaction with run?
# m_ec_lbeta_run2 <- lmer(value ~ (rewFunc +  vm_gradient17 + entropy_change_late_beta_supp + run_mc) ^2 + (1|id), df)
# summary(m_ec_lbeta_run)
# Anova(m_ec_lbeta_run, '3')

# read in entropy EV_entropy_wiz_clock
e_betas <- read_csv("Schaefer2018_200Parcels_7Networks_order_fonov_2.3mm_ants_cope_l1.csv.gz") %>% filter(l1_cope_name=="EV_entropy_wiz_clock") %>% 
  dplyr::select(id, run_number, l1_model, mask_value, value) %>% mutate(id = as.character(id),
                                                                        mask_value = as.factor(mask_value),
                                                                        run_mc  = scale(run_number, center = T, scale = F))
# merge
df <- e_betas %>% inner_join(design, by = c("id", "run_number"))  %>% inner_join(labels, by = "mask_value") %>% 
  filter(network=="Dors") %>% inner_join(dan_labels, by = "mask_value")  %>% inner_join(wbetas, by = c("id", "rewFunc"))

# inspect
e_dist <- ggplot(df, aes(rewFunc, value, color = rewFunc)) + geom_violin(draw_quantiles = .5) + ylab("Entropy fMRI beta")
e_means <- ggplot(df, aes(rewFunc, value, color = rewFunc)) + stat_summary(fun=mean, geom="point") + ylab("Entropy fMRI beta, mean")

# late beta: again, mostly DMN
m7 <- lmer(value ~ rewFunc * parcel * entropy_change_late_beta_supp + (1|id), df)
summary(m7)
Anova(m7, '3')

sm7 <- lmer(value ~ rewFunc * vm_gradient17 * entropy_change_late_beta_supp + (1|id), df)
summary(sm7)
Anova(sm7, '3')

em7 <- as_tibble(emmeans(sm7, data = df, ~vm_gradient17|entropy_change_late_beta_supp*rewFunc, at = list(entropy_change_late_beta_supp = c(-0.07553, 0.22084))))
e_ec_lbeta <- ggplot(em7, aes(x = vm_gradient17, y=emmean, ymin=asymp.LCL, ymax=asymp.UCL, color=as.factor(entropy_change_late_beta_supp))) +
  geom_point(position = position_dodge(width = .6), size=2.5) +
  geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + facet_wrap(~rewFunc) +
  theme_bw(base_size=12) +  ylab("BOLD response to entropy, clock-aligned")  +
  scale_color_manual("Entropy change\nlate beta\nsuppression", values=c("#1b3840","#4fa3b8"), labels = c("10th %ile", "90th %ile")) +
  labs(shape = "Entropy change\nlate beta\nsuppression") +
  theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) # +

setwd("~/OneDrive/collected_letters/papers/meg/plots/meg_to_fmri/")
pdf("fMRI_beta_means_by_rewFunc.pdf", height = 6, width = 6)
ggarrange(rew_means, pe_means, ec_means, e_means)
dev.off()
pdf("fMRI_beta_dist_by_rewFunc.pdf", height = 6, width = 6)
ggarrange(rew_dist, pe_dist, ec_dist, e_dist)
dev.off()


pdf("MEG_to_fMRI_DAN_betas_rewFunc.pdf", height = 12, width = 16)
ggarrange(rewom_ec_lbeta, pe_ec_lbeta, ec_lbeta, e_ec_lbeta, ncol = 2, nrow = 2)
dev.off()
