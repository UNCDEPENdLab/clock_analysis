# relating MEG to run-wise DAN fMRI betas at parcel level

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

# load meg data
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
  dplyr::select(c(id, omission_early_theta, omission_late_delta, 
                  entropy_change_early_beta_supp, entropy_change_late_beta_supp,
                  # entropy_change_early_beta_supp_ec, entropy_change_late_beta_supp_ec, 
                  vmax_late_alpha, vmax_late_beta,
                  neg_pe_early_theta, pe_late_beta_supp,
                  abspe_early_theta))
# abs_pe_late_beta_supp_ec, abs_pe_late_beta_supp))

# load fMRI betas: reward_omission, entropy, signed_pe
# reward/omission
setwd(file.path(fmri_dir, 'L1m-rew_om'))
rew_betas <- read_csv("Schaefer2018_200Parcels_7Networks_order_fonov_2.3mm_ants_cope_l1.csv.gz") %>% filter(l1_cope_name=="EV_rew_om") %>% 
  dplyr::select(id, run_number, l1_model, mask_value, value) %>% mutate(id = as.character(id),
                                                                        mask_value = as.character(mask_value),
                                                                        run_mc  = scale(run_number, center = T, scale = F))
# merge
df <- rew_betas %>% inner_join(wbetas, by = "id")

labels <- read_delim("~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/Schaefer2018_200Parcels_7Networks_order_manual.txt", 
                     delim = "\t", escape_double = FALSE, 
                     col_names = FALSE, trim_ws = TRUE) %>% dplyr::select(1:2) %>% rename(mask_value = X1, label = X2) %>%
  mutate(hemi = stringr::str_extract(label, "_([^_]+)_"), 
         hemi = stringr::str_extract(hemi, "[^_]"),
         network = substr(label, 14, 17),
         mask_value = as.factor(mask_value))
dan_labels <- read_excel("~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/MNH Dan Labels.xlsx") %>% 
  mutate(mask_value = as.character(roinum),
         parcel = word(MNHLabel, 2, sep = "_")) %>% select(mask_value, MNHLabel, parcel, Stream, Visuomotor_Gradient, Stream_Gradient, lobe)
df <- df %>% inner_join(labels, by = "mask_value") %>% filter(network=="Dors") %>% inner_join(dan_labels, by = "mask_value")


# prototype
# NB: convergence problems!!!
m1 <- lmer(value ~ run_number +  parcel * omission_early_theta + (1|id), df)
while (any(grepl("failed to converge", m1@optinfo$conv$lme4$messages) )) {
ss <- getME(m1,c("theta","fixef"))
m1 <- update(m1,start=ss,control=lmerControl(optCtrl=list(maxfun=2e4)))}
summary(m1)
Anova(m1, '3')

# em1 <- as_tibble(emmeans(m1, data = df, ~network|omission_early_theta*run_number, at = list(omission_early_theta = c(-.34, .31), run_number = c(1,8))))
# 
# rewom_rewom_theta <- ggplot(em1, aes(x=network, y=emmean, ymin=asymp.LCL, ymax=asymp.UCL, color=as.factor(omission_early_theta))) +
#   geom_point(position = position_dodge(width = .6), size=2.5) +
#   geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + facet_wrap(~run_number) +
#   theme_bw(base_size=12) +  ylab("BOLD response to reward > omission")  +
#   scale_color_manual("Reward omission\nearly theta\nresponse", values=c("#1b3840","#4fa3b8"), labels = c("10th %ile", "90th %ile")) +
#   labs(shape = "Reward omission\nearly theta\nresponse") +
#   theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
#         axis.text=element_text(size=8.5, color="grey10")) # +

# Prediction errors
setwd(file.path(fmri_dir, 'L1m-pe'))
pe_betas <- read_csv("Schaefer2018_200Parcels_7Networks_order_fonov_2.3mm_ants_cope_l1.csv.gz") %>% filter(l1_cope_name=="EV_pe") %>% 
  dplyr::select(id, run_number, l1_model, mask_value, value) %>% mutate(id = as.character(id),
                                                                        mask_value = as.factor(mask_value),
                                                                        run_mc  = scale(run_number, center = T, scale = F))
# merge
df <- pe_betas %>% inner_join(wbetas, by = "id")

labels <- read_delim("~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/Schaefer2018_200Parcels_7Networks_order_manual.txt", 
                     delim = "\t", escape_double = FALSE, 
                     col_names = FALSE, trim_ws = TRUE) %>% dplyr::select(1:2) %>% rename(mask_value = X1, label = X2) %>%
  mutate(hemi = stringr::str_extract(label, "_([^_]+)_"), 
         hemi = stringr::str_extract(hemi, "[^_]"),
         network = substr(label, 14, 17),
         mask_value = as.factor(mask_value))
df <- df %>% inner_join(labels, by = "mask_value") %>% filter(network=="Dors") %>% inner_join(dan_labels, by = "mask_value")
# PE beta prediction
m1a <- lmer(value ~ run_number * parcel * omission_early_theta + (1|id), df)
while (any(grepl("failed to converge", m1a@optinfo$conv$lme4$messages) )) {
  ss <- getME(m1a,c("theta","fixef"))
  m1a <- update(m1a,start=ss,control=lmerControl(optCtrl=list(maxfun=2e4)))}
summary(m1a)
Anova(m1a, '3')
em1a <-  as_tibble(emmeans(m1a, data = df, ~network|omission_early_theta*run_number, at = list(omission_early_theta = c(-.34, .31), run_number = c(1,8))))
pe_rewom_theta <- ggplot(em1a, aes(x=network, y=emmean, ymin=asymp.LCL, ymax=asymp.UCL, color=as.factor(omission_early_theta))) +
  geom_point(position = position_dodge(width = .6), size=2.5) +
  geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + facet_wrap(~run_number) +
  theme_bw(base_size=12) +  ylab("BOLD response to reward prediction error")  +
  scale_color_manual("Reward omission\nearly theta\nresponse", values=c("#1b3840","#4fa3b8"), labels = c("10th %ile", "90th %ile")) +
  labs(shape = "Reward omission\nearly theta\nresponse") +
  theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) # +

# nothing for PE late beta:
m2 <- lmer(value ~ run_number * parcel * pe_late_beta_supp + (1|id), df)
while (any(grepl("failed to converge", m2@optinfo$conv$lme4$messages) )) {
  ss <- getME(m2,c("theta","fixef"))
  m2 <- update(m2,start=ss,control=lmerControl(optCtrl=list(maxfun=2e4)))}

summary(m2)
Anova(m2, '3')
em2 <- as_tibble(emmeans(m2, data = df, ~network|pe_late_beta_supp*run_number, at = list(pe_late_beta_supp = c(-0.005, 0.0073), run_number = c(1, 8))))

pe_pe_beta <- ggplot(em2, aes(x=network, y=emmean, ymin=asymp.LCL, ymax=asymp.UCL, color=as.factor(pe_late_beta_supp))) +
  geom_point(position = position_dodge(width = .6), size=2.5) +
  geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + facet_wrap(~run_number) +
  theme_bw(base_size=12) +  ylab("BOLD response to reward prediction error")  +
  scale_color_manual("Signed PE\nearly theta\nresponse", values=c("#1b3840","#4fa3b8"), labels = c("10th %ile", "90th %ile")) +
  labs(shape = "Signed PE\nearly theta\nresponse") +
  theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10"))



# abs(PE) modulates control and visual
m3 <- lmer(value ~ run_number * parcel * abspe_early_theta + (1|id), df)
while (any(grepl("failed to converge", m3@optinfo$conv$lme4$messages) )) {
  ss <- getME(m3,c("theta","fixef"))
  m3 <- update(m3,start=ss,control=lmerControl(optCtrl=list(maxfun=2e4)))}

summary(m3)
Anova(m3, '3')
# n    missing   distinct       Info       Mean        Gmd        .05        .10        .25        .50        .75 
# 92800          0         58          1  0.0001946   0.004626 -0.0046630 -0.0042458 -0.0020824 -0.0005048  0.0022254 
# .90        .95 
# 0.0061380  0.0066809 
em3 <- as_tibble(emmeans(m3, data = df, ~network|abspe_early_theta*run_number, at = list(abspe_early_theta = c(-0.0042, 0.0061), run_number = c(1,8))))

pe_abspe_theta <- ggplot(em3, aes(x=network, y=emmean, ymin=asymp.LCL, ymax=asymp.UCL, color=as.factor(abspe_early_theta))) +
  geom_point(position = position_dodge(width = .6), size=2.5) +
  geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + facet_wrap(~run_number) +
  theme_bw(base_size=12) +  ylab("BOLD response to reward prediction error")  +
  scale_color_manual("Absolute PE\nearly theta\nresponse", values=c("#1b3840","#4fa3b8"), labels = c("10th %ile", "90th %ile")) +
  labs(shape = "Absolute PE\nearly theta\nresponse") +
  theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10"))

# effects of omission and abs(PE) early theta are independent
m4 <- lmer(value ~ run_number * parcel * abspe_early_theta + run_number * parcel * omission_early_theta + (1|id), df)
summary(m4)
Anova(m4, '3')
vif(m4)
em4 <- as_tibble(emmeans(m4, data = df, ~network|omission_early_theta, at = list(omission_early_theta = c(-.34, .31))))
pe_rewom_theta <- ggplot(em4, aes(x=network, y=emmean, ymin=asymp.LCL, ymax=asymp.UCL, color=as.factor(omission_early_theta))) +
  geom_point(position = position_dodge(width = .6), size=2.5) +
  geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) +
  theme_bw(base_size=12) +  ylab("BOLD response to reward prediction error")  +
  scale_color_manual("Reward omission\nearly theta\nresponse", values=c("#1b3840","#4fa3b8"), labels = c("10th %ile", "90th %ile")) +
  labs(shape = "Reward omission\nearly theta\nresponse") +
  theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) # +

# entropy change late beta: weak effects on PEs in the DMN, SM, V
m5 <- lmer(value ~ run_number * parcel * entropy_change_late_beta_supp + (1|id), df)
summary(m5)
Anova(m5, '3')
em5 <- as_tibble(emmeans(m5, data = df, ~network|entropy_change_late_beta_supp*run_number, at = list(entropy_change_late_beta_supp = c(-0.142, 0.154), run_number = c(1,8))))
pe_ec_lbeta <- ggplot(em5, aes(x=network, y=emmean, ymin=asymp.LCL, ymax=asymp.UCL, color=as.factor(entropy_change_late_beta_supp))) +
  geom_point(position = position_dodge(width = .6), size=2.5) +
  geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + facet_wrap(~run_number) +
  theme_bw(base_size=12) +  ylab("BOLD response to reward prediction error")  +
  scale_color_manual("Entropy change\nlate beta\nsuppression", values=c("#1b3840","#4fa3b8"), labels = c("10th %ile", "90th %ile")) +
  labs(shape = "Entropy change\nlate beta\nsuppression") +
  theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) # +






# read in entropy change betas
setwd(file.path(fmri_dir, 'L1m-entropy_echange'))
ec_betas <- read_csv("Schaefer2018_200Parcels_7Networks_order_fonov_2.3mm_ants_cope_l1.csv.gz") %>% filter(l1_cope_name=="EV_entropy_change_feedback") %>% 
  dplyr::select(id, run_number, l1_model, mask_value, value) %>% mutate(id = as.character(id),
                                                                        mask_value = as.factor(mask_value),
                                                                        run_mc  = scale(run_number, center = T, scale = F))
df <- ec_betas %>% inner_join(wbetas, by = "id") %>% inner_join(labels, by = "mask_value") %>% filter(network=="Dors") %>% inner_join(dan_labels, by = "mask_value")

# not much for late beta outside of DMN
m6 <- lmer(value ~ run_number * parcel * entropy_change_late_beta_supp + (1|id), df)
summary(m6)
Anova(m6, '3')

# more for early beta in DMN and somatomotor, visual
m6a <- lmer(value ~ run_number * parcel * entropy_change_early_beta_supp + (1|id), df)
summary(m6a)
Anova(m6a, '3')
# df$entropy_change_late_beta_supp 
# n  missing distinct     Info     Mean      Gmd      .05      .10      .25      .50      .75      .90      .95 
# 92800        0       58        1 -0.00563   0.1265 -0.18217 -0.14210 -0.07790 -0.01109  0.08370  0.15427  0.21376 

em6 <- as_tibble(emmeans(m6, data = df, ~network|entropy_change_late_beta_supp*run_number, at = list(entropy_change_late_beta_supp = c(-0.142, 0.154), run_number = c(1,8))))
ec_ec_lbeta <- ggplot(em6, aes(x=network, y=emmean, ymin=asymp.LCL, ymax=asymp.UCL, color=as.factor(entropy_change_late_beta_supp))) +
  geom_point(position = position_dodge(width = .6), size=2.5) +
  geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + facet_wrap(~run_number) +
  theme_bw(base_size=12) +  ylab("BOLD response to entropy change")  +
  scale_color_manual("Entropy change\nlate beta\nsuppression", values=c("#1b3840","#4fa3b8"), labels = c("10th %ile", "90th %ile")) +
  labs(shape = "Entropy change\nlate beta\nsuppression") +
  theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) # +

# read in entropy EV_entropy_wiz_clock
e_betas <- read_csv("Schaefer2018_200Parcels_7Networks_order_fonov_2.3mm_ants_cope_l1.csv.gz") %>% filter(l1_cope_name=="EV_entropy_wiz_clock") %>% 
  dplyr::select(id, run_number, l1_model, mask_value, value) %>% mutate(id = as.character(id),
                                                                        mask_value = as.factor(mask_value),
                                                                        run_mc  = scale(run_number, center = T, scale = F))
# merge
df <- e_betas %>% inner_join(wbetas, by = "id")  %>% inner_join(labels, by = "mask_value") %>% filter(network=="Dors") %>% inner_join(dan_labels, by = "mask_value")

# late beta: again, mostly DMN
m7 <- lmer(value ~ run_number * parcel * entropy_change_late_beta_supp + (1|id), df)
summary(m7)
Anova(m7, '3')
em7 <- as_tibble(emmeans(m7, data = df, ~network|entropy_change_late_beta_supp*run_number, at = list(entropy_change_late_beta_supp = c(-0.142, 0.154), run_number = c(1,8))))
e_ec_lbeta <- ggplot(em6, aes(x=network, y=emmean, ymin=asymp.LCL, ymax=asymp.UCL, color=as.factor(entropy_change_late_beta_supp))) +
  geom_point(position = position_dodge(width = .6), size=2.5) +
  geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + facet_wrap(~run_number) +
  theme_bw(base_size=12) +  ylab("BOLD response to entropy, clock-aligned")  +
  scale_color_manual("Entropy change\nlate beta\nsuppression", values=c("#1b3840","#4fa3b8"), labels = c("10th %ile", "90th %ile")) +
  labs(shape = "Entropy change\nlate beta\nsuppression") +
  theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) # +
setwd("~/OneDrive/collected_letters/papers/meg/plots/wholebrain/")
pdf("MEG_to_fMRI_runwise_betas.pdf", height = 10, width = 22)
ggarrange(rewom_rewom_theta, pe_rewom_theta, pe_pe_beta, ec_ec_lbeta, e_ec_lbeta)
dev.off()
