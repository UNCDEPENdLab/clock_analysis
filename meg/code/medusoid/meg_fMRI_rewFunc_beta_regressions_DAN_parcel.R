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

# get design
design <- get_trial_data(repo_directory = clock_folder, dataset = "mmclock_fmri", groupfixed = T) %>% select(id, run, rewFunc) %>% unique() %>%
  rename(run_number = run) %>% mutate(id = as.character(id))

design_face <- get_trial_data(repo_directory = clock_folder, dataset = "mmclock_fmri", groupfixed = T) %>% select(id, run, rewFunc, emotion) %>% unique() %>%
  rename(run_number = run) %>% mutate(id = as.character(id))


# load meg data
# wbetas <- readRDS("~/OneDrive/collected_letters/papers/meg/plots/wholebrain/betas/MEG_betas_wide_echange_vmax_reward_Nov30_2021.RDS") %>% 
wbetas <- readRDS("~/code/clock_analysis/meg/data/MEG_betas_ec_rewfunc_April_5_2022.RDS") %>% 
  mutate(entropy_change_early_beta_supp = -  entropy_change_early_beta_ec_rewfunc,
         entropy_change_late_beta_supp = - entropy_change_late_beta_ec_rewfunc) %>%
  select(id, rewFunc, entropy_change_early_beta_supp, entropy_change_late_beta_supp)

# load fMRI betas: reward_omission, entropy, signed_pe
# reward/omission
setwd(file.path(fmri_dir, 'L1m-rew_om'))
rew_betas <- read_csv("Schaefer2018_200Parcels_7Networks_order_fonov_2.3mm_ants_cope_l1.csv.gz") %>% filter(l1_cope_name=="EV_rew_om") %>% 
  dplyr::select(id, run_number, l1_model, mask_value, value) %>% mutate(id = as.character(id),
                                                                        mask_value = as.character(mask_value),
                                                                        run_mc  = scale(run_number, center = T, scale = F))
# merge
df <- rew_betas %>% inner_join(wbetas, by = "id", "rewFunc") %>% inner_join(design, by = c("id", "run_number", "rewFunc")) 

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

summary(lm(value ~ rewFunc * Stream, df))
# prototype


# streams
sm1 <- lmer(value ~ rewFunc *  Stream * entropy_change_late_beta_supp + (1|id), df)
summary(sm1)
Anova(sm1, '3')

esm1 <- as_tibble(emmeans(sm1, data = df, ~Stream|entropy_change_late_beta_supp*rewFunc, at = list(entropy_change_late_beta_supp = c(-0.07553, 0.22084), run_number = c(1,8))))
rewom_ec_lbeta <- ggplot(esm1, aes(x = Stream, y=emmean, ymin=asymp.LCL, ymax=asymp.UCL, color=as.factor(entropy_change_late_beta_supp))) +
  geom_point(position = position_dodge(width = .6), size=2.5) +
  geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + facet_wrap(~rewFunc) +
  theme_bw(base_size=12) +  ylab("BOLD response to reward > omission")  +
  scale_color_manual("Entropy change\nlate beta\nsuppression", values=c("#1b3840","#4fa3b8"), labels = c("10th %ile", "90th %ile")) +
  labs(shape = "Entropy change\nlate beta\nsuppression") +
  theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) # +



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

# PE beta prediction
# entropy change late beta
m5 <- lmer(value ~ rewFunc * Stream * entropy_change_late_beta_supp + (1|id), df)
summary(m5)
Anova(m5, '3')
em5 <- as_tibble(emmeans(m5, data = df, ~Stream|entropy_change_late_beta_supp*rewFunc, at = list(entropy_change_late_beta_supp = c(-0.07553, 0.22084), run_number = c(1,8))))
pe_ec_lbeta <- ggplot(em5, aes(x = Stream, y=emmean, ymin=asymp.LCL, ymax=asymp.UCL, color=as.factor(entropy_change_late_beta_supp))) +
  geom_point(position = position_dodge(width = .6), size=2.5) +
  geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + facet_wrap(~rewFunc) +
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
                                                                        run_mc  = scale(run_number, center = T, scale = F),
                                                                        beta_winsor = psych::winsor(value, trim = .05))
df <- ec_betas %>% inner_join(design, by = c("id", "run_number")) %>% inner_join(wbetas, by = c("id", "rewFunc")) %>% inner_join(labels, by = "mask_value") %>% filter(network=="Dors") %>% 
  inner_join(dan_labels, by = "mask_value") #%>%
#  filter(Stream!='visual-motion')

# inspect
ec_dist <-  ggplot(df, aes(rewFunc, value, color = rewFunc)) + geom_violin(draw_quantiles = .5) + ylab("Entropy change change fMRI beta")
ec_means <- ggplot(df, aes(rewFunc, value, color = rewFunc)) + stat_summary(fun=mean, geom="point") + ylab("Entropy change fMRI beta,\nmean")# not much for late beta outside of DMN


# winsorizing does not help
# m6 <- lmer(beta_winsor ~ rewFunc * Stream * entropy_change_late_beta_supp + (1|id), df)
# summary(m6)
# Anova(m6, '3')
# 


# these effects are not selective for late beta, also seen for early beta
m6a <- lmer(value ~ rewFunc * Stream * entropy_change_late_beta_supp + (1|id), df)
summary(m6a)
Anova(m6a, '3')

# em6 <- as_tibble(emmeans(m6, data = df, ~Stream|entropy_change_late_beta_supp*rewFunc, at = list(entropy_change_late_beta_supp = c(-0.07553, 0.22084), run_number = c(1,8))))
em6 <- as_tibble(emmeans(m6a, data = df, ~Stream|entropy_change_late_beta_supp*rewFunc, at = list(entropy_change_late_beta_supp = c(-0.07553, 0.22084))))
ec_ec_lbeta <- ggplot(em6, aes(x=Stream, y=emmean, ymin=asymp.LCL, ymax=asymp.UCL, color=as.factor(entropy_change_late_beta_supp))) +
  geom_point(position = position_dodge(width = .6), size=2.5) +
  geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + facet_wrap(~rewFunc) +
  theme_bw(base_size=12) +  ylab("BOLD response to entropy change")  +
  scale_color_manual("Entropy change\nlate beta\nsuppression", values=c("#1b3840","#4fa3b8"), labels = c("10th %ile", "90th %ile")) +
  labs(shape = "Entropy change\nlate beta\nsuppression") +
  theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) # +

# visuo-motor gradient
# vm6 <- lmer(value ~ rewFunc * lobe * entropy_change_late_beta_supp + (1|id), df)
# summary(vm6)
# Anova(vm6, '3')


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

sm7 <- lmer(value ~ rewFunc * Stream * entropy_change_late_beta_supp + (1|id), df)
summary(sm7)
Anova(sm7, '3')

em7 <- as_tibble(emmeans(sm7, data = df, ~Stream|entropy_change_late_beta_supp*rewFunc, at = list(entropy_change_late_beta_supp = c(-0.07553, 0.22084), run_number = c(1,8))))
e_ec_lbeta <- ggplot(em7, aes(x = Stream, y=emmean, ymin=asymp.LCL, ymax=asymp.UCL, color=as.factor(entropy_change_late_beta_supp))) +
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
ggarrange(rewom_ec_lbeta, pe_ec_lbeta, ec_ec_lbeta, e_ec_lbeta, ncol = 2, nrow = 2)
dev.off()
