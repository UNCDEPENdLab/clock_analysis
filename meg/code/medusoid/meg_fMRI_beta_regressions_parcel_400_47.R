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

# load meg data
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

# load fMRI betas: reward_omission, entropy, signed_pe
# reward/omission
setwd(file.path(fmri_dir, 'L1m-rew_om'))
rew_betas <- read_csv("Schaefer2018_200Parcels_7Networks_order_fonov_2.3mm_ants_cope_l1.csv.gz") %>% filter(l1_cope_name=="EV_rew_om") %>% 
  dplyr::select(id, run_number, l1_model, mask_value, value) %>% mutate(id = as.character(id),
                                                                        mask_value = as.character(mask_value),
                                                                        run_mc  = scale(run_number, center = T, scale = F))
# merge
df <- rew_betas %>% inner_join(design, by = c("id", "run_number")) %>% inner_join(cond_wbetas, by = c("id", "rewFunc"))

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

m1 <- lmer(value ~ rewFunc *  parcel * omission_early_theta + (1|id), df)
summary(m1)
Anova(m1, '3')
# em1 <- as_tibble(emmeans(m1, data = df, ~parcel|omission_early_theta*rewFunc, at = list(omission_early_theta = c(-.04, .61), run_number = c(1,8))))
em1 <- as_tibble(emmeans(m1, data = df, ~parcel|omission_early_theta*rewFunc, at = list(omission_early_theta = c(-.04, .61))))

rewom_rewom_theta <- ggplot(em1, aes(x=parcel, y=emmean, ymin=asymp.LCL, ymax=asymp.UCL, color=as.factor(omission_early_theta))) +
  geom_point(position = position_dodge(width = .6), size=2.5) +
  geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + facet_wrap(~rewFunc) +
  theme_bw(base_size=12) +  ylab("BOLD response to reward > omission")  +
  scale_color_manual("Reward omission\nearly theta\nresponse", values=c("#1b3840","#4fa3b8"), labels = c("10th %ile", "90th %ile")) +
  labs(shape = "Reward omission\nearly theta\nresponse") +
  theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) # +


# streams
sm1 <- lmer(value ~ rewFunc *  Stream * omission_early_theta + (1|id), df)
summary(sm1)
Anova(sm1, '3')
esm1 <- as_tibble(emmeans(sm1, data = df, ~Stream|omission_early_theta*rewFunc, at = list(omission_early_theta = c(-.04, .61))))
rewom_rewom_theta_streams <- ggplot(esm1, aes(x=Stream, y=emmean, ymin=asymp.LCL, ymax=asymp.UCL, color=as.factor(omission_early_theta))) +
  geom_point(position = position_dodge(width = .6), size=2.5) +
  geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + facet_wrap(~rewFunc) +
  theme_bw(base_size=12) +  ylab("BOLD response to reward > omission")  +
  scale_color_manual("Reward omission\nearly theta\nresponse", values=c("#1b3840","#4fa3b8"), labels = c("10th %ile", "90th %ile")) +
  labs(shape = "Reward omission\nearly theta\nresponse") +
  theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) # +

# decomposed version
dsm1 <- lmer(value ~ rewFunc *  Stream * omission_early_theta_avg +rewFunc *  Stream * om_theta_wi + (1|id), df)
summary(dsm1)
Anova(dsm1, '3')

dom1 <- as_tibble(emtrends(dsm1, data = df,  var = "omission_early_theta_avg", specs = c("Stream", "rewFunc")))
om_etheta_between_decomposed <- 
  ggplot(dom1, aes(x=Stream, y=omission_early_theta_avg.trend, ymin=asymp.LCL, ymax=asymp.UCL)) + 
  geom_point(position = position_dodge(width = .6), size=2.5) + 
  geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + 
  theme_bw(base_size=12) + facet_grid(~rewFunc) + ylab("Effect of early theta synchronization,\n between subjects") +
  theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) 

wdom1 <- as_tibble(emtrends(dsm1, data = df,  var = "om_theta_wi", specs = c("Stream", "rewFunc")))
om_etheta_within_decomposed <- 
  ggplot(wdom1, aes(x=Stream, y=om_theta_wi.trend, ymin=asymp.LCL, ymax=asymp.UCL)) + 
  geom_point(position = position_dodge(width = .6), size=2.5) + 
  geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + 
  theme_bw(base_size=12) + facet_grid(~rewFunc) + ylab("Effect of early theta synchronization,\n within-subject") +
  theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) 

# Prediction errors
setwd(file.path(fmri_dir, 'L1m-pe'))
pe_betas <- read_csv("Schaefer2018_200Parcels_7Networks_order_fonov_2.3mm_ants_cope_l1.csv.gz") %>% filter(l1_cope_name=="EV_pe") %>% 
  dplyr::select(id, run_number, l1_model, mask_value, value) %>% mutate(id = as.character(id),
                                                                        mask_value = as.factor(mask_value),
                                                                        run_mc  = scale(run_number, center = T, scale = F),
                                                                        beta_winsor = psych::winsor(value, trim = .05))
# merge
df <- pe_betas %>% inner_join(design, by = c("id", "run_number")) %>% inner_join(cond_wbetas, by = c("id", "rewFunc"))
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
m1a <- lmer(value ~ rewFunc * Stream * omission_early_theta + (rewFunc|id), df)
while (any(grepl("failed to converge", m1a@optinfo$conv$lme4$messages) )) {
  ss <- getME(m1a,c("theta","fixef"))
  m1a <- update(m1a,start=ss,control=lmerControl(optCtrl=list(maxfun=2e4)))}
summary(m1a)
Anova(m1a, '3')
em1a <-  as_tibble(emmeans(m1a, data = df, ~rewFunc|Stream*omission_early_theta, at = list(omission_early_theta = c(-.04, .61))))
# pe_rewom_theta <- ggplot(em1a, aes(x = rewFunc, y=emmean, ymin=asymp.LCL, ymax=asymp.UCL, color=as.factor(omission_early_theta))) +
pe_rewom_theta <- ggplot(em1a, aes(x = Stream, y=emmean, ymin=asymp.LCL, ymax=asymp.UCL, color=as.factor(omission_early_theta))) +
  geom_point(position = position_dodge(width = .6), size=2.5) + facet_wrap(~rewFunc) +
  geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + 
  theme_bw(base_size=12) +  ylab("BOLD response to reward prediction error")  +
  scale_color_manual("Reward omission\nearly theta\nresponse", values=c("#1b3840","#4fa3b8"), labels = c("10th %ile", "90th %ile")) +
  labs(shape = "Reward omission\nearly theta\nresponse") +
  theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) # +



# nothing for PE late beta:
m2 <- lmer(value ~ rewFunc * Stream * pe_late_beta_supp + (1|id), df)
while (any(grepl("failed to converge", m2@optinfo$conv$lme4$messages) )) {
  ss <- getME(m2,c("theta","fixef"))
  m2 <- update(m2,start=ss,control=lmerControl(optCtrl=list(maxfun=2e4)))}

summary(m2)
Anova(m2, '3')
em2 <- as_tibble(emmeans(m2, data = df, ~Stream|pe_late_beta_supp*rewFunc, at = list(pe_late_beta_supp = c(-1.750e-03, 1.057e-02))))

pe_pe_beta <- ggplot(em2, aes(x = Stream, y=emmean, ymin=asymp.LCL, ymax=asymp.UCL, color=as.factor(pe_late_beta_supp))) +
  geom_point(position = position_dodge(width = .6), size=2.5) +
  geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + facet_wrap(~rewFunc) +
  theme_bw(base_size=12) +  ylab("BOLD response to reward prediction error")  +
  scale_color_manual("Signed PE\nlate beta\nresponse", values=c("#1b3840","#4fa3b8"), labels = c("10th %ile", "90th %ile")) +
  labs(shape = "Signed PE\nlate beta\nresponse") +
  theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10"))



# abs(PE) early theta: nothing
m3 <- lmer(value ~ rewFunc * Stream * abspe_early_theta + (1|id), df)
while (any(grepl("failed to converge", m3@optinfo$conv$lme4$messages) )) {
  ss <- getME(m3,c("theta","fixef"))
  m3 <- update(m3,start=ss,control=lmerControl(optCtrl=list(maxfun=2e4)))}

summary(m3)
Anova(m3, '3')
# n    missing   distinct       Info       Mean        Gmd        .05        .10        .25        .50        .75 
# 92800          0         58          1  0.0001946   0.004626 -0.0046630 -0.0042458 -0.0020824 -0.0005048  0.0022254 
# .90        .95 
# 0.0061380  0.0066809 
em3 <- as_tibble(emmeans(m3, data = df, ~Stream|abspe_early_theta*rewFunc, at = list(abspe_early_theta = c(-0.0022299, .0078820061), run_number = c(1,8))))

pe_abspe_theta <- ggplot(em3, aes(x = Stream, y=emmean, ymin=asymp.LCL, ymax=asymp.UCL, color=as.factor(abspe_early_theta))) +
  geom_point(position = position_dodge(width = .6), size=2.5) +
  geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + facet_wrap(~rewFunc) +
  theme_bw(base_size=12) +  ylab("BOLD response to reward prediction error")  +
  scale_color_manual("Absolute PE\nearly theta\nresponse", values=c("#1b3840","#4fa3b8"), labels = c("10th %ile", "90th %ile")) +
  labs(shape = "Absolute PE\nearly theta\nresponse") +
  theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10"))

# # effects of omission and abs(PE) early theta are independent
# m4 <- lmer(value ~ rewFunc * Stream * abspe_early_theta + rewFunc * Stream * omission_early_theta + (1|id), df)
# summary(m4)
# Anova(m4, '3')
# vif(m4)
# em4 <- as_tibble(emmeans(m4, data = df, ~Stream|omission_early_theta, at = list(omission_early_theta = c(-.04, .61))))
# pe_rewom_theta <- ggplot(em4, aes(x = Stream, y=emmean, ymin=asymp.LCL, ymax=asymp.UCL, color=as.factor(omission_early_theta))) +
#   geom_point(position = position_dodge(width = .6), size=2.5) +
#   geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) +
#   theme_bw(base_size=12) +  ylab("BOLD response to reward prediction error")  +
#   scale_color_manual("Reward omission\nearly theta\nresponse", values=c("#1b3840","#4fa3b8"), labels = c("10th %ile", "90th %ile")) +
#   labs(shape = "Reward omission\nearly theta\nresponse") +
#   theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
#         axis.text=element_text(size=8.5, color="grey10")) # +

# entropy change late beta: weak effects on PEs in the DMN, SM, V
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
df <- ec_betas %>% inner_join(design, by = c("id", "run_number")) %>% inner_join(cond_wbetas, by = c("id", "rewFunc")) %>% inner_join(labels, by = "mask_value") %>% filter(network=="Dors") %>% 
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

sm6b <- lmer(value ~ rewFunc * Stream * entropy_change_late_beta_supp + rewFunc * Stream * omission_early_theta + (1|id), df)
summary(sm6b)
Anova(sm6b, '3')

esm6b <-  as_tibble(emmeans(sm6b, data = df, ~rewFunc|Stream*omission_early_theta, at = list(omission_early_theta = c(-.04, .61))))
# pe_rewom_theta <- ggplot(em1a, aes(x = rewFunc, y=emmean, ymin=asymp.LCL, ymax=asymp.UCL, color=as.factor(omission_early_theta))) +
ec_rewom_early_theta_streams <- ggplot(esm6b, aes(x = Stream, y=emmean, ymin=asymp.LCL, ymax=asymp.UCL, color=as.factor(omission_early_theta))) +
  geom_point(position = position_dodge(width = .6), size=2.5) + facet_wrap(~rewFunc) +
  geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + 
  theme_bw(base_size=12) +  ylab("BOLD response to entropy change")  +
  scale_color_manual("Reward omission\nearly theta\nresponse", values=c("#1b3840","#4fa3b8"), labels = c("10th %ile", "90th %ile")) +
  labs(shape = "Reward omission\nearly theta\nresponse") +
  theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) # +


# df$entropy_change_late_beta_supp 
# n  missing distinct     Info     Mean      Gmd      .05      .10      .25      .50      .75      .90      .95 
# 92800        0       58        1 -0.00563   0.1265 -0.18217 -0.14210 -0.07790 -0.01109  0.08370  0.15427  0.21376 

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


# decomposed version
decm1 <- lmer(value ~ rewFunc *  Stream * entropy_change_late_beta_avg + rewFunc *  Stream * ec_lbeta_wi + (1|id), df)
summary(decm1)
Anova(decm1, '3')

dec1 <- as_tibble(emtrends(decm1, data = df,  var = "entropy_change_late_beta_avg", specs = c("Stream", "rewFunc")))
ec_lbeta_between_decomposed <- 
  ggplot(dec1, aes(x=Stream, y=entropy_change_late_beta_avg.trend, ymin=asymp.LCL, ymax=asymp.UCL)) + 
  geom_point(position = position_dodge(width = .6), size=2.5) + 
  geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + 
  theme_bw(base_size=12) + facet_grid(~rewFunc) + ylab("Effect of late beta suppression,\n between subjects") +
  theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) 

wdec1 <- as_tibble(emtrends(decm1, data = df,  var = "ec_lbeta_wi", specs = c("Stream", "rewFunc")))
ec_lbeta_within_decomposed <- 
  ggplot(wdec1, aes(x=Stream, y=ec_lbeta_wi.trend, ymin=asymp.LCL, ymax=asymp.UCL)) + 
  geom_point(position = position_dodge(width = .6), size=2.5) + 
  geom_errorbar(position = position_dodge(width=0.6), width=0.4, size=0.9) + 
  theme_bw(base_size=12) + facet_grid(~rewFunc) + ylab("Effect of late beta suppression,\n within-subject") +
  theme(axis.title.x=element_blank(), panel.grid.major.x=element_blank(),
        axis.text=element_text(size=8.5, color="grey10")) 


# read in entropy EV_entropy_wiz_clock
e_betas <- read_csv("Schaefer2018_200Parcels_7Networks_order_fonov_2.3mm_ants_cope_l1.csv.gz") %>% filter(l1_cope_name=="EV_entropy_wiz_clock") %>% 
  dplyr::select(id, run_number, l1_model, mask_value, value) %>% mutate(id = as.character(id),
                                                                        mask_value = as.factor(mask_value),
                                                                        run_mc  = scale(run_number, center = T, scale = F))
# merge
df <- e_betas %>% inner_join(cond_wbetas, by = c("id", "rewFunc")) %>% inner_join(design, by = c("id", "run_number"))  %>% inner_join(labels, by = "mask_value") %>% filter(network=="Dors") %>% inner_join(dan_labels, by = "mask_value")

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
ggarrange(rewom_rewom_theta_streams, pe_rewom_theta, pe_pe_beta, pe_ec_lbeta, ec_ec_lbeta, e_ec_lbeta, ncol = 2, nrow = 3)
dev.off()
