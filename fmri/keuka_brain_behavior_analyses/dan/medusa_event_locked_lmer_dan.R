##  updated MEDUSA lmer scripts with only the newest plots 
#   for ramps (cf NCOMM Fig. 3) and post-feedback responses (cf NCOMM Fig. 4 A-D)
#   by running 'load_medusa_data.R' this also prepares data for
#   decoding analyses (cf NCOMM Fig. 4 E-G)
library(dplyr)
library(tidyverse)
library(psych)
library(ggcorrplot)
library(lme4)
library(ggpubr)
library(cowplot)
# library(sjPlot)
# library(sjmisc)
library(ggeffects)
library(viridis)
library(car)
library(data.table)
library(emmeans)
library(wesanderson)
library(directlabels)
source('~/code/Rhelpers/screen.lmerTest.R')
source('~/code/Rhelpers/vif.lme.R')

#####################

# select data
smooth_in_mask = T  # main analysis: data smoothed within mask
unsmoothed = F      # no smoothing whatsoever
newmask = F         # sensivitivy analysis: restrictive COBRA mask (default: Harvard-Oxford)
reprocess = F
if (!reprocess) {
  wide_only = T  
  tall_only = F
}

repo_directory <- "~/code/clock_analysis"
#repo_directory <- "~/Data_Analysis/clock_analysis"

#load the data, and reprocess if requested
source(file.path(repo_directory, "fmri/keuka_brain_behavior_analyses/dan/load_medusa_data_dan.R"))
setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/dan/plots/')

##### label broader regions
rt_comb <- rt_comb %>% 
  mutate(region = case_when(
    startsWith(label, "1_") ~ "prefrontal",
    startsWith(label, "1a_") ~ "fronto-opercular",
    startsWith(label, "2_") ~ "IPL",
    startsWith(label, "3_") ~ "SPL",
    startsWith(label, "4_") ~ "IPS",
    startsWith(label, "5_") ~ "temporal",
  ))

##### plot RT-aligned signals --
# p <- ggplot(rt_comb %>% filter (iti_prev>1 & iti_ideal>8 & evt_time < 9), aes(evt_time, decon_interp, color = label, group = label)) + 
#   geom_smooth(method = "gam", formula = y~splines::ns(x,3),  se = F, size = 1.2) + 
#   # scale_color_gradientn(colors = pal, guide = 'none') + 
#   xlab('Time after feedback') + ylab('DAN response (AU)') +
#   theme_bw(base_size=13) +
#   theme(legend.title = element_blank(),
#         panel.grid.major = element_line(colour = "grey45"), 
#         panel.grid.minor = element_line(colour = "grey45"), 
#         panel.background = element_rect(fill = 'grey40'),
#         axis.title.y = element_text(margin=margin(r=6)),
#         axis.title.x = element_text(margin=margin(t=6))) +
#   scale_y_continuous(breaks=c(0.3, 0.4, 0.50, 0.6)) + 
#   scale_x_continuous(expansion(add = 0.5)) + geom_vline(xintercept = 0, color = 'red', size = 1.5, lty = "dashed", alpha = .5) + 
#   scale_colour_discrete(guide = 'none') + facet_wrap(~region)
# pdf("DAN_rt_locked_signals_all_regions.pdf", width = 12, height = 12)
#  direct.label(p,list(cex = 1,"angled.boxes"))
# dev.off()


 ################
# overall DAN response to reinforcement by trial ---
fb_comb$trial_neg_inv_sc = scale(-1/fb_comb$run_trial)
# fb_comb$bin6_f = as.factor(fb_comb$bin6)
fb_comb <- fb_comb %>% mutate(trial_f = as.factor(round((run_trial + 10)/12,digits = 0))) # also a 6-bin version
tm1 <- lmer(decon_interp ~ (label + evt_time_f + trial_f)^2 + scale(rt_csv)*evt_time_f + (1 | id/run), fb_comb %>% filter (iti_prev>1 & iti_ideal>8 & evt_time < 9))
while (any(grepl("failed to converge", tm1@optinfo$conv$lme4$messages) )) {
  print(tm1@optinfo$conv$lme4$conv)
  ss <- getME(tm1,c("theta","fixef"))
  tm1 <- update(tm1, start=ss)}

summary(tm1)
vif.lme(tm1)
Anova(tm1, '3')
# emt <- as.data.frame(emmeans(tm1,specs = c("trial_neg_inv_sc", "bin_center_z"), at = list(bin_center_z = c(-2,-1, 0, 1,2), trial_neg_inv_sc = c(-2,-1, 0, 1,2))))
# emt <- as.data.frame(emmeans(tm1,specs = c("trial_neg_inv_sc", "bin6_f"), at = list(trial_neg_inv_sc = c(-2,-1, 0, 1,2))))
emt <- as.data.frame(emmeans(tm1,specs = c("trial_f", "label")))

emt <- emt %>% mutate(`DAN response` = emmean, epoch = case_when(
  trial_f == 1 ~ '0-10',
  trial_f == 2 ~ '11-20',
  trial_f == 3 ~ '21-30',
  trial_f == 4 ~ '31-40',
  trial_f == 5 ~ '41-50',
))
# Fig. 4D
#pdf("../supp/supp4d_trial_hipp_AH_PH_bin6_f.pdf", width = 3, height = 3)

# ggplot(rtvmax_comb,aes(run_trial,decon_interp, color = axis_bin, lty = reward)) + geom_smooth(method = "gam", formula = y ~ splines::ns(x,3),  se = F) + scale_color_viridis_d() + theme_dark()
fig4d <- ggplot(emt, aes(epoch, `DAN response`, color = label, group = label)) + 
  # geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),position = position_dodge(width = .5), size = .5) + geom_line(size = 1.5,position = position_dodge(width = .5)) +  theme(legend.position = "none") +
  # geom_linerange(aes(ymin = asymp.LCL, ymax = asymp.UCL),position = position_dodge(width = .5), size = .5, color = "grey80") +  geom_line(size = 1.5,position = position_dodge(width = .5)) +  theme(legend.position = "none") +
  geom_line(size = 1,position = position_dodge(width = .5)) + 
  geom_linerange(aes(ymin = emmean - SE, ymax = emmean + SE), position = position_dodge(width = .5), size = .5, color = "grey80") +
  geom_point(position = position_dodge(width = .5), size = .5, color = "black", shape=1) +
  xlab('Trial (binned)') + ylab('DAN response (AU)') +
  # scale_color_gradientn(colors = pal, guide = 'none') +  theme(legend.position = "none") +
  theme_bw(base_size=13) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_line(colour = "grey45"), 
        panel.grid.minor = element_line(colour = "grey45"), 
        panel.background = element_rect(fill = 'grey40'),
        axis.title.y = element_text(margin=margin(r=6)),
        axis.title.x = element_text(margin=margin(t=6)),
        axis.text.x = element_text(size=10)) #condense further in illustrator

ggsave("../4d_trial_DAN.pdf", fig4d, width = 8, height = 6, useDingbats=FALSE)

# post-hoc for trial * bin
em <- emmeans(tm1,specs = c("trial_f", "label"))
pwpp(em, sort = FALSE, method = "trt.vs.ctrl1", type = "response", side = ">")
CLD = emmeans::CLD(em,
          alpha=0.05,
          Letters=letters,
          adjust="tukey")



pdf("4c_fb_DAN_trial.pdf", width = 6, height = 6)
#pdf("../supp/supp4c_fb_hipp_AP_trial_anderson.pdf", width = 3, height = 3)

ggplot(fb_comb %>% filter (iti_prev>1 & iti_ideal>8 & evt_time < 9), aes(run_trial, decon_interp, color = label, group = label)) + 
  geom_smooth(method = "gam", formula = y~splines::ns(x,3),  se = F, size = 1.2) + 
  # scale_color_gradientn(colors = pal, guide = 'none') + 
  xlab('Trial') + ylab('DAN response (AU)') +
  theme_bw(base_size=13) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_line(colour = "grey45"), 
        panel.grid.minor = element_line(colour = "grey45"), 
        panel.background = element_rect(fill = 'grey40'),
        axis.title.y = element_text(margin=margin(r=6)),
        axis.title.x = element_text(margin=margin(t=6))) +
  scale_y_continuous(breaks=c(0.49, 0.50, 0.51))
dev.off()

# drill down on the effect of trial on the post-feedback timecourse

pdf("post_feedback_DAN_learning.pdf", width = 30, height = 30)
#pdf("../supp/supp4c_fb_hipp_AP_trial_anderson.pdf", width = 3, height = 3)
ggplot(fb_comb %>% filter (iti_prev>1 & iti_ideal>8 & evt_time < 9), aes(evt_time, decon_interp, color = !first10, group = !first10)) + 
  geom_smooth(method = "gam", formula = y~splines::ns(x,3),  se = F, size = 1.2) + facet_wrap(~label) +
  # scale_color_gradientn(colors = pal, guide = 'none') + 
  xlab('Time after feedback') + ylab('DAN response (AU)') +
  theme_bw(base_size=13) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_line(colour = "grey45"), 
        panel.grid.minor = element_line(colour = "grey45"), 
        panel.background = element_rect(fill = 'grey40'),
        axis.title.y = element_text(margin=margin(r=6)),
        axis.title.x = element_text(margin=margin(t=6))) +
  scale_y_continuous(breaks=c(0.49, 0.50, 0.51))
dev.off()


# # check that it's specific to entropy and not last reward
# # it is, but PH is also more active after omissions and AH, after rewards
# rm3 <- lmer(decon_interp ~ (evt_time + bin_center_z + reward_lag + side) ^2 + (evt_time_sq + bin_center_z + reward_lag + side) ^2 + decon_prev_z + reward_lag + scale(rt_csv) +  (1 | id/run), rvdf)
# summary(rm3)
# vif.lme(rm3)
# Anova(rm3, '3')
# g <- ggpredict(rm3, terms = c("evt_time","bin_center_z [-2,-1, 0, 1,2]", "reward_lag"))
# g <- plot(g, facet = F, dodge = .3)
# g + scale_color_viridis_d(option = "plasma") + theme_dark()
# 
# anova(rm1,rm2, rm3)


#### final descriptive model -- remove swing
ee6 <- lmer(decon_interp ~ label*evt_time_f*reward + scale(rt_csv)*evt_time_f + (1 | id/run) + (1 | side), fb_comb %>% filter (iti_prev>1 & iti_ideal>8 & evt_time < 9))
summary(ee6)
car::Anova(ee6, '3')
vif.lme(ee6)
em6 <- as.data.frame(emmeans(ee6, "label", by = c( "evt_time_f", "reward"), at = list( bin_center_z = c(-2,-1, 0, 1, 2)))) %>% 
  mutate(reward_text = case_when(
    reward == 'reward' ~ 'Reward',
    reward == 'omission' ~ 'Omission'
  ))
# ggplot(em6, aes(evt_time_f, emmean, color = bin_center_z)) + 
#   geom_point() + geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL)) + facet_wrap(.~reward) + scale_color_viridis() + theme_dark()

#######
# Fig. 4B (post-feedback MEDUSA)
# 
# setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/hippo/figs/final/fig_4')
# pal = wes_palette("Zissou1", 24, type = "continuous")
fig4b <- ggplot(em6, aes(as.numeric(evt_time_f), emmean, color = label, group = label)) + 
  # geom_point(position = position_dodge2(width = 1)) + geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),position = position_dodge2(width = .2)) + geom_line(position = position_dodge2(width = 1)) + facet_wrap(.~reward_text) +
  geom_line(position = position_dodge2(width = 1), size = 1) + facet_wrap(.~reward_text) +
  geom_linerange(aes(ymin = emmean - SE, ymax = emmean + SE),position = position_dodge2(width = 1), color = 'grey80') + 
  geom_point(position = position_dodge(width = 1), size = .5, color = "black", shape=1) +
  # scale_color_gradientn(colors = pal, guide = 'none') + 
  xlab("Time after feedback (seconds)") + ylab("DAN response (AU)") +
  theme_bw(base_size=13) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_line(colour = "grey45"), 
        panel.grid.minor = element_line(colour = "grey45"), 
        panel.background = element_rect(fill = 'grey40'),
        axis.title.y = element_text(margin=margin(r=6)),
        axis.title.x = element_text(margin=margin(t=6))) +
  scale_x_continuous(breaks=c(0, 2, 4, 6, 8, 10))

ggsave('4b_medusa_feedback_DAN_reward.pdf', fig4b, height = 5, width = 8, useDingbats=FALSE)

# smoothed raw data
fb_comb <- fb_comb %>%  mutate(reward_text = case_when(
  reward == 'reward' ~ 'Reward',
  reward == 'omission' ~ 'Omission'
))

pdf('4a_medusa_feedback_DAN_reward_raw_smoothed.pdf', height = 5, width = 8)
ggplot(fb_comb %>% filter (iti_prev>1 & iti_ideal>8 & evt_time < 9), aes(as.numeric(evt_time_f), decon_interp, group = label, color = label)) + 
  geom_smooth(method = "gam", formula = y~splines::ns(x,3), se = F, size = .75) +   
  # scale_color_gradientn(colors = pal, guide = 'none') + 
  xlab("Time after feedback (seconds)") + ylab("DAN response (AU)") +
  theme_bw(base_size=13) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_line(colour = "grey45"), 
        panel.grid.minor = element_line(colour = "grey45"), 
        panel.background = element_rect(fill = 'grey40'),
        axis.title.y = element_text(margin=margin(r=6)),
        axis.title.x = element_text(margin=margin(t=6))) + facet_wrap(~reward_text) +
  scale_x_continuous(breaks=c(0, 2, 4, 6, 8, 10))

dev.off()

# look at timecourses windowed by trial bins

# split by entropy - no modulation for feedback?
pdf('fb_medusa_entropy.pdf', height = 6, width = 6)
ggplot(fb_comb %>% filter (iti_prev>1 & iti_ideal>8 & evt_time < 9), aes(as.numeric(evt_time_f), decon_interp, lty = entropy)) + 
  geom_smooth(method = "gam", formula = y~splines::ns(x,3), se = T, size = .75) + facet_wrap(~label)
dev.off()
# clock:
pdf('clock_medusa_entropy.pdf', height = 6, width = 6)
ggplot(clock_comb %>% filter (iti_prev>1 & iti_ideal>1 & evt_time < 4 & !is.na(entropy)), aes(evt_time, decon_interp, lty = entropy)) + 
  geom_smooth(method = "gam", formula = y~splines::ns(x,3), se = T, size = .75) + facet_wrap(~label)
dev.off()

pdf('clock_medusa_entropy_lag.pdf', height = 6, width = 6)
ggplot(clock_comb %>% filter (iti_prev>1 & iti_ideal>1 & evt_time < 4 & !is.na(entropy_lag)), aes(evt_time, decon_interp, lty = entropy_lag)) + 
  geom_smooth(method = "gam", formula = y~splines::ns(x,3), se = T, size = .75) + facet_wrap(~label)
dev.off()

# condition
# clock:
pdf('clock_medusa_condition.pdf', height = 6, width = 6)
ggplot(clock_comb %>% filter (iti_prev>1 & iti_ideal>1 & evt_time < 4 & !is.na(entropy_lag)), aes(evt_time, decon_interp, color = rewFunc)) + 
  geom_smooth(method = "gam", formula = y~splines::ns(x,3), se = T, size = .75) + facet_wrap(~label)
dev.off()

# feedback
pdf('fb_medusa_condition.pdf', height = 6, width = 6)
ggplot(fb_comb %>% filter (iti_prev>1 & iti_ideal>6 & evt_time < 6 & !is.na(entropy_lag)), aes(evt_time, decon_interp, color = rewFunc)) + 
  geom_smooth(method = "gam", formula = y~splines::ns(x,3), se = T, size = .75) + facet_wrap(~label)
dev.off()

# reward
pdf('fb_medusa_reward.pdf', height = 8, width = 8)
ggplot(fb_comb %>% filter (iti_prev>1 & iti_ideal>8 & evt_time < 9), aes(as.numeric(evt_time_f), decon_interp, lty = rew_om)) + 
  geom_smooth(method = "gam", formula = y~splines::ns(x,3), se = T, size = .75) + facet_wrap(~label)
dev.off()


