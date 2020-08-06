##  updated MEDUSA lmer scripts with only the newest plots 
#   for ramps (Fig. 3) and post-feedback responses (Fig. 4 A-D)
#   by running 'load_medusa_data.R' this also prepares data for
#   decoding analyses (Fig. 4 E-G)
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
library(mlVAR)
library(viridis)
library(car)
library(data.table)
library(emmeans)
library(wesanderson)
#source('~/code/Rhelpers/screen.lmerTest.R')
#source('~/code/Rhelpers/vif.lme.R')

#####################

# select data
smooth_in_mask = T  # main analysis: data smoothed within mask
unsmoothed = F      # no smoothing whatsoever
newmask = F         # sensivitivy analysis: restrictive COBRA mask (default: Harvard-Oxford)
reprocess = F

repo_directory <- "~/code/clock_analysis"
#repo_directory <- "~/Data_Analysis/clock_analysis"

#load the data, and reprocess if requested
source(file.path(repo_directory, "fmri/keuka_brain_behavior_analyses/load_medusa_data_vmPFC.R"))

plots = F
#######
# RAMPS 
# filter by ITI, RT, and evt_time <3 to focus on online periods for longer-ITI trials (to minimize cross-trial contamination)
setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/vmPFC/figs/')
rvdf <- rtvmax_comb %>% filter(online == "TRUE" & iti_prev>1 & iti_ideal > 2 & rt_csv > 1 & rewFunc!="CEVR" & evt_time < 3)
rvdf <- rvdf %>% mutate(`vmPFC response` = decon_interp, entropy = case_when(
  entropy_lag == 'high' ~ 'High entropy',
  entropy_lag == 'low' ~ 'Low entropy'),
  bin6 = round((bin_num + .5)/2,digits = 0),
  bin6_f = as.factor(bin6),
  bin_num = as.factor(bin_num))


# basic descriptive model
rm1 <- lmer(decon_interp ~ evt_time*bin_center_z + evt_time_sq*bin_center_z + entropy_lag + reward_lag + (1 | id/run), rvdf)
summary(rm1)
Anova(rm1, '3')
vif(rm1)

# add entropy modulation
rm2 <- lmer(decon_interp ~ (scale(evt_time) + bin_center_z + entropy_lag) ^2 + (scale(evt_time_sq) + bin_center_z + entropy_lag) ^2 + reward_lag + scale(rt_csv) + (1 | id/run), rvdf)
summary(rm2)
vif(rm2)
Anova(rm2, '3')
anova(rm1,rm2)

em2 <- as_tibble(emmeans(rm2,specs = c("evt_time_sq", "bin_center_z", "entropy_lag", "evt_time"), at = list(bin_center_z = c(-2,-1, 0, 1,2), evt_time_sq = c(0,2,4), evt_time = c(-2,-1,0,1,2))))
em2 <- em2 %>% mutate(`vmPFC response` = emmean, `Time, squared` = as.numeric(as.character(em2$evt_time)), entropy = case_when(
  entropy_lag == 'high' ~ 'High entropy',
  entropy_lag == 'low' ~ 'Low entropy'), 
  time = case_when(
    evt_time == -2 & evt_time_sq == 4 ~ -2,
    evt_time == -1 & evt_time_sq == 2 ~ -1,
    evt_time == 0 & evt_time_sq == 0 ~ 0,
    evt_time == 1 & evt_time_sq == 2 ~ 1,
    evt_time == 2 & evt_time_sq == 4 ~ 2
  )
)

# # Supp. figs
# pdf("ramps_in_vmPFC_lin_quad.pdf", width = 6, height = 3)
# ggplot(em2, aes(time, `vmPFC response`, color = as.factor(bin_center_z))) + 
#   geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL)) + geom_line(size = 1.5) + scale_color_viridis_d() + theme_dark() + facet_wrap(~entropy) + theme(legend.position = "none") +
#   geom_vline(xintercept = 0, lty = 'dashed', color = 'red', size = 1.5) + xlab('Time') + scale_x_continuous(breaks = c(-2,-1,0,1,2)) + ylab('Hippocampal response')
# dev.off()
# 
# wesanderson palette for better visualization of the long-axis gradient
library(wesanderson)
pal = wes_palette("Zissou1", 12, type = "continuous")

if (plots) {
pdf("ramps_in_vmPFC_lin_quad_anderson.pdf", width = 6, height = 3)
ggplot(em2, aes(time, `vmPFC response`, group = bin_center_z, color = bin_center_z)) + 
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),position = position_dodge(width = .5), size = .5) + geom_line(size = 1.5, position = position_dodge(width = .5)) + facet_wrap(~entropy) + theme(legend.position = "none") +
  geom_vline(xintercept = 0, lty = 'dashed', color = 'red', size = 1.5) + xlab('Time') + scale_x_continuous(breaks = c(-2,-1,0,1,2)) + ylab('vmPFC response') +
  scale_color_gradientn(colors = pal, guide = 'none') + 
  theme(legend.title = element_blank(),
        panel.grid.major = element_line(colour = "grey45"), 
        panel.grid.minor = element_line(colour = "grey45"), 
        panel.background = element_rect(fill = 'grey40'))

dev.off()

# # tested completely general long axis bin (without assuming a linear long-axis gradient) -- BUT that worsens fit by 80 AIC points
rm2binf <- lmer(decon_interp ~ (evt_time + bin_num + entropy_lag) ^3 + (evt_time_sq + bin_num + entropy_lag) ^3 + reward_lag + scale(rt_csv) + (1 | id/run), rvdf)
summary(rm2binf)
Anova(rm2binf, '3')

# (tested the same with time as factor -- that results in a singular fit due to 0 variance for ID)
# side RE has a variance of 0


#############
# FINAL MODELS
# Fig. 3A
if (plots) {
setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/vmPFC/figs/')
pdf('3a_smoothed_ramps_vmPFC_cobra_percent_anderson.pdf', width = 4, height = 2.5)

# this bit makes supplemental figure with COBRA mask
# setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/hippo/figs/final/supp')
# pdf('supp3a_smoothed_ramps_cobra_percent_anderson.pdf', width = 4, height = 2.5)
# ggplot(rvdf[!is.na(rvdf$entropy_lag),], aes(evt_time, decon_interp, color = bin_num)) + geom_smooth(method = "loess", se = F) + scale_color_viridis_d() + theme_dark() + facet_wrap(~entropy_lag)
ggplot(rvdf[!is.na(rvdf$entropy_lag),], aes(evt_time, `vmPFC response`, color = as.numeric(bin_center_z), group = as.numeric(bin_center_z))) + 
  geom_smooth(method = 'gam', formula = y ~ splines::ns(x,3), se = F, size = .75)  + 
  facet_wrap(~entropy) + theme(legend.position = "none") + geom_vline(xintercept = 0, lty = 'dashed', color = 'white', size = 1) + xlab('Time (seconds)') + 
  scale_x_continuous(breaks = c(-2,0,2)) + ylab('vmPFC (AU)') + 
  scale_color_gradientn(colors = pal, guide = 'none') + 
  #geom_text(aes(x=-.5, y = .485, label = "RT(Vmax)"), angle = 90, color = "white", size = 2) + #this mysteriously makes the PDF load super-slowly on my ailing Mac
  theme_bw(base_size=13) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_line(colour = "grey45"), 
        panel.grid.minor = element_line(colour = "grey45"), 
        panel.background = element_rect(fill = 'grey40'),
        axis.title.y = element_text(margin=margin(r=6)),
        axis.title.x = element_text(margin=margin(t=6)))
dev.off()

# Main Fig. 3B (ramps, linear slice, general time)
rm2f <- lmer(decon_interp ~ (evt_time_f + bin_center_z + entropy_lag) ^2 + reward_lag + scale(rt_csv) + (1 | id/run), rvdf)
while (any(grepl("failed to converge", rm2f@optinfo$conv$lme4$messages) )) {
  print(rm2f@optinfo$conv$lme4$conv)
  ss <- getME(rm2f,c("theta","fixef"))
  rm2f <- update(rm2f, start=ss)}
summary(rm2f)
vif.lme(rm2f)
Anova(rm2f, '3')
em2f <- as.data.frame(emmeans(rm2f,specs = c("evt_time_f", "bin_center_z", "entropy_lag"), at = list(bin_center_z = c(-2,-1, 0, 1,2))))
em2f <- em2f %>% mutate(`Hippocampal response` = emmean, evt_time = as.numeric(as.character(em2f$evt_time_f)), entropy = case_when(
  entropy_lag == 'high' ~ 'High entropy',
  entropy_lag == 'low' ~ 'Low entropy'
))

# Fig. 3B
# setwd('../final/fig_3/')

# alternatively for supplemental figure 3:
# setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/hippo/figs/final/supp')
# pdf("supp_3b_ramps_in_AH_f_smoothed_in_mask_anderson.pdf", width = 4, height = 2.5)
setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/hippo/figs/final/fig_3/')
pdf("3b_ramps_in_vmPFC_f_smoothed_in_mask_anderson.pdf", width = 4, height = 2.5)
label1 <- expression("RT[Vmax]")
ggplot(em2f, aes(evt_time, `Hippocampal response`, color = bin_center_z, group = bin_center_z)) + 
  # geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),position = position_dodge(width = .5), size = .5) + 
  geom_line(position = position_dodge(width = .5), size = 1) +  facet_wrap(~entropy) + theme(legend.position = "none") +
  geom_linerange(aes(ymin = emmean - SE, ymax = emmean + SE),position = position_dodge(width = .5),  color = "grey80") + 
  geom_vline(xintercept = 0, lty = 'dashed', color = 'white', size = 1)+ xlab('Time (seconds)') + ylab('vmPFC response (AU)') +
  scale_color_gradientn(colors = pal, guide = 'none') + 
  #geom_text(aes(x=-.5, y = .485, label = "RT(Vmax)"), angle = 90, color = "white", size = 2) +
  theme_bw(base_size=13) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_line(colour = "grey45"), 
        panel.grid.minor = element_line(colour = "grey45"), 
        panel.background = element_rect(fill = 'grey40'),
        axis.title.y = element_text(margin=margin(r=6)),
        axis.title.x = element_text(margin=margin(t=6)))
dev.off()
}


##############################
# overall AP response to reinforcement by trial
fb_comb$trial_neg_inv_sc = scale(-1/fb_comb$run_trial)
fb_comb$bin6_f = as.factor(fb_comb$bin6)
fb_comb <- fb_comb %>% mutate(trial_f = as.factor(round((run_trial + 10)/12,digits = 0))) # also a 6-bin version
tm1 <- lmer(decon_interp ~ (bin6_f + evt_time_f + trial_f)^2 + scale(rt_csv)*evt_time_f + (1 | id/run), fb_comb %>% filter (iti_prev>1 & iti_ideal>8 & evt_time < 9))
while (any(grepl("failed to converge", tm1@optinfo$conv$lme4$messages) )) {
  print(tm1@optinfo$conv$lme4$conv)
  ss <- getME(tm1,c("theta","fixef"))
  tm1 <- update(tm1, start=ss)}

summary(tm1)
vif(tm1)
Anova(tm1, '3')
# emt <- as.data.frame(emmeans(tm1,specs = c("trial_neg_inv_sc", "bin_center_z"), at = list(bin_center_z = c(-2,-1, 0, 1,2), trial_neg_inv_sc = c(-2,-1, 0, 1,2))))
# emt <- as.data.frame(emmeans(tm1,specs = c("trial_neg_inv_sc", "bin6_f"), at = list(trial_neg_inv_sc = c(-2,-1, 0, 1,2))))
emt <- as.data.frame(emmeans(tm1,specs = c("trial_f", "bin6_f")))

emt <- emt %>% mutate(`vmPFC response` = emmean, epoch = case_when(
  trial_f == 1 ~ '0-10',
  trial_f == 2 ~ '11-20',
  trial_f == 3 ~ '21-30',
  trial_f == 4 ~ '31-40',
  trial_f == 5 ~ '41-50',
))
# Fig. 4D
#pdf("../supp/supp4d_trial_hipp_AH_PH_bin6_f.pdf", width = 3, height = 3)

# ggplot(rtvmax_comb,aes(run_trial,decon_interp, color = axis_bin, lty = reward)) + geom_smooth(method = "gam", formula = y ~ splines::ns(x,3),  se = F) + scale_color_viridis_d() + theme_dark()
fig4d <- ggplot(emt, aes(epoch, `vmPFC response`, color = as.numeric(bin6_f), group = as.numeric(bin6_f))) + 
  # geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),position = position_dodge(width = .5), size = .5) + geom_line(size = 1.5,position = position_dodge(width = .5)) +  theme(legend.position = "none") +
  # geom_linerange(aes(ymin = asymp.LCL, ymax = asymp.UCL),position = position_dodge(width = .5), size = .5, color = "grey80") +  geom_line(size = 1.5,position = position_dodge(width = .5)) +  theme(legend.position = "none") +
  geom_line(size = 1,position = position_dodge(width = .5)) +  theme(legend.position = "none") +
  geom_linerange(aes(ymin = emmean - SE, ymax = emmean + SE), position = position_dodge(width = .5), size = .5, color = "grey80") +
  geom_point(position = position_dodge(width = .5), size = .5, color = "black", shape=1) +
  xlab('Trial (binned)') + ylab('vmPFC response (AU)') +
  scale_color_gradientn(colors = pal, guide = 'none') + 
  theme_bw(base_size=13) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_line(colour = "grey45"), 
        panel.grid.minor = element_line(colour = "grey45"), 
        panel.background = element_rect(fill = 'grey40'),
        axis.title.y = element_text(margin=margin(r=6)),
        axis.title.x = element_text(margin=margin(t=6)),
        axis.text.x = element_text(size=10)) #condense further in illustrator

# ggsave("../fig_4/4d_trial_hipp_AH_PH_bin6_f.pdf", fig4d, width = 3, height = 3, useDingbats=FALSE)

# post-hoc for trial * bin
em <- emmeans(tm1,specs = c("trial_f", "bin6_f"))
pwpp(em, sort = FALSE, method = "trt.vs.ctrl1", type = "response", side = ">")
CLD = emmeans::CLD(em,
          alpha=0.05,
          Letters=letters,
          adjust="tukey")

# Fig. 4C

# first inspect all data that go into this analysis -- nothing too worrisome.  Some subjects have constricted evt_time ranges (responded mostly early).
# also, more variability in AH than in PH
pdf('ind_ramps_cobra_percent.pdf', width = 30, height = 30)
ggplot(rvdf[!is.na(rvdf$entropy_lag),], aes(evt_time, decon_interp, color = bin_num)) + geom_jitter() + scale_color_viridis_d() + theme_dark() + facet_wrap(id~entropy_lag)
dev.off()



pdf("4c_fb_vmPFC_AP_trial_anderson.pdf", width = 3, height = 3)
#pdf("../supp/supp4c_fb_hipp_AP_trial_anderson.pdf", width = 3, height = 3)

ggplot(fb_comb %>% filter (iti_prev>1 & iti_ideal>8 & evt_time < 9), aes(run_trial, decon_interp, color = bin6, group = bin6)) + 
  geom_smooth(method = "gam", formula = y~splines::ns(x,3),  se = T, size = 1.2) + 
  scale_color_gradientn(colors = pal, guide = 'none') + 
  xlab('Trial') + ylab('vmPFC response (AU)') +
  theme_bw(base_size=13) +
  # theme(legend.title = element_blank(),
  theme(
        panel.grid.major = element_line(colour = "grey45"), 
        panel.grid.minor = element_line(colour = "grey45"), 
        panel.background = element_rect(fill = 'grey40'),
        axis.title.y = element_text(margin=margin(r=6)),
        axis.title.x = element_text(margin=margin(t=6))) +
  scale_y_continuous(breaks=c(0.49, 0.50, 0.51))
dev.off()



# check that it's specific to entropy and not last reward
# it is, but PH is also more active after omissions and AH, after rewards
rm3 <- lmer(decon_interp ~ (evt_time + bin_center_z + reward_lag) ^2 + (evt_time_sq + bin_center_z + reward_lag) ^2 + decon_prev_z + reward_lag + scale(rt_csv) +  (1 | id/run), rvdf)
summary(rm3)
vif.lme(rm3)
Anova(rm3, '3')
g <- ggpredict(rm3, terms = c("evt_time","bin_center_z [-2,-1, 0, 1,2]", "reward_lag"))
g <- plot(g, facet = F, dodge = .3)
g + scale_color_viridis_d(option = "plasma") + theme_dark()

anova(rm1,rm2, rm3)


#### final descriptive model -- remove swing
ee6 <- lmer(decon_interp ~ bin_center_z*evt_time_f*reward + scale(rt_csv)*evt_time_f + (1 | id/run), fb_comb %>% filter (iti_prev>1 & iti_ideal>8 & evt_time < 9))
summary(ee6)
car::Anova(ee6, '3')
vif.lme(ee6)
em6 <- as.data.frame(emmeans(ee6, "bin_center_z", by = c( "evt_time_f", "reward"), at = list( bin_center_z = c(-2,-1, 0, 1, 2)))) %>% 
  mutate(reward_text = case_when(
    reward == 'reward' ~ 'Reward',
    reward == 'omission' ~ 'Omission'
  ))
# ggplot(em6, aes(evt_time_f, emmean, color = bin_center_z)) + 
#   geom_point() + geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL)) + facet_wrap(.~reward) + scale_color_viridis() + theme_dark()

#######
# Fig. 4B (post-feedback MEDUSA)
# 
setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/vmPFC/figs')
pal = wes_palette("Zissou1", 24, type = "continuous")
fig4b <- ggplot(em6, aes(as.numeric(evt_time_f), emmean, color = bin_center_z, group = bin_center_z)) + 
  # geom_point(position = position_dodge2(width = 1)) + geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),position = position_dodge2(width = .2)) + geom_line(position = position_dodge2(width = 1)) + facet_wrap(.~reward_text) +
  geom_line(position = position_dodge2(width = 1), size = 1) + facet_wrap(.~reward_text) +
  geom_linerange(aes(ymin = emmean - SE, ymax = emmean + SE),position = position_dodge2(width = 1), color = 'grey80') + 
  geom_point(position = position_dodge(width = 1), size = .5, color = "black", shape=1) +
  scale_color_gradientn(colors = pal, guide = 'none') + xlab("Time after feedback (seconds)") + ylab("vmPFC response (AU)") +
  theme_bw(base_size=13) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_line(colour = "grey45"), 
        panel.grid.minor = element_line(colour = "grey45"), 
        panel.background = element_rect(fill = 'grey40'),
        axis.title.y = element_text(margin=margin(r=6)),
        axis.title.x = element_text(margin=margin(t=6))) +
  scale_x_continuous(breaks=c(0, 2, 4, 6, 8, 10))

ggsave('4b_medusa_feedback_vmPFC_reward_anderson.pdf', fig4b, height = 3, width = 5, useDingbats=FALSE)

# smoothed raw data
fb_comb <- fb_comb %>%  mutate(reward_text = case_when(
  reward == 'reward' ~ 'Reward',
  reward == 'omission' ~ 'Omission'
))

# more reward modulation caudally
pdf('4a_medusa_feedback_vmPFC_reward_anderson_raw_smoothed.pdf', height = 3, width = 5)
ggplot(fb_comb %>% filter (iti_prev>1 & iti_ideal>8 & evt_time < 9), aes(as.numeric(evt_time_f), decon_interp, group = bin_center_z, color = bin_center_z)) + 
  geom_smooth(method = "gam", formula = y~splines::ns(x,3), se = F, size = .75) +   
  scale_color_gradientn(colors = pal) + xlab("Time after feedback (seconds)") + ylab("Hippocampal response (AU)") +
  theme_bw(base_size=13) +
  theme( panel.grid.major = element_line(colour = "grey45"), 
        panel.grid.minor = element_line(colour = "grey45"), 
        panel.background = element_rect(fill = 'grey40'),
        axis.title.y = element_text(margin=margin(r=6)),
        axis.title.x = element_text(margin=margin(t=6))) + facet_wrap(~reward_text) +
  scale_x_continuous(breaks=c(0, 2, 4, 6, 8, 10))

dev.off()
}

