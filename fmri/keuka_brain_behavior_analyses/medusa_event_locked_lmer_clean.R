##  updated MEDUSA lmer scripts with only the newest plots 
#   for ramps and post-feedback responses 

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
source('~/code/Rhelpers/screen.lmerTest.R')
source('~/code/Rhelpers/vif.lme.R')

#####################
reprocess = F
unsmoothed = F
newmask = F
smooth_in_mask = T

repo_directory <- "~/code/clock_analysis"
# repo_directory <- "~/Data_Analysis/clock_analysis"

#load the data, and reprocess if requested
source(file.path(repo_directory, "fmri/keuka_brain_behavior_analyses/load_medusa_data.R"))


#######
# RAMPS 

# preliminary analyses
# filter by ITI, RT, and evt_time <3
setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/hippo/figs/ramps/')
rvdf <- rtvmax_comb %>% filter(online == "TRUE" & iti_prev>1 & iti_ideal > 2 & rt_csv > 1 & rewFunc!="CEVR" & evt_time < 3)
rvdf <- rvdf %>% mutate(`Hippocampal response` = decon_interp, entropy = case_when(
  entropy_lag == 'high' ~ 'High entropy',
  entropy_lag == 'low' ~ 'Low entropy'),
  bin6 = round((bin_num + .5)/2,digits = 0),
  bin6_f = as.factor(bin6),
  bin_num = as.factor(bin_num))


# strangely we only see ramps on long-ITI trials  
# rm* models do not converge with cobra percentchange
rm1 <- lmer(decon_interp ~ evt_time*bin_center_z + evt_time_sq*bin_center_z + entropy_lag + reward_lag + (1 | id/run), rvdf)
summary(rm1)
Anova(rm1, '3')
vif.lme(rm1)

# add entropy modulation
rm2 <- lmer(decon_interp ~ (evt_time + bin_center_z + entropy_lag) ^2 + (evt_time_sq + bin_center_z + entropy_lag) ^2 + reward_lag + scale(rt_csv) + (1 | id/run), rvdf)
summary(rm2)
vif.lme(rm2)
Anova(rm2, '3')
anova(rm1,rm2)

em2 <- as_tibble(emmeans(rm2,specs = c("evt_time_sq", "bin_center_z", "entropy_lag", "evt_time"), at = list(bin_center_z = c(-2,-1, 0, 1,2), evt_time_sq = c(0,2,4), evt_time = c(-2,-1,0,1,2))))
em2 <- em2 %>% mutate(`Hippocampal response` = emmean, `Time, squared` = as.numeric(as.character(em2$evt_time)), entropy = case_when(
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
# Supp. figs
pdf("ramps_in_AH_lin_quad.pdf", width = 6, height = 3)
ggplot(em2, aes(time, `Hippocampal response`, color = as.factor(bin_center_z))) + 
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL)) + geom_line(size = 1.5) + scale_color_viridis_d() + theme_dark() + facet_wrap(~entropy) + theme(legend.position = "none") +
  geom_vline(xintercept = 0, lty = 'dashed', color = 'red', size = 1.5) + xlab('Time') + scale_x_continuous(breaks = c(-2,-1,0,1,2)) + ylab('Hippocampal response')
dev.off()

#wesanderson version 
library(wesanderson)
pal = wes_palette("Zissou1", 12, type = "discrete")
pdf("ramps_in_AH_lin_quad_anderson.pdf", width = 6, height = 3)
ggplot(em2, aes(time, `Hippocampal response`, group = bin_center_z, color = bin_center_z)) + 
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),position = position_dodge(width = .5), size = .5) + geom_line(size = 1.5, position = position_dodge(width = .5)) + facet_wrap(~entropy) + theme(legend.position = "none") +
  geom_vline(xintercept = 0, lty = 'dashed', color = 'red', size = 1.5) + xlab('Time') + scale_x_continuous(breaks = c(-2,-1,0,1,2)) + ylab('Hippocampal response') +
  scale_color_gradientn(colors = pal, guide = 'none') + 
  theme(legend.title = element_blank(),
        panel.grid.major = element_line(colour = "grey45"), 
        panel.grid.minor = element_line(colour = "grey45"), 
        panel.background = element_rect(fill = 'grey40'))

dev.off()

# pdf("ramps_in_AH_cobra.pdf", width = 6, height = 3)
# ggplot(em2, aes(evt_time_sq, `Hippocampal response`, color = as.factor(bin_center_z))) + 
#   geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL)) + geom_line(size = 1.5) + scale_color_viridis_d() + theme_dark() + facet_wrap(~entropy) + theme(legend.position = "none") +
#   geom_vline(xintercept = 0, lty = 'dashed', color = 'red', size = 1.5) + xlab('Time, squared') + scale_x_continuous(breaks = c(0,2,4)) + ylab('Hippocampal response')
# dev.off()
#   
# # add completely general bin -- BUT that worsens fit by 80 AIC points
rm2binf <- lmer(decon_interp ~ (evt_time + bin_num + entropy_lag) ^3 + (evt_time_sq + bin_num + entropy_lag) ^3 + reward_lag + scale(rt_csv) + (1 | id/run), rvdf)
summary(rm2binf)
Anova(rm2binf, '3')


# test the same with time as factor -- that results in a singular fit due to 0 variance for ID
# side RE has a variance of 0
# THE MOST CONVINCING MODEL

#############
# Main Fig. 4
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
# pdf("ramps_in_AH_f_cobra.pdf", width = 6, height = 3)
# ggplot(em2f, aes(evt_time, `Hippocampal response`, color = as.factor(bin_center_z))) + 
#   geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL)) + geom_line(size = 1.5) + scale_color_viridis_d() + theme_dark() + facet_wrap(~entropy) + theme(legend.position = "none") +
#   geom_vline(xintercept = 0, lty = 'dashed', color = 'red', size = 1.5)+ xlab('Time') + ylab('Hippocampal response')
# dev.off()
# anova(rm1,rm2,rm2f, rm2binf)
# 
# wesanderson
pdf("ramps_in_AH_f_cobra_anderson.pdf", width = 5, height = 3)
ggplot(em2f, aes(evt_time, `Hippocampal response`, color = bin_center_z, group = bin_center_z)) + 
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),position = position_dodge(width = .5), size = .5) + geom_line(size = 1.5,position = position_dodge(width = .5)) +  facet_wrap(~entropy) + theme(legend.position = "none") +
  geom_vline(xintercept = 0, lty = 'dashed', color = 'red', size = 1.5)+ xlab('Time') + ylab('Hippocampal response') +
  scale_color_gradientn(colors = pal, guide = 'none') + 
  theme(legend.title = element_blank(),
        panel.grid.major = element_line(colour = "grey45"), 
        panel.grid.minor = element_line(colour = "grey45"), 
        panel.background = element_rect(fill = 'grey40'))
dev.off()


# 
# # make bin a factor
#   # 3-way interaction is NS with any decon
#   rm2ff <- lmer(decon_interp ~ (evt_time_f + bin6_f + entropy_lag) ^3 + scale(rt_csv)*evt_time_f + (1 | id/run), rvdf)
#   summary(rm2ff)
#   vif.lme(rm2f)
#   anova(rm1,rm2,rm2ff)
#   Anova(rm2ff, '3')
#   em2ff <- as.data.frame(emmeans(rm2ff,specs = c("evt_time_f", "bin6_f", "entropy_lag")))
#   em2ff$hipp_response <- em2ff$emmean
#   anova(rm2f,rm2ff)
#   pdf("ramps_in_AH_ff.pdf", width = 8, height = 6)
#   ggplot(em2ff, aes(evt_time_f, hipp_response, group = bin_num, color = bin_num)) + geom_point() +
#     geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL)) + geom_line() + facet_wrap(~entropy_lag) + scale_color_viridis_d() + theme_dark()
#   dev.off()

# also plot smoothed raw ramps data

pdf('smoothed_ramps_cobra_percent.pdf', width = 6, height = 3)
# ggplot(rvdf[!is.na(rvdf$entropy_lag),], aes(evt_time, decon_interp, color = bin_num)) + geom_smooth(method = "loess", se = F) + scale_color_viridis_d() + theme_dark() + facet_wrap(~entropy_lag)
ggplot() + stat_smooth(data = rvdf[!is.na(rvdf$entropy_lag),], aes(evt_time, `Hippocampal response`, color = as.factor(bin_center_z)), geom = 'line', method = "loess", se = T)  + 
  scale_color_viridis_d() + theme_dark() + facet_wrap(~entropy) + theme(legend.position = "none") + geom_vline(xintercept = 0, lty = 'dashed', color = 'red', size = 1.5) + xlab('Time') + 
  scale_x_continuous(breaks = c(-2,0,2)) + ylab('Hippocampal response')
dev.off()

# wesanderson
pdf('smoothed_ramps_cobra_percent_anderson.pdf', width = 5, height = 3)
# ggplot(rvdf[!is.na(rvdf$entropy_lag),], aes(evt_time, decon_interp, color = bin_num)) + geom_smooth(method = "loess", se = F) + scale_color_viridis_d() + theme_dark() + facet_wrap(~entropy_lag)
ggplot(rvdf[!is.na(rvdf$entropy_lag),], aes(evt_time, `Hippocampal response`, color = as.numeric(bin_center_z), group = as.numeric(bin_center_z))) + 
  geom_smooth(method = 'gam', formula = y ~ splines::ns(x,3), se = F)  + 
  facet_wrap(~entropy) + theme(legend.position = "none") + geom_vline(xintercept = 0, lty = 'dashed', color = 'red', size = 1.5) + xlab('Time') + 
  scale_x_continuous(breaks = c(-2,0,2)) + ylab('Hippocampal response') + 
  scale_color_gradientn(colors = pal, guide = 'none') + 
  theme(legend.title = element_blank(),
        panel.grid.major = element_line(colour = "grey45"), 
        panel.grid.minor = element_line(colour = "grey45"), 
        panel.background = element_rect(fill = 'grey40'))
dev.off()


# try and combine two plots
pdf("ramps_in_AH_combined_cobra_percent.pdf", width = 6, height = 3)
ggplot(em2f, aes(evt_time, `Hippocampal response`, group = bin_center_z, color = bin_center_z)) + geom_point() + geom_line() +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL)) + scale_color_viridis_c() + theme_dark() + facet_wrap(~entropy) +
  stat_smooth(data = rvdf[!is.na(rvdf$entropy_lag),], aes(evt_time, decon_interp, group = bin_center_z, color = bin_center_z), geom = 'line', alpha = .2, method = "loess", se = F)  + theme(legend.position = "none") + ylab('Hippocampal response')
dev.off()

# RTvmax-aligned
pal = wes_palette("Zissou1", 24, type = "continuous")
pdf("trial_rtvmax_hipp_AH_PH_gam.pdf", width = 11, height = 8)
# ggplot(rtvmax_comb,aes(run_trial,decon_interp, color = axis_bin, lty = reward)) + geom_smooth(method = "gam", formula = y ~ splines::ns(x,3),  se = F) + scale_color_viridis_d() + theme_dark()
ggplot(rtvmax_comb,aes(run_trial,decon_interp, color = axis_bin)) + geom_smooth(method = "gam", formula = y~splines::ns(x,4)) + 
  scale_color_gradientn(colors = pal, guide = 'none') + 
  theme(legend.title = element_blank(),
        panel.grid.major = element_line(colour = "grey45"), 
        panel.grid.minor = element_line(colour = "grey45"), 
        panel.background = element_rect(fill = 'grey40'))

dev.off()


# overall AP response by trial
fb_comb$trial_neg_inv_sc = scale(-1/fb_comb$run_trial)
fb_comb$bin6_f = as.factor(fb_comb$bin6)
fb_comb <- fb_comb %>% mutate(trial_f = as.factor(round((run_trial + 10)/12,digits = 0))) # also a 6-bin version
tm1 <- lmer(decon_interp ~ (bin6_f + evt_time_f + trial_f)^2 + scale(rt_csv)*evt_time_f + (1 | id/run) + (1 | side), fb_comb %>% filter (iti_prev>1 & iti_ideal>8 & evt_time < 9))
while (any(grepl("failed to converge", tm1@optinfo$conv$lme4$messages) )) {
  print(tm1@optinfo$conv$lme4$conv)
  ss <- getME(tm1,c("theta","fixef"))
  tm1 <- update(tm1, start=ss)}

summary(tm1)
vif.lme(tm1)
Anova(tm1, '3')
# emt <- as.data.frame(emmeans(tm1,specs = c("trial_neg_inv_sc", "bin_center_z"), at = list(bin_center_z = c(-2,-1, 0, 1,2), trial_neg_inv_sc = c(-2,-1, 0, 1,2))))
# emt <- as.data.frame(emmeans(tm1,specs = c("trial_neg_inv_sc", "bin6_f"), at = list(trial_neg_inv_sc = c(-2,-1, 0, 1,2))))
emt <- as.data.frame(emmeans(tm1,specs = c("trial_f", "bin6_f")))

emt <- emt %>% mutate(`Hippocampal response` = emmean, epoch = case_when(
  trial_f == 1 ~ '0-10',
  trial_f == 2 ~ '11-20',
  trial_f == 3 ~ '21-30',
  trial_f == 4 ~ '31-40',
  trial_f == 5 ~ '41-50',
))
pdf("../early_late/trial_hipp_AH_PH_bin6_f.pdf", width = 3, height = 3)
# ggplot(rtvmax_comb,aes(run_trial,decon_interp, color = axis_bin, lty = reward)) + geom_smooth(method = "gam", formula = y ~ splines::ns(x,3),  se = F) + scale_color_viridis_d() + theme_dark()
ggplot(emt, aes(epoch, `Hippocampal response`, color = as.numeric(bin6_f), group = as.numeric(bin6_f))) + 
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),position = position_dodge(width = .5), size = .5) + geom_line(size = 1.5,position = position_dodge(width = .5)) +  theme(legend.position = "none") +
  xlab('Trial') + ylab('Hippocampal response') +
  scale_color_gradientn(colors = pal, guide = 'none') + 
  theme(legend.title = element_blank(),
        panel.grid.major = element_line(colour = "grey45"), 
        panel.grid.minor = element_line(colour = "grey45"), 
        panel.background = element_rect(fill = 'grey40'))
dev.off()





pdf("../early_late/fb_hipp_AP_trial_anderson.pdf", width = 3, height = 3)
ggplot(fb_comb %>% filter (iti_prev>1 & iti_ideal>8 & evt_time < 9), aes(run_trial, decon_interp, color = bin6, group = bin6)) + geom_smooth(method = "gam", formula = y~splines::ns(x,3),  se = T) + 
  scale_color_gradientn(colors = pal, guide = 'none') + 
  xlab('Trial') + ylab('Hippocampal response') +
  theme(legend.title = element_blank(),
        panel.grid.major = element_line(colour = "grey45"), 
        panel.grid.minor = element_line(colour = "grey45"), 
        panel.background = element_rect(fill = 'grey40'))

dev.off()
# inspect all data that go into this analysis -- nothing too worrisome, but hard to read.  Some subjects have constricted evt_time ranges (responded mostly early).
# also, more variability in AH than in PH
pdf('ind_ramps_cobra_percent.pdf', width = 30, height = 30)
ggplot(rvdf[!is.na(rvdf$entropy_lag),], aes(evt_time, decon_interp, color = bin_num)) + geom_jitter() + scale_color_viridis_d() + theme_dark() + facet_wrap(id~entropy_lag)
dev.off()


# check that it's specific to entropy and not last reward
# it is, but PH is also more active after omissions and AH, after rewards
rm3 <- lmer(decon_interp ~ (evt_time + bin_center_z + reward_lag + side) ^2 + (evt_time_sq + bin_center_z + reward_lag + side) ^2 + decon_prev_z + reward_lag + scale(rt_csv) +  (1 | id/run), rvdf)
summary(rm3)
vif.lme(rm3)
Anova(rm3, '3')
g <- ggpredict(rm3, terms = c("evt_time","bin_center_z [-2,-1, 0, 1,2]", "reward_lag"))
g <- plot(g, facet = F, dodge = .3)
g + scale_color_viridis_d(option = "plasma") + theme_dark()

anova(rm1,rm2, rm3)

########
# offline activity in PH vs AH
# commented out historic exploratory analyses

# # comparison model without time * bin
#   om0 <- lmer(decon_interp ~ evt_time_f + bin_center_z*reward + decon_prev_z + scale(rt_csv)*evt_time_f + (1 | id/run) + (1 | side), fb_comb %>% filter (iti_prev>1 & iti_ideal>8 & evt_time < 9))
# # without time * reward * bin
#   om0a <- lmer(decon_interp ~ evt_time_f*bin_center_z + bin_center_z*reward + decon_prev_z + scale(rt_csv)*evt_time_f + (1 | id/run) + (1 | side), fb_comb %>% filter (iti_prev>1 & iti_ideal>8 & evt_time < 9))
#   
#   om1 <- lmer(decon_interp ~ evt_time_f*bin_center_z*reward + decon_prev_z + scale(rt_csv)*evt_time_f + (1 | id/run) + (1 | side), fb_comb %>% filter (iti_prev>1 & iti_ideal>8 & evt_time < 9))
#   summary(om1)
#   vif.lme(om1)
#   car::Anova(om1)
#   
#   anova(om0, om0a, om1)
#   g <- ggpredict(om1, terms = c("evt_time_f", "bin_center_z [-2,-1, 0, 1,2]", "reward"))
#   g <- plot(g, facet = F, dodge = .4)
#   pdf("offline_reward_AH_PH.pdf", width = 6, height = 6)
#   g + scale_color_viridis_d(option = "plasma") + theme_dark()
#   dev.off()
#   
#   om2 <- lmer(decon_interp ~ evt_time_f*bin_center_z*abs_pe_f + decon_prev_z + scale(rt_csv)*evt_time_f + (1 | id/run) + (1 | side), fb_comb %>% filter (iti_prev>1 & iti_ideal>8 & evt_time < 9))
#   summary(om2)
#   vif.lme(om2)
#   car::Anova(om2)
#   g <- ggpredict(om2, terms = c("evt_time_f", "bin_center_z [-2,-1, 0, 1,2]", "abs_pe_f"))
#   g <- plot(g, facet = F, dodge = .4)
#   pdf("offline_abs_pe_AH_PH.pdf", width = 6, height = 6)
#   g + scale_color_viridis_d(option = "plasma") + theme_dark()
#   dev.off()
#   anova(om1, om2)
#   
# # moderation by gamma?
#   
#   om3 <- lmer(decon_interp ~ evt_time_f*bin_center_z*abs_pe_f*reward + decon_prev_z + scale(rt_csv)*evt_time_f + (1 | id/run) + (1 | side), fb_comb %>% filter (iti_prev>1 & iti_ideal>8 & evt_time < 9))
#   summary(om3)
#   vif.lme(om3)
#   car::Anova(om3)
#   g <- ggpredict(om3, terms = c("evt_time_f", "abs_pe_f", "reward"))
#   g <- plot(g, facet = F, dodge = .4)
#   pdf("offline_abs_pe_reward_AH_PH.pdf", width = 6, height = 6)
#   g + scale_color_viridis_d() + theme_dark()
#   dev.off()
#   
# # reversion in MEDUSA
# # need a measure of closeness to RTvmax: let's start with EV
#   e1 = 38
#   xtabs(~bin_center + swing_above_median, fb_comb %>% filter(is.na(decon_interp)))
#   
#   pdf("lose_switch_reversal_to_vmax.pdf", height = 20, width = 20)
# #lty = swing_above_median, 
#   ggplot(fb_comb, aes(evt_time, decon_interp, color = factor(bin_center), lty = next_swing_above_median)) +
#       stat_smooth(method = 'gam',se = F) + scale_color_viridis_d("TEST") + facet_grid
#   dev.off()
# # Do shorter ITIs disrupt processing?  No, they just seem to fall asleep during ITIs
#   itim1 <- lmer(scale(rt_csv) ~ scale(rt_lag)*scale(rt_vmax)*scale(iti_prev) + (1 | id/run), trial_df )
#   summary(itim1)
#   
#   ################
# # responses as a function of trial and explore/exploit
#   
# # preliminary plots: trial
# # trial by condition
#   pdf("trial_fb_hipp_AH_PH_gam_rew.pdf", width = 11, height = 8)
#   ggplot(fb_comb,aes(run_trial,decon_interp, color = axis_bin, lty = reward)) + geom_smooth(method = "gam", se = F) + facet_wrap(~rewFunc) + scale_color_viridis_d() + theme_dark()
#   dev.off()
# # ...and reward
#   pdf("trial_fb_hipp_AH_PH_loess_rew.pdf", width = 11, height = 8)
#   ggplot(fb_comb,aes(run_trial,decon_interp, color = axis_bin, lty = reward)) + geom_smooth(method = "loess", se = F) + facet_wrap(~rewFunc) + scale_color_viridis_d() + theme_dark()
#   dev.off()
#   
#   
# # wave form by trial
#   pdf("fb_hipp_AP_trial_rew.pdf", width = 11, height = 8)
#   ggplot(fb_comb, aes(evt_time, decon_interp, color = axis_bin, lty = reward)) + geom_smooth(method = "gam", se = F) + facet_grid(first10 ~ rewFunc) + scale_color_viridis_d() + theme_dark()
#   dev.off()
#   
# # wave form as a function of next swing?
#   fb_comb <- fb_comb %>% group_by(id,run) %>% mutate(swing_above_median_lead = lead(swing_above_median)) %>% ungroup()
#   
#   pdf("fb_hipp_AP_next_swing_loess_rewFunc.pdf", width = 11, height = 8)
#   ggplot(fb_comb[!is.na(fb_comb$swing_above_median_lead) & fb_comb$evt_time < 8,], aes(evt_time, decon_interp, color = axis_bin, lty = swing_above_median_lead)) + geom_smooth(method = "loess", se = F) + facet_grid(rewFunc~reward)+ scale_color_viridis_d() + theme_dark()
#   dev.off()
#   fb_comb$trial_neg_inv <- -1/fb_comb$run_trial
#   
#   # final plot without swing:
#   pdf("fb_hipp_AP_next_swing_loess_rewFunc.pdf", width = 11, height = 8)
#   ggplot(fb_comb[!is.na(fb_comb$swing_above_median_lead) & fb_comb$evt_time < 8,], aes(evt_time, decon_interp, color = axis_bin, lty = swing_above_median_lead)) + geom_smooth(method = "loess", se = F) + facet_grid(rewFunc~reward)+ scale_color_viridis_d() + theme_dark()
#   dev.off()
#   fb_comb$trial_neg_inv <- -1/fb_comb$run_trial
#   
#   
#   ee1 <- lmer(decon_interp ~ bin_center_z*trial_neg_inv + scale(rt_csv)*evt_time + (1 | id/run) + (1 | side), fb_comb %>% filter (iti_prev>1 & iti_ideal>8 & evt_time < 9))
#   summary(ee1)
#   vif.lme(ee1)
#   car::Anova(ee1)
#   g <- ggpredict(ee1, terms = c(""))
#   g <- plot(g, facet = F, dodge = .4)
#   pdf("offline_abs_pe_reward_AH_PH.pdf", width = 6, height = 6)
#   g + scale_color_viridis_d() + theme_dark()
#   dev.off()
#   
#   ee2 <- lmer(decon_interp ~ bin_center_z*trial_neg_inv*rewFunc + scale(rt_csv)*evt_time + (1 | id/run) + (1 | side), fb_comb %>% filter (iti_prev>1 & iti_ideal>8 & evt_time < 9))
#   summary(ee2)
#   car::Anova(ee2)
#   anova(ee1,ee2)
#   library(emmeans)
#   em2 <- emmeans(ee2, "bin_center_z", by = c( "rewFunc","trial_neg_inv"), at = list( bin_center_z = c(-2,2), trial_neg_inv = c(-1, -.1,-.05, -.02)))
#   plot(em2, horiz = F)
#   df2 <- as.data.frame(em2)
#   
#   ee3 <- lmer(decon_interp ~ bin_center_z*evt_time*swing_above_median_lead + scale(rt_csv)*evt_time + (1 | id/run) + (1 | side), fb_comb %>% filter (iti_prev>1 & iti_ideal>8 & evt_time < 9))
#   summary(ee3)
#   car::Anova(ee3)
#   
# # + completely general time
#   ee4 <- lmer(decon_interp ~ bin_center_z*evt_time_f*swing_above_median_lead + scale(rt_csv)*evt_time_f + (1 | id/run) + (1 | side), fb_comb %>% filter (iti_prev>1 & iti_ideal>8 & evt_time < 9))
#   summary(ee4)
#   car::Anova(ee4, '3')
#   em4 <- as.data.frame(emmeans(ee4, "bin_center_z", by = c( "evt_time_f","swing_above_median_lead"), at = list( bin_center_z = c(-2,2))))
#   ggplot(em4, aes(evt_time_f, emmean, color = bin_center_z, lty = swing_above_median_lead)) + geom_point() + geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL))
#   
# # + reward
#   ee5 <- lmer(decon_interp ~ bin_center_z*evt_time_f*swing_above_median_lead*reward + scale(rt_csv)*evt_time_f + (1 | id/run) + (1 | side), fb_comb %>% filter (iti_prev>1 & iti_ideal>8 & evt_time < 9))
#   summary(ee5)
#   car::Anova(ee5, '3')
#   vif.lme(ee5)
#   anova(ee1,ee2,ee3,ee4,ee5)

#### final descriptive model -- remove swing
ee6 <- lmer(decon_interp ~ bin_center_z*evt_time_f*reward + scale(rt_csv)*evt_time_f + (1 | id/run) + (1 | side), fb_comb %>% filter (iti_prev>1 & iti_ideal>8 & evt_time < 9))
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
# Fig. 5
# 
setwd('../early_late')
pal = wes_palette("Zissou1", 24, type = "continuous")

pdf('medusa_feedback_ph_ah_reward_anderson.pdf', height = 3, width = 5)
ggplot(em6, aes(as.numeric(evt_time_f), emmean, color = bin_center_z, group = bin_center_z)) + 
  geom_point(position = position_dodge2(width = 1)) + geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),position = position_dodge2(width = .2)) + geom_line(position = position_dodge2(width = 1)) + facet_wrap(.~reward_text) +
  scale_color_gradientn(colors = pal, guide = 'none') + xlab("Time after feedback, seconds") + ylab("Hippocampal response") +
  theme(legend.title = element_blank(),
        panel.grid.major = element_line(colour = "grey45"), 
        panel.grid.minor = element_line(colour = "grey45"), 
        panel.background = element_rect(fill = 'grey40'))

dev.off()

# smoothed raw data
fb_comb <- fb_comb %>%  mutate(reward_text = case_when(
  reward == 'reward' ~ 'Reward',
  reward == 'omission' ~ 'Omission'
))

pdf('medusa_feedback_ph_ah_reward_anderson_raw_smoothed.pdf', height = 3, width = 5)
ggplot(fb_comb %>% filter (iti_prev>1 & iti_ideal>8 & evt_time < 9), aes(as.numeric(evt_time_f), decon_interp, group = bin_center_z, color = bin_center_z)) + geom_smooth(method = "gam", formula = y~splines::ns(x,3), se = F) +   
  scale_color_gradientn(colors = pal, guide = 'none') + xlab("Time after feedback, seconds") + ylab("Hippocampal response") +
  theme(legend.title = element_blank(),
        panel.grid.major = element_line(colour = "grey45"), 
        panel.grid.minor = element_line(colour = "grey45"), 
        panel.background = element_rect(fill = 'grey40')) + facet_wrap(~reward_text)

dev.off()


#models: each one has 1 3-way interaction e.g., bin_num_f*evt_time_f*scale(rt_csv)
# 2 3-way interactions: 

# replicate lm decoding analyses
# scale(-1/run_trial)*rewFunc + reward + scale(rt_csv) + scale(rt_vmax_lag) + scale(rt_vmax_change) + v_entropy_wi
fb_comb$bin_num_f <- as.factor(fb_comb$bin_num)

# just because we can
dm1 <- lmer(decon_interp ~ 
              bin_num_f*evt_time_f*scale(rt_csv) + 
              bin_num_f*evt_time_f*scale(rt_vmax_lag) +
              (1 | id/run) + (1 | side), fb_comb %>% filter (evt_time < 5))
summary(dm1)
car::Anova(dm1, '3')
vif(dm1)
library(emmeans)
r1 <- emtrends(dm1, var = 'rt_vmax_lag', specs = c('bin_num_f','evt_time_f'), data = fb_comb %>% filter (iti_prev>1 & iti_ideal>8 & evt_time < 9))
r1 <- as.data.frame(r1)
ggplot(r1, aes(evt_time_f, bin_num_f, fill = rt_vmax_lag.trend)) + 
  geom_tile() + scale_fill_viridis_c(option = "plasma")

# reduce this monstrosity to just one effect of interest
dm2 <- lmer(decon_interp ~ 
              bin_num_f*evt_time_f*scale(rt_vmax_lag) +
              (1 | id/run) + (1 | side), fb_comb %>% filter (iti_prev>1 & iti_ideal>8 & evt_time < 9))
summary(dm2)
car::Anova(dm2, '3')
r2 <- emtrends(dm2, var = 'rt_vmax_lag', specs = c('bin_num_f','evt_time_f'), data = fb_comb %>% filter (iti_prev>1 & iti_ideal>8 & evt_time < 9))
r2 <- as.data.frame(r2)
ggplot(r2, aes(evt_time_f, bin_num_f, color = rt_vmax_lag.trend)) + geom_tile()

dm3 <- lmer(decon_interp ~ 
              bin_num_f*evt_time_f*scale(v_entropy_wi_change) +
              (1 | id/run) + (1 | side), fb_comb %>% filter (iti_prev>1 & iti_ideal>8 & evt_time < 9))
summary(dm3)
car::Anova(dm3, '3')
r3 <- emtrends(dm3, var = 'rt_vmax_lag', specs = c('bin_num_f','evt_time_f'), data = fb_comb %>% filter (iti_prev>1 & iti_ideal>8 & evt_time < 9))
r3 <- as.data.frame(r3)
ggplot(r2, aes(evt_time_f, bin_num_f, color = rt_vmax_lag.trend)) + geom_tile()

# not even reward??
dm4 <- lmer(decon_interp ~ 
              bin_num_f*evt_time_f*reward + side +
              (1 | id/run) , fb_comb %>% filter (evt_time < 8))
summary(dm4)
car::Anova(dm4, '3')
em <- emmeans(dm4, )

r4 <- as.data.frame(emmeans(dm4, specs = c('bin_num_f','evt_time_f', 'reward')))
ggplot(r4, aes(evt_time_f, bin_num_f,  fill = emmean)) + geom_tile() + facet_wrap(~reward)



# # collinearity checks: it's all fine, only non-essential collinearity
# fb_comb$evt_time_f_sum <- fb_comb$evt_time_f
# contrasts(fb_comb$evt_time_f_sum) <- contr.sum(12)
# fb_comb$swing_above_median_lead_sum <- fb_comb$swing_above_median_lead
# contrasts(fb_comb$swing_above_median_lead_sum) <- contr.sum(2)
# fb_comb$reward_sum <- fb_comb$reward
# contrasts(fb_comb$reward_sum) <- contr.sum(2)
# ee5sum <- lmer(decon_interp ~ bin_center_z*evt_time_f_sum*swing_above_median_lead_sum*reward_sum + scale(rt_csv)*evt_time_f_sum + (1 | id/run) + (1 | side), fb_comb %>% filter (iti_prev>1 & iti_ideal>8 & evt_time < 9))
# summary(ee5sum)
# vif.lme(ee5sum)

em5 <- as.data.frame(emmeans(ee5, "bin_center_z", by = c( "evt_time_f","swing_above_median_lead", "reward"), at = list( bin_center_z = c(-2,2))))
ggplot(em5, aes(evt_time_f, emmean, color = bin_center_z, lty = swing_above_median_lead)) + 
  geom_point() + geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL)) + facet_wrap(.~reward) + scale_color_viridis() + theme_dark()

# + trial?
ee7 <- lmer(decon_interp ~ bin_center_z*evt_time_f*swing_above_median_lead*reward*run_trial + scale(rt_csv)*evt_time_f + (1 | id/run) + (1 | side), fb_comb %>% filter (iti_prev>1 & iti_ideal>8 & evt_time < 9))
summary(ee7)
car::Anova(ee7, '3')


ggplot(df2, aes(-1/trial_neg_inv, emmean, color = bin_center_z)) + geom_point() + geom_linerange(aes(ymin = asymp.LCL, ymax = asymp.UCL)) + facet_wrap(~rewFunc) + xlab("Trial")

#############
## VAR of PH -> AH
##############
# library(mlVAR)
# # reduce to bivariate AH PH
# # dayvar = run_trial
# # beepvar = evt_time
# # mlVAR(fb_var, c("AH", "PH"), id, lags = 1, dayvar, beepvar,
# #       estimator = c("default", "lmer", "lm","Mplus"),
# #       contemporaneous = c("default", "correlated",
# #                           "orthogonal", "fixed", "unique"), temporal =
# #         c("default", "correlated", "orthogonal", "fixed",
# #           "unique"), nCores = 1, verbose = TRUE, compareToLags,
# #       scale = TRUE, scaleWithin = FALSE, AR = FALSE,
# #       MplusSave = TRUE, MplusName = "mlVAR", iterations = "(2000)",
# #       chains = nCores, signs, orthogonal
# # )
# 
# # cor <- corr.test(fb_wide[,45:68], use = "complete")
# 
# # v0 <- mlVAR(fb_var, vars = c("AH", "PH"), idvar = "id", lags = 3, dayvar = "run_trial", beepvar = "evt_time",
# #             estimator = "lmer",
# #             contemporaneous = "correlated", temporal = "fixed",
# #             nCores = 8, verbose = TRUE, compareToLags = 2,
# #             scale = TRUE, scaleWithin = FALSE, AR = FALSE,
# #             iterations = "(2000)",
# #             chains = nCores
# # )
# 
# # just the left HIPP
# load('~/Box/SCEPTIC_fMRI/var/feedback_hipp_wide_ts.Rdata')
# # vl1 <- mlVAR(fb_wide, vars = names(fb_wide[grep('_l', names(fb_wide))]), idvar = "id", lags = 1, dayvar = "run_trial", beepvar = "evt_time",
# #             estimator = "lmer",
# #             contemporaneous = "correlated", temporal = "fixed",
# #             nCores = 8, verbose = TRUE, compareToLags = 1,
# #             scale = TRUE, scaleWithin = FALSE, AR = FALSE,
# #             iterations = "(2000)",
# #             chains = nCores
# # )
# # save(vl1,"vl1.Rdata")
# # vl1$output
# 
# vr1 <- mlVAR(fb_wide, vars = names(fb_wide[grep('_r', names(fb_wide))]), idvar = "id", lags = 1, dayvar = "run_trial", beepvar = "evt_time",
#             estimator = "lmer",
#             contemporaneous = "correlated", temporal = "fixed",
#             nCores = 8, verbose = TRUE, compareToLags = 1,
#             scale = TRUE, scaleWithin = FALSE, AR = FALSE,
#             iterations = "(2000)",
#             chains = nCores
# )
# save(vr1,"vr1.Rdata")
# layout(t(1:2))
# plot(v0, "temporal", title = "True temporal relationships", layout = "circle")
# PH and online exploration

# recode_vec <- bin_centers
# names(bin_centers) <- bin_levels
# 
# vv <- recode(clock_comb$axis_bin, !!!bin_centers)
# 
# clock_comb <- clock_comb %>% mutate(axis_cont=recode(axis_bin, !!!bin_centers))


# test the interaction with contingency
# mc1 <- lmer(decon_interp ~ rewFunc*bin_center*rt_csv + 
#               side + evt_time*rt_csv + reward * evt_time * bin_center + 
#               reward_lag * evt_time * bin_center + (1|id/run), clock_comb %>% filter(iti_prev>7 & iti_ideal>4))
# summary(mc1)
# vif.lme(mc1)
# g <- ggpredict(mc1, terms = c("rewFunc", "bin_center [.1,.3,.5,.7]", "rt_csv [1,2,3]"))
# g <- plot(g, facet = F, dodge = .4)
# g + scale_color_viridis_d(option = "plasma") + theme_dark()
# 
# # including shorter preceding ITIs does not help
# 
# # include time course
# clock_comb$evt_time_f <- as.factor(clock_comb$evt_time + 1)
# mc2 <- lmer(decon_interp ~ rewFunc*bin_center*rt_csv*evt_time_f + 
#               side + (1|id/run), clock_comb %>% filter(iti_prev>7))
# summary(mc2)
# car::Anova(mc2)
# # can't plot a four-way interaction, but let's look at the 3-way
# library(emmeans)
# g <- emmip(mc2, bin_center ~  rt_csv | rewFunc, at = list(bin_center = c(.1, .3, .5, .7), rt_csv = c(1,2,3)) )
# g + scale_color_viridis_d(option = "plasma") + theme_dark()
# 
# 
# # control for current reinforcement
# mc3 <- lmer(decon_interp ~ rewFunc*bin_center*rt_csv*evt_time_f*reward + 
#               side + (1|id/run), clock_comb %>% filter(iti_prev>7))
# summary(mc3)
# car::Anova(mc3)
# 
# # is anterior signal more auto-correlated?
# ma1 <- lmer(decon_interp ~ scale(decon_prev)*scale(bin_center) + side + (1|id/run) + (1 | side), clock_comb)
# 
# # is there a lag posterior -> anterior?
# ml1 <- lmer(decon_interp ~ scale(ah_mean_lag)*scale(bin_center) + (1|id/run) + (1|side), clock_comb %>% filter(bin_center<.2))
# summary(ml1)
# ml2 <- lmer(decon_interp ~ scale(ph_mean_lag)*scale(bin_center) + (1|id/run) + (1|side), clock_comb %>% filter(bin_center>.2))
# summary(ml2)
# anova(ml1,ml2)
# 
# g <- ggpredict(mc2, terms = c("evt_time_f", "rewFunc", "bin_center [.1,.3,.5,.7]", "rt_csv [1,2,3]"))
# g <- plot(g, facet = F, dodge = .4)
# g + scale_color_viridis_d(option = "plasma") + theme_dark()
# 
# 
# 
# # how does pre-response activity predict rt swings?
# ms1 <- lmer(decon_interp ~ decon_prev*bin_center + swing_above_median*bin_center + (1 | id/run) + (1 | side), clock_comb %>% filter(iti_prev > 7 & evt_time < 0))
# summary(ms1)
# vif.lme(ms1)
# g <- ggpredict(ms1, terms=c("swing_above_median", "bin_center", "side"), pretty = FALSE)
# plot(g, facet = TRUE)
# ggplot(g, aes(as.factor(x), predicted, color = group)) + geom_line()
# plot_model(ms1)
# 
# # try the right direction decon -> RT swing -- that does not work
# ms2 <- lmer(decon_interp ~ bin_center*evt_time*entropy_lag +  (1 | id/run), clock_comb %>% filter(iti_prev > 7))
# summary(ms2)
# vif.lme(ms2)
# g <- plot_model(ms2, show.values = TRUE)
# g + ylim(-.02,.02)
# 
# mm2 <- lmer(decon_interp ~ scale(decon_prev)*scale(iti_ideal) + (1 | id/run), clock_comb %>% filter(side=="l" & iti_ideal > 1))
# summary(mm2)
# 
# 
# 
# clock_comb <- clock_comb %>% group_by(id, run, run_trial) %>% 
#   mutate(decon_prev = dplyr::lag(decon_interp, 1, by="run_trial")) %>% ungroup()
# 

############################
# Loose plots, will organize
# clock_sum <- clock %>% group_by(id, run, evt_time, axis_bin,side) %>% summarise(mdecon_interp = mean(decon_interp))

# pdf('clock_locked_means.pdf', width = 14, height = 8)
# #ggplot(clock_sum, aes(evt_time, mdecon_interp, color = axis_bin)) + stat_smooth(method = 'gam') + facet_wrap(~side)
# ggplot(clock_sum, aes(factor(evt_time), mdecon_interp, color = axis_bin)) + 
#   stat_summary(fun.data=mean_cl_boot, geom="pointrange", position=position_dodge(width=0.4)) + facet_wrap(~side)
# dev.off()

#hist(clock_comb$iti_ideal)
# 
# clock_sum <- clock_comb %>% filter(iti_ideal < 3) %>% group_by(id, run, evt_time, axis_bin,side, reward, reward_lag, rt_above_1s) %>% summarise(mdecon_interp = mean(decon_interp))
# clock_sum_first10 <- clock_comb %>% filter(iti_ideal < 3) %>% group_by(id, run, evt_time, axis_bin,side, first10, rt_above_1s) %>% summarise(mdecon_interp = mean(decon_interp))
# 
# setwd(file.path(repo_directory, "fmri/keuka_brain_behavior_analyses/plots"))
# pdf('clock_locked_means_filter_lt_3s_iti.pdf', width = 14, height = 8)
# #ggplot(clock_sum, aes(evt_time, mdecon_interp, color = axis_bin)) + stat_smooth(method = 'gam') + facet_wrap(~side)
# cl <- ggplot(clock_sum, aes(evt_time, mdecon_interp, color = axis_bin)) + 
#   # stat_summary(fun.data=mean_cl_boot, geom="pointrange", position=position_dodge(width=0.4)) + facet_wrap(~side)
#   stat_summary(fun.y=mean, geom="line") + facet_wrap(~side)
# dev.off()
# 
# 
# fb_sum <- fb %>% group_by(id, run, evt_time, axis_bin,side) %>% summarise(mdecon_interp = mean(decon_interp))
# 
# # pdf('fb_locked_means.pdf', width = 14, height = 8)
# # #ggplot(fb_sum, aes(evt_time, mdecon_interp, color = axis_bin)) + stat_smooth(method = 'gam') + facet_wrap(~side)
# # ggplot(fb_sum, aes(factor(evt_time), mdecon_interp, color = axis_bin)) + 
# #   stat_summary(fun.data=mean_cl_boot, geom="pointrange", position=position_dodge(width=0.4)) + facet_wrap(~side)
# # dev.off()
# 
# #hist(fb_comb$iti_ideal)
# 
# fb_sum <- fb_comb %>% filter(iti_ideal < 3) %>% group_by(id, run, evt_time, axis_bin,side) %>% summarise(mdecon_interp = mean(decon_interp))
# fb_sum_first10 <- fb_comb %>% filter(iti_ideal < 3) %>% group_by(id, run, evt_time, axis_bin,side, first10) %>% summarise(mdecon_interp = mean(decon_interp))
# 
# pdf('fb_locked_means_filter_lt_3s_iti.pdf', width = 14, height = 8)
# #ggplot(fb_sum, aes(evt_time, mdecon_interp, color = axis_bin)) + stat_smooth(method = 'gam') + facet_wrap(~side)
# feed <- ggplot(fb_sum, aes(evt_time, mdecon_interp, color = axis_bin)) + 
#   # stat_summary(fun.data=mean_cl_boot, geom="pointrange", position=position_dodge(width=0.4)) + facet_wrap(~side)
#   stat_summary(fun.y=mean, geom="line") + facet_wrap(~side)
# 
# dev.off()
# 
# 
# rtvmax_sum <- rtvmax %>% group_by(id, run, evt_time, axis_bin,side) %>% summarise(mdecon_interp = mean(decon_interp))
# 
# pdf('rtvmax_locked_means.pdf', width = 14, height = 8)
# #ggplot(rtvmax_sum, aes(evt_time, mdecon_interp, color = axis_bin)) + stat_smooth(method = 'gam') + facet_wrap(~side)
# ggplot(rtvmax_sum, aes(factor(evt_time), mdecon_interp, color = axis_bin)) + 
#   stat_summary(fun.data=mean_cl_boot, geom="pointrange", position=position_dodge(width=0.4)) + facet_wrap(~side)
# dev.off()
# 
# #hist(rtvmax_comb$iti_ideal)
# 
# rtvmax_sum <- rtvmax_comb %>% filter(iti_ideal < 3) %>% group_by(id, run, evt_time, axis_bin,side) %>% summarise(mdecon_interp = mean(decon_interp))
# rtvmax_sum_vmax <- rtvmax_comb %>% filter(iti_ideal < 3) %>% group_by(id, run, evt_time, axis_bin,side, v_max_above_median) %>% summarise(mdecon_interp = mean(decon_interp))
# rtvmax_sum_vmax <- rtvmax_comb %>% filter(iti_ideal < 3) %>% group_by(id, run, evt_time, axis_bin,side, v_max_above_median, first10) %>% summarise(mdecon_interp = mean(decon_interp))
# 
# pdf('rtvmax_locked_means_filter_lt_3s_iti.pdf', width = 14, height = 8)
# #ggplot(rtvmax_sum, aes(evt_time, mdecon_interp, color = axis_bin)) + stat_smooth(method = 'gam') + facet_wrap(~side)
# # ggplot(rtvmax_sum, aes(factor(evt_time), mdecon_interp, color = axis_bin)) + 
# rtv <- ggplot(rtvmax_sum, aes(evt_time, mdecon_interp, color = axis_bin)) + 
#   # stat_summary(fun.data=mean_cl_boot, geom="pointrange", position=position_dodge(width=0.4)) + facet_wrap(~side)
#   stat_summary(fun.y=mean, geom="line") + facet_wrap(~side)
# 
# dev.off()
# 
# pdf('events_locked_means_filter_lt_3s_iti.pdf', width = 8, height = 16)
# ggarrange(cl,feed,rtv,labels = c("clock", "feedback", "RT_Vmax"), ncol = 1, nrow = 3)
# dev.off()
# # they look surprisingly similar across events
# # examine modulation by RT swing for clock, reward for feedback, Vmax for RTvmax
# clock_sum_swing <- clock_comb %>% filter(iti_ideal < 3 & !is.na(swing_above_median)) %>% group_by(id, run, evt_time, axis_bin,side, swing_above_median, rt_above_1s) %>% summarise(mdecon_interp = mean(decon_interp))
# clock_sum_swing_early <- clock_comb %>% filter(iti_ideal < 3 & !is.na(swing_above_median)) %>% group_by(id, run, evt_time, axis_bin,side, swing_above_median, first10, rt_above_1s) %>% summarise(mdecon_interp = mean(decon_interp))
# 
# # rt swing for clock-aligned
# pdf('clock_locked_by_rtswing_filter_3s_iti.pdf', width = 14, height = 8)
# #ggplot(clock_sum, aes(evt_time, mdecon_interp, color = axis_bin)) + stat_smooth(method = 'gam') + facet_wrap(~side)
# ggplot(clock_sum_swing, aes(evt_time, mdecon_interp, color = axis_bin, lty = swing_above_median)) + 
#   # stat_summary(fun.data=mean_cl_boot, geom="pointrange", position=position_dodge(width=0.4)) + facet_wrap(~side)
#   stat_summary(fun.y=mean, geom="line") + facet_wrap(~side)
# dev.off()
# 
# pdf('clock_locked_early_by_rtswing_filter_3s_iti.pdf', width = 14, height = 8)
# #ggplot(clock_sum, aes(evt_time, mdecon_interp, color = axis_bin)) + stat_smooth(method = 'gam') + facet_wrap(~side)
# ce <- ggplot(clock_sum_swing_early, aes(evt_time, mdecon_interp, color = axis_bin, lty = swing_above_median)) + 
#   # stat_summary(fun.data=mean_cl_boot, geom="pointrange", position=position_dodge(width=0.4)) + facet_wrap(~side)
#   stat_summary(fun.y=mean, geom="line") + facet_grid(first10~side)
# dev.off()
# 
# # split by current RT
# pdf('clock_locked_early_by_rt_filter_3s_iti.pdf', width = 14, height = 8)
# #ggplot(clock_sum, aes(evt_time, mdecon_interp, color = axis_bin)) + stat_smooth(method = 'gam') + facet_wrap(~side)
# ggplot(clock_sum_swing_early, aes(evt_time, mdecon_interp, color = axis_bin, lty = rt_above_1s)) + 
#   # stat_summary(fun.data=mean_cl_boot, geom="pointrange", position=position_dodge(width=0.4)) + facet_wrap(~side)
#   stat_summary(fun.y=mean, geom="line") + facet_grid(first10~side)
# dev.off()
# 
# 
# #reward for feedback-aligned
# fb_sum_rew <- fb_comb %>% filter(iti_ideal < 3) %>% group_by(id, run, evt_time, axis_bin,side, reward) %>% summarise(mdecon_interp = mean(decon_interp))
# fb_sum_rew_early <- fb_comb %>% filter(iti_ideal < 3) %>% group_by(id, run, evt_time, axis_bin,side, reward, first10,reward_lag) %>% summarise(mdecon_interp = mean(decon_interp))
# 
# pdf('fb_locked_by_reward_filter_lt_3s_iti.pdf', width = 10, height = 8)
# #ggplot(fb_sum, aes(evt_time, mdecon_interp, color = axis_bin)) + stat_smooth(method = 'gam') + facet_wrap(~side)
# ggplot(fb_sum_rew, aes(evt_time, mdecon_interp, color = axis_bin, lty = reward)) + 
#   # stat_summary(fun.data=mean_cl_boot, geom="pointrange", position=position_dodge(width=0.4)) + facet_wrap(~side)
#   stat_summary(fun.y=mean, geom="line") + facet_wrap(~side)
# dev.off()
# 
# pdf('fb_locked_early_by_reward_filter_lt_3s_iti.pdf', width = 10, height = 8)
# #ggplot(fb_sum, aes(evt_time, mdecon_interp, color = axis_bin)) + stat_smooth(method = 'gam') + facet_wrap(~side)
# fbe <- ggplot(fb_sum_rew_early, aes(evt_time, mdecon_interp, color = axis_bin, lty = reward)) + 
#   # stat_summary(fun.data=mean_cl_boot, geom="pointrange", position=position_dodge(width=0.4)) + facet_wrap(~side)
#   stat_summary(fun.y=mean, geom="line") + facet_wrap(first10~side)
# dev.off()
# 
# pdf('fb_locked_early_by_reward_lag_filter_lt_3s_iti.pdf', width = 10, height = 8)
# #ggplot(fb_sum, aes(evt_time, mdecon_interp, color = axis_bin)) + stat_smooth(method = 'gam') + facet_wrap(~side)
# ggplot(fb_sum_rew_early, aes(evt_time, mdecon_interp, color = axis_bin, lty = reward_lag)) + 
#   # stat_summary(fun.data=mean_cl_boot, geom="pointrange", position=position_dodge(width=0.4)) + facet_wrap(~side)
#   stat_summary(fun.y=mean, geom="line") + facet_wrap(~first10)
# dev.off()
# 
# pdf('modulated_clock_fb_early_late.pdf', width = 16, height = 8)
# ggarrange(ce,fbe,labels = c("clock", "feedback"), ncol = 2, nrow = 1)
# dev.off()
# 
# 
# # v_max for rt_vmax
# pdf('rtvmax_by_vmax_locked_means_filter_lt_3s_iti.pdf', width = 14, height = 8)
# #ggplot(rtvmax_sum, aes(evt_time, mdecon_interp, color = axis_bin)) + stat_smooth(method = 'gam') + facet_wrap(~side)
# # ggplot(rtvmax_sum, aes(factor(evt_time), mdecon_interp, color = axis_bin)) + 
# ggplot(rtvmax_sum_vmax, aes(evt_time, mdecon_interp, color = axis_bin, lty = v_max_above_median)) + 
#   # stat_summary(fun.data=mean_cl_boot, geom="pointrange", position=position_dodge(width=0.4)) + facet_wrap(~side)
#   stat_summary(fun.y=mean, geom="line") + facet_wrap(~side)
# 
# dev.off()
# 
# # combine plots for event type modulation by trial type
# swing <- ggplot(clock_sum_swing, aes(evt_time, mdecon_interp, color = axis_bin, lty = swing_above_median)) + 
#   # stat_summary(fun.data=mean_cl_boot, geom="pointrange", position=position_dodge(width=0.4)) + facet_wrap(~side)
#   stat_summary(fun.y=mean, geom="line") + facet_wrap(~side)
# rew <- ggplot(fb_sum_rew, aes(evt_time, mdecon_interp, color = axis_bin, lty = reward)) + 
#   # stat_summary(fun.data=mean_cl_boot, geom="pointrange", position=position_dodge(width=0.4)) + facet_wrap(~side)
#   stat_summary(fun.y=mean, geom="line") + facet_wrap(~side)
# vmax <- ggplot(rtvmax_sum_vmax, aes(evt_time, mdecon_interp, color = axis_bin, lty = v_max_above_median)) + 
#   # stat_summary(fun.data=mean_cl_boot, geom="pointrange", position=position_dodge(width=0.4)) + facet_wrap(~side)
#   stat_summary(fun.y=mean, geom="line") + facet_wrap(~side)
# pdf('modulated_events_locked_means_filter_lt_3s_iti.pdf', width = 16, height = 8)
# ggarrange(swing,rew,vmax,labels = c("clock", "feedback", "RT_Vmax"), ncol = 3, nrow = 1)
# dev.off()
