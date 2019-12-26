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
source('~/code/Rhelpers/screen.lmerTest.R')
source('~/code/Rhelpers/vif.lme.R')

#####################
reprocess = F
unsmoothed = F
newmask = T
smooth_in_mask = F

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
pal = wes_palette("Zissou1", 12, type = "continuous")
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
# Fig. 3A
setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/hippo/figs/final/fig_3/')
pdf('3a_smoothed_ramps_cobra_percent_anderson.pdf', width = 4, height = 2.5)

# this bit makes supplemental figure with COBRA mask
# setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/hippo/figs/final/supp')
# pdf('supp3a_smoothed_ramps_cobra_percent_anderson.pdf', width = 4, height = 2.5)
# ggplot(rvdf[!is.na(rvdf$entropy_lag),], aes(evt_time, decon_interp, color = bin_num)) + geom_smooth(method = "loess", se = F) + scale_color_viridis_d() + theme_dark() + facet_wrap(~entropy_lag)
ggplot(rvdf[!is.na(rvdf$entropy_lag),], aes(evt_time, `Hippocampal response`, color = as.numeric(bin_center_z), group = as.numeric(bin_center_z))) + 
  geom_smooth(method = 'gam', formula = y ~ splines::ns(x,3), se = F, size = .75)  + 
  facet_wrap(~entropy) + theme(legend.position = "none") + geom_vline(xintercept = 0, lty = 'dashed', color = 'white', size = 1) + xlab('Time') + 
  scale_x_continuous(breaks = c(-2,0,2)) + ylab('Hippocampal response') + 
  scale_color_gradientn(colors = pal, guide = 'none') +   scale_color_gradientn(colors = pal, guide = 'none') + 
  geom_text(aes(x=-.5, y = .485, label = "RT(Vmax)"), angle = 90, color = "white", size = 2) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_line(colour = "grey45"), 
        panel.grid.minor = element_line(colour = "grey45"), 
        panel.background = element_rect(fill = 'grey40'))
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
# pdf("ramps_in_AH_f_cobra.pdf", width = 6, height = 3)
# ggplot(em2f, aes(evt_time, `Hippocampal response`, color = as.factor(bin_center_z))) + 
#   geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL)) + geom_line(size = 1.5) + scale_color_viridis_d() + theme_dark() + facet_wrap(~entropy) + theme(legend.position = "none") +
#   geom_vline(xintercept = 0, lty = 'dashed', color = 'red', size = 1.5)+ xlab('Time') + ylab('Hippocampal response')
# dev.off()
# anova(rm1,rm2,rm2f, rm2binf)
# 
# Fig. 3B
# setwd('../final/fig_3/')

# alternatively for supplemental figure 3:
# setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/hippo/figs/final/supp')
# pdf("supp_3b_ramps_in_AH_f_smoothed_in_mask_anderson.pdf", width = 4, height = 2.5)
setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/hippo/figs/final/fig_3/')
pdf("3b_ramps_in_AH_f_smoothed_in_mask_anderson.pdf", width = 4, height = 2.5)
label1 <- expression("RT[Vmax]")
ggplot(em2f, aes(evt_time, `Hippocampal response`, color = bin_center_z, group = bin_center_z)) + 
  # geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),position = position_dodge(width = .5), size = .5) + 
  geom_line(position = position_dodge(width = .5), size = 1) +  facet_wrap(~entropy) + theme(legend.position = "none") +
  geom_linerange(aes(ymin = emmean - SE, ymax = emmean + SE),position = position_dodge(width = .5),  color = "grey80") + 
  geom_vline(xintercept = 0, lty = 'dashed', color = 'white', size = 1)+ xlab('Time, seconds') + ylab('Hippocampal response') +
  scale_color_gradientn(colors = pal, guide = 'none') + geom_text(aes(x=-.5, y = .485, label = "RT(Vmax)"), angle = 90, color = "white", size = 2) +
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

# pdf('smoothed_ramps_cobra_percent.pdf', width = 6, height = 3)
# # ggplot(rvdf[!is.na(rvdf$entropy_lag),], aes(evt_time, decon_interp, color = bin_num)) + geom_smooth(method = "loess", se = F) + scale_color_viridis_d() + theme_dark() + facet_wrap(~entropy_lag)
# ggplot() + stat_smooth(data = rvdf[!is.na(rvdf$entropy_lag),], aes(evt_time, `Hippocampal response`, color = as.factor(bin_center_z)), geom = 'line', method = "loess", se = T)  + 
#   scale_color_viridis_d() + theme_dark() + facet_wrap(~entropy) + theme(legend.position = "none") + geom_vline(xintercept = 0, lty = 'dashed', color = 'red', size = 1.5) + xlab('Time') + 
#   scale_x_continuous(breaks = c(-2,0,2)) + ylab('Hippocampal response')
# dev.off()

# 
# 
# # try and combine two plots
# pdf("ramps_in_AH_combined_cobra_percent.pdf", width = 6, height = 3)
# ggplot(em2f, aes(evt_time, `Hippocampal response`, group = bin_center_z, color = bin_center_z)) + geom_point() + geom_line() +
#   geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL)) + scale_color_viridis_c() + theme_dark() + facet_wrap(~entropy) +
#   stat_smooth(data = rvdf[!is.na(rvdf$entropy_lag),], aes(evt_time, decon_interp, group = bin_center_z, color = bin_center_z), geom = 'line', alpha = .2, method = "loess", se = F)  + theme(legend.position = "none") + ylab('Hippocampal response')
# dev.off()
# 
# # RTvmax-aligned
# pal = wes_palette("Zissou1", 24, type = "continuous")
# pdf("trial_rtvmax_hipp_AH_PH_gam.pdf", width = 11, height = 8)
# # ggplot(rtvmax_comb,aes(run_trial,decon_interp, color = axis_bin, lty = reward)) + geom_smooth(method = "gam", formula = y ~ splines::ns(x,3),  se = F) + scale_color_viridis_d() + theme_dark()
# ggplot(rtvmax_comb,aes(run_trial,decon_interp, color = axis_bin)) + geom_smooth(method = "gam", formula = y~splines::ns(x,4)) + 
#   scale_color_gradientn(colors = pal, guide = 'none') + 
#   theme(legend.title = element_blank(),
#         panel.grid.major = element_line(colour = "grey45"), 
#         panel.grid.minor = element_line(colour = "grey45"), 
#         panel.background = element_rect(fill = 'grey40'))
# 
# dev.off()


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
# Fig. 4D
# pdf("../fig_4/4d_trial_hipp_AH_PH_bin6_f.pdf", width = 3, height = 3)
pdf("../supp/supp4d_trial_hipp_AH_PH_bin6_f.pdf", width = 3, height = 3)

# ggplot(rtvmax_comb,aes(run_trial,decon_interp, color = axis_bin, lty = reward)) + geom_smooth(method = "gam", formula = y ~ splines::ns(x,3),  se = F) + scale_color_viridis_d() + theme_dark()
ggplot(emt, aes(epoch, `Hippocampal response`, color = as.numeric(bin6_f), group = as.numeric(bin6_f))) + 
  # geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),position = position_dodge(width = .5), size = .5) + geom_line(size = 1.5,position = position_dodge(width = .5)) +  theme(legend.position = "none") +
  # geom_linerange(aes(ymin = asymp.LCL, ymax = asymp.UCL),position = position_dodge(width = .5), size = .5, color = "grey80") +  geom_line(size = 1.5,position = position_dodge(width = .5)) +  theme(legend.position = "none") +
  geom_linerange(aes(ymin = emmean - SE, ymax = emmean + SE),position = position_dodge(width = .5), size = .5, color = "grey80")  +  geom_line(size = 1,position = position_dodge(width = .5)) +  theme(legend.position = "none") +
    xlab('Trial') + ylab('Hippocampal response') +
  scale_color_gradientn(colors = pal, guide = 'none') + 
  theme(legend.title = element_blank(),
        panel.grid.major = element_line(colour = "grey45"), 
        panel.grid.minor = element_line(colour = "grey45"), 
        panel.background = element_rect(fill = 'grey40'))
dev.off()


# post-hoc for trial * bin
em <- emmeans(tm1,specs = c("trial_f", "bin6_f"))
pwpp(em, sort = FALSE, method = "trt.vs.ctrl1", type = "response", side = ">")
CLD = emmeans::CLD(em,
          alpha=0.05,
          Letters=letters,
          adjust="tukey")

# Fig. 4C
# pdf("../fig_4/4c_fb_hipp_AP_trial_anderson.pdf", width = 3, height = 3)
pdf("../supp/supp4c_fb_hipp_AP_trial_anderson.pdf", width = 3, height = 3)

ggplot(fb_comb %>% filter (iti_prev>1 & iti_ideal>8 & evt_time < 9), aes(run_trial, decon_interp, color = bin6, group = bin6)) + 
  geom_smooth(method = "gam", formula = y~splines::ns(x,3),  se = T, size = 1) + 
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
# Fig. 4B (post-feedback MEDUSA)
# 
setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/hippo/figs/final/fig_4')
pal = wes_palette("Zissou1", 24, type = "continuous")
pdf('4b_ medusa_feedback_ph_ah_reward_anderson.pdf', height = 2.5, width = 4)
ggplot(em6, aes(as.numeric(evt_time_f), emmean, color = bin_center_z, group = bin_center_z)) + 
  # geom_point(position = position_dodge2(width = 1)) + geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),position = position_dodge2(width = .2)) + geom_line(position = position_dodge2(width = 1)) + facet_wrap(.~reward_text) +
  geom_linerange(aes(ymin = emmean - SE, ymax = emmean + SE),position = position_dodge2(width = 1), color = 'grey80') + geom_line(position = position_dodge2(width = 1), size = 1) + facet_wrap(.~reward_text) +
  
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

pdf('4a_medusa_feedback_ph_ah_reward_anderson_raw_smoothed.pdf', height = 3, width = 5)
ggplot(fb_comb %>% filter (iti_prev>1 & iti_ideal>8 & evt_time < 9), aes(as.numeric(evt_time_f), decon_interp, group = bin_center_z, color = bin_center_z)) + 
  geom_smooth(method = "gam", formula = y~splines::ns(x,3), se = F, size = .75) +   
  scale_color_gradientn(colors = pal, guide = 'none') + xlab("Time after feedback, seconds") + ylab("Hippocampal response") +
  theme(legend.title = element_blank(),
        panel.grid.major = element_line(colour = "grey45"), 
        panel.grid.minor = element_line(colour = "grey45"), 
        panel.background = element_rect(fill = 'grey40')) + facet_wrap(~reward_text)

dev.off()


# #models: each one has 1 3-way interaction e.g., bin_num_f*evt_time_f*scale(rt_csv)
# # 2 3-way interactions: 
# 
# # replicate lm decoding analyses
# # scale(-1/run_trial)*rewFunc + reward + scale(rt_csv) + scale(rt_vmax_lag) + scale(rt_vmax_change) + v_entropy_wi
# fb_comb$bin_num_f <- as.factor(fb_comb$bin_num)
# 
# # just because we can
# dm1 <- lmer(decon_interp ~ 
#               bin_num_f*evt_time_f*scale(rt_csv) + 
#               bin_num_f*evt_time_f*scale(rt_vmax_lag) +
#               (1 | id/run) + (1 | side), fb_comb %>% filter (evt_time < 5))
# summary(dm1)
# car::Anova(dm1, '3')
# vif(dm1)
# library(emmeans)
# r1 <- emtrends(dm1, var = 'rt_vmax_lag', specs = c('bin_num_f','evt_time_f'), data = fb_comb %>% filter (iti_prev>1 & iti_ideal>8 & evt_time < 9))
# r1 <- as.data.frame(r1)
# ggplot(r1, aes(evt_time_f, bin_num_f, fill = rt_vmax_lag.trend)) + 
#   geom_tile() + scale_fill_viridis_c(option = "plasma")
# 
# # reduce this monstrosity to just one effect of interest
# dm2 <- lmer(decon_interp ~ 
#               bin_num_f*evt_time_f*scale(rt_vmax_lag) +
#               (1 | id/run) + (1 | side), fb_comb %>% filter (iti_prev>1 & iti_ideal>8 & evt_time < 9))
# summary(dm2)
# car::Anova(dm2, '3')
# r2 <- emtrends(dm2, var = 'rt_vmax_lag', specs = c('bin_num_f','evt_time_f'), data = fb_comb %>% filter (iti_prev>1 & iti_ideal>8 & evt_time < 9))
# r2 <- as.data.frame(r2)
# ggplot(r2, aes(evt_time_f, bin_num_f, color = rt_vmax_lag.trend)) + geom_tile()
# 
# dm3 <- lmer(decon_interp ~ 
#               bin_num_f*evt_time_f*scale(v_entropy_wi_change) +
#               (1 | id/run) + (1 | side), fb_comb %>% filter (iti_prev>1 & iti_ideal>8 & evt_time < 9))
# summary(dm3)
# car::Anova(dm3, '3')
# r3 <- emtrends(dm3, var = 'rt_vmax_lag', specs = c('bin_num_f','evt_time_f'), data = fb_comb %>% filter (iti_prev>1 & iti_ideal>8 & evt_time < 9))
# r3 <- as.data.frame(r3)
# ggplot(r2, aes(evt_time_f, bin_num_f, color = rt_vmax_lag.trend)) + geom_tile()
# 
# # not even reward??
# dm4 <- lmer(decon_interp ~ 
#               bin_num_f*evt_time_f*reward + side +
#               (1 | id/run) , fb_comb %>% filter (evt_time < 8))
# summary(dm4)
# car::Anova(dm4, '3')
# em <- emmeans(dm4, )
# 
# r4 <- as.data.frame(emmeans(dm4, specs = c('bin_num_f','evt_time_f', 'reward')))
# ggplot(r4, aes(evt_time_f, bin_num_f,  fill = emmean)) + geom_tile() + facet_wrap(~reward)
# 


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
# 
# em5 <- as.data.frame(emmeans(ee5, "bin_center_z", by = c( "evt_time_f","swing_above_median_lead", "reward"), at = list( bin_center_z = c(-2,2))))
# ggplot(em5, aes(evt_time_f, emmean, color = bin_center_z, lty = swing_above_median_lead)) + 
#   geom_point() + geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL)) + facet_wrap(.~reward) + scale_color_viridis() + theme_dark()
# 
# # + trial?
# ee7 <- lmer(decon_interp ~ bin_center_z*evt_time_f*swing_above_median_lead*reward*run_trial + scale(rt_csv)*evt_time_f + (1 | id/run) + (1 | side), fb_comb %>% filter (iti_prev>1 & iti_ideal>8 & evt_time < 9))
# summary(ee7)
# car::Anova(ee7, '3')
# 
# 
# ggplot(df2, aes(-1/trial_neg_inv, emmean, color = bin_center_z)) + geom_point() + geom_linerange(aes(ymin = asymp.LCL, ymax = asymp.UCL)) + facet_wrap(~rewFunc) + xlab("Trial")
