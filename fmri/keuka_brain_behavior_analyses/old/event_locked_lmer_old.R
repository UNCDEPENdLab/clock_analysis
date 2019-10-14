library(dplyr)
library(tidyverse)
library(psych)
library(ggcorrplot)
library(lme4)
library(ggpubr)

trial_df <- read_csv(file.path("~/code/clock_analysis/fmri/data/mmclock_fmri_decay_factorize_selective_psequate_mfx_trial_statistics.csv.gz")) %>%
  mutate(trial=as.numeric(trial))
trial_df <- trial_df %>%
  group_by(id, run) %>%  
  dplyr::mutate(rt_swing = abs(c(NA, diff(rt_csv))), #compute rt_swing within run and subject
                rt_swing_lr = abs(log(rt_csv/lag(rt_csv))),
                clock_onset_prev=dplyr::lag(clock_onset, 1, by="run"),
                rt_lag = lag(rt_csv) ,
                reward = score_csv>0,
                reward_lag = lag(reward),
                omission_lag = lag(score_csv==0),
                rt_vmax_lag = lag(rt_vmax),
                v_entropy_wi = scale(v_entropy),
                entropy = case_when(
                  v_entropy_wi > mean(v_entropy_wi) ~ "high",
                  v_entropy_wi < mean(v_entropy_wi) ~ "low",
                  TRUE ~ NA_character_),
                rt_change = rt_csv - rt_lag,
                rt_above_1s = rt_csv > 1000,
                swing_above_median = abs(rt_change) > median(abs(na.omit(rt_change))),
                pe_max_lag = lag(pe_max), 
                pe_max_lag2 = lag(pe_max_lag),
                pe_max_lag3 = lag(pe_max_lag2),
                abs_pe_max_lag = abs(pe_max_lag), 
                rt_vmax_change = rt_vmax - rt_vmax_lag,
                feedback_onset_prev = lag(feedback_onset),
                v_max_above_median = v_max > median(na.omit(v_max)),
                rt_bin=case_when(
                  rt_csv >= 0 & rt_csv <= 1000 ~ '0-1s',
                  rt_csv >= 1001 & rt_csv <= 2000 ~ '1-2s',
                  rt_csv >= 2001 & rt_csv <= 3000 ~ '2-3s',
                  rt_csv >= 3001 & rt_csv <= 4000 ~ '3-4s',
                  TRUE ~ NA_character_),
                run_trial=case_when(
                  trial >= 1 & trial <= 50 ~ trial,
                  trial >= 51 & trial <= 100 ~ trial - 50, #dplyr/rlang has gotten awfully picky about data types!!
                  trial >= 101 & trial <= 150 ~ trial - 100,
                  trial >= 151 & trial <= 200 ~ trial - 150,
                  trial >= 201 & trial <= 250 ~ trial - 200,
                  trial >= 251 & trial <= 300 ~ trial - 250,
                  trial >= 301 & trial <= 350 ~ trial - 300,
                  trial >= 351 & trial <= 400 ~ trial - 350,
                  TRUE ~ NA_real_),
                first10  = run_trial<11) %>% ungroup() %>%
  dplyr::mutate(rt_csv=rt_csv/1000, rt_vmax=rt_vmax/10) %>% 
  mutate(rt_vmax_cum=clock_onset + rt_vmax, rt_vmax_cum_lag = lag(rt_vmax_cum))



setwd('~/Box Sync/SCEPTIC_fMRI/')
l <- read_csv("long_axis_l_2.3mm_clock_onset_decon_locked.csv.gz") %>% mutate(side = 'l')
r <- read_csv("long_axis_r_2.3mm_clock_onset_decon_locked.csv.gz") %>% mutate(side = 'r')
clock <- rbind(l,r)
# load feedback
l <- read_csv("long_axis_l_2.3mm_feedback_onset_decon_locked.csv.gz") %>% mutate(side = 'l')
r <- read_csv("long_axis_r_2.3mm_feedback_onset_decon_locked.csv.gz") %>% mutate(side = 'r')
fb <- rbind(l,r)
# load RT(vmax)
l <- read_csv("long_axis_l_2.3mm_rt_vmax_cum_decon_locked.csv.gz") %>% mutate(side = 'l')
r <- read_csv("long_axis_r_2.3mm_rt_vmax_cum_decon_locked.csv.gz") %>% mutate(side = 'r')
rtvmax <- rbind(l,r)


clock_comb <- trial_df %>% select(id, run, run_trial, iti_ideal, score_csv, clock_onset, clock_onset_prev, 
                                  swing_above_median, first10,reward, reward_lag, rt_above_1s, rt_bin, rt_csv, entropy) %>% mutate(rewom=if_else(score_csv > 0, "rew", "om")) %>%
  group_by(id, run) %>% mutate(iti_prev=dplyr::lag(iti_ideal, by="run_trial")) %>% ungroup() %>%
  inner_join(clock)
fb_comb <- trial_df %>% select(id, run, run_trial, iti_ideal, score_csv, feedback_onset, feedback_onset_prev, 
                               reward,reward_lag, first10, rt_above_1s, rt_bin, rt_csv, entropy) %>% mutate(rewom=if_else(score_csv > 0, "rew", "om")) %>%
  group_by(id, run) %>% mutate(iti_prev=dplyr::lag(iti_ideal, by="run_trial")) %>% ungroup() %>%
  inner_join(fb)
rtvmax_comb <- trial_df %>% select(id, run, run_trial, iti_ideal, score_csv, rt_vmax_cum, rt_vmax_cum_lag, 
                                   v_max_above_median, first10, rt_bin, rt_csv, entropy) %>% mutate(rewom=if_else(score_csv > 0, "rew", "om")) %>%
  group_by(id, run) %>% mutate(iti_prev=dplyr::lag(iti_ideal, by="run_trial")) %>% ungroup() %>%
  inner_join(rtvmax)


#####################################
# Michael's proper plotting function
plot_by_summary <- function(alignment, trial_split=NULL, facet_by, filter_expr=NULL) {
  if (alignment == "clock") {
    df <- clock_comb
  } else if (alignment == "feedback") {
    df <- fb_comb
  } else if (alignment == "rtvmax") {
    df <- rtvmax_comb
  } else { stop("what is: ", alignment) }
  
  gg <- enquo(trial_split)
  df <- df %>% filter(!is.na(!!gg))

  if (!is.null(filter_expr)) {
    #fe <- enquo(filter_expr)
    #df <- df %>% filter(!!fe)
    df <- df %>% filter_(filter_expr)
  }
  
  # if (!missing(facet_by)) {
  #   fb <- enquo(facet_by)
  #   df_sum <- df %>% mutate(evt_time=evt_time+1) %>% group_by(id, run, evt_time, axis_bin, side, !!gg, !!fb) %>%
  #     summarise(mdecon_interp = mean(decon_interp)) %>% ungroup()
  # } else {
    df_sum <- df %>% mutate(evt_time=evt_time+1) %>% group_by(id, run, evt_time, axis_bin, side, !!gg) %>% #, !!fb) %>%
      summarise(mdecon_interp = mean(decon_interp)) %>% ungroup()
  # }
  
  g <- ggplot(df_sum, aes(x=evt_time, y=mdecon_interp, color = axis_bin, lty = !!gg)) +
    stat_summary(fun.y=mean, geom="line")
  
  # browser()
  # if (!missing(facet_by)) {
  #   g <- g + facet_grid(side ~ vars(facet_by))
  # } else {
    g <- g + facet_wrap(~side) + scale_colour_viridis_d() + theme_bw()
  # }
  
  return(g)
}

#######################
# Plots made with function
xx <- plot_by_summary("clock", trial_split=swing_above_median, facet_by=NULL, filter_expr=NULL)
xxx <- plot_by_summary("feedback", trial_split=reward_lag, facet_by=NULL, filter_expr=NULL)
#plot_by_summary("clock", trial_split=rt_above_1s, facet_by=NULL, filter_expr="iti_prev>1")


c1 <- plot_by_summary("clock", trial_split=swing_above_median, filter_expr = 'rt_csv<1')
c2 <- plot_by_summary("clock", trial_split=swing_above_median, filter_expr = 'rt_csv>1 & rt_csv<=2')
c3 <- plot_by_summary("clock", trial_split=swing_above_median, filter_expr = 'rt_csv>2 & rt_csv<=3')
c4 <- plot_by_summary("clock", trial_split=swing_above_median, filter_expr = 'rt_csv>3')

# pdf("clock_by_swing_rt_bin.pdf", width = 20, height = 8)
# ggarrange(c1,c2,c3,c4, ncol = 4)
# dev.off()

f1 <- plot_by_summary("feedback", trial_split=reward, filter_expr = 'rt_csv<1')
f2 <- plot_by_summary("feedback", trial_split=reward, filter_expr = 'rt_csv>1 & rt_csv<=2')
f3 <- plot_by_summary("feedback", trial_split=reward, filter_expr = 'rt_csv>2 & rt_csv<=3')
f4 <- plot_by_summary("feedback", trial_split=reward, filter_expr = 'rt_csv>3')

fe1 <- plot_by_summary("feedback", trial_split=entropy, filter_expr = 'rt_csv<1')
fe2 <- plot_by_summary("feedback", trial_split=entropy, filter_expr = 'rt_csv>1 & rt_csv<=2')
fe3 <- plot_by_summary("feedback", trial_split=entropy, filter_expr = 'rt_csv>2 & rt_csv<=3')
fe4 <- plot_by_summary("feedback", trial_split=entropy, filter_expr = 'rt_csv>3')


pdf("clock_by_swing_and_feedback_by_reward_and_entropy_by_rt_bin.pdf", width = 26, height = 21)
ggarrange(c1,c2,c3,c4,f1,f2,f3,f4,fe1,fe2,fe3,fe4,ncol = 4, nrow = 3, labels = c("RT: 0-1s", "1-2s", "2-3s","3-4s"))
dev.off()


############################
# Loose plots, will organize
# clock_sum <- clock %>% group_by(id, run, evt_time, axis_bin,side) %>% summarise(mdecon_interp = mean(decon_interp))

# pdf('clock_locked_means.pdf', width = 14, height = 8)
# #ggplot(clock_sum, aes(evt_time+1, mdecon_interp, color = axis_bin)) + stat_smooth(method = 'gam') + facet_wrap(~side)
# ggplot(clock_sum, aes(factor(evt_time), mdecon_interp, color = axis_bin)) + 
#   stat_summary(fun.data=mean_cl_boot, geom="pointrange", position=position_dodge(width=0.4)) + facet_wrap(~side)
# dev.off()

#hist(clock_comb$iti_ideal)

clock_sum <- clock_comb %>% filter(iti_ideal < 3) %>% group_by(id, run, evt_time, axis_bin,side, reward, reward_lag, rt_above_1s) %>% summarise(mdecon_interp = mean(decon_interp))
clock_sum_first10 <- clock_comb %>% filter(iti_ideal < 3) %>% group_by(id, run, evt_time, axis_bin,side, first10, rt_above_1s) %>% summarise(mdecon_interp = mean(decon_interp))

setwd('~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/plots/')
pdf('clock_locked_means_filter_lt_3s_iti.pdf', width = 14, height = 8)
#ggplot(clock_sum, aes(evt_time+1, mdecon_interp, color = axis_bin)) + stat_smooth(method = 'gam') + facet_wrap(~side)
cl <- ggplot(clock_sum, aes(evt_time+1, mdecon_interp, color = axis_bin)) + 
  # stat_summary(fun.data=mean_cl_boot, geom="pointrange", position=position_dodge(width=0.4)) + facet_wrap(~side)
  stat_summary(fun.y=mean, geom="line") + facet_wrap(~side)
dev.off()


fb_sum <- fb %>% group_by(id, run, evt_time, axis_bin,side) %>% summarise(mdecon_interp = mean(decon_interp))

# pdf('fb_locked_means.pdf', width = 14, height = 8)
# #ggplot(fb_sum, aes(evt_time+1, mdecon_interp, color = axis_bin)) + stat_smooth(method = 'gam') + facet_wrap(~side)
# ggplot(fb_sum, aes(factor(evt_time), mdecon_interp, color = axis_bin)) + 
#   stat_summary(fun.data=mean_cl_boot, geom="pointrange", position=position_dodge(width=0.4)) + facet_wrap(~side)
# dev.off()

#hist(fb_comb$iti_ideal)

fb_sum <- fb_comb %>% filter(iti_ideal < 3) %>% group_by(id, run, evt_time, axis_bin,side) %>% summarise(mdecon_interp = mean(decon_interp))
fb_sum_first10 <- fb_comb %>% filter(iti_ideal < 3) %>% group_by(id, run, evt_time, axis_bin,side, first10) %>% summarise(mdecon_interp = mean(decon_interp))

pdf('fb_locked_means_filter_lt_3s_iti.pdf', width = 14, height = 8)
#ggplot(fb_sum, aes(evt_time+1, mdecon_interp, color = axis_bin)) + stat_smooth(method = 'gam') + facet_wrap(~side)
feed <- ggplot(fb_sum, aes(evt_time+1, mdecon_interp, color = axis_bin)) + 
  # stat_summary(fun.data=mean_cl_boot, geom="pointrange", position=position_dodge(width=0.4)) + facet_wrap(~side)
  stat_summary(fun.y=mean, geom="line") + facet_wrap(~side)

dev.off()


rtvmax_sum <- rtvmax %>% group_by(id, run, evt_time, axis_bin,side) %>% summarise(mdecon_interp = mean(decon_interp))

pdf('rtvmax_locked_means.pdf', width = 14, height = 8)
#ggplot(rtvmax_sum, aes(evt_time+1, mdecon_interp, color = axis_bin)) + stat_smooth(method = 'gam') + facet_wrap(~side)
ggplot(rtvmax_sum, aes(factor(evt_time), mdecon_interp, color = axis_bin)) + 
  stat_summary(fun.data=mean_cl_boot, geom="pointrange", position=position_dodge(width=0.4)) + facet_wrap(~side)
dev.off()

#hist(rtvmax_comb$iti_ideal)

rtvmax_sum <- rtvmax_comb %>% filter(iti_ideal < 3) %>% group_by(id, run, evt_time, axis_bin,side) %>% summarise(mdecon_interp = mean(decon_interp))
rtvmax_sum_vmax <- rtvmax_comb %>% filter(iti_ideal < 3) %>% group_by(id, run, evt_time, axis_bin,side, v_max_above_median) %>% summarise(mdecon_interp = mean(decon_interp))
rtvmax_sum_vmax <- rtvmax_comb %>% filter(iti_ideal < 3) %>% group_by(id, run, evt_time, axis_bin,side, v_max_above_median, first10) %>% summarise(mdecon_interp = mean(decon_interp))

pdf('rtvmax_locked_means_filter_lt_3s_iti.pdf', width = 14, height = 8)
#ggplot(rtvmax_sum, aes(evt_time+1, mdecon_interp, color = axis_bin)) + stat_smooth(method = 'gam') + facet_wrap(~side)
# ggplot(rtvmax_sum, aes(factor(evt_time), mdecon_interp, color = axis_bin)) + 
rtv <- ggplot(rtvmax_sum, aes(evt_time+1, mdecon_interp, color = axis_bin)) + 
  # stat_summary(fun.data=mean_cl_boot, geom="pointrange", position=position_dodge(width=0.4)) + facet_wrap(~side)
  stat_summary(fun.y=mean, geom="line") + facet_wrap(~side)

dev.off()

pdf('events_locked_means_filter_lt_3s_iti.pdf', width = 8, height = 16)
ggarrange(cl,feed,rtv,labels = c("clock", "feedback", "RT_Vmax"), ncol = 1, nrow = 3)
dev.off()
# they look surprisingly similar across events
# examine modulation by RT swing for clock, reward for feedback, Vmax for RTvmax
clock_sum_swing <- clock_comb %>% filter(iti_ideal < 3 & !is.na(swing_above_median)) %>% group_by(id, run, evt_time, axis_bin,side, swing_above_median, rt_above_1s) %>% summarise(mdecon_interp = mean(decon_interp))
clock_sum_swing_early <- clock_comb %>% filter(iti_ideal < 3 & !is.na(swing_above_median)) %>% group_by(id, run, evt_time, axis_bin,side, swing_above_median, first10, rt_above_1s) %>% summarise(mdecon_interp = mean(decon_interp))

# rt swing for clock-aligned
pdf('clock_locked_by_rtswing_filter_3s_iti.pdf', width = 14, height = 8)
#ggplot(clock_sum, aes(evt_time+1, mdecon_interp, color = axis_bin)) + stat_smooth(method = 'gam') + facet_wrap(~side)
ggplot(clock_sum_swing, aes(evt_time+1, mdecon_interp, color = axis_bin, lty = swing_above_median)) + 
  # stat_summary(fun.data=mean_cl_boot, geom="pointrange", position=position_dodge(width=0.4)) + facet_wrap(~side)
  stat_summary(fun.y=mean, geom="line") + facet_wrap(~side)
dev.off()

pdf('clock_locked_early_by_rtswing_filter_3s_iti.pdf', width = 14, height = 8)
#ggplot(clock_sum, aes(evt_time+1, mdecon_interp, color = axis_bin)) + stat_smooth(method = 'gam') + facet_wrap(~side)
ce <- ggplot(clock_sum_swing_early, aes(evt_time+1, mdecon_interp, color = axis_bin, lty = swing_above_median)) + 
  # stat_summary(fun.data=mean_cl_boot, geom="pointrange", position=position_dodge(width=0.4)) + facet_wrap(~side)
  stat_summary(fun.y=mean, geom="line") + facet_grid(first10~side)
dev.off()

# split by current RT
pdf('clock_locked_early_by_rt_filter_3s_iti.pdf', width = 14, height = 8)
#ggplot(clock_sum, aes(evt_time+1, mdecon_interp, color = axis_bin)) + stat_smooth(method = 'gam') + facet_wrap(~side)
ggplot(clock_sum_swing_early, aes(evt_time+1, mdecon_interp, color = axis_bin, lty = rt_above_1s)) + 
  # stat_summary(fun.data=mean_cl_boot, geom="pointrange", position=position_dodge(width=0.4)) + facet_wrap(~side)
  stat_summary(fun.y=mean, geom="line") + facet_grid(first10~side)
dev.off()


#reward for feedback-aligned
fb_sum_rew <- fb_comb %>% filter(iti_ideal < 3) %>% group_by(id, run, evt_time, axis_bin,side, reward) %>% summarise(mdecon_interp = mean(decon_interp))
fb_sum_rew_early <- fb_comb %>% filter(iti_ideal < 3) %>% group_by(id, run, evt_time, axis_bin,side, reward, first10,reward_lag) %>% summarise(mdecon_interp = mean(decon_interp))

pdf('fb_locked_by_reward_filter_lt_3s_iti.pdf', width = 10, height = 8)
#ggplot(fb_sum, aes(evt_time+1, mdecon_interp, color = axis_bin)) + stat_smooth(method = 'gam') + facet_wrap(~side)
ggplot(fb_sum_rew, aes(evt_time+1, mdecon_interp, color = axis_bin, lty = reward)) + 
  # stat_summary(fun.data=mean_cl_boot, geom="pointrange", position=position_dodge(width=0.4)) + facet_wrap(~side)
  stat_summary(fun.y=mean, geom="line") + facet_wrap(~side)
dev.off()

pdf('fb_locked_early_by_reward_filter_lt_3s_iti.pdf', width = 10, height = 8)
#ggplot(fb_sum, aes(evt_time+1, mdecon_interp, color = axis_bin)) + stat_smooth(method = 'gam') + facet_wrap(~side)
fbe <- ggplot(fb_sum_rew_early, aes(evt_time+1, mdecon_interp, color = axis_bin, lty = reward)) + 
  # stat_summary(fun.data=mean_cl_boot, geom="pointrange", position=position_dodge(width=0.4)) + facet_wrap(~side)
  stat_summary(fun.y=mean, geom="line") + facet_wrap(first10~side)
dev.off()

pdf('fb_locked_early_by_reward_lag_filter_lt_3s_iti.pdf', width = 10, height = 8)
#ggplot(fb_sum, aes(evt_time+1, mdecon_interp, color = axis_bin)) + stat_smooth(method = 'gam') + facet_wrap(~side)
ggplot(fb_sum_rew_early, aes(evt_time+1, mdecon_interp, color = axis_bin, lty = reward_lag)) + 
  # stat_summary(fun.data=mean_cl_boot, geom="pointrange", position=position_dodge(width=0.4)) + facet_wrap(~side)
  stat_summary(fun.y=mean, geom="line") + facet_wrap(~first10)
dev.off()

pdf('modulated_clock_fb_early_late.pdf', width = 16, height = 8)
ggarrange(ce,fbe,labels = c("clock", "feedback"), ncol = 2, nrow = 1)
dev.off()


# v_max for rt_vmax
pdf('rtvmax_by_vmax_locked_means_filter_lt_3s_iti.pdf', width = 14, height = 8)
#ggplot(rtvmax_sum, aes(evt_time+1, mdecon_interp, color = axis_bin)) + stat_smooth(method = 'gam') + facet_wrap(~side)
# ggplot(rtvmax_sum, aes(factor(evt_time), mdecon_interp, color = axis_bin)) + 
ggplot(rtvmax_sum_vmax, aes(evt_time+1, mdecon_interp, color = axis_bin, lty = v_max_above_median)) + 
  # stat_summary(fun.data=mean_cl_boot, geom="pointrange", position=position_dodge(width=0.4)) + facet_wrap(~side)
  stat_summary(fun.y=mean, geom="line") + facet_wrap(~side)

dev.off()

# combine plots for event type modulation by trial type
swing <- ggplot(clock_sum_swing, aes(evt_time+1, mdecon_interp, color = axis_bin, lty = swing_above_median)) + 
  # stat_summary(fun.data=mean_cl_boot, geom="pointrange", position=position_dodge(width=0.4)) + facet_wrap(~side)
  stat_summary(fun.y=mean, geom="line") + facet_wrap(~side)
rew <- ggplot(fb_sum_rew, aes(evt_time+1, mdecon_interp, color = axis_bin, lty = reward)) + 
  # stat_summary(fun.data=mean_cl_boot, geom="pointrange", position=position_dodge(width=0.4)) + facet_wrap(~side)
  stat_summary(fun.y=mean, geom="line") + facet_wrap(~side)
vmax <- ggplot(rtvmax_sum_vmax, aes(evt_time+1, mdecon_interp, color = axis_bin, lty = v_max_above_median)) + 
  # stat_summary(fun.data=mean_cl_boot, geom="pointrange", position=position_dodge(width=0.4)) + facet_wrap(~side)
  stat_summary(fun.y=mean, geom="line") + facet_wrap(~side)
pdf('modulated_events_locked_means_filter_lt_3s_iti.pdf', width = 16, height = 8)
ggarrange(swing,rew,vmax,labels = c("clock", "feedback", "RT_Vmax"), ncol = 3, nrow = 1)
dev.off()


####### 
#models

#temporal dependency of decon estimates
clock_comb <- clock_comb %>% group_by(id, run, run_trial) %>% 
  mutate(decon_prev = dplyr::lag(decon_interp, 1, by="run_trial"), 
         telapsed=clock_onset - clock_onset_prev) %>%
  ungroup() %>%
  mutate(decon_prev_z=as.vector(scale(decon_prev)), iti_ideal_z=as.vector(scale(iti_ideal)))


fb_comb <- fb_comb %>% group_by(id, run, run_trial) %>% 
  mutate(decon_prev = dplyr::lag(decon_interp, 1, by="run_trial"), 
         telapsed=feedback_onset - feedback_onset_prev) %>%
  ungroup() %>%
  mutate(decon_prev_z=as.vector(scale(decon_prev)), iti_ideal_z=as.vector(scale(iti_ideal)))


hist(clock_comb$telapsed)
hist(fb_comb$telapsed)

 mm <- lmer(decon_interp ~ decon_prev*iti_prev + (1 | id/run), clock_comb %>% filter(side=="r" & evt_time > 0))
 summary(mm)
# 
# mm2 <- lmer(decon_interp ~ decon_prev*iti_ideal_z*evt_time + (1 | id/run), clock_comb %>% filter(side=="l" & evt_time > 0))
# summary(mm2)

mm2 <- lmer(decon_interp ~ decon_prev*telapsed*evt_time + (1 | id/run), clock_comb %>% filter(side=="l" & evt_time > 0))
summary(mm2)

# modulation by PE
m1 <- lmer(decon_interp ~ decon_prev + reward * axis_bin + entropy * axis_bin + (1 | id/run), fb_comb )
summary(m1)

clock_comb %>% filter()

#which(is.na(clock_comb$axis_bin))
clock_hack <- clock_comb %>%
  mutate(bin_low = as.numeric(sub("[^\\d]+([\\d+\\.]+),.*", "\\1", axis_bin, perl=TRUE)),
         bin_high =as.numeric(sub("[^\\d]+[\\d+\\.]+,([\\d+\\.]+)\\]", "\\1", axis_bin, perl=TRUE)))

clock_hack$bin_center <- rowMeans(clock_hack[, c("bin_low", "bin_high")])

fb_hack <- fb_comb %>%
  mutate(bin_low = as.numeric(sub("[^\\d]+([\\d+\\.]+),.*", "\\1", axis_bin, perl=TRUE)),
         bin_high =as.numeric(sub("[^\\d]+[\\d+\\.]+,([\\d+\\.]+)\\]", "\\1", axis_bin, perl=TRUE)))

fb_hack$bin_center <- rowMeans(fb_hack[, c("bin_low", "bin_high")])

# recode_vec <- bin_centers
# names(bin_centers) <- bin_levels
# 
# vv <- recode(clock_comb$axis_bin, !!!bin_centers)
# 
# clock_comb <- clock_comb %>% mutate(axis_cont=recode(axis_bin, !!!bin_centers))

# simple entropy model
m1e <- lmer(decon_interp ~ decon_prev + entropy * bin_center * evt_time + (1 | id/run) + (1 | side), fb_hack)
summary(m1e)
car::Anova(m1e)

mm2 <- lmer(decon_interp ~ decon_prev*telapsed*evt_time + (1 | id/run), clock_comb %>% filter(side=="l" & evt_time > 0))
summary(mm2)



table(clock_comb$telapsed)

head(clock_comb)

library(sjPlot)
plot_model(mm, type="pred", terms=c("decon_prev", "iti_ideal_z"))

clock_comb <- clock_comb %>% group_by(id, run, run_trial) %>% 
  mutate(decon_prev = dplyr::lag(decon_interp, 1, by="run_trial")) %>% ungroup()


