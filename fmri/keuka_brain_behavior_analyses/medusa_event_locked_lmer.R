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
# read in, process; go with "long" [-1:10] clock windows for now, will censor later
#####################
plots = F
reprocess = F
analyze = F
unsmoothed = F
newmask = F
smooth_in_mask = T

repo_directory <- "~/code/clock_analysis"
# repo_directory <- "~/Data_Analysis/clock_analysis"

#load the data, and reprocess if requested
source(file.path(repo_directory, "fmri/keuka_brain_behavior_analyses/load_medusa_data.R"))


# Michael's plotting function
#####################################
plot_by_summary <- function(alignment, trial_split=NULL, facet_by, filter_expr=NULL, change = FALSE) {
  if (alignment == "clock") {
    df <- clock_comb
  } else if (alignment == "feedback") {
    df <- fb_comb
  } else if (alignment == "rtvmax") {
    df <- rtvmax_comb
  } else if (alignment == "m1_clock") {
    df <- m1L_clock_comb
  } else if (alignment == "v1_clock") {
    df <- v1_clock_comb
  } else if (alignment == "m1_feedback") {
    df <- m1L_feedback_comb
  } else if (alignment == "v1_feedback") {
    df <- v1_feedback_comb
  } else { stop("what is: ", alignment) }
  
  # if (is.null(change)) {change = FALSE}
  
  gg <- enquo(trial_split)
  df <- df %>% filter(!is.na(!!gg))
  
  if (!is.null(filter_expr)) {
    #fe <- enquo(filter_expr)
    #df <- df %>% filter(!!fe)
    df <- df %>% filter_(filter_expr)
  }
  # if (!missing(facet_by)) {
  #   fb <- enquo(facet_by)
  #   df_sum <- df %>% mutate(evt_time=evt_time) %>% group_by(id, run, evt_time, axis_bin, side, !!gg, !!fb) %>%
  #     summarise(mdecon_interp = mean(decon_interp)) %>% ungroup()
  # } else {
  if (alignment == "clock" | alignment == "feedback" | alignment == "rtvmax") {
    if (change) {
      df_sum <- df %>% group_by(id, run, evt_time, axis_bin, side, !!gg) %>% #, !!fb) %>%
        summarise(mdecon_interp = mean(abs(decon_change))) %>% ungroup() # version with abs signal change
    } else {
      df_sum <- df %>% group_by(id, run, evt_time, axis_bin, side, !!gg) %>% #, !!fb) %>%
        summarise(mdecon_interp = mean(decon_interp)) %>% ungroup()  
    }
    # }
    g <- ggplot(df_sum, aes(x=evt_time, y=mdecon_interp, color = axis_bin, lty = !!gg)) +
      stat_summary(fun.y=mean, geom="line")
    # stat_summary(fun.data=mean_cl_boot, geom="pointrange", position=position_dodge(width=0.4))
    g <- g + facet_wrap(~side) + theme_dark()
  } else if (alignment == "v1_clock" | alignment == "v1_feedback") {
    if (change) {
      df_sum <- df %>% group_by(id, run, evt_time, side, !!gg) %>% #, !!fb) %>%
        summarise(mdecon_interp = mean(abs(decon_change))) %>% ungroup() # version with abs signal change
    } else {
      df_sum <- df %>% group_by(id, run, evt_time, side, !!gg) %>% #, !!fb) %>%
        summarise(mdecon_interp = mean(decon_interp)) %>% ungroup()  
    }
    # }
    g <- ggplot(df_sum, aes(x=evt_time, y=mdecon_interp, lty = !!gg)) +
      stat_summary(fun.y=mean, geom="line")
    g <- g + facet_wrap(~side) + theme_gray()
  } else if (alignment == "m1_clock" | alignment == 'm1_feedback') {
    if (change) {
      df_sum <- df %>% group_by(id, run, evt_time, axis_bin, !!gg) %>% #, !!fb) %>%
        summarise(mdecon_interp = mean(abs(decon_change))) %>% ungroup() # version with abs signal change
    } else {
      df_sum <- df %>% group_by(id, run, evt_time, axis_bin, !!gg) %>% #, !!fb) %>%
        summarise(mdecon_interp = mean(decon_interp)) %>% ungroup()  
    }
    # }
    g <- ggplot(df_sum, aes(x=evt_time, y=mdecon_interp, lty = !!gg)) +
      stat_summary(fun.y=mean, geom="line") + theme_gray()
  }
  # browser()
  # if (!missing(facet_by)) {
  #   g <- g + facet_grid(side ~ vars(facet_by))
  # } else {
  g <- g + scale_colour_viridis_d(option = "C") + geom_vline(xintercept = 0, linetype = 'dotted')
  # }
  
  return(g)
}

if (plots) {
  
  ##########################
  # PLOTS
  ##########################
  # Plots made with function
  # Tactics: make line plots and test significance in lmer
  
  setwd(file.path(repo_directory, "/fmri/keuka_brain_behavior_analyses/plots"))
  
  # offline activity post-feedback
  ch <- plot_by_summary("feedback", trial_split = reward, filter_expr = "iti_ideal>7 & evt_time>=0", change = FALSE)
  pdf("feedback_by_reward_across_long_axis.pdf", width = 10, height = 8)
  ch
  dev.off()
  
  # online activity during clock
  ch <- plot_by_summary("clock", trial_split = rewFunc, filter_expr = 'online=="TRUE"', change = FALSE)
  pdf("clock_online_only_by_rewFunc_across_long_axis.pdf", width = 10, height = 8)
  ch
  dev.off()
  
  
  # ramping
  ##############
  # check if ramping is robust to previous/current ITIs
  # r <- plot_by_summary("rtvmax", trial_split = reward_lag, filter_expr = 'online == "TRUE" & iti_prev>2 & iti_ideal>3', change = FALSE)
  r0 <- plot_by_summary("rtvmax", trial_split = entropy_lag, filter_expr = 'online == "TRUE"', change = FALSE)
  r1 <- plot_by_summary("rtvmax", trial_split = entropy_lag, filter_expr = 'online == "TRUE" & iti_prev>1', change = FALSE)
  r2 <- plot_by_summary("rtvmax", trial_split = entropy_lag, filter_expr = 'online == "TRUE" & iti_prev>2', change = FALSE)
  r3 <- plot_by_summary("rtvmax", trial_split = entropy_lag, filter_expr = 'online == "TRUE" & iti_prev>3', change = FALSE)
  r4 <- plot_by_summary("rtvmax", trial_split = entropy_lag, filter_expr = 'online == "TRUE" & iti_prev>4', change = FALSE)
  r5 <- plot_by_summary("rtvmax", trial_split = entropy_lag, filter_expr = 'online == "TRUE" & iti_prev>5', change = FALSE)
  pdf("rtvmax_online_by_entropy_lag_by_iti_prev_0_5s.pdf", width = 24, height = 16)
  ggarrange(r0,r1,r2,r3,r4,r5, nrow = 3, ncol = 2, labels = c("all", "iti_prev>1s", "iti_prev>2s", "iti_prev>3s", "iti_prev>4s", "iti_prev>5s"), font.label = list(color = "red"))
  dev.off()
  
  # I don't think this current ITI should matter, but let's check -- YES, but iti_ideal>2s looks best
  r20 <- plot_by_summary("rtvmax", trial_split = entropy_lag, filter_expr = 'online == "TRUE" & iti_prev>2', change = FALSE)
  r22 <- plot_by_summary("rtvmax", trial_split = entropy_lag, filter_expr = 'online == "TRUE" & iti_prev>2 & iti_ideal>2', change = FALSE)
  r24 <- plot_by_summary("rtvmax", trial_split = entropy_lag, filter_expr = 'online == "TRUE" & iti_prev>2 & iti_ideal>4', change = FALSE)
  pdf("rtvmax_online_by_entropy_lag_by_iti_ideal_0_4s.pdf", width = 24, height = 6)
  ggarrange(r20,r22,r24, nrow = 1, ncol = 3, labels = c("iti_prev>2", "iti_prev>2 & iti > 2s", "iti_prev>2 & iti > 4s"), font.label = list(color = "red"))
  dev.off()
  
  # is ramping seen across the RT range? Great -- we see it with the longest RTs!!!  Should definitely take >1s...
  r22_0 <- plot_by_summary("rtvmax", trial_split = entropy_lag, filter_expr = 'online == "TRUE" & iti_prev>2 & iti_ideal>2 & rt_csv <1', change = FALSE)
  r22_1 <- plot_by_summary("rtvmax", trial_split = entropy_lag, filter_expr = 'online == "TRUE" & iti_prev>2 & iti_ideal>2 & rt_csv >1 & rt_csv<2', change = FALSE)
  r22_2 <- plot_by_summary("rtvmax", trial_split = entropy_lag, filter_expr = 'online == "TRUE" & iti_prev>2 & iti_ideal>2 & rt_csv >2 & rt_csv<3', change = FALSE)
  r22_3 <- plot_by_summary("rtvmax", trial_split = entropy_lag, filter_expr = 'online == "TRUE" & iti_prev>2 & iti_ideal>2 & rt_csv >3 & rt_csv<4', change = FALSE)
  pdf("rtvmax_online_by_entropy_lag_by_rt.pdf", width = 18, height = 16)
  ggarrange(r22_0,r22_1,r22_2,r22_3, nrow = 2, ncol = 2, labels = c("RT < 1s", "1s > RT > 2s", "2s > RT > 3s", "3s > RT > 4s"), font.label = list(color = "red"))
  dev.off()
  
  # does it vary by contingency? Excellent: looks best in learnable, but effect of 
  ri <- plot_by_summary("rtvmax", trial_split = entropy_lag, filter_expr = 'online == "TRUE" & iti_prev>2 & iti_ideal>2 & rt_csv >1 & rewFunc=="IEV"', change = FALSE)
  rd <- plot_by_summary("rtvmax", trial_split = entropy_lag, filter_expr = 'online == "TRUE" & iti_prev>2 & iti_ideal>2 & rt_csv >1 & rewFunc=="DEV"', change = FALSE)
  rc <- plot_by_summary("rtvmax", trial_split = entropy_lag, filter_expr = 'online == "TRUE" & iti_prev>2 & iti_ideal>2 & rt_csv >1 & rewFunc=="CEV"', change = FALSE)
  rr <- plot_by_summary("rtvmax", trial_split = entropy_lag, filter_expr = 'online == "TRUE" & iti_prev>2 & iti_ideal>2 & rt_csv >1 & rewFunc=="CEVR"', change = FALSE)
  pdf("rtvmax_online_by_entropy_lag_by_rewFunc.pdf", width = 18, height = 12)
  ggarrange(rd, ri, rc,rr, nrow = 2, ncol = 2, labels = c("DEV", "IEV", "CEV", "CEVR"), font.label = list(color = "red"))
  dev.off()
  
  
  # is it the same with reward/omission?  Good: no, it is not.  But look at how CEVR is strangely reversed!
  ri <- plot_by_summary("rtvmax", trial_split = reward, filter_expr = 'online == "TRUE" & iti_prev>2 & iti_ideal>2 & rt_csv >1 & rewFunc=="IEV"', change = FALSE)
  rd <- plot_by_summary("rtvmax", trial_split = reward, filter_expr = 'online == "TRUE" & iti_prev>2 & iti_ideal>2 & rt_csv >1 & rewFunc=="DEV"', change = FALSE)
  rc <- plot_by_summary("rtvmax", trial_split = reward, filter_expr = 'online == "TRUE" & iti_prev>2 & iti_ideal>2 & rt_csv >1 & rewFunc=="CEV"', change = FALSE)
  rr <- plot_by_summary("rtvmax", trial_split = reward, filter_expr = 'online == "TRUE" & iti_prev>2 & iti_ideal>2 & rt_csv >1 & rewFunc=="CEVR"', change = FALSE)
  pdf("rtvmax_online_by_reward_by_rewFunc.pdf", width = 18, height = 12)
  ggarrange(rd, ri, rc,rr, nrow = 2, ncol = 2, labels = c("DEV", "IEV", "CEV", "CEVR"), font.label = list(color = "red"))
  dev.off()
  
  # ramps are more rampant late in learning
  re <- plot_by_summary("rtvmax", trial_split = entropy_lag, filter_expr = 'online == "TRUE" & iti_prev>2 & iti_ideal>2 & rt_csv >1 & first10 & (rewFunc=="DEV" | rewFunc=="IEV")', change = FALSE)
  rl <- plot_by_summary("rtvmax", trial_split = entropy_lag, filter_expr = 'online == "TRUE" & iti_prev>2 & iti_ideal>2 & rt_csv >1 & !first10  & (rewFunc=="DEV" | rewFunc=="IEV")', change = FALSE)
  pdf("rtvmax_online_by_entropy_by_early.pdf", width = 10, height = 10)
  ggarrange(re, rl, nrow = 2, ncol = 1, labels = c("First 10 trials", "Last 40 trials"), font.label = list(color = "red"))
  dev.off()
  
  # PH and online exploration
  # what are the predictions? Higher PH on RT swing trials
  plot_by_summary("clock", trial_split = swing_above_median, filter_expr = 'online == "TRUE" & iti_prev>2 & iti_ideal>2 & rt_csv >1 & first10 & (rewFunc=="DEV" | rewFunc=="IEV")', change = FALSE)
  plot_by_summary("clock", trial_split = swing_above_median, filter_expr = 'online == "TRUE" & iti_prev>2 & iti_ideal>2 & rt_csv >1 & !first10 & (rewFunc=="DEV" | rewFunc=="IEV")', change = FALSE)
  
  rd <- plot_by_summary("rtvmax", trial_split = reward, filter_expr = 'online == "TRUE" & iti_prev>2 & iti_ideal>2 & rt_csv >1 & rewFunc=="DEV"', change = FALSE)
  rc <- plot_by_summary("rtvmax", trial_split = reward, filter_expr = 'online == "TRUE" & iti_prev>2 & iti_ideal>2 & rt_csv >1 & rewFunc=="CEV"', change = FALSE)
  rr <- plot_by_summary("rtvmax", trial_split = reward, filter_expr = 'online == "TRUE" & iti_prev>2 & iti_ideal>2 & rt_csv >1 & rewFunc=="CEVR"', change = FALSE)
  pdf("rtvmax_online_by_reward_by_rewFunc.pdf", width = 18, height = 12)
  ggarrange(rd, ri, rc,rr, nrow = 2, ncol = 2, labels = c("DEV", "IEV", "CEV", "CEVR"), font.label = list(color = "red"))
  dev.off()
  
  
  # HIPP vs m1 vs v1
  ######################
  
  # feedback, offline processing
  # by reward
  fr <- plot_by_summary("feedback", trial_split = reward, filter_expr = 'iti_prev > 1 & iti_ideal > 8 & evt_time < 9')
  fe <- plot_by_summary("feedback", trial_split = entropy, filter_expr = 'iti_prev > 1 & iti_ideal > 8 & evt_time < 9')
  
  fp <- plot_by_summary("feedback", trial_split = abs_pe_f, filter_expr = 'iti_prev > 1 & iti_ideal > 8 & evt_time < 9')
  
  pdf("hipp_feedback_offline_rew_entropy.pdf", width = 10, height = 9)
  ggarrange(fr,fe ,ncol = 1, nrow = 2, labels = c("By last reward", "By entropy"), font.label = list(color = "red"))
  dev.off()
  
  
  pdf("hipp_feedback_offline_rew_entropy.pdf", width = 10, height = 9)
  ggarrange(fr,fe ,ncol = 1, nrow = 2, labels = c("By last reward", "By entropy"), font.label = list(color = "red"))
  dev.off()
  pdf("hipp_feedback_offline_abs_pe.pdf", width = 10, height = 5)
  fr
  dev.off()
  
  # They do show prescience that may be explained by the choice (early/late depending on contingency)
  
  h1 <- plot_by_summary("feedback", trial_split = reward, filter_expr = 'online == "FALSE" & iti_prev > 1 & iti_ideal > 8 & rt_csv<1 & evt_time < 9')
  h2 <- plot_by_summary("feedback", trial_split = reward, filter_expr = 'online == "FALSE" & iti_prev > 1 & iti_ideal > 8 & rt_csv>1 & rt_csv<=2 & evt_time < 9')
  h3 <- plot_by_summary("feedback", trial_split = reward, filter_expr = 'online == "FALSE" & iti_prev > 1 & iti_ideal > 8 & rt_csv>2 & rt_csv<=3 & evt_time < 9')
  h4 <- plot_by_summary("feedback", trial_split = reward, filter_expr = 'online == "FALSE" & iti_prev > 1 & iti_ideal > 8 & rt_csv>3 & rt_csv<=4 & evt_time < 9')
  # h1 <- h1 + geom_vline(xintercept = 0:1)
  # h2 <- h2 + geom_vline(xintercept = 1:2)
  # h3 <- h3 + geom_vline(xintercept = 2:3)
  # h4 <- h4 + geom_vline(xintercept = 3:4) 
  
  pdf("hipp_feedback_offline_by_reward.pdf", width = 16, height = 10)
  ggarrange(h1, h2, h3, h4 ,ncol = 2, nrow = 2, labels = c("HIPP <1", "HIPP 1-2", "HIPP 2-3", "HIPP 3-4"), font.label = list(color = "red"))
  dev.off()
  
  m1 <- plot_by_summary("m1_feedback", trial_split = rewFunc, filter_expr = 'iti_prev > 4 & evt_time < 6 & rt_csv<1')
  m2 <- plot_by_summary("m1_feedback", trial_split = rewFunc, filter_expr = 'iti_prev > 4 & evt_time < 6 & rt_csv>1 & rt_csv<=2')
  m3 <- plot_by_summary("m1_feedback", trial_split = rewFunc, filter_expr = 'iti_prev > 4 & evt_time < 6 & rt_csv>2 & rt_csv<=3')
  m4 <- plot_by_summary("m1_feedback", trial_split = rewFunc, filter_expr = 'iti_prev > 4 & evt_time < 6 & rt_csv>3')
  # m1 <- m1 + geom_vline(xintercept = 0:1)
  # m2 <- m2 + geom_vline(xintercept = 1:2)
  # m3 <- m3 + geom_vline(xintercept = 2:3)
  # m4 <- m4 + geom_vline(xintercept = 3:4) 
  m <- ggarrange(m1, m2, m3, m4, ncol = 2, nrow = 2, labels = c("Left m1 <1", "Left m1 1-2", "Left m1 2-3", "Left m1 3-4"), font.label = list(color = "red"))
  
  v1 <- plot_by_summary("v1_feedback", trial_split = rewFunc, filter_expr = 'iti_prev > 4 & evt_time < 6 & rt_csv<1')
  v2 <- plot_by_summary("v1_feedback", trial_split = rewFunc, filter_expr = 'iti_prev > 4 & evt_time < 6 & rt_csv>1 & rt_csv<=2')
  v3 <- plot_by_summary("v1_feedback", trial_split = rewFunc, filter_expr = 'iti_prev > 4 & evt_time < 6 & rt_csv>2 & rt_csv<=3')
  v4 <- plot_by_summary("v1_feedback", trial_split = rewFunc, filter_expr = 'iti_prev > 4 & evt_time < 6 & rt_csv>3')
  # v1 <- v1 + geom_vline(xintercept = 0:1)
  # v2 <- v2 + geom_vline(xintercept = 1:2)
  # v3 <- v3 + geom_vline(xintercept = 2:3)
  # v4 <- v4 + geom_vline(xintercept = 3:4) 
  v <- ggarrange(v1, v2, v3, v4, ncol = 2, nrow = 2, labels = c("V1 <1", "V1 1-2", "V1 2-3", "V1 3-4"), font.label = list(color = "red"))
  
  pdf("feedback_hipp_with_m1_v1_by_rt_rewFunc.pdf", width = 16, height = 16)
  ggarrange(h,v,m, ncol = 1, nrow = 3)
  dev.off()
  
  #plot_by_summary("clock", trial_split=rt_above_1s, facet_by=NULL, filter_expr="iti_prev>1")
  c1 <- plot_by_summary("clock", trial_split=swing_above_median, filter_expr = 'rt_csv<1')
  c2 <- plot_by_summary("clock", trial_split=swing_above_median, filter_expr = 'rt_csv>1 & rt_csv<=2')
  c3 <- plot_by_summary("clock", trial_split=swing_above_median, filter_expr = 'rt_csv>2 & rt_csv<=3')
  c4 <- plot_by_summary("clock", trial_split=swing_above_median, filter_expr = 'rt_csv>3')
  
  
  #clock, online
  
  h1 <- plot_by_summary("clock", trial_split = rewFunc, filter_expr = 'iti_prev > 2 & online == "TRUE" & rt_csv<1')
  h2 <- plot_by_summary("clock", trial_split = rewFunc, filter_expr = 'iti_prev > 2 & online == "TRUE" & rt_csv>1 & rt_csv<=2')
  h3 <- plot_by_summary("clock", trial_split = rewFunc, filter_expr = 'iti_prev > 2 & online == "TRUE" & rt_csv>2 & rt_csv<=3')
  h4 <- plot_by_summary("clock", trial_split = rewFunc, filter_expr = 'iti_prev > 2 & online == "TRUE" & rt_csv>3 & rt_csv<=4')
  # h1 <- h1 + geom_vline(xintercept = 0:1)
  # h2 <- h2 + geom_vline(xintercept = 1:2)
  # h3 <- h3 + geom_vline(xintercept = 2:3)
  # h4 <- h4 + geom_vline(xintercept = 3:4) 
  h <- ggarrange(h1, h2, h3, h4 ,ncol = 2, nrow = 2, labels = c("HIPP <1", "HIPP 1-2", "HIPP 2-3", "HIPP 3-4"), font.label = list(color = "red"))
  
  m1 <- plot_by_summary("m1_clock", trial_split = rewFunc, filter_expr = 'iti_prev > 2 & online == "TRUE" & rt_csv<1')
  m2 <- plot_by_summary("m1_clock", trial_split = rewFunc, filter_expr = 'iti_prev > 2 & online == "TRUE" & rt_csv>1 & rt_csv<=2')
  m3 <- plot_by_summary("m1_clock", trial_split = rewFunc, filter_expr = 'iti_prev > 2 & online == "TRUE" & rt_csv>2 & rt_csv<=3')
  m4 <- plot_by_summary("m1_clock", trial_split = rewFunc, filter_expr = 'iti_prev > 2 & online == "TRUE" & rt_csv>3')
  # m1 <- m1 + geom_vline(xintercept = 0:1)
  # m2 <- m2 + geom_vline(xintercept = 1:2)
  # m3 <- m3 + geom_vline(xintercept = 2:3)
  # m4 <- m4 + geom_vline(xintercept = 3:4) 
  m <- ggarrange(m1, m2, m3, m4, ncol = 2, nrow = 2, labels = c("Left m1 <1", "Left m1 1-2", "Left m1 2-3", "Left m1 3-4"), font.label = list(color = "red"))
  
  v1 <- plot_by_summary("v1_clock", trial_split = rewFunc, filter_expr = 'iti_prev > 2 & online == "TRUE" & rt_csv<1')
  v2 <- plot_by_summary("v1_clock", trial_split = rewFunc, filter_expr = 'iti_prev > 2 & online == "TRUE" & rt_csv>1 & rt_csv<=2')
  v3 <- plot_by_summary("v1_clock", trial_split = rewFunc, filter_expr = 'iti_prev > 2 & online == "TRUE" & rt_csv>2 & rt_csv<=3')
  v4 <- plot_by_summary("v1_clock", trial_split = rewFunc, filter_expr = 'iti_prev > 2 & online == "TRUE" & rt_csv>3')
  # v1 <- v1 + geom_vline(xintercept = 0:1)
  # v2 <- v2 + geom_vline(xintercept = 1:2)
  # v3 <- v3 + geom_vline(xintercept = 2:3)
  # v4 <- v4 + geom_vline(xintercept = 3:4) 
  v <- ggarrange(v1, v2, v3, v4, ncol = 2, nrow = 2, labels = c("V1 <1", "V1 1-2", "V1 2-3", "V1 3-4"), font.label = list(color = "red"))
  
  pdf("clock_hipp_online_with_m1_v1_by_rt_rewFunc.pdf", width = 16, height = 16)
  ggarrange(h,v,m, ncol = 1, nrow = 3)
  dev.off()
  
  
  
  # pdf("clock_by_swing_rt_bin.pdf", width = 20, height = 8)
  # ggarrange(c1,c2,c3,c4, ncol = 4)
  # dev.off()
  
  # f1 <- plot_by_summary("feedback", trial_split=reward, filter_expr = 'rt_csv<1')
  # f2 <- plot_by_summary("feedback", trial_split=reward, filter_expr = 'rt_csv>1 & rt_csv<=2')
  # f3 <- plot_by_summary("feedback", trial_split=reward, filter_expr = 'rt_csv>2 & rt_csv<=3')
  # f4 <- plot_by_summary("feedback", trial_split=reward, filter_expr = 'rt_csv>3')
  # 
  # fe1 <- plot_by_summary("feedback", trial_split=entropy, filter_expr = 'rt_csv<1')
  # fe2 <- plot_by_summary("feedback", trial_split=entropy, filter_expr = 'rt_csv>1 & rt_csv<=2')
  # fe3 <- plot_by_summary("feedback", trial_split=entropy, filter_expr = 'rt_csv>2 & rt_csv<=3')
  # fe4 <- plot_by_summary("feedback", trial_split=entropy, filter_expr = 'rt_csv>3')
  # 
  # 
  # pdf("clock_by_swing_and_feedback_by_reward_and_entropy_by_rt_bin.pdf", width = 26, height = 21)
  # ggarrange(c1,c2,c3,c4,f1,f2,f3,f4,fe1,fe2,fe3,fe4,ncol = 4, nrow = 3, labels = c("RT: 0-1s", "1-2s", "2-3s","3-4s"))
  # dev.off()
  # 
  # # effect of past reward, filtering short iti trials
  # cl1 <- plot_by_summary("clock", trial_split=reward_lag, filter_expr = 'bin_center <.80 & rt_csv<1 & iti_ideal>3 & iti_prev>7 & reward == TRUE') 
  # cl1 <- cl1 + geom_vline(xintercept = 0:1)
  # cl2 <- plot_by_summary("clock", trial_split=reward_lag, filter_expr = 'bin_center <.80 & rt_csv>1 & rt_csv<=2 & iti_ideal>3 & iti_prev>7 & reward == TRUE')
  # cl2 <- cl2 + geom_vline(xintercept = 1:2)
  # cl3 <- plot_by_summary("clock", trial_split=reward_lag, filter_expr = 'bin_center <.80 & rt_csv>2 & rt_csv<=3 & iti_ideal>3 & iti_prev>7 & reward == TRUE')
  # cl3 <- cl3 + geom_vline(xintercept = 2:3)
  # cl4 <- plot_by_summary("clock", trial_split=reward_lag, filter_expr = 'bin_center <.80 & rt_csv>3 & iti_ideal>3 & iti_prev>7 & reward == TRUE')
  # cl4 <- cl4 + geom_vline(xintercept = 3:4) 
  # pdf("clock_by_lagged_reward_by_rt_bin.pdf", width = 16, height = 10)
  # # ggarrange(cl1,cl2,cl3,cl4,ncol = 4, nrow = 1, labels = c("RT: 0-1s", "1-2s", "2-3s","3-4s"))
  # plot_grid(cl1,cl2,cl3,cl4,ncol = 2, align = 'hv',  labels = c("0-1s", "1-2s", "2-3s","3-4s"))
  # dev.off()
  # 
  # ce1 <- plot_by_summary("clock", trial_split=entropy_lag, filter_expr = 'bin_center <.80 & rt_csv<1 & iti_ideal>3 & iti_prev>7 & reward == TRUE') 
  # ce1 <- ce1 + geom_vline(xintercept = 0:1)
  # ce2 <- plot_by_summary("clock", trial_split=entropy_lag, filter_expr = 'bin_center <.80 & rt_csv>1 & rt_csv<=2 & iti_ideal>3 & iti_prev>7 & reward == TRUE')
  # ce2 <- ce2 + geom_vline(xintercept = 1:2)
  # ce3 <- plot_by_summary("clock", trial_split=entropy_lag, filter_expr = 'bin_center <.80 & rt_csv>2 & rt_csv<=3 & iti_ideal>3 & iti_prev>7 & reward == TRUE')
  # ce3 <- ce3 + geom_vline(xintercept = 2:3)
  # ce4 <- plot_by_summary("clock", trial_split=entropy_lag, filter_expr = 'bin_center <.80 & rt_csv>3 & iti_ideal>3 & iti_prev>7 & reward == TRUE')
  # ce4 <- ce4 + geom_vline(xintercept = 3:4) 
  # pdf("clock_by_lagged_entropy_by_rt_bin)_rew_only.pdf", width = 16, height = 10)
  # # ggarrange(ce1,ce2,ce3,ce4,ncol = 4, nrow = 1, labels = c("RT: 0-1s", "1-2s", "2-3s","3-4s"))
  # plot_grid(ce1,ce2,ce3,ce4,ncol = 2, align = 'hv',  labels = c("0-1s", "1-2s", "2-3s","3-4s"))
  # dev.off()
  # 
  # cs1 <- plot_by_summary("clock", trial_split=swing_above_median, filter_expr = 'bin_center <.80 & rt_csv<1 & iti_ideal>4 & iti_prev>8 & reward == TRUE') 
  # cs1 <- cs1 + geom_vline(xintercept = 0:1)
  # cs2 <- plot_by_summary("clock", trial_split=swing_above_median, filter_expr = 'bin_center <.80 & rt_csv>1 & rt_csv<=2 & iti_ideal>4 & iti_prev>8 & reward == TRUE')
  # cs2 <- cs2 + geom_vline(xintercept = 1:2)
  # cs3 <- plot_by_summary("clock", trial_split=swing_above_median, filter_expr = 'bin_center <.80 & rt_csv>2 & rt_csv<=3 & iti_ideal>4 & iti_prev>8 & reward == TRUE')
  # cs3 <- cs3 + geom_vline(xintercept = 2:3)
  # cs4 <- plot_by_summary("clock", trial_split=swing_above_median, filter_expr = 'bin_center <.80 & rt_csv>3 & iti_ideal>4 & iti_prev>8 & reward == TRUE')
  # cs4 <- cs4 + geom_vline(xintercept = 3:4) 
  # pdf("clock_by_rt_swing_by_rt_bin)_rew_only.pdf", width = 16, height = 10)
  # # ggarrange(ce1,ce2,ce3,ce4,ncol = 4, nrow = 1, labels = c("RT: 0-1s", "1-2s", "2-3s","3-4s"))
  # ggarrange(cs1,cs2,cs3,cs4,ncol = 2,nrow = 2, align = 'hv',  labels = c("0-1s", "1-2s", "2-3s","3-4s"))
  # dev.off()
  # 
  # # by contingency -- learnable, rewarded only
  # cc1 <- plot_by_summary("clock", trial_split=rewFunc, filter_expr = 'bin_center <.80 & rt_csv<1 & iti_ideal>3 & iti_prev>6 & reward == TRUE & rewFunc !="CEVR" & rewFunc !="CEV"' ) 
  # cc1 <- cc1 + geom_vline(xintercept = 0:1)
  # cc2 <- plot_by_summary("clock", trial_split=rewFunc, filter_expr = 'bin_center <.80 & rt_csv>1 & rt_csv<=2 & iti_ideal>3 & iti_prev>6 & reward == TRUE & rewFunc !="CEVR" & rewFunc !="CEV"')
  # cc2 <- cc2 + geom_vline(xintercept = 1:2)
  # cc3 <- plot_by_summary("clock", trial_split=rewFunc, filter_expr = 'bin_center <.80 & rt_csv>2 & rt_csv<=3 & iti_ideal>3 & iti_prev>6 & reward == TRUE & rewFunc !="CEVR" & rewFunc !="CEV"')
  # cc3 <- cc3 + geom_vline(xintercept = 2:3)
  # cc4 <- plot_by_summary("clock", trial_split=rewFunc, filter_expr = 'bin_center <.80 & rt_csv>3 & iti_ideal>3 & iti_prev>6 & reward == TRUE & rewFunc !="CEVR" & rewFunc !="CEV"')
  # cc4 <- cc4 + geom_vline(xintercept = 3:4) 
  # pdf("clock_by_learnable_rewFunc_by_rt_bin_rew_only.pdf", width = 16, height = 12)
  # # ggarrange(ce1,ce2,ce3,ce4,ncol = 4, nrow = 1, labels = c("RT: 0-1s", "1-2s", "2-3s","3-4s"))
  # ggarrange(cc1,cc2,cc3,cc4,ncol = 2,nrow = 2, align = 'hv',  labels = c("0-1s", "1-2s", "2-3s","3-4s"))
  # dev.off()
  # 
  # # omission only
  # occ1 <- plot_by_summary("clock", trial_split=rewFunc, filter_expr = 'bin_center <.80 & rt_csv<1 & iti_ideal>3 & iti_prev>6 & reward == FALSE & rewFunc !="CEVR" & rewFunc !="CEV"' ) 
  # occ1 <- occ1 + geom_vline(xintercept = 0:1)
  # occ2 <- plot_by_summary("clock", trial_split=rewFunc, filter_expr = 'bin_center <.80 & rt_csv>1 & rt_csv<=2 & iti_ideal>3 & iti_prev>6 & reward == FALSE & rewFunc !="CEVR" & rewFunc !="CEV"')
  # occ2 <- occ2 + geom_vline(xintercept = 1:2)
  # occ3 <- plot_by_summary("clock", trial_split=rewFunc, filter_expr = 'bin_center <.80 & rt_csv>2 & rt_csv<=3 & iti_ideal>3 & iti_prev>6 & reward == FALSE & rewFunc !="CEVR" & rewFunc !="CEV"')
  # occ3 <- occ3 + geom_vline(xintercept = 2:3)
  # occ4 <- plot_by_summary("clock", trial_split=rewFunc, filter_expr = 'bin_center <.80 & rt_csv>3 & iti_ideal>3 & iti_prev>6 & reward == FALSE & rewFunc !="CEVR" & rewFunc !="CEV"')
  # occ4 <- occ4 + geom_vline(xintercept = 3:4) 
  # pdf("clock_by_learnable_rewFunc_by_rt_bin_omission_only.pdf", width = 16, height = 12)
  # # ggarrange(ce1,ce2,ce3,ce4,ncol = 4, nrow = 1, labels = c("RT: 0-1s", "1-2s", "2-3s","3-4s"))
  # ggarrange(occ1,occ2,occ3,occ4,ncol = 2,nrow = 2, align = 'hv',  labels = c("0-1s", "1-2s", "2-3s","3-4s"))
  # dev.off()
  # 
  # # these effects should be stronger on late trials, check
  # # by contingency -- learnable, rewarded only; not enough long preceding ITIs, reduced to >4
  # lcc1 <- plot_by_summary("clock", trial_split=rewFunc, filter_expr = 'bin_center <.80 & rt_csv<1 & iti_ideal>3 & iti_prev>5 & reward == TRUE & rewFunc !="CEVR" & rewFunc !="CEV" & first10 == FALSE' ) 
  # lcc1 <- lcc1 + geom_vline(xintercept = 0:1)
  # lcc2 <- plot_by_summary("clock", trial_split=rewFunc, filter_expr = 'bin_center <.80 & rt_csv>1 & rt_csv<=2 & iti_ideal>3 & iti_prev>5 & reward == TRUE & rewFunc !="CEVR" & rewFunc !="CEV" & first10 == FALSE')
  # lcc2 <- lcc2 + geom_vline(xintercept = 1:2)
  # lcc3 <- plot_by_summary("clock", trial_split=rewFunc, filter_expr = 'bin_center <.80 & rt_csv>2 & rt_csv<=3 & iti_ideal>3 & iti_prev>5 & reward == TRUE & rewFunc !="CEVR" & rewFunc !="CEV" & first10 == FALSE')
  # lcc3 <- lcc3 + geom_vline(xintercept = 2:3)
  # lcc4 <- plot_by_summary("clock", trial_split=rewFunc, filter_expr = 'bin_center <.80 & rt_csv>3 & iti_ideal>3 & iti_prev>5 & reward == TRUE & rewFunc !="CEVR" & rewFunc !="CEV" & first10 == FALSE')
  # lcc4 <- lcc4 + geom_vline(xintercept = 3:4) 
  # pdf("clock_by_learnable_rewFunc_by_rt_bin_rew_only_last40.pdf", width = 16, height = 12)
  # # ggarrange(ce1,ce2,ce3,ce4,ncol = 4, nrow = 1, labels = c("RT: 0-1s", "1-2s", "2-3s","3-4s"))
  # ggarrange(lcc1,lcc2,lcc3,lcc4,ncol = 2,nrow = 2, align = 'hv',  labels = c("0-1s", "1-2s", "2-3s","3-4s"))
  # dev.off()
  
  # omission only
  occ1 <- plot_by_summary("clock", trial_split=rewFunc, filter_expr = 'bin_center <.80 & rt_csv<1 & iti_ideal>3 & iti_prev>6 & reward == FALSE & rewFunc !="CEVR" & rewFunc !="CEV"' ) 
  occ1 <- occ1 + geom_vline(xintercept = 0:1)
  occ2 <- plot_by_summary("clock", trial_split=rewFunc, filter_expr = 'bin_center <.80 & rt_csv>1 & rt_csv<=2 & iti_ideal>3 & iti_prev>6 & reward == FALSE & rewFunc !="CEVR" & rewFunc !="CEV"')
  occ2 <- occ2 + geom_vline(xintercept = 1:2)
  occ3 <- plot_by_summary("clock", trial_split=rewFunc, filter_expr = 'bin_center <.80 & rt_csv>2 & rt_csv<=3 & iti_ideal>3 & iti_prev>6 & reward == FALSE & rewFunc !="CEVR" & rewFunc !="CEV"')
  occ3 <- occ3 + geom_vline(xintercept = 2:3)
  occ4 <- plot_by_summary("clock", trial_split=rewFunc, filter_expr = 'bin_center <.80 & rt_csv>3 & iti_ideal>3 & iti_prev>6 & reward == FALSE & rewFunc !="CEVR" & rewFunc !="CEV"')
  occ4 <- occ4 + geom_vline(xintercept = 3:4) 
  pdf("clock_by_learnable_rewFunc_by_rt_bin_omission_only.pdf", width = 16, height = 12)
  # ggarrange(ce1,ce2,ce3,ce4,ncol = 4, nrow = 1, labels = c("RT: 0-1s", "1-2s", "2-3s","3-4s"))
  ggarrange(occ1,occ2,occ3,occ4,ncol = 2,nrow = 2, align = 'hv',  labels = c("0-1s", "1-2s", "2-3s","3-4s"))
  dev.off()
  
  # do they show prescience of the reward once you control for RT and contingency?
  devcc1 <- plot_by_summary("clock", trial_split=reward, filter_expr = 'bin_center <.80 & rt_csv<1 & iti_ideal>3 & iti_prev>6 & rewFunc =="DEV"') 
  devcc1 <- devcc1 + geom_vline(xintercept = 0:1)
  devcc2 <- plot_by_summary("clock", trial_split=reward, filter_expr = 'bin_center <.80 & rt_csv>1 & rt_csv<=2 & iti_ideal>3 & iti_prev>6  & rewFunc =="DEV"')
  devcc2 <- devcc2 + geom_vline(xintercept = 1:2)
  devcc3 <- plot_by_summary("clock", trial_split=reward, filter_expr = 'bin_center <.80 & rt_csv>2 & rt_csv<=3 & iti_ideal>3 & iti_prev>6  & rewFunc =="DEV"')
  devcc3 <- devcc3 + geom_vline(xintercept = 2:3)
  devcc4 <- plot_by_summary("clock", trial_split=reward, filter_expr = 'bin_center <.80 & rt_csv>3 & iti_ideal>3 & iti_prev>6 & rewFunc =="DEV"')
  devcc4 <- devcc4 + geom_vline(xintercept = 3:4) 
  
  pdf("clock_by_reward_by_rt_bin_DEV_only.pdf", width = 16, height = 12)
  # ggarrange(ce1,ce2,ce3,ce4,ncol = 4, nrow = 1, labels = c("RT: 0-1s", "1-2s", "2-3s","3-4s"))
  ggarrange(devcc1,devcc2,devcc3,devcc4,ncol = 2,nrow = 2, align = 'hv',  labels = c("0-1s", "1-2s", "2-3s","3-4s"))
  dev.off()
  
  ievcc1 <- plot_by_summary("clock", trial_split=reward, filter_expr = 'bin_center <.80 & rt_csv<1 & iti_ideal>3 & iti_prev>6 & rewFunc =="IEV"') 
  ievcc1 <- ievcc1 + geom_vline(xintercept = 0:1)
  ievcc2 <- plot_by_summary("clock", trial_split=reward, filter_expr = 'bin_center <.80 & rt_csv>1 & rt_csv<=2 & iti_ideal>3 & iti_prev>6  & rewFunc =="IEV"')
  ievcc2 <- ievcc2 + geom_vline(xintercept = 1:2)
  ievcc3 <- plot_by_summary("clock", trial_split=reward, filter_expr = 'bin_center <.80 & rt_csv>2 & rt_csv<=3 & iti_ideal>3 & iti_prev>6  & rewFunc =="IEV"')
  ievcc3 <- ievcc3 + geom_vline(xintercept = 2:3)
  ievcc4 <- plot_by_summary("clock", trial_split=reward, filter_expr = 'bin_center <.80 & rt_csv>3 & iti_ideal>3 & iti_prev>6 & rewFunc =="IEV"')
  ievcc4 <- ievcc4 + geom_vline(xintercept = 3:4) 
  
  pdf("clock_by_reward_by_rt_bin_iev_only.pdf", width = 16, height = 12)
  # ggarrange(ce1,ce2,ce3,ce4,ncol = 4, nrow = 1, labels = c("RT: 0-1s", "1-2s", "2-3s","3-4s"))
  ggarrange(ievcc1,ievcc2,ievcc3,ievcc4,ncol = 2,nrow = 2, align = 'hv',  labels = c("0-1s", "1-2s", "2-3s","3-4s"))
  dev.off()
  
  # by contingency -- unlearnable, reward
  ucc1 <- plot_by_summary("clock", trial_split=rewFunc, filter_expr = 'bin_center <.80 & rt_csv<1 & iti_ideal>3 & iti_prev>6 & reward == TRUE & rewFunc !="IEV" & rewFunc !="DEV"' ) 
  ucc1 <- ucc1 + geom_vline(xintercept = 0:1)
  ucc2 <- plot_by_summary("clock", trial_split=rewFunc, filter_expr = 'bin_center <.80 & rt_csv>1 & rt_csv<=2 & iti_ideal>3 & iti_prev>6 & reward == TRUE & rewFunc !="IEV" & rewFunc !="DEV"')
  ucc2 <- ucc2 + geom_vline(xintercept = 1:2)
  ucc3 <- plot_by_summary("clock", trial_split=rewFunc, filter_expr = 'bin_center <.80 & rt_csv>2 & rt_csv<=3 & iti_ideal>3 & iti_prev>6 & reward == TRUE & rewFunc !="IEV" & rewFunc !="DEV"')
  ucc3 <- ucc3 + geom_vline(xintercept = 2:3)
  ucc4 <- plot_by_summary("clock", trial_split=rewFunc, filter_expr = 'bin_center <.80 & rt_csv>3 & iti_ideal>3 & iti_prev>6 & reward == TRUE & rewFunc !="IEV" & rewFunc !="DEV"')
  ucc4 <- ucc4 + geom_vline(xintercept = 3:4) 
  pdf("clock_by_unlearnable_rewFunc_by_rt_bin_rew_only.pdf", width = 16, height = 12)
  # ggarrange(ce1,ce2,ce3,ce4,ncol = 4, nrow = 1, labels = c("RT: 0-1s", "1-2s", "2-3s","3-4s"))
  ggarrange(ucc1,ucc2,ucc3,ucc4,ncol = 2,nrow = 2, align = 'hv',  labels = c("0-1s", "1-2s", "2-3s","3-4s"))
  dev.off()
  
  # greater autocorrelation in anteror -- SIC: run abs change version
  ar1 <- plot_by_summary("clock", trial_split = reward, filter_expr = 'bin_center <.80')
  ar1 <- ar1 + ylab('Absolute signal change from last trial')
  pdf("clock_abs_signal_change_by_bin_center.pdf", height = 8, width = 10)
  ar1
  dev.off()
  
}


# LMER
####### 
#LMER models
#########

vif.lme <- function (fit) {
  ## adapted from rms::vif
  v <- vcov(fit)
  nam <- names(fixef(fit))
  ## exclude intercepts
  ns <- sum(1 * (nam == "Intercept" | nam == "(Intercept)"))
  if (ns > 0) {
    v <- v[-(1:ns), -(1:ns), drop = FALSE]
    nam <- nam[-(1:ns)] }
  d <- diag(v)^0.5
  v <- diag(solve(v/(d %o% d)))
  names(v) <- nam
  v }


#temporal dependency of decon estimates
#####
# hist(clock_comb$telapsed)
# hist(fb_comb$telapsed)


if (analyze) {
  #   mm <- lmer(decon_interp ~ decon_prev*iti_prev + (1 | id/run), clock_comb %>% filter(evt_time == 1))
  #   summary(mm)
  #   mm2 <- lmer(decon_interp ~ decon_prev*telapsed*evt_time + (1 | id/run), clock_comb %>% filter(side=="l" & evt_time > 0))
  #   summary(mm2)
  #   m1 <- lmer(decon_interp ~ (scale(decon_prev) + scale(bin_center) + evt_time)^3 + (1 | id/run) + (1 | side), fb_comb %>% filter(evt_time>0))
  #   summary(m1)
  #   
  # # does entropy differentially modulate signal autocorrelation along the axis?
  #   cm2 <- lmer(decon_interp ~ (decon_prev_z + bin_center_z + online + entropy)^3 + (1 | id/run) + (1 | side), clock_comb %>% filter (iti_prev>4 & iti_ideal >2))
  #   summary(cm2)
  #   vif.lme(cm2)
  #   g <- ggpredict(cm2, terms = c("entropy_lag", "bin_center_z [-2,0,2]"))
  #   g <- ggpredict(cm2, terms = c("online", "bin_center_z [-2,0,2]"))
  # # g <- ggpredict(cm2, terms = c("entropy_lag", "decon_prev [.1,.5,.9]"))
  #   g <- plot(g, facet = F, dodge = .4)
  #   g + scale_color_viridis_d(option = "plasma") + theme_dark()
  #   
  #   fm2 <- lmer(decon_interp ~ (decon_prev_z + bin_center_z + online)^2 + (1 | id/run) + (1 | side), clock_comb %>% filter (iti_prev>4 & iti_ideal >3))
  #   summary(fm2)
  #   vif.lme(fm2)
  # # g <- ggpredict(fm2, terms = c("entropy", "bin_center [-2,0,2]"))
  #   g <- ggpredict(fm2, terms = c("online", "bin_center_z [-2,0,2]"))
  #   g <- plot(g, facet = F, dodge = .2)
  #   pdf("ah_vs_ph_online_lmer.pdf", width = 6, height = 6)
  #   g + scale_color_viridis_d(option = "plasma") + theme_dark()
  #   dev.off()
  
  
  # test for ramps in AH in rtvmax-aligned data: quadratic term
  #######
  # Ramps figures
  
    # filter by ITI, RT, and evt_time <3
  setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/hippo/figs/ramps/')
  rvdf <- rtvmax_comb %>% filter(online == "TRUE" & iti_prev>1 & iti_ideal > 2 & rt_csv > 1 & rewFunc!="CEVR" & evt_time < 3)
  rvdf$bin_num <- as.factor(rvdf$bin_num)
  rvdf <- rvdf %>% mutate(`Hippocampal response` = decon_interp, entropy = case_when(
    entropy_lag == 'high' ~ 'High entropy',
    entropy_lag == 'low' ~ 'Low entropy'
  ))
  
  # strangely we only see ramps on long-ITI trials  
  # rm* models do not converge with cobra percentchange
  rm1 <- lmer(decon_interp ~ evt_time*bin_center_z + evt_time_sq*bin_center_z + entropy_lag + reward_lag + (1 | id/run), rvdf)
  summary(rm1)
  Anova(rm1, '3')
  vif.lme(rm1)
  
  # add entropy modulation
  # no longer a 3-way interaction on Cobra-masked data
  rm2 <- lmer(decon_interp ~ (evt_time + bin_center_z + entropy_lag) ^2 + (evt_time_sq + bin_center_z + entropy_lag) ^2 + reward_lag + scale(rt_csv) + (1 | id/run), rvdf)
  summary(rm2)
  vif.lme(rm2)
  Anova(rm2, '3')
  anova(rm1,rm2,rm2f)
  
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
  pdf("ramps_in_AH_lin_quad_cobra_demean.pdf", width = 6, height = 3)
  ggplot(em2, aes(time, `Hippocampal response`, color = as.factor(bin_center_z))) + 
    geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL)) + geom_line(size = 1.5) + scale_color_viridis_d() + theme_dark() + facet_wrap(~entropy) + theme(legend.position = "none") +
    geom_vline(xintercept = 0, lty = 'dashed', color = 'red', size = 1.5) + xlab('Time') + scale_x_continuous(breaks = c(-2,-1,0,1,2)) + ylab('Hippocampal response')
  dev.off()
  
  #wesanderson version 
  library(wesanderson)
  pal = wes_palette("Zissou1", 12, type = "discrete")
  pdf("ramps_in_AH_lin_quad_cobra_demean_anderson.pdf", width = 6, height = 3)
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
  vif.lme(rm2)
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
  
  
  
  # # make bin a factor
  #   # 3-way interaction is NS with any decon
  #   rm2ff <- lmer(decon_interp ~ (evt_time_f + bin_num + entropy_lag) ^2 + reward_lag + scale(rt_csv)*evt_time_f + (1 | id/run), rvdf)
  #   summary(rm2ff)
  #   vif.lme(rm2f)
  #   Anova(rm2ff, '3')
  #   em2ff <- as.data.frame(emmeans(rm2ff,specs = c("evt_time_f", "bin_num", "entropy_lag")))
  #   em2ff$hipp_response <- em2ff$emmean
  #   anova(rm2f,rm2ff)
  #   pdf("ramps_in_AH_ff.pdf", width = 8, height = 6)
  #   ggplot(em2ff, aes(evt_time_f, hipp_response, group = bin_num, color = bin_num)) + geom_point() + 
  #     geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL)) + geom_line() + facet_wrap(~entropy_lag) + scale_color_viridis_d() + theme_dark()
  #   dev.off()
  
  # also plot smoothed raw data
  rvdf <- rvdf %>% mutate(entropy = case_when(
    entropy_lag == 'high' ~ 'High entropy',
    entropy_lag == 'low' ~ 'Low entropy'))
  pdf('smoothed_ramps_cobra_percent.pdf', width = 6, height = 3)
  # ggplot(rvdf[!is.na(rvdf$entropy_lag),], aes(evt_time, decon_interp, color = bin_num)) + geom_smooth(method = "loess", se = F) + scale_color_viridis_d() + theme_dark() + facet_wrap(~entropy_lag)
  ggplot() + stat_smooth(data = rvdf[!is.na(rvdf$entropy_lag),], aes(evt_time, `Hippocampal response`, color = as.factor(bin_center_z)), geom = 'line', method = "loess", se = F)  + 
    scale_color_viridis_d() + theme_dark() + facet_wrap(~entropy) + theme(legend.position = "none") + geom_vline(xintercept = 0, lty = 'dashed', color = 'red', size = 1.5) + xlab('Time') + 
    scale_x_continuous(breaks = c(-2,0,2)) + ylab('Hippocampal response')
  dev.off()
  
  # wesanderson
  pdf('smoothed_ramps_cobra_percent_anderson.pdf', width = 5, height = 3)
  # ggplot(rvdf[!is.na(rvdf$entropy_lag),], aes(evt_time, decon_interp, color = bin_num)) + geom_smooth(method = "loess", se = F) + scale_color_viridis_d() + theme_dark() + facet_wrap(~entropy_lag)
  ggplot() + stat_smooth(data = rvdf[!is.na(rvdf$entropy_lag),], aes(evt_time, `Hippocampal response`, color = bin_center_z, group = bin_center_z), geom = 'line', method = "loess", se = F)  + 
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
}
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
