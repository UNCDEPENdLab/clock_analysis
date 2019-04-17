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

#####################
# read in, process

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


trial_df <- read_csv(file.path("~/code/clock_analysis/fmri/data/mmclock_fmri_decay_factorize_selective_psequate_mfx_trial_statistics.csv.gz")) %>%
  mutate(trial=as.numeric(trial))
trial_df <- trial_df %>%
  group_by(id, run) %>%  
  dplyr::mutate(rt_swing = abs(c(NA, diff(rt_csv))), #compute rt_swing within run and subject
                rt_swing_lr = abs(log(rt_csv/lag(rt_csv))),
                clock_onset_prev=dplyr::lag(clock_onset, 1, by="run"),
                rt_lag = lag(rt_csv) ,
                reward = as.factor(score_csv>0),
                reward_lag = lag(reward),
                omission_lag = lag(score_csv==0),
                rt_vmax_lag = lag(rt_vmax),
                v_entropy_wi = scale(v_entropy),
                entropy = case_when(
                  v_entropy_wi > mean(v_entropy_wi) ~ "high",
                  v_entropy_wi < mean(v_entropy_wi) ~ "low",
                  TRUE ~ NA_character_),
                entropy_lag = lag(entropy),
                rt_change = rt_csv - rt_lag,
                rt_above_1s = rt_csv > 1000,
                swing_above_median = as.factor(abs(rt_change) > median(abs(na.omit(rt_change)))),
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
                  rewFunc <- as.factor(rewFunc),
                  first10  = run_trial<11) %>% ungroup() %>%
  dplyr::mutate(rt_csv=rt_csv/1000, rt_vmax=rt_vmax/10) %>% 
  mutate(rt_vmax_cum=clock_onset + rt_vmax, rt_vmax_cum_lag = lag(rt_vmax_cum))

trial_df$rewFunc <- as.factor(trial_df$rewFunc)
levels(trial_df$rewFunc) <- c("DEV", "IEV", "CEV", "CEVR")


clock_comb <- trial_df %>% select(id, run, run_trial, iti_ideal, score_csv, clock_onset, clock_onset_prev, rt_lag, rewFunc,
                                  swing_above_median, first10,reward, reward_lag, rt_above_1s, rt_bin, rt_csv, entropy, entropy_lag) %>% mutate(rewom=if_else(score_csv > 0, "rew", "om")) %>%
  group_by(id, run) %>% mutate(iti_prev=dplyr::lag(iti_ideal, by="run_trial")) %>% ungroup() %>%
  inner_join(clock)
fb_comb <- trial_df %>% select(id, run, run_trial, iti_ideal, score_csv, feedback_onset, feedback_onset_prev,  rt_lag,rewFunc,
                               reward,reward_lag, first10, rt_above_1s, rt_bin, rt_csv, entropy, entropy_lag) %>% mutate(rewom=if_else(score_csv > 0, "rew", "om")) %>%
  group_by(id, run) %>% mutate(iti_prev=dplyr::lag(iti_ideal, by="run_trial")) %>% ungroup() %>%
  inner_join(fb)
rtvmax_comb <- trial_df %>% select(id, run, run_trial, iti_ideal, score_csv, rt_vmax_cum, rt_vmax_cum_lag,  rt_lag,rewFunc,
                                   v_max_above_median, first10, rt_bin, rt_csv, entropy, entropy_lag) %>% mutate(rewom=if_else(score_csv > 0, "rew", "om")) %>%
  group_by(id, run) %>% mutate(iti_prev=dplyr::lag(iti_ideal, by="run_trial")) %>% ungroup() %>%
  inner_join(rtvmax)

# 20% of clock- and 32% of feedback-aligned timepoints are from the next trial: censor
clock_comb$decon_interp[clock_comb$evt_time+1 > clock_comb$rt_csv + clock_comb$iti_ideal] <- NA
fb_comb$decon_interp[fb_comb$evt_time+1 > clock_comb$iti_ideal] <- NA


# numeric axis slice positions
clock_comb <- clock_comb %>%
  mutate(bin_low = as.numeric(sub("[^\\d]+([\\d+\\.]+),.*", "\\1", axis_bin, perl=TRUE)),
         bin_high =as.numeric(sub("[^\\d]+[\\d+\\.]+,([\\d+\\.]+)\\]", "\\1", axis_bin, perl=TRUE)))

clock_comb$bin_center <- rowMeans(clock_comb[, c("bin_low", "bin_high")])

fb_comb <- fb_comb %>%
  mutate(bin_low = as.numeric(sub("[^\\d]+([\\d+\\.]+),.*", "\\1", axis_bin, perl=TRUE)),
         bin_high =as.numeric(sub("[^\\d]+[\\d+\\.]+,([\\d+\\.]+)\\]", "\\1", axis_bin, perl=TRUE)))

fb_comb$bin_center <- rowMeans(fb_comb[, c("bin_low", "bin_high")])

# lags
# # dplyr just freezes on this one, omit for now
clock_comb <- clock_comb %>% group_by(id, run, side, axis_bin, evt_time) %>%
  mutate(decon_prev = dplyr::lag(decon_interp, order_by = run_trial),
         telapsed=clock_onset - clock_onset_prev
         ) %>%
  ungroup() %>%
  mutate(decon_prev_z=as.vector(scale(decon_prev)), iti_ideal_z=as.vector(scale(iti_ideal)))


fb_comb <- fb_comb %>% group_by(id, run, axis_bin, side, evt_time) %>%
  mutate(decon_prev = dplyr::lag(decon_interp, order_by = run_trial),
         telapsed=feedback_onset - feedback_onset_prev) %>%
  ungroup() %>%
  mutate(decon_prev_z=as.vector(scale(decon_prev)), iti_ideal_z=as.vector(scale(iti_ideal)))

# # I wonder whether we should subtract the previous trial's signal to isolate unique variation
clock_comb$decon_change <- clock_comb$decon_interp - clock_comb$decon_prev
fb_comb$decon_change <- fb_comb$decon_interp - fb_comb$decon_prev

# lagged AH and PH mean signal by evt_time -- spot-checked these, correct
clock_comb <- clock_comb %>% group_by(id, run, run_trial, side, evt_time) %>%
  mutate(ah_mean = mean(decon_interp[bin_center>.2 & bin_center<.8]),
         ph_mean = mean(decon_interp[bin_center<.2])) %>%
ungroup()
clock_comb <- clock_comb %>% group_by(id, run, run_trial, side, bin_center) %>%
  mutate(ah_mean_lag = dplyr::lag(ah_mean, order_by = evt_time),
         ph_mean_lag = dplyr::lag(ph_mean, order_by = evt_time)) %>%
  ungroup()


fb_comb <- fb_comb %>% group_by(id, run, run_trial, evt_time, side) %>%
  mutate(ah_mean = mean(decon_interp[bin_center>.2 & bin_center<.8]),
         ph_mean = mean(decon_interp[bin_center<.2]), 
         ah_mean_lag = dplyr::lag(ah_mean, order_by = evt_time),
         ph_mean_lag = dplyr::lag(ph_mean, order_by = evt_time)) %>%
  ungroup()
View(clock_comb)

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
    # summarise(mdecon_interp = mean(decon_interp)) %>% ungroup()
    summarise(mdecon_interp = mean(abs(decon_change))) %>% ungroup() # version with abs signal change
  # }
  
  g <- ggplot(df_sum, aes(x=evt_time, y=mdecon_interp, color = axis_bin, lty = !!gg)) +
    stat_summary(fun.y=mean, geom="line")
  
  # browser()
  # if (!missing(facet_by)) {
  #   g <- g + facet_grid(side ~ vars(facet_by))
  # } else {
  g <- g + facet_wrap(~side) + scale_colour_viridis_d(option = "C") + theme_dark() + geom_vline(xintercept = 0, linetype = 'dotted')
  # }
  
  return(g)
}

##########################
# preliminaries stop here
##########################
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

# effect of past reward, filtering short iti trials
cl1 <- plot_by_summary("clock", trial_split=reward_lag, filter_expr = 'bin_center <.80 & rt_csv<1 & iti_ideal>3 & iti_prev>7 & reward == TRUE') 
cl1 <- cl1 + geom_vline(xintercept = 0:1)
cl2 <- plot_by_summary("clock", trial_split=reward_lag, filter_expr = 'bin_center <.80 & rt_csv>1 & rt_csv<=2 & iti_ideal>3 & iti_prev>7 & reward == TRUE')
cl2 <- cl2 + geom_vline(xintercept = 1:2)
cl3 <- plot_by_summary("clock", trial_split=reward_lag, filter_expr = 'bin_center <.80 & rt_csv>2 & rt_csv<=3 & iti_ideal>3 & iti_prev>7 & reward == TRUE')
cl3 <- cl3 + geom_vline(xintercept = 2:3)
cl4 <- plot_by_summary("clock", trial_split=reward_lag, filter_expr = 'bin_center <.80 & rt_csv>3 & iti_ideal>3 & iti_prev>7 & reward == TRUE')
cl4 <- cl4 + geom_vline(xintercept = 3:4) 
pdf("clock_by_lagged_reward_by_rt_bin.pdf", width = 16, height = 10)
# ggarrange(cl1,cl2,cl3,cl4,ncol = 4, nrow = 1, labels = c("RT: 0-1s", "1-2s", "2-3s","3-4s"))
plot_grid(cl1,cl2,cl3,cl4,ncol = 2, align = 'hv',  labels = c("0-1s", "1-2s", "2-3s","3-4s"))
dev.off()

ce1 <- plot_by_summary("clock", trial_split=entropy_lag, filter_expr = 'bin_center <.80 & rt_csv<1 & iti_ideal>3 & iti_prev>7 & reward == TRUE') 
ce1 <- ce1 + geom_vline(xintercept = 0:1)
ce2 <- plot_by_summary("clock", trial_split=entropy_lag, filter_expr = 'bin_center <.80 & rt_csv>1 & rt_csv<=2 & iti_ideal>3 & iti_prev>7 & reward == TRUE')
ce2 <- ce2 + geom_vline(xintercept = 1:2)
ce3 <- plot_by_summary("clock", trial_split=entropy_lag, filter_expr = 'bin_center <.80 & rt_csv>2 & rt_csv<=3 & iti_ideal>3 & iti_prev>7 & reward == TRUE')
ce3 <- ce3 + geom_vline(xintercept = 2:3)
ce4 <- plot_by_summary("clock", trial_split=entropy_lag, filter_expr = 'bin_center <.80 & rt_csv>3 & iti_ideal>3 & iti_prev>7 & reward == TRUE')
ce4 <- ce4 + geom_vline(xintercept = 3:4) 
pdf("clock_by_lagged_entropy_by_rt_bin)_rew_only.pdf", width = 16, height = 10)
# ggarrange(ce1,ce2,ce3,ce4,ncol = 4, nrow = 1, labels = c("RT: 0-1s", "1-2s", "2-3s","3-4s"))
plot_grid(ce1,ce2,ce3,ce4,ncol = 2, align = 'hv',  labels = c("0-1s", "1-2s", "2-3s","3-4s"))
dev.off()

cs1 <- plot_by_summary("clock", trial_split=swing_above_median, filter_expr = 'bin_center <.80 & rt_csv<1 & iti_ideal>4 & iti_prev>8 & reward == TRUE') 
cs1 <- cs1 + geom_vline(xintercept = 0:1)
cs2 <- plot_by_summary("clock", trial_split=swing_above_median, filter_expr = 'bin_center <.80 & rt_csv>1 & rt_csv<=2 & iti_ideal>4 & iti_prev>8 & reward == TRUE')
cs2 <- cs2 + geom_vline(xintercept = 1:2)
cs3 <- plot_by_summary("clock", trial_split=swing_above_median, filter_expr = 'bin_center <.80 & rt_csv>2 & rt_csv<=3 & iti_ideal>4 & iti_prev>8 & reward == TRUE')
cs3 <- cs3 + geom_vline(xintercept = 2:3)
cs4 <- plot_by_summary("clock", trial_split=swing_above_median, filter_expr = 'bin_center <.80 & rt_csv>3 & iti_ideal>4 & iti_prev>8 & reward == TRUE')
cs4 <- cs4 + geom_vline(xintercept = 3:4) 
pdf("clock_by_rt_swing_by_rt_bin)_rew_only.pdf", width = 16, height = 10)
# ggarrange(ce1,ce2,ce3,ce4,ncol = 4, nrow = 1, labels = c("RT: 0-1s", "1-2s", "2-3s","3-4s"))
ggarrange(cs1,cs2,cs3,cs4,ncol = 2,nrow = 2, align = 'hv',  labels = c("0-1s", "1-2s", "2-3s","3-4s"))
dev.off()

# by contingency -- learnable, rewarded only
cc1 <- plot_by_summary("clock", trial_split=rewFunc, filter_expr = 'bin_center <.80 & rt_csv<1 & iti_ideal>3 & iti_prev>6 & reward == TRUE & rewFunc !="CEVR" & rewFunc !="CEV"' ) 
cc1 <- cc1 + geom_vline(xintercept = 0:1)
cc2 <- plot_by_summary("clock", trial_split=rewFunc, filter_expr = 'bin_center <.80 & rt_csv>1 & rt_csv<=2 & iti_ideal>3 & iti_prev>6 & reward == TRUE & rewFunc !="CEVR" & rewFunc !="CEV"')
cc2 <- cc2 + geom_vline(xintercept = 1:2)
cc3 <- plot_by_summary("clock", trial_split=rewFunc, filter_expr = 'bin_center <.80 & rt_csv>2 & rt_csv<=3 & iti_ideal>3 & iti_prev>6 & reward == TRUE & rewFunc !="CEVR" & rewFunc !="CEV"')
cc3 <- cc3 + geom_vline(xintercept = 2:3)
cc4 <- plot_by_summary("clock", trial_split=rewFunc, filter_expr = 'bin_center <.80 & rt_csv>3 & iti_ideal>3 & iti_prev>6 & reward == TRUE & rewFunc !="CEVR" & rewFunc !="CEV"')
cc4 <- cc4 + geom_vline(xintercept = 3:4) 
pdf("clock_by_learnable_rewFunc_by_rt_bin_rew_only.pdf", width = 16, height = 12)
# ggarrange(ce1,ce2,ce3,ce4,ncol = 4, nrow = 1, labels = c("RT: 0-1s", "1-2s", "2-3s","3-4s"))
ggarrange(cc1,cc2,cc3,cc4,ncol = 2,nrow = 2, align = 'hv',  labels = c("0-1s", "1-2s", "2-3s","3-4s"))
dev.off()

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


####### 
#models

#temporal dependency of decon estimates


hist(clock_comb$telapsed)
hist(fb_comb$telapsed)

mm <- lmer(decon_interp ~ decon_prev*iti_prev + (1 | id/run), clock_comb %>% filter(evt_time == 1))
summary(mm)
# 
# mm2 <- lmer(decon_interp ~ decon_prev*iti_ideal_z*evt_time + (1 | id/run), clock_comb %>% filter(side=="l" & evt_time > 0))
# summary(mm2)

mm2 <- lmer(decon_interp ~ decon_prev*telapsed*evt_time + (1 | id/run), clock_comb %>% filter(side=="l" & evt_time > 0))
summary(mm2)

# modulation by PE

#which(is.na(clock_comb$axis_bin))


m1 <- lmer(decon_interp ~ scale(decon_prev) * reward * scale(bin_center) * evt_time + (1 | id/run) + (1 | side), fb_comb)
summary(m1)



# recode_vec <- bin_centers
# names(bin_centers) <- bin_levels
# 
# vv <- recode(clock_comb$axis_bin, !!!bin_centers)
# 
# clock_comb <- clock_comb %>% mutate(axis_cont=recode(axis_bin, !!!bin_centers))


# test the interaction with contingency
mc1 <- lmer(decon_interp ~ rewFunc*bin_center*rt_csv + 
              side + evt_time*rt_csv + reward * evt_time * bin_center + 
              reward_lag * evt_time * bin_center + (1|id/run), clock_comb %>% filter(iti_prev>7))
summary(mc1)
vif.lme(mc1)
g <- ggpredict(mc1, terms = c("rewFunc", "bin_center [.1,.3,.5,.7]", "rt_csv [1,2,3]"))
g <- plot(g, facet = F, dodge = .4)
g + scale_color_viridis_d(option = "plasma") + theme_dark()

# including shorter preceding ITIs does not help

# include time course
clock_comb$evt_time_f <- as.factor(clock_comb$evt_time + 1)
mc2 <- lmer(decon_interp ~ rewFunc*bin_center*rt_csv*evt_time_f + 
              side + (1|id/run), clock_comb %>% filter(iti_prev>7))
summary(mc2)
car::Anova(mc2)
# can't plot a four-way interaction, but let's look at the 3-way
library(emmeans)
g <- emmip(mc2, bin_center ~  rt_csv | rewFunc, at = list(bin_center = c(.1, .3, .5, .7), rt_csv = c(1,2,3)) )
g + scale_color_viridis_d(option = "plasma") + theme_dark()


# control for current reinforcement
mc3 <- lmer(decon_interp ~ rewFunc*bin_center*rt_csv*evt_time_f*reward + 
              side + (1|id/run), clock_comb %>% filter(iti_prev>7))
summary(mc3)
car::Anova(mc3)

# is anterior signal more auto-correlated?
ma1 <- lmer(decon_interp ~ scale(decon_prev)*scale(bin_center) + side + (1|id/run) + (1 | side), clock_comb)

# is there a lag posterior -> anterior?
ml1 <- lmer(decon_interp ~ scale(ah_mean_lag)*scale(bin_center) + (1|id/run) + (1|side), clock_comb %>% filter(iti_prev > 7))

g <- ggpredict(mc2, terms = c("evt_time_f", "rewFunc", "bin_center [.1,.3,.5,.7]", "rt_csv [1,2,3]"))
g <- plot(g, facet = F, dodge = .4)
g + scale_color_viridis_d(option = "plasma") + theme_dark()



# how does pre-response activity predict rt swings?
ms1 <- lmer(decon_interp ~ decon_prev*bin_center + swing_above_median*bin_center + (1 | id/run) + (1 | side), clock_comb %>% filter(iti_prev > 7 & evt_time < 0))
summary(ms1)
vif.lme(ms1)
g <- ggpredict(ms1, terms=c("swing_above_median", "bin_center", "side"), pretty = FALSE)
plot(g, facet = TRUE)
ggplot(g, aes(as.factor(x), predicted, color = group)) + geom_line()
plot_model(ms1)

# try the right direction decon -> RT swing -- that does not work
ms2 <- lmer(decon_interp ~ bin_center*evt_time*entropy_lag +  (1 | id/run), clock_comb %>% filter(iti_prev > 7))
summary(ms2)
vif.lme(ms2)
g <- plot_model(ms2, show.values = TRUE)
g + ylim(-.02,.02)

mm2 <- lmer(decon_interp ~ scale(decon_prev)*scale(iti_ideal) + (1 | id/run), clock_comb %>% filter(side=="l" & iti_ideal > 1))
summary(mm2)



clock_comb <- clock_comb %>% group_by(id, run, run_trial) %>% 
  mutate(decon_prev = dplyr::lag(decon_interp, 1, by="run_trial")) %>% ungroup()


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
