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
# read in, process; go with "long" [-1:10] clock windows for now, will censor later
#####################
plots = F
reprocess = T
analyze = F

setwd('~/Box/SCEPTIC_fMRI/deconvolved_evt_locked/')

if (reprocess) {
l <- read_csv("long_axis_l_2.3mm_clock_long_decon_locked.csv.gz") %>% mutate(side = 'l')
r <- read_csv("long_axis_r_2.3mm_clock_long_decon_locked.csv.gz") %>% mutate(side = 'r')
clock <- rbind(l,r) 
clock <- clock[clock$id!=11347,]
# load feedback, SWITCH TO LONG
l <- read_csv("long_axis_l_2.3mm_feedback_long_decon_locked.csv.gz") %>% mutate(side = 'l')
r <- read_csv("long_axis_r_2.3mm_feedback_long_decon_locked.csv.gz") %>% mutate(side = 'r')
fb <- rbind(l,r)
fb <- fb[fb$id!=11347,]
# load RT(vmax)
l <- read_csv("long_axis_l_2.3mm_rt_vmax_cum_decon_locked.csv.gz") %>% mutate(side = 'l')
r <- read_csv("long_axis_r_2.3mm_rt_vmax_cum_decon_locked.csv.gz") %>% mutate(side = 'r')
rtvmax <- rbind(l,r)
rtvmax <- rtvmax[rtvmax$id!=11347,]
# load V1 clock
l <- read_csv("l_v1_2.3mm_clock_long_decon_locked.csv.gz") %>% mutate(side = 'l')
r <- read_csv("r_v1_2.3mm_clock_long_decon_locked.csv.gz") %>% mutate(side = 'r')
v1clock <- rbind(l,r)
v1clock <- v1clock[v1clock$id!=11347,]
# load M1
m1Lclock <- read_csv("l_motor_2.3mm_clock_long_decon_locked.csv.gz")
m1Lclock <- m1Lclock[m1Lclock$id!=11347,]
# V1 feedback
l <- read_csv("l_v1_2.3mm_feedback_onset_decon_locked.csv.gz") %>% mutate(side = 'l')
r <- read_csv("r_v1_2.3mm_feedback_onset_decon_locked.csv.gz") %>% mutate(side = 'r')
v1_fb <- rbind(l,r)
v1_fb <- v1_fb[v1_fb$id!=11347,]
# load M1 fb
m1L_fb <- read_csv("l_motor_2.3mm_feedback_onset_decon_locked.csv.gz")
m1L_fb <- m1L_fb[m1L_fb$id!=11347,]


trial_df <- read_csv(file.path("~/code/clock_analysis/fmri/data/mmclock_fmri_decay_factorize_selective_psequate_mfx_trial_statistics.csv.gz")) %>%
  mutate(trial=as.numeric(trial)) %>%
  group_by(id, run) %>%  
  dplyr::mutate(rt_swing = abs(c(NA, diff(rt_csv))), #compute rt_swing within run and subject
                rt_swing_lr = abs(log(rt_csv/lag(rt_csv))),
                clock_onset_prev=dplyr::lag(clock_onset, 1, by="run"),
                rt_lag = lag(rt_csv) ,
                reward = case_when(
                  score_csv>0 ~ "reward",
                  score_csv==0 ~ "omission",
                  TRUE ~ NA_character_),
                reward = as.factor(reward),
                reward_lag = lag(reward),
                iti_prev = lag(iti_ideal),
                omission_lag = lag(score_csv==0),
                rt_vmax_lag = lag(rt_vmax),
                v_entropy_wi = scale(v_entropy),
                v_entropy_wi_lead = lead(v_entropy_wi),
                v_entropy_wi_change = v_entropy_wi_lead - v_entropy_wi,
                entropy = case_when(
                  v_entropy_wi > mean(v_entropy_wi) ~ "high",
                  v_entropy_wi < mean(v_entropy_wi) ~ "low",
                  TRUE ~ NA_character_),
                entropy_lag = lag(entropy),
                rt_change = rt_csv - rt_lag,
                rt_above_1s = rt_csv > 1000,
                swing_above_median = as.factor(abs(rt_change) > median(abs(na.omit(rt_change)))),
                next_swing_above_median = lead(swing_above_median),
                pe_max_lag = lag(pe_max), 
                pe_max_lag2 = lag(pe_max_lag),
                pe_max_lag3 = lag(pe_max_lag2),
                abs_pe_max_lag = abs(pe_max_lag), 
                abs_pe_f  = case_when(
                  abs(pe_max) > mean(abs(pe_max)) ~ 'high abs. PE',
                  abs(pe_max) < mean(abs(pe_max)) ~ 'low abs. PE',
                  TRUE ~ NA_character_),
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
  dplyr::mutate(rt_csv=rt_csv/1000, rt_lag = rt_lag/1000, rt_vmax=rt_vmax/10, rt_vmax_lag = rt_vmax_lag/10) %>% # careful not to do this multiple times
  mutate(rt_vmax_cum=clock_onset + rt_vmax, rt_vmax_cum_lag = lag(rt_vmax_cum))
trial_df <- trial_df %>% group_by(id) %>% mutate(total_earnings = sum(score_csv)) %>% ungroup()
trial_df$rewFunc <- as.factor(trial_df$rewFunc)
levels(trial_df$rewFunc) <- c("DEV", "IEV", "CEV", "CEVR")

# add parameters
params <- read_csv("~/code/clock_analysis/fmri/data/mmclock_fmri_decay_factorize_selective_psequate_mfx_sceptic_global_statistics.csv")

trial_df <- inner_join(trial_df, params, by = "id")

clock_comb <- trial_df %>% select(id, run, run_trial, iti_ideal, score_csv, clock_onset, clock_onset_prev, rt_lag, rewFunc,
                                  swing_above_median, first10,reward, reward_lag, rt_above_1s, rt_bin, rt_csv, entropy, entropy_lag, gamma, total_earnings) %>% mutate(rewom=if_else(score_csv > 0, "rew", "om")) %>%
  group_by(id, run) %>% mutate(iti_prev=dplyr::lag(iti_ideal, by="run_trial")) %>% ungroup() %>%
  inner_join(clock)
fb_comb <- trial_df %>% select(id, run, run_trial, iti_ideal, score_csv, feedback_onset, feedback_onset_prev, rt_lag, rewFunc,
                               swing_above_median, first10,reward, reward_lag, rt_above_1s, rt_bin, rt_csv, entropy, entropy_lag, 
                               abs_pe_f, gamma, total_earnings, ev,next_swing_above_median, rt_vmax_change, rt_vmax_lag, v_entropy_wi, v_entropy_wi_change
                               ) %>% mutate(rewom=if_else(score_csv > 0, "rew", "om")) %>%
  group_by(id, run) %>% mutate(iti_prev=dplyr::lag(iti_ideal, by="run_trial")) %>% ungroup() %>%
  inner_join(fb)
rtvmax_comb <- trial_df %>% select(id, run, run_trial, iti_ideal, score_csv, clock_onset, clock_onset_prev, rt_lag, rewFunc,
                                   swing_above_median, first10,reward, reward_lag, rt_above_1s, rt_bin, rt_csv, entropy, entropy_lag, rt_vmax, gamma, total_earnings) %>% mutate(rewom=if_else(score_csv > 0, "rew", "om")) %>%
  group_by(id, run) %>% mutate(iti_prev=dplyr::lag(iti_ideal, by="run_trial")) %>% ungroup() %>%
  inner_join(rtvmax)
v1_clock_comb <- trial_df %>% select(id, run, run_trial, iti_ideal, score_csv, clock_onset, clock_onset_prev, rt_lag, rewFunc,
                                     swing_above_median, first10,reward, reward_lag, rt_above_1s, rt_bin, rt_csv, entropy, entropy_lag) %>% mutate(rewom=if_else(score_csv > 0, "rew", "om")) %>%
  group_by(id, run) %>% mutate(iti_prev=dplyr::lag(iti_ideal, by="run_trial")) %>% ungroup() %>%
  inner_join(v1clock)
m1L_clock_comb <- trial_df %>% select(id, run, run_trial, iti_ideal, score_csv, clock_onset, clock_onset_prev, rt_lag, rewFunc,
                                      swing_above_median, first10,reward, reward_lag, rt_above_1s, rt_bin, rt_csv, entropy, entropy_lag) %>% mutate(rewom=if_else(score_csv > 0, "rew", "om")) %>%
  group_by(id, run) %>% mutate(iti_prev=dplyr::lag(iti_ideal, by="run_trial")) %>% ungroup() %>%
  inner_join(m1Lclock)

v1_feedback_comb <- trial_df %>% select(id, run, run_trial, iti_ideal, score_csv, feedback_onset, feedback_onset_prev, rt_lag, rewFunc,
                                     swing_above_median, first10,reward, reward_lag, rt_above_1s, rt_bin, rt_csv, entropy, entropy_lag) %>% mutate(rewom=if_else(score_csv > 0, "rew", "om")) %>%
  group_by(id, run) %>% mutate(iti_prev=dplyr::lag(iti_ideal, by="run_trial")) %>% ungroup() %>%
  inner_join(v1_fb)
m1L_feedback_comb <- trial_df %>% select(id, run, run_trial, iti_ideal, score_csv, feedback_onset, feedback_onset_prev, rt_lag, rewFunc,
                                      swing_above_median, first10,reward, reward_lag, rt_above_1s, rt_bin, rt_csv, entropy, entropy_lag) %>% mutate(rewom=if_else(score_csv > 0, "rew", "om")) %>%
  group_by(id, run) %>% mutate(iti_prev=dplyr::lag(iti_ideal, by="run_trial")) %>% ungroup() %>%
  inner_join(m1L_fb)


# 20% of clock- and 32% of feedback-aligned timepoints are from the next trial: censor
clock_comb$decon_interp[clock_comb$evt_time > clock_comb$rt_csv + clock_comb$iti_ideal] <- NA
fb_comb$decon_interp[fb_comb$evt_time > fb_comb$iti_ideal] <- NA
m1L_clock_comb$decon_interp[m1L_clock_comb$evt_time > m1L_clock_comb$iti_ideal] <- NA
v1_clock_comb$decon_interp[v1_clock_comb$evt_time > v1_clock_comb$iti_ideal] <- NA
m1L_feedback_comb$decon_interp[m1L_feedback_comb$evt_time > m1L_feedback_comb$iti_ideal] <- NA
v1_feedback_comb$decon_interp[v1_feedback_comb$evt_time > v1_feedback_comb$iti_ideal] <- NA

# code on- and offline periods and evt_time^2 for plotting models
clock_comb <- clock_comb %>% mutate(online = evt_time > -1 & evt_time < rt_csv)
clock_comb$online <- as.factor(clock_comb$online)
v1_clock_comb <- v1_clock_comb %>% mutate(online = evt_time > -1 & evt_time <= rt_csv + 1)
v1_clock_comb$online <- as.factor(v1_clock_comb$online)
m1L_clock_comb <- m1L_clock_comb %>% mutate(online = evt_time > -1 & evt_time <= rt_csv + 1)
m1L_clock_comb$online <- as.factor(m1L_clock_comb$online)

# try more stringent feedback online window
fb_comb <- fb_comb %>% mutate(online = evt_time > -rt_csv & evt_time < 0)
fb_comb$online <- as.factor(fb_comb$online)
fb_comb$evt_time_sq <- fb_comb$evt_time^2
# from clock onset (-rtvmax) until feedback (rt_csv-rt_vmax)
rtvmax_comb <- rtvmax_comb %>% mutate(online = evt_time > -rt_vmax & evt_time < (rt_csv-rt_vmax))
# try more restrictive window
# rtvmax_comb <- rtvmax_comb %>% mutate(online = evt_time > -rt_vmax & evt_time < (rt_csv-rt_vmax-1))
rtvmax_comb$online <- as.factor(rtvmax_comb$online)
rtvmax_comb$evt_time_sq <- rtvmax_comb$evt_time^2

# add evt_time as factor
clock_comb$evt_time_f <- as.factor(clock_comb$evt_time)
fb_comb$evt_time_f <- as.factor(fb_comb$evt_time)
rtvmax_comb$evt_time_f <- as.factor(rtvmax_comb$evt_time)
# numeric axis slice positions
clock_comb <- clock_comb %>%
  mutate(bin_low = as.numeric(sub("[^\\d]+([\\d+\\.]+),.*", "\\1", axis_bin, perl=TRUE)),
         bin_high =as.numeric(sub("[^\\d]+[\\d+\\.]+,([\\d+\\.]+)\\]", "\\1", axis_bin, perl=TRUE)))

clock_comb$bin_center <- rowMeans(clock_comb[, c("bin_low", "bin_high")])
clock_comb$bin_center_z <- scale(clock_comb$bin_center)

fb_comb <- fb_comb %>%
  mutate(bin_low = as.numeric(sub("[^\\d]+([\\d+\\.]+),.*", "\\1", axis_bin, perl=TRUE)),
         bin_high =as.numeric(sub("[^\\d]+[\\d+\\.]+,([\\d+\\.]+)\\]", "\\1", axis_bin, perl=TRUE)))

fb_comb$bin_center <- rowMeans(fb_comb[, c("bin_low", "bin_high")])
fb_comb$bin_center_z <- scale(fb_comb$bin_center)

rtvmax_comb <- rtvmax_comb %>%
  mutate(bin_low = as.numeric(sub("[^\\d]+([\\d+\\.]+),.*", "\\1", axis_bin, perl=TRUE)),
         bin_high =as.numeric(sub("[^\\d]+[\\d+\\.]+,([\\d+\\.]+)\\]", "\\1", axis_bin, perl=TRUE)))

rtvmax_comb$bin_center <- rowMeans(rtvmax_comb[, c("bin_low", "bin_high")])
rtvmax_comb$bin_center_z <- scale(rtvmax_comb$bin_center)

# lags -- these take very long

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

rtvmax_comb <- rtvmax_comb %>% group_by(id, run, axis_bin, side, evt_time) %>%
  mutate(decon_prev = dplyr::lag(decon_interp, order_by = run_trial),
         telapsed=clock_onset - clock_onset_prev) %>%
  ungroup() %>%
  mutate(decon_prev_z=as.vector(scale(decon_prev)), iti_ideal_z=as.vector(scale(iti_ideal)))


v1_clock_comb <- v1_clock_comb %>% group_by(id, run, axis_bin, side, evt_time) %>%
  mutate(decon_prev = dplyr::lag(decon_interp, order_by = run_trial),
         telapsed=clock_onset - clock_onset_prev) %>%
  ungroup() %>%
  mutate(decon_prev_z=as.vector(scale(decon_prev)), iti_ideal_z=as.vector(scale(iti_ideal)))

m1L_clock_comb <- m1L_clock_comb %>% group_by(id, run, axis_bin, evt_time) %>%
  mutate(decon_prev = dplyr::lag(decon_interp, order_by = run_trial),
         telapsed=clock_onset - clock_onset_prev) %>%
  ungroup() %>%
  mutate(decon_prev_z=as.vector(scale(decon_prev)), iti_ideal_z=as.vector(scale(iti_ideal)))

v1_feedback_comb <- v1_feedback_comb %>% group_by(id, run, axis_bin, side, evt_time) %>%
  mutate(decon_prev = dplyr::lag(decon_interp, order_by = run_trial),
         telapsed=feedback_onset - feedback_onset_prev) %>%
  ungroup() %>%
  mutate(decon_prev_z=as.vector(scale(decon_prev)), iti_ideal_z=as.vector(scale(iti_ideal)))

m1L_feedback_comb <- m1L_feedback_comb %>% group_by(id, run, axis_bin, evt_time) %>%
  mutate(decon_prev = dplyr::lag(decon_interp, order_by = run_trial),
         telapsed=feedback_onset - feedback_onset_prev) %>%
  ungroup() %>%
  mutate(decon_prev_z=as.vector(scale(decon_prev)), iti_ideal_z=as.vector(scale(iti_ideal)))

# # I wonder whether we should subtract the previous trial's signal to isolate unique variation
clock_comb$decon_change <- clock_comb$decon_interp - clock_comb$decon_prev
fb_comb$decon_change <- fb_comb$decon_interp - fb_comb$decon_prev
m1L_clock_comb$decon_change <- m1L_clock_comb$decon_interp - m1L_clock_comb$decon_prev
v1_clock_comb$decon_change <- v1_clock_comb$decon_interp - v1_clock_comb$decon_prev
m1L_feedback_comb$decon_change <- m1L_feedback_comb$decon_interp - m1L_feedback_comb$decon_prev
v1_feedback_comb$decon_change <- v1_feedback_comb$decon_interp - v1_feedback_comb$decon_prev

# lagged AH and PH mean signal by evt_time -- spot-checked these, correct
clock_comb <- clock_comb %>% group_by(id, run, run_trial, side, evt_time) %>%
  mutate(ah_mean = mean(decon_interp[bin_center>.2 & bin_center<.8]),
         ph_mean = mean(decon_interp[bin_center<.2])) %>%
  ungroup() %>% group_by(id, run, run_trial, side, bin_center) %>%
  mutate(ah_mean_lag = dplyr::lag(ah_mean, order_by = evt_time),
         ph_mean_lag = dplyr::lag(ph_mean, order_by = evt_time)) %>%
  ungroup()


fb_comb <- fb_comb %>% group_by(id, run, run_trial, evt_time, side) %>%
  mutate(ah_mean = mean(decon_interp[bin_center>.2 & bin_center<.8]),
         ph_mean = mean(decon_interp[bin_center<.2])) %>%
  ungroup() %>% group_by(id, run, run_trial, side, bin_center) %>%
  mutate(ah_mean_lag = dplyr::lag(ah_mean, order_by = evt_time),
         ph_mean_lag = dplyr::lag(ph_mean, order_by = evt_time)) %>%
  ungroup()
# View(clock_comb)
# summarize for VAR models
fb_var <- fb_comb %>% group_by(id, run, run_trial, evt_time, side, online, entropy, reward, reward_lag, abs_pe_f) %>% 
  summarise(AH = mean(ah_mean), PH = mean(ph_mean)) %>% ungroup()
View(fb_var)

myspread <- function(df, key, value) {
  # quote key
  keyq <- rlang::enquo(key)
  # break value vector into quotes
  valueq <- rlang::enquo(value)
  s <- rlang::quos(!!valueq)
  df %>% gather(variable, value, !!!s) %>%
    unite(temp, !!keyq, variable) %>%
    spread(temp, value)
}
fb_comb <- fb_comb %>% group_by(id,run,run_trial,evt_time,side) %>% mutate(bin_num = rank(bin_center)) %>% ungroup()
fb_comb <- fb_comb %>% mutate(bin6 = round((bin_num + .5)/2,digits = 0)) # also a 6-bin version
fb_wide <- fb_comb %>% select(id, run, run_trial, evt_time, side, bin_num, decon_interp) %>% spread(key = side, decon_interp) %>% myspread(bin_num, c("l", "r"))
names(fb_wide)[5:28] <- paste("hipp", names(fb_wide)[5:28], sep = "_")
fb_wide_ex <- inner_join(fb_wide, trial_df[,c("id", "run", "run_trial", "pe_max", "reward", "v_entropy_wi")], by = c("id", "run", "run_trial"))
fb_wide6 <- fb_comb %>% select(id, run, run_trial, evt_time, side, bin6, decon_interp) %>% group_by(id, run, run_trial, evt_time, side, bin6) %>% summarise(decon6  = mean(decon_interp)) %>% spread(key = side, decon6) %>% myspread(bin6, c("l", "r"))
names(fb_wide6)[5:length(names(fb_wide6))] <- paste("hipp", names(fb_wide6)[5:length(names(fb_wide6))], sep = "_")
fb_wide6_ex <- inner_join(fb_wide6, trial_df[,c("id", "run", "run_trial", "pe_max", "reward", "v_entropy_wi")], by = c("id", "run", "run_trial"))

# save fb responses for trial-wise prediction analyses

library(data.table)
slices <- names(fb_wide)[5:28]
fb_wide_t <- dcast(setDT(fb_wide), id + run + run_trial ~ evt_time, value.var = slices)


setwd("~/Box/SCEPTIC_fMRI/var/")
save(fb_wide, fb_wide_ex, fb_wide6, fb_wide6_ex, file = "feedback_hipp_wide_ts.Rdata")
save(fb_comb, file = "feedback_hipp_tall_ts.Rdata")
save(fb_wide_t, file = 'feedback_hipp_tallest_by_timepoint_decon.Rdata')


clock_comb <- clock_comb %>% group_by(id,run,run_trial,evt_time,side) %>% mutate(bin_num = rank(bin_center)) %>% ungroup()
# take only online event times
clock_wide <- clock_comb %>% filter(online==T) %>% select(id, run, run_trial, evt_time, side, bin_num, decon_interp) %>% spread(key = side, decon_interp) %>% myspread(bin_num, c("l", "r"))
names(clock_wide)[5:28] <- paste("hipp", names(clock_wide)[5:28], sep = "_")
clock_wide_ex <- inner_join(clock_wide, trial_df[,c("id", "run", "run_trial", "pe_max", "reward", "v_entropy_wi", "swing_above_median")], by = c("id", "run", "run_trial"))
setwd("~/Box/SCEPTIC_fMRI/var/")

# setwd("~/Box/SCEPTIC_fMRI/var/")
save(clock_wide, clock_wide_ex, file = "clock_hipp_wide_ts.Rdata")


}
if (!reprocess) {
  setwd("~/Box/SCEPTIC_fMRI/var/")
  load('clock_hipp_wide_ts.Rdata')
  load('feedback_hipp_wide_ts.Rdata')
  # load("feedback_hipp_tall.Rdata")
  load('feedback_hipp_tall_ts.Rdata')
  
}
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

setwd('~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/plots/')

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
mm <- lmer(decon_interp ~ decon_prev*iti_prev + (1 | id/run), clock_comb %>% filter(evt_time == 1))
summary(mm)
mm2 <- lmer(decon_interp ~ decon_prev*telapsed*evt_time + (1 | id/run), clock_comb %>% filter(side=="l" & evt_time > 0))
summary(mm2)
m1 <- lmer(decon_interp ~ (scale(decon_prev) + scale(bin_center) + evt_time)^3 + (1 | id/run) + (1 | side), fb_comb %>% filter(evt_time>0))
summary(m1)

# does entropy differentially modulate signal autocorrelation along the axis?
cm2 <- lmer(decon_interp ~ (decon_prev_z + bin_center_z + online + entropy)^3 + (1 | id/run) + (1 | side), clock_comb %>% filter (iti_prev>4 & iti_ideal >2))
summary(cm2)
vif.lme(cm2)
g <- ggpredict(cm2, terms = c("entropy_lag", "bin_center_z [-2,0,2]"))
g <- ggpredict(cm2, terms = c("online", "bin_center_z [-2,0,2]"))
# g <- ggpredict(cm2, terms = c("entropy_lag", "decon_prev [.1,.5,.9]"))
g <- plot(g, facet = F, dodge = .4)
g + scale_color_viridis_d(option = "plasma") + theme_dark()

fm2 <- lmer(decon_interp ~ (decon_prev_z + bin_center_z + online)^2 + (1 | id/run) + (1 | side), clock_comb %>% filter (iti_prev>4 & iti_ideal >3))
summary(fm2)
vif.lme(fm2)
# g <- ggpredict(fm2, terms = c("entropy", "bin_center [-2,0,2]"))
g <- ggpredict(fm2, terms = c("online", "bin_center_z [-2,0,2]"))
g <- plot(g, facet = F, dodge = .2)
pdf("ah_vs_ph_online_lmer.pdf", width = 6, height = 6)
g + scale_color_viridis_d(option = "plasma") + theme_dark()
dev.off()

# test for ramps in AH in rtvmax-aligned data: quadratic term

rm1 <- lmer(decon_interp ~ evt_time*bin_center_z + evt_time_sq*bin_center_z + decon_prev_z + entropy_lag + reward_lag + (1 | id/run) + (1 | side), rtvmax_comb %>% filter (online == "TRUE" & iti_prev>2 & iti_ideal>2 & rt_csv >1 & (rewFunc!="CEVR")))
summary(rm1)
vif.lme(rm1)

# add entropy modulation
rm2 <- lmer(decon_interp ~ evt_time*bin_center_z*entropy_lag + evt_time_sq*bin_center_z*entropy_lag + decon_prev_z + reward_lag + scale(rt_csv)*evt_time + (1 | id/run) + (1 | side), rtvmax_comb %>% filter (online == "TRUE" & iti_prev>2 & iti_ideal>2 & rt_csv >1 & (rewFunc!="CEVR")))
summary(rm2)
vif.lme(rm2)
car::Anova(rm2)
anova(rm1,rm2)

# plot out the ramp and entropy modulation
g <- ggpredict(rm2, terms = c("evt_time_sq [0,1,2,3,4]", "bin_center_z [-2,-1, 0, 1,2]", "entropy_lag"))
g <- plot(g, facet = F, dodge = .2)
pdf("ramps_in_AH.pdf", width = 6, height = 6)
g + scale_color_viridis_d(option = "plasma") + theme_dark()
dev.off()

# nesting within trial
rm2t <- lmer(decon_interp ~ evt_time*bin_center_z*entropy_lag + evt_time_sq*bin_center_z*entropy_lag + decon_prev_z + reward_lag + scale(rt_csv)*evt_time + (1 | id/run/run_trial) + (1 | side), rtvmax_comb %>% filter (online == "TRUE" & iti_prev>2 & iti_ideal>2 & rt_csv >1 & (rewFunc!="CEVR")))
summary(rm2t)
vif.lme(rm2)
car::Anova(rm2, '3')


# test the same with time as factor
rm2f <- lmer(decon_interp ~ evt_time_f*bin_center_z*entropy_lag + decon_prev_z + reward_lag + scale(rt_csv)*evt_time_f + (1 | id/run) + (1 | side), rtvmax_comb %>% filter (online == "TRUE" & iti_prev>2 & iti_ideal>2 & rt_csv >1 & (rewFunc!="CEVR") ))
summary(rm2f)
vif.lme(rm2f)
car::Anova(rm2f)
g <- ggpredict(rm2f, terms = c("evt_time_f", "bin_center_z [-2,-1, 0, 1,2]", "entropy_lag"))
g <- plot(g, facet = F, dodge = .4)
pdf("ramps_in_AH_f.pdf", width = 6, height = 6)
g + scale_color_viridis_d(option = "plasma") + theme_dark()
dev.off()

# check that it's specific to entropy and not last reward
# it is, but PH is also more active after omissions and AH, after rewards
rm3 <- lmer(decon_interp ~ evt_time*bin_center_z + evt_time_sq*bin_center_z*reward_lag + decon_prev_z + entropy_lag + (1 | id/run) + (1 | side), rtvmax_comb %>% filter (online == "TRUE" & iti_prev>2 & iti_ideal>2 & rt_csv >1 & (rewFunc!="CEVR") ))
summary(rm3)
vif.lme(rm3)
car::Anova(rm3)
g <- ggpredict(rm3, terms = c("bin_center_z [-2,-1, 0, 1,2]", "reward_lag"))
g <- plot(g, facet = F, dodge = .3)
g + scale_color_viridis_d(option = "plasma") + theme_dark()

anova(rm1,rm2, rm3)

# offline activity in PH vs AH
# comparison model without time * bin
om0 <- lmer(decon_interp ~ evt_time_f + bin_center_z*reward + decon_prev_z + scale(rt_csv)*evt_time_f + (1 | id/run) + (1 | side), fb_comb %>% filter (iti_prev>1 & iti_ideal>8 & evt_time < 9))
# without time * reward * bin
om0a <- lmer(decon_interp ~ evt_time_f*bin_center_z + bin_center_z*reward + decon_prev_z + scale(rt_csv)*evt_time_f + (1 | id/run) + (1 | side), fb_comb %>% filter (iti_prev>1 & iti_ideal>8 & evt_time < 9))

om1 <- lmer(decon_interp ~ evt_time_f*bin_center_z*reward + decon_prev_z + scale(rt_csv)*evt_time_f + (1 | id/run) + (1 | side), fb_comb %>% filter (iti_prev>1 & iti_ideal>8 & evt_time < 9))
summary(om1)
vif.lme(om1)
car::Anova(om1)

anova(om0, om0a, om1)
g <- ggpredict(om1, terms = c("evt_time_f", "bin_center_z [-2,-1, 0, 1,2]", "reward"))
g <- plot(g, facet = F, dodge = .4)
pdf("offline_reward_AH_PH.pdf", width = 6, height = 6)
g + scale_color_viridis_d(option = "plasma") + theme_dark()
dev.off()

om2 <- lmer(decon_interp ~ evt_time_f*bin_center_z*abs_pe_f + decon_prev_z + scale(rt_csv)*evt_time_f + (1 | id/run) + (1 | side), fb_comb %>% filter (iti_prev>1 & iti_ideal>8 & evt_time < 9))
summary(om2)
vif.lme(om2)
car::Anova(om2)
g <- ggpredict(om2, terms = c("evt_time_f", "bin_center_z [-2,-1, 0, 1,2]", "abs_pe_f"))
g <- plot(g, facet = F, dodge = .4)
pdf("offline_abs_pe_AH_PH.pdf", width = 6, height = 6)
g + scale_color_viridis_d(option = "plasma") + theme_dark()
dev.off()
anova(om1, om2)

# moderation by gamma?

om3 <- lmer(decon_interp ~ evt_time_f*bin_center_z*abs_pe_f*reward + decon_prev_z + scale(rt_csv)*evt_time_f + (1 | id/run) + (1 | side), fb_comb %>% filter (iti_prev>1 & iti_ideal>8 & evt_time < 9))
summary(om3)
vif.lme(om3)
car::Anova(om3)
g <- ggpredict(om3, terms = c("evt_time_f", "abs_pe_f", "reward"))
g <- plot(g, facet = F, dodge = .4)
pdf("offline_abs_pe_reward_AH_PH.pdf", width = 6, height = 6)
g + scale_color_viridis_d() + theme_dark()
dev.off()

# reversion in MEDUSA
# need a measure of closeness to RTvmax: let's start with EV
e1 = 38
xtabs(~bin_center + swing_above_median, fb_comb %>% filter(is.na(decon_interp)))

pdf("lose_switch_reversal_to_vmax.pdf", height = 20, width = 20)
#lty = swing_above_median, 
ggplot(fb_comb, aes(evt_time, decon_interp, color = factor(bin_center), lty = next_swing_above_median)) +
  stat_smooth(method = 'gam',se = F) + scale_color_viridis_d("TEST") + facet_grid
dev.off()
# Do shorter ITIs disrupt processing?  No, they just seem to fall asleep during ITIs
itim1 <- lmer(scale(rt_csv) ~ scale(rt_lag)*scale(rt_vmax)*scale(iti_prev) + (1 | id/run), trial_df )
summary(itim1)

################
# responses as a function of trial and explore/exploit

# preliminary plots: trial
# trial by condition
pdf("trial_fb_hipp_AH_PH_gam_rew.pdf", width = 11, height = 8)
ggplot(fb_comb,aes(run_trial,decon_interp, color = axis_bin, lty = reward)) + geom_smooth(method = "gam", se = F) + facet_wrap(~rewFunc) + scale_color_viridis_d() + theme_dark()
dev.off()
# ...and reward
pdf("trial_fb_hipp_AH_PH_loess_rew.pdf", width = 11, height = 8)
ggplot(fb_comb,aes(run_trial,decon_interp, color = axis_bin, lty = reward)) + geom_smooth(method = "loess", se = F) + facet_wrap(~rewFunc) + scale_color_viridis_d() + theme_dark()
dev.off()


# wave form by trial
pdf("fb_hipp_AP_trial_rew.pdf", width = 11, height = 8)
ggplot(fb_comb, aes(evt_time, decon_interp, color = axis_bin, lty = reward)) + geom_smooth(method = "gam", se = F) + facet_grid(first10 ~ rewFunc) + scale_color_viridis_d() + theme_dark()
dev.off()

# wave form as a function of next swing?
fb_comb <- fb_comb %>% group_by(id,run) %>% mutate(swing_above_median_lead = lead(swing_above_median)) %>% ungroup()

pdf("fb_hipp_AP_next_swing_loess_rewFunc.pdf", width = 11, height = 8)
ggplot(fb_comb[!is.na(fb_comb$swing_above_median_lead) & fb_comb$evt_time < 8,], aes(evt_time, decon_interp, color = axis_bin, lty = swing_above_median_lead)) + geom_smooth(method = "loess", se = F) + facet_grid(rewFunc~reward)+ scale_color_viridis_d() + theme_dark()
dev.off()
fb_comb$trial_neg_inv <- -1/fb_comb$run_trial

ee1 <- lmer(decon_interp ~ bin_center_z*trial_neg_inv + scale(rt_csv)*evt_time + (1 | id/run) + (1 | side), fb_comb %>% filter (iti_prev>1 & iti_ideal>8 & evt_time < 9))
summary(ee1)
vif.lme(ee1)
car::Anova(ee1)
g <- ggpredict(ee1, terms = c(""))
g <- plot(g, facet = F, dodge = .4)
pdf("offline_abs_pe_reward_AH_PH.pdf", width = 6, height = 6)
g + scale_color_viridis_d() + theme_dark()
dev.off()

ee2 <- lmer(decon_interp ~ bin_center_z*trial_neg_inv*rewFunc + scale(rt_csv)*evt_time + (1 | id/run) + (1 | side), fb_comb %>% filter (iti_prev>1 & iti_ideal>8 & evt_time < 9))
summary(ee2)
car::Anova(ee2)
anova(ee1,ee2)
library(emmeans)
em2 <- emmeans(ee2, "bin_center_z", by = c( "rewFunc","trial_neg_inv"), at = list( bin_center_z = c(-2,2), trial_neg_inv = c(-1, -.1,-.05, -.02)))
plot(em2, horiz = F)
df2 <- as.data.frame(em2)

ee3 <- lmer(decon_interp ~ bin_center_z*evt_time*swing_above_median_lead + scale(rt_csv)*evt_time + (1 | id/run) + (1 | side), fb_comb %>% filter (iti_prev>1 & iti_ideal>8 & evt_time < 9))
summary(ee3)
car::Anova(ee3)

# + completely general time
ee4 <- lmer(decon_interp ~ bin_center_z*evt_time_f*swing_above_median_lead + scale(rt_csv)*evt_time_f + (1 | id/run) + (1 | side), fb_comb %>% filter (iti_prev>1 & iti_ideal>8 & evt_time < 9))
summary(ee4)
car::Anova(ee4, '3')
em4 <- as.data.frame(emmeans(ee4, "bin_center_z", by = c( "evt_time_f","swing_above_median_lead"), at = list( bin_center_z = c(-2,2))))
ggplot(em4, aes(evt_time_f, emmean, color = bin_center_z, lty = swing_above_median_lead)) + geom_point() + geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL))

# + reward
ee5 <- lmer(decon_interp ~ bin_center_z*evt_time_f*swing_above_median_lead*reward + scale(rt_csv)*evt_time_f + (1 | id/run) + (1 | side), fb_comb %>% filter (iti_prev>1 & iti_ideal>8 & evt_time < 9))
summary(ee5)
car::Anova(ee5, '3')
vif.lme(ee5)
anova(ee1,ee2,ee3,ee4,ee5)

# replicate lm decoding analyses
# scale(-1/run_trial)*rewFunc + reward + scale(rt_csv) + scale(rt_vmax_lag) + scale(rt_vmax_change) + v_entropy_wi
fb_comb$bin_num_f <- as.factor(fb_comb$bin_num)
dm1 <- lmer(decon_interp ~ bin_num_f*evt_time_f*scale(-1/run_trial)*rewFunc + 
              bin_num_f*evt_time_f*scale(rt_csv) + 
              bin_num_f*evt_time_f*scale(rt_vmax_lag) +
              bin_num_f*evt_time_f*scale(rt_vmax_change) + 
              bin_num_f*evt_time_f*scale(v_entropy_wi) +
              bin_num_f*evt_time_f*scale(v_entropy_wi_change) +
              (1 | id/run) + (1 | side), fb_comb %>% filter (iti_prev>1 & iti_ideal>8 & evt_time < 9))
summary(dm1)
car::Anova(dm1, '3')
vif(dm1)
library(emmeans)
r1 <- emtrends(dm1, var = 'rt_vmax_change', specs = c('bin_num_f','evt_time_f'), data = fb_comb %>% filter (iti_prev>1 & iti_ideal>8 & evt_time < 9))
r1 <- as.data.frame(r1)
ggplot(r1, aes(evt_time_f, bin_num_f, color = rt_vmax_change.trend)) + geom_tile()

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
              (1 | id/run) , fb_comb %>% filter (iti_prev>1 & iti_ideal>8 & evt_time < 9))
summary(dm4)
car::Anova(dm4, '3')
em <- emmeans(dm4, )
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
ee6 <- lmer(decon_interp ~ bin_center_z*evt_time_f*swing_above_median_lead*reward*run_trial + scale(rt_csv)*evt_time_f + (1 | id/run) + (1 | side), fb_comb %>% filter (iti_prev>1 & iti_ideal>8 & evt_time < 9))
summary(ee6)
car::Anova(ee6, '3')


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
# setwd('~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/plots/')
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
