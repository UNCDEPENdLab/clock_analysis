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

clock_comb <- clock_comb %>%
  mutate(bin_low = as.numeric(sub("[^\\d]+([\\d+\\.]+),.*", "\\1", axis_bin, perl=TRUE)),
         bin_high =as.numeric(sub("[^\\d]+[\\d+\\.]+,([\\d+\\.]+)\\]", "\\1", axis_bin, perl=TRUE)))

clock_comb$bin_center <- rowMeans(clock_comb[, c("bin_low", "bin_high")])

fb_comb <- fb_comb %>%
  mutate(bin_low = as.numeric(sub("[^\\d]+([\\d+\\.]+),.*", "\\1", axis_bin, perl=TRUE)),
         bin_high =as.numeric(sub("[^\\d]+[\\d+\\.]+,([\\d+\\.]+)\\]", "\\1", axis_bin, perl=TRUE)))

fb_comb$bin_center <- rowMeans(fb_comb[, c("bin_low", "bin_high")])

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
################
ma1 <- lmer(scale(decon_interp) ~ scale(decon_prev)*scale(bin_center) + side + (1|id/run), clock_comb)
summary(ma1)
# cut out superimposed trials -- same results with <10% of the data
ma1_7 <- lmer(scale(decon_interp) ~ scale(decon_prev)*scale(bin_center) + side + (1|id/run), clock_comb %>% filter(iti_prev>7))
summary(ma1_7)


