library(dplyr)
library(tidyverse)
library(psych)
library(ggcorrplot)
library(lme4)


trial_df <- read_csv(file.path("~/code/clock_analysis/fmri/data/mmclock_fmri_decay_factorize_selective_psequate_mfx_trial_statistics.csv.gz")) %>%
  mutate(trial=as.numeric(trial))
trial_df <- trial_df %>%
  group_by(id, run) %>%  
  dplyr::mutate(rt_swing = abs(c(NA, diff(rt_csv))), #compute rt_swing within run and subject
                rt_swing_lr = abs(log(rt_csv/lag(rt_csv))),
                clock_onset_prev=dplyr::lag(clock_onset, 1, by="run"),
                rt_lag = lag(rt_csv) ,
                omission_lag = lag(score_csv==0),
                rt_vmax_lag = lag(rt_vmax),
                v_entropy_wi = scale(v_entropy),
                run_trial=case_when(
                  trial >= 1 & trial <= 50 ~ trial,
                  trial >= 51 & trial <= 100 ~ trial - 50, #dplyr/rlang has gotten awfully picky about data types!!
                  trial >= 101 & trial <= 150 ~ trial - 100,
                  trial >= 151 & trial <= 200 ~ trial - 150,
                  trial >= 201 & trial <= 250 ~ trial - 200,
                  trial >= 251 & trial <= 300 ~ trial - 250,
                  trial >= 301 & trial <= 350 ~ trial - 300,
                  trial >= 351 & trial <= 400 ~ trial - 350,
                  TRUE ~ NA_real_)) %>% ungroup() %>%
  dplyr::mutate(rt_csv=rt_csv/1000, rt_vmax=rt_vmax/10) %>% 
  mutate(rt_vmax_cum=clock_onset + rt_vmax)



setwd('~/Box Sync/SCEPTIC_fMRI/')
l <- read_csv("long_axis_l_2.3mm_clock_onset_decon_locked.csv.gz") %>% mutate(side = 'l')
r <- read_csv("long_axis_r_2.3mm_clock_onset_decon_locked.csv.gz") %>% mutate(side = 'r')
clock <- rbind(l,r)

clock_sum <- clock %>% group_by(id, run, evt_time, axis_bin,side) %>% summarise(mdecon_interp = mean(decon_interp))

pdf('clock_locked_means.pdf', width = 14, height = 8)
#ggplot(clock_sum, aes(evt_time, mdecon_interp, color = axis_bin)) + stat_smooth(method = 'gam') + facet_wrap(~side)
ggplot(clock_sum, aes(factor(evt_time), mdecon_interp, color = axis_bin)) + 
  stat_summary(fun.data=mean_cl_boot, geom="pointrange", position=position_dodge(width=0.4)) + facet_wrap(~side)
dev.off()

clock_comb <- trial_df %>% select(id, run, run_trial, iti_ideal, score_csv, clock_onset, clock_onset_prev) %>% mutate(rewom=if_else(score_csv > 0, "rew", "om")) %>%
  group_by(id, run) %>% mutate(iti_prev=dplyr::lag(iti_ideal, by="run_trial")) %>% ungroup() %>%
  inner_join(clock)
#hist(clock_comb$iti_ideal)

clock_sum <- clock_comb %>% filter(iti_ideal < 3) %>% group_by(id, run, evt_time, axis_bin,side) %>% summarise(mdecon_interp = mean(decon_interp))

pdf('clock_locked_means_filter_lt_3s_iti.pdf', width = 14, height = 8)
#ggplot(clock_sum, aes(evt_time, mdecon_interp, color = axis_bin)) + stat_smooth(method = 'gam') + facet_wrap(~side)
ggplot(clock_sum, aes(factor(evt_time), mdecon_interp, color = axis_bin)) + 
  stat_summary(fun.data=mean_cl_boot, geom="pointrange", position=position_dodge(width=0.4)) + facet_wrap(~side)
dev.off()

#temporal dependency of decon estimates
clock_comb <- clock_comb %>% group_by(id, run, run_trial) %>% 
  mutate(decon_prev = dplyr::lag(decon_interp, 1, by="run_trial"), 
         telapsed=clock_onset - clock_onset_prev) %>%
  ungroup() %>%
  mutate(decon_prev_z=as.vector(scale(decon_prev)), iti_ideal_z=as.vector(scale(iti_ideal)))


hist(clock_comb$telapsed)
mm <- lmer(decon_interp ~ decon_prev*iti_ideal_z + (1 | id/run), clock_comb %>% filter(side=="l" & evt_time > 0))
summary(mm)

mm2 <- lmer(decon_interp ~ decon_prev*iti_ideal_z*evt_time + (1 | id/run), clock_comb %>% filter(side=="l" & evt_time > 0))
summary(mm2)

mm2 <- lmer(decon_interp ~ decon_prev*telapsed*evt_time + (1 | id/run), clock_comb %>% filter(side=="l" & evt_time > 0))
summary(mm2)

clock_comb %>% filter()

#which(is.na(clock_comb$axis_bin))
clock_hack <- clock_comb %>%
  mutate(bin_low = as.numeric(sub("[^\\d]+([\\d+\\.]+),.*", "\\1", axis_bin, perl=TRUE)),
         bin_high =as.numeric(sub("[^\\d]+[\\d+\\.]+,([\\d+\\.]+)\\]", "\\1", axis_bin, perl=TRUE)))

clock_hack$bin_center <- rowMeans(clock_hack[, c("bin_low", "bin_high")])

# recode_vec <- bin_centers
# names(bin_centers) <- bin_levels
# 
# vv <- recode(clock_comb$axis_bin, !!!bin_centers)
# 
# clock_comb <- clock_comb %>% mutate(axis_cont=recode(axis_bin, !!!bin_centers))


mm2 <- lmer(decon_interp ~ decon_prev*telapsed*evt_time + (1 | id/run), clock_comb %>% filter(side=="l" & evt_time > 0))
summary(mm2)



table(clock_comb$telapsed)

head(clock_comb)

library(sjPlot)
plot_model(mm, type="pred", terms=c("decon_prev", "iti_ideal_z"))

mm2 <- lmer(decon_interp ~ scale(decon_prev)*scale(iti_ideal) + (1 | id/run), clock_comb %>% filter(side=="l" & iti_ideal > 1))
summary(mm2)

clock_comb <- clock_comb %>% group_by(id, run, run_trial) %>% 
  mutate(decon_prev = dplyr::lag(decon_interp, 1, by="run_trial")) %>% ungroup()


