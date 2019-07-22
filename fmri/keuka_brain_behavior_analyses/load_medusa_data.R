medusa_dir="~/Box/SCEPTIC_fMRI/deconvolved_evt_locked"
cache_dir="~/Box/SCEPTIC_fMRI/var"
if (!exists("reprocess") || !is.logical(reprocess)) { reprocess=FALSE } #default

stopifnot(dir.exists(medusa_dir))  
stopifnot(dir.exists(cache_dir))  

cwd <- getwd()
setwd(medusa_dir)

if (!reprocess) {
  message("Loading medusa data from cache: ", cache_dir)
  
  load(file.path(cache_dir, 'clock_hipp_wide_ts.Rdata'))
  load(file.path(cache_dir, 'feedback_hipp_wide_ts.Rdata'))
  # load(file.path(cache_dir, "feedback_hipp_tall.Rdata"))
  load(file.path(cache_dir, 'feedback_hipp_tall_ts.Rdata'))
  load(file.path(cache_dir, "sceptic_trial_df_for_medusa.RData"))
  
} else {
  
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
  
  
  trial_df <- read_csv(file.path(repo_directory, "fmri/data/mmclock_fmri_decay_factorize_selective_psequate_mfx_trial_statistics.csv.gz")) %>%
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
  params <- read_csv(file.path(repo_directory, "fmri/data/mmclock_fmri_decay_factorize_selective_psequate_mfx_sceptic_global_statistics.csv"))
  
  trial_df <- inner_join(trial_df, params, by = "id")
  
  clock_comb <- trial_df %>% select(id, run, run_trial, iti_ideal, score_csv, clock_onset, clock_onset_prev, rt_lag, rewFunc,
          swing_above_median, first10,reward, reward_lag, rt_above_1s, rt_bin, rt_csv, entropy, entropy_lag, gamma, total_earnings) %>% mutate(rewom=if_else(score_csv > 0, "rew", "om")) %>%
      group_by(id, run) %>% mutate(iti_prev=dplyr::lag(iti_ideal, by="run_trial")) %>% ungroup() %>%
      inner_join(clock)
  fb_comb <- trial_df %>% select(id, run, run_trial, iti_ideal, score_csv, feedback_onset, feedback_onset_prev, rt_lag, rewFunc,
          swing_above_median, first10,reward, reward_lag, rt_above_1s, rt_bin, rt_csv, entropy, entropy_lag, abs_pe_f, gamma, total_earnings, ev,next_swing_above_median) %>% mutate(rewom=if_else(score_csv > 0, "rew", "om")) %>%
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
  
  # I wonder whether we should subtract the previous trial's signal to isolate unique variation
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
  
  #View(fb_var)
  
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
  
  slices <- names(fb_wide)[5:28]
  fb_wide_t <- dcast(setDT(fb_wide), id + run + run_trial ~ evt_time, value.var = slices)
  
  save(fb_wide, fb_wide_ex, fb_wide6, fb_wide6_ex, file = file.path(cache_dir, "feedback_hipp_wide_ts.Rdata"))
  save(fb_comb, file = file.path(cache_dir, "feedback_hipp_tall_ts.Rdata"))
  save(fb_wide_t, file = file.path(cache_dir, 'feedback_hipp_tallest_by_timepoint_decon.Rdata'))
  
  clock_comb <- clock_comb %>% group_by(id,run,run_trial,evt_time,side) %>% mutate(bin_num = rank(bin_center)) %>% ungroup()
  
  # take only online event times
  clock_wide <- clock_comb %>% filter(online==T) %>% select(id, run, run_trial, evt_time, side, bin_num, decon_interp) %>% spread(key = side, decon_interp) %>% myspread(bin_num, c("l", "r"))
  names(clock_wide)[5:28] <- paste("hipp", names(clock_wide)[5:28], sep = "_")
  clock_wide_ex <- inner_join(clock_wide, trial_df[,c("id", "run", "run_trial", "pe_max", "reward", "v_entropy_wi", "swing_above_median")], by = c("id", "run", "run_trial"))
  
  save(clock_wide, clock_wide_ex, file = file.path(cache_dir, "clock_hipp_wide_ts.Rdata"))
  
  save(trial_df, file=file.path(cache_dir, "sceptic_trial_df_for_medusa.RData"))
}

#reset working directory
setwd(cwd)

