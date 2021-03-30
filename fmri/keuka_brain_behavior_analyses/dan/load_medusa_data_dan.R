# reads in data
# normally called from medusa_event_locked_lmer_dan.R

# if running separately, uncomment lines 4-8
# reprocess = T
# if (!reprocess) {
#   wide_only = F  
#   tall_only = T
# }

medusa_dir = "~/Box/SCEPTIC_fMRI/dan_medusa/"
if (replicate_compression) {medusa_dir = "~/Box/SCEPTIC_fMRI/dan_medusa/new_locking_compression"}
cache_dir = "~/Box/SCEPTIC_fMRI/dan_medusa/cache"
repo_directory <- "~/code/clock_analysis"

# reprocess = T # for troubleshooting only
if (!exists("reprocess") || !is.logical(reprocess)) { reprocess=FALSE } #default

stopifnot(dir.exists(medusa_dir))  
stopifnot(dir.exists(cache_dir))  

cwd <- getwd()
setwd(medusa_dir)

if (!reprocess) {
  message("Loading medusa data from cache: ", cache_dir)
  if (wide_only) {
    load(file.path(cache_dir, 'clock_dan_wide_ts.Rdata'))
    # load(file.path(cache_dir, 'feedback_dan_wide_ts.Rdata'))
    load(file.path(cache_dir, 'rt_dan_wide_ts.Rdata'))
  } else if (tall_only) {
    load(file.path(cache_dir, 'clock_dan_tall_ts.Rdata'))
    load(file.path(cache_dir, 'rt_dan_tall_ts.Rdata'))
    load(file.path(cache_dir, "sceptic_trial_df_for_medusa.RData"))
  } else if (streams) {
    load(file.path(cache_dir, 'clock_dan_streams.Rdata'))
    load(file.path(cache_dir, 'rt_dan_streams.Rdata'))
  } else if (visuomotor) {
    load(file.path(cache_dir, 'clock_dan_visuomotor.Rdata'))
    load(file.path(cache_dir, 'rt_dan_visuomotor.Rdata'))
  } else {
    load(file.path(cache_dir, 'clock_dan_wide_ts.Rdata'))
    load(file.path(cache_dir, 'clock_dan_tall_ts.Rdata'))
    # load(file.path(cache_dir, 'feedback_dan_wide_ts.Rdata'))
    load(file.path(cache_dir, 'rt_dan_wide_ts.Rdata'))
    load(file.path(cache_dir, 'rt_dan_tall_ts.Rdata'))
    load(file.path(cache_dir, "sceptic_trial_df_for_medusa.RData"))
  }
} else { 
  # load clock ----
  clock <- as_tibble(read_csv("Schaefer_DorsAttn_2.3mm_clock_long_decon_locked.csv.gz"))
  clock$atlas_value <- as.character(clock$atlas_value)
  
  if (replicate_compression) {clock <- clock %>%
    mutate(run_trial = trial - (run - 1)*50, 
           decon_interp = decon_mean, 
           sd_interp = decon_sd) %>% ungroup()
  }
  clock <- clock %>% arrange(id, run, run_trial, evt_time) 
  # load RT ----
  if (!replicate_compression) {
    
    rt <- as_tibble(read_csv("Schaefer_DorsAttn_2.3mm_rt_long_decon_locked.csv.gz"))
    rt$atlas_value <- as.character(rt$atlas_value)
    if (replicate_compression) {rt <- rt %>%
      mutate(run_trial = trial - (run - 1)*50, 
             decon_interp = decon_mean) %>% ungroup()
    }
    rt <- rt %>% arrange(id, run, run_trial, evt_time)
  }
  message("Labeling regions")
  # add manual labels for Schaefer areas ----
  labels <- as_tibble(read_excel("~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/MNH Dan Labels.xlsx")) %>%
    select(c("roinum", "plot_label", "Stream", "Visuomotor_Gradient", "Stream_Gradient"))
  
  names(labels) <- c("atlas_value","label_short", "stream", "visuomotor_grad", "stream_grad")
  labels$stream_grad <- as.numeric(labels$stream_grad)
  labels <- labels %>% arrange(visuomotor_grad, stream_grad) %>% mutate(
    side  = case_when(
      grepl("L_", label_short) ~ "L",
      grepl("R_", label_short) ~ "R"),
    label_short = substr(label_short, 3, length(label_short)),
    label = paste(visuomotor_grad, stream_grad, label_short, side, sep = "_"),
    stream_side = paste0(stream, "_", side),
    visuomotor_side = paste0(visuomotor_grad, "_", side)) %>% 
    select(c(label, label_short, side, atlas_value, stream, visuomotor_grad, stream_grad, stream_side, visuomotor_side))
  
  # fb <- merge(fb, labels)
  clock <- merge(clock, labels)
  if (!replicate_compression) {
    rt <- merge(rt, labels)
  }
  # load trial-level data ----
  message("Loading trial-level data")
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
                  rt_vmax_change = rt_vmax - rt_vmax_lag,
                  rt_vmax_next = lead(rt_vmax),
                  rt_vmax_change_next = rt_vmax_next - rt_vmax,
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
                  abs_pe = abs(pe_max),
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
  
  
  u_df <- read_csv("~/Box/SCEPTIC_fMRI/sceptic_model_fits/mmclock_fmri_fixed_uv_ureset_mfx_trial_statistics.csv.gz")
  u_df <- u_df %>% select(id, run, trial, u_chosen, u_chosen_lag, u_chosen_change)
  
  trial_df <- inner_join(trial_df,u_df)
  
  # load fixed entropy and RT_vmax
  fixed <-  read_csv(file.path(repo_directory, "fmri/data/mmclock_fmri_fixed_fixedparams_fmri_ffx_trial_statistics.csv.gz"))
  fixed <- as_tibble(fixed) %>% select(c(id, run, trial, v_entropy, rt_vmax, pe_max)) %>%
    rename(v_entropy_fixed = v_entropy, rt_vmax_fixed = rt_vmax, pe_max_fixed = pe_max) %>% group_by(id, run) %>%
    mutate(
      rt_vmax_lag_fixed = lag(rt_vmax_fixed),
      rt_vmax_change_fixed = rt_vmax_fixed - rt_vmax_lag_fixed,
      rt_vmax_next_fixed = lead(rt_vmax_fixed),
      rt_vmax_change_next_fixed = rt_vmax_next_fixed - rt_vmax_fixed,
      v_entropy_wi_fixed = scale(v_entropy_fixed),
      v_entropy_wi_lead_fixed = lead(v_entropy_wi_fixed),
      v_entropy_wi_change_fixed = v_entropy_wi_lead_fixed - v_entropy_wi_fixed,
    )
  
  trial_df <- inner_join(trial_df, fixed, by = c("id", "trial", "run"))
  # add parameters ----
  params <- read_csv(file.path(repo_directory, "fmri/data/mmclock_fmri_decay_factorize_selective_psequate_mfx_sceptic_global_statistics.csv"))
  
  trial_df <- inner_join(trial_df, params, by = "id")
  
  # transform for MEDUSA ----
  message("Transforming for MEDUSA")
  
  clock_comb <- trial_df %>% select(id, run, run_trial, iti_ideal, score_csv, clock_onset, clock_onset_prev, rt_lag, rewFunc,
                                    swing_above_median, first10,reward, reward_lag, rt_above_1s, rt_bin, rt_csv, v_entropy_wi, entropy, entropy_lag, 
                                    gamma, total_earnings, u_chosen, u_chosen_lag, u_chosen_change) %>% mutate(rewom=if_else(score_csv > 0, "rew", "om")) %>%
    group_by(id, run) %>% mutate(iti_prev=dplyr::lag(iti_ideal, by="run_trial")) %>% ungroup() %>%
    inner_join(clock) %>% arrange(id, run, run_trial, evt_time)
  if (!replicate_compression) {
    rt_comb <- trial_df %>% select(id, run, run_trial, iti_ideal, score_csv, feedback_onset, feedback_onset_prev, rt_lag, rewFunc,
                                   swing_above_median, first10,reward, reward_lag, rt_above_1s, rt_bin, rt_csv, entropy, entropy_lag, abs_pe_f, rt_vmax_lag, rt_vmax_change,
                                   gamma, total_earnings, ev,next_swing_above_median, u_chosen, u_chosen_lag, u_chosen_change, rt_vmax_change_next) %>% mutate(rewom=if_else(score_csv > 0, "rew", "om")) %>%
      group_by(id, run) %>% mutate(iti_prev=dplyr::lag(iti_ideal, by="run_trial")) %>% ungroup() %>%
      inner_join(rt) %>% arrange(id, run, run_trial, evt_time)
  }
  
  # 20% of clock- and 32% of feedback-aligned timepoints are from the next trial: censor
  clock_comb$decon_interp[clock_comb$evt_time > clock_comb$rt_csv + clock_comb$iti_ideal] <- NA
  clock_comb$sd_interp[clock_comb$evt_time > clock_comb$rt_csv + clock_comb$iti_ideal] <- NA
  # code on- and offline periods and evt_time^2 for plotting models
  clock_comb <- clock_comb %>% mutate(online = evt_time > -1 & evt_time < rt_csv & evt_time<4)
  clock_comb$online <- as.factor(clock_comb$online)
  
  # fb_comb$decon_interp[fb_comb$evt_time > fb_comb$iti_ideal] <- NA
  # fb_comb$sd_interp[fb_comb$evt_time > fb_comb$iti_ideal] <- NA
  if (!replicate_compression) {
    rt_comb$decon_interp[rt_comb$evt_time > rt_comb$iti_ideal] <- NA
    rt_comb$sd_interp[rt_comb$evt_time > rt_comb$iti_ideal] <- NA
    rt_comb <- rt_comb %>% mutate(online = evt_time > -1 & evt_time < rt_csv & evt_time<4)
    rt_comb$online <- as.factor(rt_comb$online)
  }
  
  # add evt_time as factor
  clock_comb$evt_time_f <- as.factor(clock_comb$evt_time)
  
  # both online and offline event times
  clock_wide <- clock_comb %>% select(id, run, run_trial, evt_time, label, decon_interp) %>% 
    group_by(id, run, run_trial) %>%
    pivot_wider(names_from = c(label, evt_time), values_from = decon_interp) %>% ungroup()
  # clock_wide_cens <- clock_comb %>% select(id, run, run_trial, evt_time, label, decon_interp, online) %>% 
  #   group_by(id, run, run_trial) %>% filter(evt_time < 1 | online == "TRUE") %>% select(!online) %>%
  #   pivot_wider(names_from = c(label, evt_time), values_from = decon_interp) %>% ungroup()
  clock_streams <- clock_comb %>% select(id, run, run_trial, evt_time, label, decon_interp, stream_side) %>% 
    group_by(id, run, run_trial) %>%
    pivot_wider(names_from = c(stream_side, evt_time), values_from = decon_interp) %>% ungroup()
  clock_visuomotor <- clock_comb %>% select(id, run, run_trial, evt_time, label, decon_interp, visuomotor_side) %>% 
    group_by(id, run, run_trial) %>%
    pivot_wider(names_from = c(visuomotor_side, evt_time), values_from = decon_interp) %>% ungroup()
  
  
  # for coxme -- don't filter online times so that we can later interpolate
  clock_cox <- clock_comb %>% select(id, run, run_trial, evt_time, online, label, decon_interp) %>% 
    group_by(id, run, run_trial) %>%
    pivot_wider(names_from = c(label), values_from = decon_interp) %>% ungroup()
  # save clock ----
  if (!replicate_compression) {
    
    message("Saving to cache")
    save(clock_wide, file = file.path(cache_dir, "clock_dan_wide_ts.Rdata"))
    save(clock_comb, file = file.path(cache_dir, "clock_dan_tall_ts.Rdata"))
    save(clock_cox, file = file.path(cache_dir, "clock_dan_medusa_for_coxme.Rdata"))
    save(clock_streams,  file = file.path(cache_dir, "clock_dan_streams.Rdata"))
    save(clock_visuomotor,  file = file.path(cache_dir, "clock_dan_visuomotor.Rdata"))
    
    
    # save RT ----
    # take all preceding timepoints for RT_wide
    rt_comb <- rt_comb %>% arrange(id, run, run_trial, evt_time)
    rt_wide <- rt_comb %>%  select(id, run, run_trial, evt_time, label, decon_interp) %>% 
      group_by(id, run, run_trial) %>%
      pivot_wider(names_from = c(label, evt_time), values_from = decon_interp) %>% ungroup()
    rt_streams <- rt_comb %>% select(id, run, run_trial, evt_time, label, decon_interp, stream_side) %>% 
      group_by(id, run, run_trial) %>%
      pivot_wider(names_from = c(stream_side, evt_time), values_from = decon_interp) %>% ungroup()
    rt_visuomotor <- rt_comb %>% select(id, run, run_trial, evt_time, label, decon_interp, visuomotor_side) %>% 
      group_by(id, run, run_trial) %>%
      pivot_wider(names_from = c(visuomotor_side, evt_time), values_from = decon_interp) %>% ungroup()
    
    
    
    
    save(rt_wide, file = file.path(cache_dir, "rt_dan_wide_ts.Rdata"))
    save(rt_comb, file = file.path(cache_dir, "rt_dan_tall_ts.Rdata"))
    save(trial_df, file=file.path(cache_dir, "sceptic_trial_df_for_medusa.RData"))
    save(rt_streams,  file = file.path(cache_dir, "rt_dan_streams.Rdata"))
    save(rt_visuomotor,  file = file.path(cache_dir, "rt_dan_visuomotor.Rdata"))
    
  }
}
#reset working directory
setwd(cwd)

