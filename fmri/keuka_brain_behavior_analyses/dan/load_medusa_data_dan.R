library(dplyr)
library(tidyr)
library(haven)
library(readxl)
library(readr)
# reads in data
# normally called from medusa_event_locked_lmer_dan.R
# if running separately, uncomment lines 4-8
use_new_pipeline <- TRUE
# reprocess <- T
# if (!reprocess) {
#   wide_only <- F
#   tall_only <- T
# }

# medusa_dir = "~/Box/SCEPTIC_fMRI/dan_medusa/"
medusa_dir <- "/Volumes/GoogleDrive/My Drive/SCEPTIC_fMRI/dan_medusa"
if (use_new_pipeline) {
  #medusa_dir <- "~/Box/SCEPTIC_fMRI/dan_medusa/new_locking_compression"
  medusa_dir <- "/Volumes/GoogleDrive/My Drive/SCEPTIC_fMRI/dan_medusa/new_locking_compression"
}
# cache_dir = "~/Box/SCEPTIC_fMRI/dan_medusa/cache"
cache_dir <- file.path(medusa_dir, "cache")
repo_directory <- "~/code/clock_analysis"
# repo_directory <- "~/Data_Analysis/clock_analysis"

source(file.path(repo_directory, "fmri/keuka_brain_behavior_analyses/dan/get_trial_data.R"))

#reprocess = TRUE # for troubleshooting only
if (!exists("reprocess") || !is.logical(reprocess)) {
  reprocess <- FALSE
} # default

stopifnot(dir.exists(medusa_dir))
stopifnot(dir.exists(cache_dir))

if (!exists("wide_only") || !is.logical(wide_only)) {
  wide_only <- FALSE
} # default

if (!exists("tall_only") || !is.logical(tall_only)) {
  tall_only <- FALSE
} # default

if (!exists("streams") || !is.logical(streams)) {
  streams <- FALSE
} # default

if (!exists("visuomotor") || !is.logical(visuomotor)) {
  visuomotor <- FALSE
} # default

if (!exists("visuomotor_long") || !is.logical(visuomotor_long)) {
  visuomotor_long <- FALSE
} # default


if (!reprocess) {
  message("Loading medusa data from cache: ", cache_dir)
  load(file.path(cache_dir, "sceptic_trial_df_for_medusa.RData")) # always load trial_df
  if (wide_only) {
    load(file.path(cache_dir, "clock_dan_wide_ts.Rdata"))
    # load(file.path(cache_dir, 'feedback_dan_wide_ts.Rdata'))
    load(file.path(cache_dir, "rt_dan_wide_ts.Rdata"))
  } else if (tall_only) {
    load(file.path(cache_dir, "clock_dan_tall_ts.Rdata"))
    load(file.path(cache_dir, "rt_dan_tall_ts.Rdata"))
    load(file.path(cache_dir, "sceptic_trial_df_for_medusa.RData"))
  } else if (streams) {
    load(file.path(cache_dir, "clock_dan_streams.Rdata"))
    load(file.path(cache_dir, "rt_dan_streams.Rdata"))
  } else if (visuomotor) {
    load(file.path(cache_dir, "clock_dan_visuomotor.Rdata"))
    load(file.path(cache_dir, "rt_dan_visuomotor.Rdata"))
  } else if (visuomotor_long) {
    load(file.path(cache_dir, "clock_dan_visuomotor_long.Rdata"))
    load(file.path(cache_dir, "clock_dan_visuomotor_long_online.Rdata"))
    load(file.path(cache_dir, "rt_dan_visuomotor_long.Rdata"))
  } else {
    load(file.path(cache_dir, "clock_dan_wide_ts.Rdata"))
    load(file.path(cache_dir, "clock_dan_tall_ts.Rdata"))
    # load(file.path(cache_dir, 'feedback_dan_wide_ts.Rdata'))
    load(file.path(cache_dir, "rt_dan_wide_ts.Rdata"))
    load(file.path(cache_dir, "rt_dan_tall_ts.Rdata"))
    load(file.path(cache_dir, "sceptic_trial_df_for_medusa.RData"))
  }
} else {

  # load clock ----
  clock <- as_tibble(read_csv(file.path(medusa_dir, "Schaefer_DorsAttn_2.3mm_clock_long_decon_locked.csv.gz")))
  clock$atlas_value <- as.character(clock$atlas_value)

  if (use_new_pipeline) {
    clock <- clock %>%
      mutate(run_trial = trial - (run - 1) * 50) %>%
      ungroup() %>%
      dplyr::rename(
        decon_interp = decon_mean,
        sd_interp = decon_sd
      ) %>%
      dplyr::select(-decon_median)
  }
  clock <- clock %>% arrange(id, run, run_trial, evt_time)

  # load RT ----
  if (use_new_pipeline) {
    rt <- as_tibble(read_csv(file.path(medusa_dir, "Schaefer_DorsAttn_2.3mm_rt8_decon_locked.csv.gz")))
    rt$atlas_value <- as.character(rt$atlas_value)
  } else {
    rt <- as_tibble(read_csv(file.path(medusa_dir, "Schaefer_DorsAttn_2.3mm_rt_long_decon_locked.csv.gz")))
  }
  
  if (use_new_pipeline) {
    rt <- rt %>%
      mutate(run_trial = trial - (run - 1) * 50) %>%
      ungroup() %>%
      dplyr::rename(
        decon_interp = decon_mean,
        sd_interp = decon_sd
      ) %>%
      dplyr::select(-decon_median)
  }
  rt <- rt %>% arrange(id, run, run_trial, evt_time)

  # ---- add manual labels for Schaefer areas ----
  message("Labeling regions")

  labels <- as_tibble(read_excel(file.path(repo_directory, "fmri/keuka_brain_behavior_analyses/dan/MNH Dan Labels.xlsx"))) %>%
    select(c("roinum", "plot_label", "Stream", "Visuomotor_Gradient", "Stream_Gradient"))

  names(labels) <- c("atlas_value", "label_short", "stream", "visuomotor_grad", "stream_grad")
  labels$stream_grad <- as.numeric(labels$stream_grad)
  labels <- labels %>%
    arrange(visuomotor_grad, stream_grad) %>%
    mutate(
      atlas_value = as.character(atlas_value), # to match modifications above
      side = case_when(
        grepl("L_", label_short) ~ "L",
        grepl("R_", label_short) ~ "R"
      ),
      label_short = substr(label_short, 3, length(label_short)), # drop off side from label
      label = paste(visuomotor_grad, stream_grad, label_short, side, sep = "_"),
      stream_side = paste0(stream, "_", side),
      visuomotor_side = paste0(visuomotor_grad, "_", side)
    ) %>%
    select(c(label, label_short, side, atlas_value, stream, visuomotor_grad, stream_grad, stream_side, visuomotor_side))

  clock <- dplyr::full_join(clock, labels, by = c("atlas_value"))
  rt <- dplyr::full_join(rt, labels, by="atlas_value")
  
  # load trial-level data ----
  message("Loading trial-level data")

  # add parameters ----
  trial_df <- get_trial_data(repo_directory)
  
  # transform for MEDUSA ----

  message("Transforming for MEDUSA")

  clock_comb <- trial_df %>%
    select(
      id, run, run_trial, iti_ideal, score_csv, clock_onset, clock_onset_prev, rt_lag, rewFunc,
      swing_above_median, first10, reward, reward_lag, rt_above_1s, rt_bin, rt_csv, v_entropy_wi, entropy_split, entropy_split_lag,
      total_earnings, u_chosen, u_chosen_lag, u_chosen_change
    ) %>%
    mutate(rew_om = if_else(score_csv > 0, "rew", "om")) %>%
    group_by(id, run) %>%
    mutate(iti_prev = dplyr::lag(iti_ideal, by = "run_trial")) %>%
    ungroup() %>%
    inner_join(clock) %>%
    arrange(id, run, run_trial, evt_time)
  
  rm(clock) # not used below

  rt_comb <- trial_df %>%
    select(
      id, run, run_trial, iti_ideal, score_csv, feedback_onset, feedback_onset_prev, rt_lag, rewFunc,
      swing_above_median, first10, reward, reward_lag, rt_above_1s, rt_bin, rt_csv, entropy_split, entropy_split_lag, abs_pe_f, rt_vmax_lag, rt_vmax_change,
      total_earnings, ev, next_swing_above_median, u_chosen, u_chosen_lag, u_chosen_change, rt_vmax_change_next
    ) %>%
    mutate(rew_om = if_else(score_csv > 0, "rew", "om")) %>%
    group_by(id, run) %>%
    mutate(iti_prev = dplyr::lag(iti_ideal, by = "run_trial")) %>%
    ungroup() %>%
    inner_join(rt) %>%
    arrange(id, run, run_trial, evt_time)
  
  rm(rt)
  
  # testing
  # clock_comb$cens <- clock_comb$evt_time > clock_comb$rt_csv + clock_comb$iti_ideal
  # clock_comb %>%
  #   dplyr::filter(atlas_value == 135) %>%
  #   select(id, run, evt_time, rt_csv, iti_ideal, clock_onset, decon_interp, cens) %>%
  #   print(n = 100)

  fb_time <- 0.91 # Feedback is 0.9 seconds. Lengthen slightly to allow for small rounding/duration error (keep more evt_time samples)

  # LEFT CENSOR: About 19% of evt_time < 0 trials bleed into previous trial
  #prop.table(table(clock_comb$evt_time < -1*clock_comb$iti_prev))
  #summary(clock_comb$evt_time[clock_comb$evt_time < -1*clock_comb$iti_prev])
  clock_comb$decon_interp[clock_comb$evt_time < -1*clock_comb$iti_prev] <- NA
  clock_comb$sd_interp[clock_comb$evt_time < -1*clock_comb$iti_prev] <- NA
  
  # CLOCK VARIANT 1: right censor at onset of feedback with 100ms tolerance (so that a 900ms RT still gets evt_time 1.0)
  clock_comb_online <- clock_comb
  clock_comb_online$decon_interp[clock_comb$evt_time > clock_comb$rt_csv + .1] <- NA
  clock_comb_online$sd_interp[clock_comb$evt_time > clock_comb$rt_csv + .1] <- NA
  
  # CLOCK VARIANT 2: right censor at onset of next trial
  #prop.table(table(clock_comb$evt_time > clock_comb$rt_csv + clock_comb$iti_ideal + fb_time))
  clock_comb$decon_interp[clock_comb$evt_time > clock_comb$rt_csv + clock_comb$iti_ideal + fb_time] <- NA
  clock_comb$sd_interp[clock_comb$evt_time > clock_comb$rt_csv + clock_comb$iti_ideal + fb_time] <- NA

  # ff <- clock_comb %>% dplyr::filter(is.na(decon_interp))
  # hh <- clock_comb %>% dplyr::filter(!is.na(decon_interp))
  # 
  # xtabs(~evt_time, ff)
  # xtabs(~evt_time, hh)
  
  # code on- and offline periods and evt_time^2 for plotting models
  clock_comb <- clock_comb %>% mutate(online = evt_time > -1 & evt_time < rt_csv & evt_time < 4)
  clock_comb$online <- as.factor(clock_comb$online)

  # RIGHT CENSOR: any times after feedback + ITI period
  #prop.table(table(rt_comb$evt_time > rt_comb$iti_ideal + fb_time))
  rt_comb$decon_interp[rt_comb$evt_time > rt_comb$iti_ideal + fb_time] <- NA
  rt_comb$sd_interp[rt_comb$evt_time > rt_comb$iti_ideal + fb_time] <- NA
  
  # LEFT CENSOR: any time before the RT + previous ITI
  #prop.table(table(rt_comb$evt_time < -1*(rt_comb$iti_prev + rt_comb$rt_csv)))
  rt_comb$decon_interp[rt_comb$evt_time < -1*(rt_comb$iti_prev + rt_comb$rt_csv)] <- NA
  rt_comb$sd_interp[rt_comb$evt_time < -1*(rt_comb$iti_prev + rt_comb$rt_csv)] <- NA
  
  rt_comb <- rt_comb %>% mutate(online = evt_time > -1 & evt_time < 1.1) # treat the 0 and 1 second bins as online for RT
  rt_comb$online <- as.factor(rt_comb$online)

  # add evt_time as factor
  clock_comb$evt_time_f <- as.factor(clock_comb$evt_time)

  # both online and offline event times
  clock_wide <- clock_comb %>%
    select(id, run, run_trial, evt_time, label, decon_interp) %>%
    group_by(id, run, run_trial) %>%
    pivot_wider(names_from = c(label, evt_time), values_from = decon_interp) %>%
    ungroup()
  
  # clock_wide_cens <- clock_comb %>% select(id, run, run_trial, evt_time, label, decon_interp, online) %>%
  #   group_by(id, run, run_trial) %>% filter(evt_time < 1 | online == "TRUE") %>% select(!online) %>%
  #   pivot_wider(names_from = c(label, evt_time), values_from = decon_interp) %>% ungroup()

  clock_streams <- clock_comb %>%
    select(id, run, run_trial, evt_time, decon_interp, stream_side) %>%
    group_by(id, run, run_trial) %>%
    pivot_wider(names_from = c(stream_side, evt_time), values_from = decon_interp, values_fn=mean) %>%
    ungroup()
  
  # Note: rather than keeping parcel-level details, values_fn=mean will take the average of decon_interp for
  # each combination of visuomotor_side and evt_time, thereby generating a (subject x trials) x (visuomotor side x evt_time) matrix
  clock_visuomotor <- clock_comb %>%
    select(id, run, run_trial, evt_time, visuomotor_side, decon_interp) %>%
    pivot_wider(names_from = c(visuomotor_side, evt_time), values_from = decon_interp, values_fn = mean) %>%
    ungroup()
  
  # format for mixed_by (long)
  clock_visuomotor_long <- clock_comb %>%
    select(id, run, run_trial, evt_time, visuomotor_side, decon_interp) %>%
    group_by(id, run, run_trial, evt_time, visuomotor_side) %>%
    dplyr::summarise(decon_interp = mean(decon_interp, na.rm=TRUE)) %>%
    ungroup()

  clock_visuomotor_long_online <- clock_comb_online %>%
    select(id, run, run_trial, evt_time, visuomotor_side, decon_interp) %>%
    group_by(id, run, run_trial, evt_time, visuomotor_side) %>%
    dplyr::summarise(decon_interp = mean(decon_interp, na.rm=TRUE)) %>%
    ungroup()

  # for coxme -- don't filter online times so that we can later interpolate
  clock_cox <- clock_comb %>%
    select(id, run, run_trial, evt_time, online, label, decon_interp) %>%
    group_by(id, run, run_trial) %>%
    pivot_wider(names_from = c(label), values_from = decon_interp) %>%
    ungroup()

  # save clock ----
  message("Saving to cache")
  save(clock_wide, file = file.path(cache_dir, "clock_dan_wide_ts.Rdata"))
  save(clock_comb, file = file.path(cache_dir, "clock_dan_tall_ts.Rdata"))
  save(clock_cox, file = file.path(cache_dir, "clock_dan_medusa_for_coxme.Rdata"))
  save(clock_streams, file = file.path(cache_dir, "clock_dan_streams.Rdata"))
  save(clock_visuomotor, file = file.path(cache_dir, "clock_dan_visuomotor.Rdata"))
  save(clock_visuomotor_long, file = file.path(cache_dir, "clock_dan_visuomotor_long.Rdata"))
  save(clock_visuomotor_long_online, file = file.path(cache_dir, "clock_dan_visuomotor_long_online.Rdata"))

  # save RT ----
  # take all preceding timepoints for RT_wide
  rt_comb <- rt_comb %>% arrange(id, run, run_trial, evt_time)
  rt_wide <- rt_comb %>%
    select(id, run, run_trial, evt_time, label, decon_interp) %>%
    group_by(id, run, run_trial) %>%
    pivot_wider(names_from = c(label, evt_time), values_from = decon_interp) %>%
    ungroup()
  rt_streams <- rt_comb %>%
    select(id, run, run_trial, evt_time, decon_interp, stream_side) %>%
    group_by(id, run, run_trial) %>%
    pivot_wider(names_from = c(stream_side, evt_time), values_from = decon_interp, values_fn=mean) %>%
    ungroup()
  rt_visuomotor <- rt_comb %>%
    select(id, run, run_trial, evt_time, decon_interp, visuomotor_side) %>%
    group_by(id, run, run_trial) %>%
    pivot_wider(names_from = c(visuomotor_side, evt_time), values_from = decon_interp, values_fn=mean) %>%
    ungroup()

  # format for mixed_by (long)
  rt_visuomotor_long <- rt_comb %>%
    select(id, run, run_trial, evt_time, visuomotor_side, decon_interp) %>%
    group_by(id, run, run_trial, evt_time, visuomotor_side) %>%
    dplyr::summarise(decon_interp = mean(decon_interp, na.rm = TRUE)) %>%
    ungroup()

  save(rt_wide, file = file.path(cache_dir, "rt_dan_wide_ts.Rdata"))
  save(rt_comb, file = file.path(cache_dir, "rt_dan_tall_ts.Rdata"))
  save(trial_df, file = file.path(cache_dir, "sceptic_trial_df_for_medusa.RData"))
  save(rt_streams, file = file.path(cache_dir, "rt_dan_streams.Rdata"))
  save(rt_visuomotor, file = file.path(cache_dir, "rt_dan_visuomotor.Rdata"))
  save(rt_visuomotor_long, file = file.path(cache_dir, "rt_dan_visuomotor_long.Rdata"))
}
