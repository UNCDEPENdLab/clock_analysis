# reads in data
# normally called from medusa_event_locked_lmer_dan.R
# if running separately, uncomment lines 4-8
replicate_compression <- FALSE
# reprocess <- T
# if (!reprocess) {
#   wide_only <- F
#   tall_only <- T
# }

# medusa_dir = "~/Box/SCEPTIC_fMRI/dan_medusa/"
medusa_dir <- "/Volumes/GoogleDrive/My Drive/SCEPTIC_fMRI/dan_medusa"
if (replicate_compression) {
  medusa_dir <- "~/Box/SCEPTIC_fMRI/dan_medusa/new_locking_compression"
}
# cache_dir = "~/Box/SCEPTIC_fMRI/dan_medusa/cache"
cache_dir <- file.path(medusa_dir, "cache")
# repo_directory <- "~/code/clock_analysis"
repo_directory <- "~/Data_Analysis/clock_analysis"

source(file.path(repo_directory, "fmri/keuka_brain_behavior_analyses/dan/get_trial_data.R"))

# reprocess = T # for troubleshooting only
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

  if (replicate_compression) {
    clock <- clock %>%
      mutate(
        run_trial = trial - (run - 1) * 50,
        decon_interp = decon_mean,
        sd_interp = decon_sd
      ) %>%
      ungroup()
  }
  clock <- clock %>% arrange(id, run, run_trial, evt_time)

  # load RT ----
  if (!replicate_compression) {
    rt <- as_tibble(read_csv(file.path(medusa_dir, "Schaefer_DorsAttn_2.3mm_rt_long_decon_locked.csv.gz")))
    rt$atlas_value <- as.character(rt$atlas_value)
    if (replicate_compression) {
      rt <- rt %>%
        mutate(
          run_trial = trial - (run - 1) * 50,
          decon_interp = decon_mean
        ) %>%
        ungroup()
    }
    rt <- rt %>% arrange(id, run, run_trial, evt_time)
  }

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
  if (!replicate_compression) {
    rt <- merge(rt, labels)
  }
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
      gamma, total_earnings, u_chosen, u_chosen_lag, u_chosen_change
    ) %>%
    mutate(rew_om = if_else(score_csv > 0, "rew", "om")) %>%
    group_by(id, run) %>%
    mutate(iti_prev = dplyr::lag(iti_ideal, by = "run_trial")) %>%
    ungroup() %>%
    inner_join(clock) %>%
    arrange(id, run, run_trial, evt_time)

  if (!replicate_compression) {
    rt_comb <- trial_df %>%
      select(
        id, run, run_trial, iti_ideal, score_csv, feedback_onset, feedback_onset_prev, rt_lag, rewFunc,
        swing_above_median, first10, reward, reward_lag, rt_above_1s, rt_bin, rt_csv, entropy_split, entropy_split_lag, abs_pe_f, rt_vmax_lag, rt_vmax_change,
        gamma, total_earnings, ev, next_swing_above_median, u_chosen, u_chosen_lag, u_chosen_change, rt_vmax_change_next
      ) %>%
      mutate(rew_om = if_else(score_csv > 0, "rew", "om")) %>%
      group_by(id, run) %>%
      mutate(iti_prev = dplyr::lag(iti_ideal, by = "run_trial")) %>%
      ungroup() %>%
      inner_join(rt) %>%
      arrange(id, run, run_trial, evt_time)
  }

  # testing
  # clock_comb$cens <- clock_comb$evt_time > clock_comb$rt_csv + clock_comb$iti_ideal
  # clock_comb %>%
  #   dplyr::filter(atlas_value == 135) %>%
  #   select(id, run, evt_time, rt_csv, iti_ideal, clock_onset, decon_interp, cens) %>%
  #   print(n = 100)

  # 20% of clock- and 32% of feedback-aligned timepoints are from the next trial: censor
  clock_comb$decon_interp[clock_comb$evt_time > clock_comb$rt_csv + clock_comb$iti_ideal] <- NA
  clock_comb$sd_interp[clock_comb$evt_time > clock_comb$rt_csv + clock_comb$iti_ideal] <- NA

  # code on- and offline periods and evt_time^2 for plotting models
  clock_comb <- clock_comb %>% mutate(online = evt_time > -1 & evt_time < rt_csv & evt_time < 4)
  clock_comb$online <- as.factor(clock_comb$online)

  # fb_comb$decon_interp[fb_comb$evt_time > fb_comb$iti_ideal] <- NA
  # fb_comb$sd_interp[fb_comb$evt_time > fb_comb$iti_ideal] <- NA
  if (!replicate_compression) {
    rt_comb$decon_interp[rt_comb$evt_time + 0.9 > rt_comb$iti_ideal] <- NA # feedback is 0.9s fixed duration
    rt_comb$sd_interp[rt_comb$evt_time + 0.9 > rt_comb$iti_ideal] <- NA
    rt_comb <- rt_comb %>% mutate(online = evt_time > -1 & evt_time < 1.1) # treat the 0 and 1 second bins as online for RT
    rt_comb$online <- as.factor(rt_comb$online)
  }

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

  # for coxme -- don't filter online times so that we can later interpolate
  clock_cox <- clock_comb %>%
    select(id, run, run_trial, evt_time, online, label, decon_interp) %>%
    group_by(id, run, run_trial) %>%
    pivot_wider(names_from = c(label), values_from = decon_interp) %>%
    ungroup()

  # save clock ----
  if (!replicate_compression) {
    message("Saving to cache")
    save(clock_wide, file = file.path(cache_dir, "clock_dan_wide_ts.Rdata"))
    save(clock_comb, file = file.path(cache_dir, "clock_dan_tall_ts.Rdata"))
    save(clock_cox, file = file.path(cache_dir, "clock_dan_medusa_for_coxme.Rdata"))
    save(clock_streams, file = file.path(cache_dir, "clock_dan_streams.Rdata"))
    save(clock_visuomotor, file = file.path(cache_dir, "clock_dan_visuomotor.Rdata"))
    save(clock_visuomotor_long, file = file.path(cache_dir, "clock_dan_visuomotor_long.Rdata"))

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
}
