# centralize this for readability
# note that the VBA outputs contain a variable called 'rt_next' that is actually the *lagged* RT -- amend this right off the bat.

get_trial_data <- function(repo_directory=NULL, trials_per_run=50) {
  checkmate::assert_directory_exists(repo_directory)
  #trial_df <- read_csv(file.path(repo_directory, "fmri/data/mmclock_fmri_decay_factorize_selective_psequate_mfx_trial_statistics.csv.gz")) %>%
  
  trial_df <- read.csv(file.path(repo_directory, "fmri/data/mmclock_fmri_decay_factorize_selective_psequate_fixedparams_ffx_trial_statistics.csv.gz")) %>%
    arrange(id, run, trial) %>% # make sure trials are consecutive
    
    # group-level mutations
    mutate(
      trial = as.numeric(trial),
      run_trial = trial - (run - 1) * !!trials_per_run,
      trial_neg_inv = -1000 / run_trial,
      trial_neg_inv_sc = as.vector(scale(trial_neg_inv)),
      rt_csv = rt_csv / 1000, # respecify in terms of seconds
      rt_csv_sc = as.vector(scale(rt_csv)),
      rt_vmax_sc = as.vector(scale(rt_vmax)),
      rt_vmax = rt_vmax / 10, # put into seconds
      rewFunc = factor(rewFunc),
      rewFunc = relevel(rewFunc, ref = "DEV"), # switch DEV to reference
      outcome = case_when(
        score_csv > 0 ~ 'Reward',
        score_csv == 0 ~ "Omission")
    ) %>%

    #run-level mutations
    group_by(id, run) %>%
    dplyr::mutate(
      rt_swing = abs(c(NA, diff(rt_csv))), # compute rt_swing within run and subject
      rt_swing_lr = abs(log(rt_csv / lag(rt_csv))),
      clock_onset_prev = dplyr::lag(clock_onset, 1, order_by = "run_trial"),
      rt_next = lead(rt_csv),
      rt_next_sc = lead(rt_csv_sc),
      rt_lag = dplyr::lag(rt_csv),
      rt_lag_sc = dplyr::lag(rt_csv_sc),
      rt_lag2_sc = dplyr::lag(rt_csv_sc, 2),
      rt_lag3_sc = dplyr::lag(rt_csv_sc, 3),
      reward = factor(case_when( # redundant with outcome, but leaving here for now
        score_csv > 0 ~ "reward",
        score_csv == 0 ~ "omission",
        TRUE ~ NA_character_
      )),
      last_outcome = dplyr::lag(outcome),
      reward_lag = dplyr::lag(reward),
      iti_prev = lag(iti_ideal),
      omission_lag = lag(score_csv == 0),
      rt_vmax_lag = lag(rt_vmax),
      rt_vmax_lag_sc = lag(rt_vmax_sc),
      rt_vmax_change = rt_vmax - rt_vmax_lag,
      rt_vmax_next = lead(rt_vmax),
      rt_vmax_change_next = rt_vmax_next - rt_vmax,
      v_entropy_wi = as.vector(scale(v_entropy)),
      v_entropy_wi_lead = lead(v_entropy_wi),
      v_entropy_wi_change = v_entropy_wi_lead - v_entropy_wi,
      entropy_split = case_when(
        v_entropy_wi > mean(v_entropy_wi, na.rm=TRUE) ~ "high",
        v_entropy_wi < mean(v_entropy_wi, na.rm=TRUE) ~ "low",
        TRUE ~ NA_character_
      ),
      entropy_split_lag = lag(entropy_split),
      rt_change = rt_csv - rt_lag,
      rt_above_1s = rt_csv > 1,
      swing_above_median = as.factor(abs(rt_change) > median(abs(na.omit(rt_change)))),
      next_swing_above_median = lead(swing_above_median),
      pe_max_lag = lag(pe_max),
      pe_max_lag2 = lag(pe_max_lag),
      pe_max_lag3 = lag(pe_max_lag2),
      abs_pe_max_lag = abs(pe_max_lag),
      abs_pe_f = case_when(
        abs(pe_max) > mean(abs(pe_max)) ~ "high abs. PE",
        abs(pe_max) < mean(abs(pe_max)) ~ "low abs. PE",
        TRUE ~ NA_character_
      ),
      abs_pe = abs(pe_max),
      rt_vmax_change = rt_vmax - rt_vmax_lag,
      feedback_onset_prev = lag(feedback_onset),
      v_max_above_median = v_max > median(na.omit(v_max)),
      v_max_wi = as.vector(scale(v_max)),
      rt_bin = case_when(
        rt_csv >= 0 & rt_csv <= 1 ~ "0-1s",
        rt_csv > 1 & rt_csv <= 2 ~ "1-2s",
        rt_csv > 2 & rt_csv <= 3 ~ "2-3s",
        rt_csv > 3 & rt_csv <= 4 ~ "3-4s",
        TRUE ~ NA_character_
      ),
      first10 = run_trial < 11,
      rt_vmax_cum = clock_onset + rt_vmax,
      rt_vmax_cum_lag <- lag(rt_vmax_cum)
    ) %>%
    ungroup() %>%
    group_by(id) %>%
    mutate(total_earnings = sum(score_csv)) %>%
    ungroup()

  
  u_df <- read_csv(file.path(repo_directory, "fmri/data/mmclock_fmri_fixed_uv_ureset_fixedparams_fmri_ffx_trial_statistics.csv.gz")) %>%
    dplyr::select(
      id, run, trial, u_chosen, u_chosen_quantile, u_chosen_lag,
      u_chosen_quantile_lag, u_chosen_change, u_chosen_quantile_change
    )

  trial_df <- inner_join(trial_df, u_df, by=c("id", "run", "trial"))

  # load full entropy and RT_vmax
  full <- read_csv(file.path(repo_directory, "fmri/data/mmclock_fmri_fixed_fixedparams_fmri_ffx_trial_statistics.csv.gz"))
  full <- as_tibble(full) %>%
    dplyr::select(c(id, run, trial, v_entropy, rt_vmax, pe_max)) %>%
    mutate(rt_vmax = rt_vmax / 10) %>% # put into seconds
    dplyr::rename(v_entropy_full = v_entropy, rt_vmax_full = rt_vmax, pe_max_full = pe_max) %>%
    group_by(id, run) %>%
    mutate(
      rt_vmax_lag_full = lag(rt_vmax_full),
      rt_vmax_change_full = rt_vmax_full - rt_vmax_lag_full,
      rt_vmax_next_full = lead(rt_vmax_full),
      rt_vmax_change_next_full = rt_vmax_next_full - rt_vmax_full,
      v_entropy_wi_full = as.vector(scale(v_entropy_full)),
      v_entropy_wi_lead_full = lead(v_entropy_wi_full),
      v_entropy_wi_change_full = v_entropy_wi_lead_full - v_entropy_wi_full,
    ) %>%
    ungroup()

  trial_df <- inner_join(trial_df, full, by = c("id", "run", "trial"))

  params <- read_csv(file.path(repo_directory, "fmri/data/mmclock_fmri_decay_factorize_selective_psequate_mfx_sceptic_global_statistics.csv"))

  trial_df <- inner_join(trial_df, params, by = c("dataset", "id"))

  return(trial_df)
}
