# centralize this for readability
# note that the VBA outputs contain a variable called 'rt_next' that is actually the *lagged* RT -- amend this right off the bat.


get_kldsum <- function(v1, v2) {
  require(LaplacesDemon)
  stopifnot(length(v1) == length(v2))
  if (any(is.na(v1)) || any(is.na(v2))) { return(NA_real_) }
  kk <- KLD(v1, v2)
  return(kk$sum.KLD.px.py)
}

# get_kldsum(1:10, 2:11)
# get_kldsum(c(1.2, 1.5, 2.5), c(0.6, 1.2, 1.5))
# get_kldsum(c(1.2, 1.5, 2.5), c(2.5, 1.2, 1.5))
# get_kldsum(c(1.2, 1.5, 2.5), c(0, 1.5, 2.5))


get_trial_data <- function(repo_directory=NULL, dataset="mmclock_fmri", groupfixed=TRUE) {
  checkmate::assert_directory_exists(repo_directory)
  
  if (dataset=="mmclock_fmri") {
    trials_per_run <- 50
    full <- read_csv(file.path(repo_directory, "fmri/data/mmclock_fmri_fixed_fixedparams_fmri_ffx_trial_statistics.csv.gz"))
    u_df <- read_csv(file.path(repo_directory, "fmri/data/mmclock_fmri_fixed_uv_ureset_fixedparams_fmri_ffx_trial_statistics.csv.gz"))
    
    if (isTRUE(groupfixed)) {
      trial_df <- read.csv(file.path(repo_directory, "fmri/data/mmclock_fmri_decay_factorize_selective_psequate_fixedparams_ffx_trial_statistics.csv.gz"))
    } else {
      trial_df <- read_csv(file.path(repo_directory, "fmri/data/mmclock_fmri_decay_factorize_selective_psequate_mfx_trial_statistics.csv.gz"))
    }
  } else if (dataset == "mmclock_meg") {
    trials_per_run <- 63
    full <- read_csv(file.path(repo_directory, "meg/data/mmclock_meg_fixed_fixedparams_meg_ffx_trial_statistics.csv.gz"))
    u_df <- read_csv(file.path(repo_directory, "meg/data/mmclock_meg_fixed_uv_ureset_fixedparams_meg_ffx_trial_statistics.csv.gz"))
    
    if (isTRUE(groupfixed)) {
      trial_df <- read.csv(file.path(repo_directory, "meg/data/mmclock_meg_decay_factorize_selective_psequate_fixedparams_meg_ffx_trial_statistics.csv.gz"))
    } else {
      trial_df <- read_csv(file.path(repo_directory, "meg/data/mmclock_meg_decay_factorize_selective_psequate_mfx_trial_statistics.csv.gz"))
    }
    
    trial_df <- trial_df %>% dplyr::rename(clock_onset = starttime) %>%
      mutate(iti_ideal = 0, feedback_onset = 0)
  } else {
    stop("Don't know how to interpret dataset")
  }
  
  trial_df <- trial_df %>%
  
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
        score_csv == 0 ~ "Omission"),
      rew_om = as.numeric(score_csv > 0), # 1/0 representation for fMRI GLMs
      rew_om_c = rew_om - 0.5 # effect coding variant for interactions
    ) %>%

    #run-level mutations
    group_by(id, run) %>%
    dplyr::mutate(
      rt_swing = abs(c(NA, diff(rt_csv))), # compute rt_swing within run and subject
      rt_swing_lr = abs(log(rt_csv / lag(rt_csv))),
      clock_onset_prev = dplyr::lag(clock_onset, 1, order_by = run_trial),
      rt_next = lead(rt_csv),
      rt_next_sc = lead(rt_csv_sc),
      rt_lag = dplyr::lag(rt_csv, order_by=run_trial),
      rt_lag2 = dplyr::lag(rt_lag, order_by=run_trial),
      rt_lag3 = dplyr::lag(rt_lag2, order_by=run_trial),
      rt_lag4 = dplyr::lag(rt_lag3, order_by=run_trial),
      rt_lag5 = dplyr::lag(rt_lag4, order_by=run_trial),
      rt_lag_sc = dplyr::lag(rt_csv_sc),
      rt_lag2_sc = dplyr::lag(rt_csv_sc, 2),
      rt_lag3_sc = dplyr::lag(rt_csv_sc, 3),
      reward = factor(case_when( # redundant with outcome, but leaving here for now
        score_csv > 0 ~ "reward",
        score_csv == 0 ~ "omission",
        TRUE ~ NA_character_
      )),
      last_outcome = dplyr::lag(outcome, order_by=run_trial),
      outcome_lag = dplyr::lag(outcome, order_by=run_trial), # synonym, but last_outcome is used in some code...
      outcome_lag2 = dplyr::lag(outcome_lag),
      reward_lag = dplyr::lag(reward),
      iti_prev = lag(iti_ideal),
      outcome_lag = dplyr::lag(outcome, order_by=run_trial),
      omission_lag = dplyr::lag(score_csv == 0, order_by=run_trial),
      omission_lag2 = dplyr::lag(omission_lag),
      omission_lag3 = dplyr::lag(omission_lag2),
      omission_lag4 = dplyr::lag(omission_lag3),
      omission_lag5 = dplyr::lag(omission_lag4),
      rt_vmax_lag = lag(rt_vmax),
      rt_vmax_lag_sc = lag(rt_vmax_sc),
      rt_vmax_change = rt_vmax - rt_vmax_lag,
      rt_vmax_next = lead(rt_vmax),
      rt_vmax_change_next = rt_vmax_next - rt_vmax,
      v_entropy_wi = as.vector(scale(v_entropy)),
      v_entropy_wi_lead = lead(v_entropy_wi),
      v_entropy_wi_change = v_entropy_wi_lead - v_entropy_wi, # change in entropy after feedback on this trial (good for RT alignment)
      v_entropy_wi_change_lag = lag(v_entropy_wi_change), # change in entropy following trial - 1 update -- good for clock alignment?
      v_entropy_wi_change_lag2 = lag(v_entropy_wi_change_lag), # change in entropy following tria - 2 update
      
      v_entropy_lag = dplyr::lag(v_entropy, 1, order_by=run_trial), # lagged unscaled entropy
      #v_entropy_change_old = v_entropy - v_entropy_lag, #change in entropy calculation used in get_mmy3_trial_df -- reflects change update after outcome on t-1
      v_entropy_lead = lead(v_entropy),
      v_entropy_change = v_entropy_lead - v_entropy, # no run z-scoring
      v_entropy_change_lag = lag(v_entropy_change), # SAME AS v_entropy - v_entropy_lag
      v_entropy_change_pos = v_entropy_change*(v_entropy_change >= 0),
      v_entropy_change_neg = abs(v_entropy_change*(v_entropy_change < 0)),
      v_entropy_wi_change_pos = v_entropy_wi_change*(v_entropy_wi_change >= 0),
      v_entropy_wi_change_neg = abs(v_entropy_wi_change*(v_entropy_wi_change < 0)), # always scale upright so that 'hot' blobs indicate sensitivity to decrease
      
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
      abs_pe = abs(pe_max),
      abs_pe_c = abs_pe - mean(abs_pe, na.rm=TRUE), # run-centered absolute PE
      abspexrew = abs_pe_c*rew_om_c, # interaction of abs PE and centered reward (0.5/-0.5 coding)
      abs_pe_lag = lag(abs_pe),
      abs_pe_f = case_when(
        abs_pe > mean(abs_pe) ~ "high abs. PE",
        abs_pe < mean(abs_pe) ~ "low abs. PE",
        TRUE ~ NA_character_
      ),
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
      rt_vmax_cum_lag = lag(rt_vmax_cum)
    ) %>%
    ungroup() %>%
    group_by(id) %>%
    mutate(total_earnings = sum(score_csv)) %>%
    ungroup()
  
  # handle KLD calculations
  trial_df <- trial_df %>% rowwise() %>% 
    mutate(
      kld4 = get_kldsum(c(rt_lag4, rt_lag3, rt_lag2, rt_lag), c(rt_lag5, rt_lag4, rt_lag3, rt_lag2)),
      kld3 = get_kldsum(c(rt_lag3, rt_lag2, rt_lag), c(rt_lag4, rt_lag3, rt_lag2))
    ) %>%
    ungroup() %>% 
    group_by(id, run) %>%
    mutate(
      kld3_lag = dplyr::lag(kld3, order_by=run_trial),
      kld4_lag = dplyr::lag(kld4, order_by=run_trial),
      kld3_cum2 = kld3 + kld3_lag,
      kld4_cum2 = kld4 + kld4_lag,
      log_kld3 = log(kld3 + .00001), # to handle large positive skew
      log_kld3_cum2 = log(kld3_cum2 + .00001)
    ) %>% ungroup()
  
  u_df <- u_df %>%
    dplyr::select(
      id, run, trial, u_chosen, u_chosen_quantile, u_chosen_lag,
      u_chosen_quantile_lag, u_chosen_change, u_chosen_quantile_change
    )

  trial_df <- inner_join(trial_df, u_df, by=c("id", "run", "trial"))

  # load full entropy and RT_vmax
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

  if (dataset=="mmclock_meg") {
    trial_df <- trial_df %>% 
      tidyr::separate(id, sep="_", into=c("id", "date")) %>%
      mutate(Subject=as.integer(id)) %>%
      dplyr::select(-feedback_onset, -iti_ideal)
  }
  
  # params <- read_csv(file.path(repo_directory, "fmri/data/mmclock_fmri_decay_factorize_selective_psequate_mfx_sceptic_global_statistics.csv")) %>%
  #   dplyr::select(-model)
  # 
  # trial_df <- inner_join(trial_df, params, by = c("dataset", "id"))

  return(trial_df)
}
