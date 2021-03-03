# reads in meg data
# normally called from ...
# if running separately, uncomment lines 4-9
# reprocess = T
# if (!reprocess) {
#   wide_only = F
#   tall_only = T
# }

medusa_dir = "~/Box/SCEPTIC_fMRI/MEG_20Hz/"
cache_dir = "~/Box/SCEPTIC_fMRI/MEG_20Hz/cache"
repo_directory <- "~/code/clock_analysis"

# Kai’s guidance on sensors is: ‘So for FEF, I say focus on 612/613, 543/542, 1022/1023, 
# For IPS, 1823, 1822, 2222,2223.’
fef_sensors <- c("0612","0613", "0542", "0543","1022")
ips_sensors <- c("1823", "1822", "2222","2223")
all_sensors <- c(fef_sensors,ips_sensors)

# reprocess = T # for troubleshooting only
if (!exists("reprocess") || !is.logical(reprocess)) { reprocess=FALSE } #default

stopifnot(dir.exists(medusa_dir))  
stopifnot(dir.exists(cache_dir))  

cwd <- getwd()
setwd(medusa_dir)

if (!reprocess) {
  message("Loading MEG medusa data from cache: ", cache_dir)
  load(file.path(cache_dir, "meg_sceptic_trial_df_for_medusa.RData"))
  if (wide_only) {
    # load(file.path(cache_dir, 'clock_dan_wide_ts.Rdata'))
    # load(file.path(cache_dir, 'feedback_dan_wide_ts.Rdata'))
    load(file.path(cache_dir, "meg_rt_dan_wide_ts.Rdata"))  } else if (tall_only) {
    # load(file.path(cache_dir, 'clock_dan_tall_ts.Rdata'))
    load(file.path(cache_dir, 'meg_rt_dan_tall_ts.Rdata'))
  } else {
    # load(file.path(cache_dir, 'meg_clock_dan_wide_ts.Rdata'))
    # load(file.path(cache_dir, 'meg_clock_dan_tall_ts.Rdata'))
    # load(file.path(cache_dir, 'feedback_dan_wide_ts.Rdata'))
    load(file.path(cache_dir, 'rt_dan_wide_ts.Rdata'))
    load(file.path(cache_dir, 'rt_dan_tall_ts.Rdata'))
  } 
} else {
  # load rt-aligned "ips" files ----
  rt <- data.frame(matrix(ncol = 5, nrow = 0))
  colnames(rt) <- colnames(as_tibble(readRDS(paste0("MEG", ips_sensors[1], "_20Hz.rds"))))
  for (sensor in all_sensors) {
    c <- as_tibble(readRDS(paste0("MEG", sensor, "_20Hz.rds")))
    c$sensor <- as.character(sensor)
    rt <- rbind(rt, c)
  }
  rt <- rt %>% rename(id = Subject, trial = Trial, run = Run, evt_time = Time, signal = Signal)
  
  # # add manual labels at some point ----
  # load trial-level data ----
  message("Loading trial-level data")
  trial_df <-  read_csv("~/Box/SCEPTIC_fMRI/sceptic_model_fits/mmclock_meg_decay_factorize_selective_psequate_fixedparams_meg_ffx_trial_statistics.csv.gz") %>%
    mutate(trial=as.numeric(trial)) %>%
    group_by(id, run) %>%  
    dplyr::mutate(rt_swing = abs(c(NA, diff(rt_csv))), #compute rt_swing within run and subject
                  rt_swing_lr = abs(log(rt_csv/lag(rt_csv))),
                  # clock_onset_prev=dplyr::lag(clock_onset, 1, by="run"),
                  rt_lag = lag(rt_csv) ,
                  reward = case_when(
                    score_csv>0 ~ "reward",
                    score_csv==0 ~ "omission",
                    TRUE ~ NA_character_),
                  reward = as.factor(reward),
                  reward_lag = lag(reward),
                  # iti_prev = lag(iti_ideal),
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
                  # feedback_onset_prev = lag(feedback_onset),
                  v_max_above_median = v_max > median(na.omit(v_max)),
                  last_outcome = reward_lag,
                  run_trial=1:63, 
                  first10  = run_trial<11) %>% ungroup() 
  trial_df <- trial_df %>% dplyr::mutate(rt_csv=rt_csv/1000, rt_lag = rt_lag/1000, rt_vmax=rt_vmax/10, rt_vmax_lag = rt_vmax_lag/10) # careful not to do this multiple times
  trial_df <- trial_df %>% group_by(id) %>% mutate(total_earnings = sum(score_csv)) %>% ungroup() %>% mutate(id = as.integer(substr(id, 1, 5)))
  trial_df$rewFunc <- as.factor(trial_df$rewFunc)
  levels(trial_df$rewFunc) <- c("DEV", "IEV", "CEV", "CEVR")
  
  # don't really trust uncertainty
  # u_df <- read_csv("~/Box/SCEPTIC_fMRI/sceptic_model_fits/mmclock_fmri_fixed_uv_ureset_mfx_trial_statistics.csv.gz")
  # u_df <- u_df %>% select(id, run, trial, u_chosen, u_chosen_lag, u_chosen_change)
  # 
  # trial_df <- inner_join(trial_df,u_df)
  
  # # load fixed entropy and RT_vmax
  # fixed <-  read_csv(file.path(repo_directory, "fmri/data/mmclock_fmri_fixed_fixedparams_fmri_ffx_trial_statistics.csv.gz"))
  # fixed <- as_tibble(fixed) %>% select(c(id, run, trial, v_entropy, rt_vmax, pe_max)) %>%
  #   rename(v_entropy_fixed = v_entropy, rt_vmax_fixed = rt_vmax, pe_max_fixed = pe_max) %>% group_by(id, run) %>%
  #   mutate(
  #     rt_vmax_lag_fixed = lag(rt_vmax_fixed),
  #     rt_vmax_change_fixed = rt_vmax_fixed - rt_vmax_lag_fixed,
  #     rt_vmax_next_fixed = lead(rt_vmax_fixed),
  #     rt_vmax_change_next_fixed = rt_vmax_next_fixed - rt_vmax_fixed,
  #     v_entropy_wi_fixed = scale(v_entropy_fixed),
  #     v_entropy_wi_lead_fixed = lead(v_entropy_wi_fixed),
  #     v_entropy_wi_change_fixed = v_entropy_wi_lead_fixed - v_entropy_wi_fixed,
  #   )
  # 
  # trial_df <- inner_join(trial_df, fixed, by = c("id", "trial", "run"))
  
  # transform for MEDUSA ----
  message("Transforming for MEDUSA")
  
  rt_comb <- trial_df %>% select(id, run, trial, run_trial, score_csv,  rt_lag, rewFunc,
                                 swing_above_median, first10,reward, reward_lag, rt_above_1s, rt_csv, entropy, entropy_lag, abs_pe_f, rt_vmax_lag, rt_vmax_change,
                                 ev,rt_vmax_change_next) %>% mutate(rewom=if_else(score_csv > 0, "rew", "om")) %>%
    group_by(id, run) %>% ungroup() %>%
    inner_join(rt) %>% arrange(id, run, run_trial, evt_time)
  
  rt_comb$evt_time_f <- as.factor(rt_comb$evt_time)
  
  message("Saving to cache")
  
  # save RT ----
  rt_comb <- rt_comb %>% arrange(id, run, run_trial, evt_time)
  rt_wide <- rt_comb %>%  select(id, run, run_trial, evt_time, sensor, signal) %>% 
    group_by(id, run, run_trial) %>%
    pivot_wider(names_from = c(sensor, evt_time), values_from = signal)
  
  save(rt_wide, file = file.path(cache_dir, "meg_rt_dan_wide_ts.Rdata"))
  save(rt_comb, file = file.path(cache_dir, "meg_rt_dan_tall_ts.Rdata"))
  save(trial_df, file=file.path(cache_dir, "meg_sceptic_trial_df_for_medusa.RData"))
}

#reset working directory
setwd(cwd)

