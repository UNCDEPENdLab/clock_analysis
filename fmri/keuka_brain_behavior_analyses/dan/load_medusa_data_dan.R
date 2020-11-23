# reads in data

medusa_dir = "~/Box/SCEPTIC_fMRI/dan_medusa"
cache_dir = "~/Box/SCEPTIC_fMRI/dan_medusa/cache"
repo_directory <- "~/code/clock_analysis"

# reprocess = T
if (!exists("reprocess") || !is.logical(reprocess)) { reprocess=FALSE } #default

stopifnot(dir.exists(medusa_dir))  
stopifnot(dir.exists(cache_dir))  

cwd <- getwd()
setwd(medusa_dir)

if (!reprocess) {
  message("Loading medusa data from cache: ", cache_dir)
  
  load(file.path(cache_dir, 'clock_dan_wide_ts.Rdata'))
  load(file.path(cache_dir, 'feedback_dan_wide_ts.Rdata'))
  # load(file.path(cache_dir, "feedback_hipp_tall.Rdata"))
  load(file.path(cache_dir, 'feedback_dan_tall_ts.Rdata'))
  load(file.path(cache_dir, "sceptic_trial_df_for_medusa.RData"))
} else { 
  # load clock
  clock <- as_tibble(read_csv("Schaefer_DorsAttn_2.3mm_clock_long_decon_locked.csv.gz"))
  clock$atlas_value <- as.character(clock$atlas_value)
  clock <- clock %>% arrange(id, run, run_trial, evt_time)
  # load feedback
  fb <- as_tibble(read_csv("Schaefer_DorsAttn_2.3mm_feedback_long_decon_locked.csv.gz"))
  fb$atlas_value <- as.character(fb$atlas_value)
  fb <- fb %>% arrange(id, run, run_trial, evt_time)
  
}

# add manual labels for Schaeffer areas
labels <- as_tibble(read_delim("~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/Schaefer2018_200Parcels_7Networks_order_manual.txt", 
                               "\t", escape_double = FALSE, col_names = FALSE, 
                               trim_ws = TRUE)) %>% select(1:4)
names(labels) <- c("atlas_value", "label_long", "label_short", "entropy_signal")
labels <- labels %>% filter(grepl("DorsAtt", label_long))  %>% mutate(
  side  = case_when(
    grepl("LH", label_long) ~ "L",
    grepl("RH", label_long) ~ "R"
  ),
  label_short_side = paste(label_short, side, sep ="_"),
  label_long1 = substr(label_long, 23, 100),
  label = case_when(
    label_short!="0" ~ label_short_side,
    label_short=="0" ~ label_long1
  )
) %>% select(c(label, side, atlas_value))

fb <- merge(fb, labels)
clock <- merge(clock, labels)
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


# add parameters
params <- read_csv(file.path(repo_directory, "fmri/data/mmclock_fmri_decay_factorize_selective_psequate_mfx_sceptic_global_statistics.csv"))

trial_df <- inner_join(trial_df, params, by = "id")

message("Transforming for MEDUSA")

clock_comb <- trial_df %>% select(id, run, run_trial, iti_ideal, score_csv, clock_onset, clock_onset_prev, rt_lag, rewFunc,
                                  swing_above_median, first10,reward, reward_lag, rt_above_1s, rt_bin, rt_csv, v_entropy_wi, entropy, entropy_lag, 
                                  gamma, total_earnings, u_chosen, u_chosen_lag, u_chosen_change) %>% mutate(rewom=if_else(score_csv > 0, "rew", "om")) %>%
  group_by(id, run) %>% mutate(iti_prev=dplyr::lag(iti_ideal, by="run_trial")) %>% ungroup() %>%
  inner_join(clock)
fb_comb <- trial_df %>% select(id, run, run_trial, iti_ideal, score_csv, feedback_onset, feedback_onset_prev, rt_lag, rewFunc,
                               swing_above_median, first10,reward, reward_lag, rt_above_1s, rt_bin, rt_csv, entropy, entropy_lag, abs_pe_f, rt_vmax_lag, rt_vmax_change,
                               gamma, total_earnings, ev,next_swing_above_median, u_chosen, u_chosen_lag, u_chosen_change) %>% mutate(rewom=if_else(score_csv > 0, "rew", "om")) %>%
  group_by(id, run) %>% mutate(iti_prev=dplyr::lag(iti_ideal, by="run_trial")) %>% ungroup() %>%
  inner_join(fb)

# 20% of clock- and 32% of feedback-aligned timepoints are from the next trial: censor
clock_comb$decon_interp[clock_comb$evt_time > clock_comb$rt_csv + clock_comb$iti_ideal] <- NA
clock_comb$sd_interp[clock_comb$evt_time > clock_comb$rt_csv + clock_comb$iti_ideal] <- NA
fb_comb$decon_interp[fb_comb$evt_time > fb_comb$iti_ideal] <- NA
fb_comb$sd_interp[fb_comb$evt_time > fb_comb$iti_ideal] <- NA

# code on- and offline periods and evt_time^2 for plotting models
clock_comb <- clock_comb %>% mutate(online = evt_time > -1 & evt_time < rt_csv & evt_time<4)
clock_comb$online <- as.factor(clock_comb$online)

# use more stringent feedback online window
fb_comb <- fb_comb %>% mutate(online = evt_time > -rt_csv & evt_time < 0)
fb_comb$online <- as.factor(fb_comb$online)
fb_comb$evt_time_sq <- fb_comb$evt_time^2


# add evt_time as factor
clock_comb$evt_time_f <- as.factor(clock_comb$evt_time)
fb_comb$evt_time_f <- as.factor(fb_comb$evt_time)
# rtvmax_comb$evt_time_f <- as.factor(rtvmax_comb$evt_time)

# # lags -- these take very long
# clock_comb <- clock_comb %>% group_by(id, run, side, axis_bin, evt_time) %>%
#   mutate(decon_prev = dplyr::lag(decon_interp, order_by = run_trial),
#          telapsed=clock_onset - clock_onset_prev
#   ) %>%
#   ungroup() %>%
#   mutate(decon_prev_z=as.vector(scale(decon_prev)), iti_ideal_z=as.vector(scale(iti_ideal)))
# 
# 
# fb_comb <- fb_comb %>% group_by(id, run, axis_bin, side, evt_time) %>%
#   mutate(decon_prev = dplyr::lag(decon_interp, order_by = run_trial),
#          telapsed=feedback_onset - feedback_onset_prev) %>%
#   ungroup() %>%
#   mutate(decon_prev_z=as.vector(scale(decon_prev)), iti_ideal_z=as.vector(scale(iti_ideal)))
# 
# rtvmax_comb <- rtvmax_comb %>% group_by(id, run, axis_bin, side, evt_time) %>%
#   mutate(decon_prev = dplyr::lag(decon_interp, order_by = run_trial),
#          telapsed=clock_onset - clock_onset_prev) %>%
#   ungroup() %>%
#   mutate(decon_prev_z=as.vector(scale(decon_prev)), iti_ideal_z=as.vector(scale(iti_ideal)))

# myspread <- function(df, key, value) {
#   # quote key
#   keyq <- rlang::enquo(key)
#   # break value vector into quotes
#   valueq <- rlang::enquo(value)
#   s <- rlang::quos(!!valueq)
#   df %>% gather(variable, value, !!!s) %>%
#     unite(temp, !!keyq, variable) %>%
#     spread(temp, value)
# }

# fb_comb <- fb_comb %>% group_by(id,run,run_trial,evt_time,side) %>% mutate(bin_num = rank(bin_center)) %>% ungroup()
# fb_comb <- fb_comb %>% mutate(bin6 = round((bin_num + .5)/2,digits = 0)) # also a 6-bin version
# fb_wide <- fb_comb %>% select(id, run, run_trial, evt_time, label, decon_interp) %>% spread(key = label, decon_interp) #%>% myspread(bin_num, c("l", "r"))
# names(fb_wide)[5:28] <- paste("hipp", names(fb_wide)[5:28], sep = "_")
# fb_wide_ex <- inner_join(fb_wide, trial_df[,c("id", "run", "run_trial", "pe_max", "reward", "v_entropy_wi")], by = c("id", "run", "run_trial"))
# fb_wide6 <- fb_comb %>% select(id, run, run_trial, evt_time, side, bin6, decon_interp) %>% group_by(id, run, run_trial, evt_time, side, bin6) %>% summarise(decon6  = mean(decon_interp)) %>% spread(key = side, decon6) %>% myspread(bin6, c("l", "r"))
# names(fb_wide6)[5:length(names(fb_wide6))] <- paste("hipp", names(fb_wide6)[5:length(names(fb_wide6))], sep = "_")
# fb_wide6_ex <- inner_join(fb_wide6, trial_df[,c("id", "run", "run_trial", "pe_max", "reward", "v_entropy_wi")], by = c("id", "run", "run_trial"))

# make wide fb df with side as observation to examine laterality

fb_wide <- fb_comb %>% select(id, run, run_trial, evt_time, label, decon_interp) %>%  
  pivot_wider(id_cols = c(id, run, run_trial), names_from = c(label, evt_time), values_from = decon_interp) 

# names(fb_wide)[4:148] <- paste("dan", names(fb_wide)[4:148], sep = "_")


# save fb responses for trial-wise prediction analyses

# slices <- names(fb_wide)[4:28]
# # fb_wide_t <- dcast(setDT(fb_wide), id + run + run_trial ~ evt_time, value.var = slices)
# fb_wide_t <- fb_wide %>% pivot_wider(names_from = evt_time, values_from = slices)

message("Saving to cache")
save(fb_wide, file = file.path(cache_dir, "feedback_dan_wide_ts.Rdata"))
save(fb_comb, file = file.path(cache_dir, "feedback_dan_tall_ts.Rdata"))
# save(fb_wide_t, fb_wide_bl, file = file.path(cache_dir, 'feedback_hipp_widest_by_timepoint_decon.Rdata'))
# rtvmax_comb <- rtvmax_comb %>% group_by(id,run,run_trial,evt_time,side) %>% mutate(bin_num = rank(bin_center)) %>% ungroup()
# 
# save(rtvmax_comb, file = file.path(cache_dir, "rtvmax_hipp_tall_ts.Rdata"))

# clock_comb <- clock_comb %>% group_by(id,run,run_trial,evt_time,side) %>% mutate(bin_num = rank(bin_center)) %>% ungroup()

# take only online event times
clock_wide <- clock_comb %>% filter(online==T) %>% select(id, run, run_trial, evt_time, label, decon_interp) %>% 
  group_by(id, run, run_trial) %>%
  pivot_wider(names_from = c(label, evt_time), values_from = decon_interp)

clock_cox <- clock_comb %>% filter(online==T) %>% select(id, run, run_trial, evt_time, label, decon_interp) %>% 
  group_by(id, run, run_trial) %>%
  pivot_wider(names_from = c(label), values_from = decon_interp)
# names(clock_wide)[4:length(names(clock_wide))] <- paste("dan", names(clock_wide)[4:length(names(clock_wide))], sep = "_")
# clock_wide_t <- clock_wide %>% pivot_wider(names_from = evt_time, values_from = slices)
# clock_wide_ex <- inner_join(clock_wide, trial_df[,c("id", "run", "run_trial", "pe_max", "reward", "v_entropy_wi", "swing_above_median")], by = c("id", "run", "run_trial"))
# 
save(clock_wide,  file = file.path(cache_dir, "clock_dan_wide_ts.Rdata"))
save(clock_comb, file = file.path(cache_dir, "clock_dan_tall_ts.Rdata"))
save(clock_cox, file = file.path(cache_dir, "clock_dan_medusa_for_coxme.Rdata"))
save(trial_df, file=file.path(cache_dir, "sceptic_trial_df_for_medusa.RData"))


#reset working directory
setwd(cwd)

