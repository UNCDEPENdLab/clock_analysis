library(modelr)
library(tidyverse)
library(lme4)
library(afex)
library(broom)
library(broom.mixed) #plays will with afex p-values in lmer wrapper
library(ggpubr)
library(car)
library(viridis)
library(psych)
library(corrplot)
library(foreach)
library(doParallel)
# library(data.table)
repo_directory <- "~/code/clock_analysis"
medusa_dir = "~/Box/SCEPTIC_fMRI/MEG_TimeFreq/"
decode_plot_dir = "~/OneDrive/collected_letters/papers/meg/plots/rt_decode/"
rt_plot_dir = "~/OneDrive/collected_letters/papers/meg/plots/rt_rt//"

stopifnot(dir.exists(medusa_dir))  

setwd(medusa_dir)

# options, files ----
parallel = T
decode = T  # main analysis analogous to Fig. 4 E-G in NComm 2020
rt_predict = T # predicts next response based on signal and behavioral variables
online = F # whether to analyze clock-aligned ("online") or RT-aligned ("offline") responses
exclude_first_run = F
reg_diagnostics = F
start_time = -3
end_time = .94
log = T # whether to log-transform power

# # Kai’s guidance on sensors is: ‘So for FEF, I say focus on 612/613, 543/542, 1022/1023, 
# # For IPS, 1823, 1822, 2222,2223.’
# fef_sensors <- c("0612","0613", "0542", "0543","1022")
# ips_sensors <- c("1823", "1822", "2222","2223")
# dan_sensors <- c(fef_sensors,ips_sensors)
files <- list.files(medusa_dir)
files <- files[grepl("MEG", files)]
all_sensors <- substr(files, 4,7)
all_sensors <- all_sensors[1]

# # take first few for testing
# all_sensors <- all_sensors[1:2] # TEST ONLY
scale2 <- function(x, na.rm = FALSE) (x - mean(x, na.rm = na.rm)) / sd(x, na.rm)


# get behavioral data ----
message("Loading behavioral data")
trial_df <-  read_csv("~/Box/SCEPTIC_fMRI/sceptic_model_fits/mmclock_meg_decay_factorize_selective_psequate_fixedparams_meg_ffx_trial_statistics.csv.gz") %>%
  mutate(trial=as.numeric(trial)) %>% ungroup() %>%   mutate(rt_csv_sc = scale(rt_csv),
                                                             rt_vmax_sc = scale(rt_vmax)) %>% group_by(id, run) %>%  
  dplyr::mutate(rt_swing = abs(c(NA, diff(rt_csv))), #compute rt_swing within run and subject
                rt_swing_lr = abs(log(rt_csv/lag(rt_csv))),
                rt_csv=rt_csv/1000, # careful not to do this multiple times
                rt_vmax=rt_vmax/10, 
                rt_next = lead(rt_csv),
                rt_lag = lag(rt_csv),
                rt_lag_sc = lag(rt_csv_sc),
                reward = case_when(
                  score_csv>0 ~ "Reward",
                  score_csv==0 ~ "Omission",
                  TRUE ~ NA_character_),
                reward = as.factor(reward),
                reward_lag = lag(reward),
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
                abs_pe_lag = abs(pe_max_lag), 
                abs_pe_f  = case_when(
                  abs(pe_max) > mean(abs(pe_max)) ~ 'high abs. PE',
                  abs(pe_max) < mean(abs(pe_max)) ~ 'low abs. PE',
                  TRUE ~ NA_character_),
                abs_pe = abs(pe_max),
                rt_vmax_change = rt_vmax - rt_vmax_lag,
                v_max_above_median = v_max > median(na.omit(v_max)),
                last_outcome = reward_lag,
                run_trial=1:63, 
                first10  = run_trial<11,
                rt_change = rt_next - rt_csv_sc,
                rt_vmax_lead = lead(rt_vmax),
                v_entropy_wi_lead = lead(v_entropy_wi),
                v_entropy_wi_change = v_entropy_wi_lead-v_entropy_wi,
                v_entropy_wi_change_lag = lag(v_entropy_wi_change),
                outcome = reward,
                abs_pe = abs(pe_max),
                abs_pe_lag = lag(abs_pe),
                v_max_wi = scale(v_max),
                v_max_wi_lag = lag(v_max_wi),
                v_entropy_wi = scale(v_entropy),
                v_max_b = mean(na.omit(v_max)),
                v_entropy_b = mean(na.omit(v_entropy)),
                rt_change = rt_csv - rt_lag,
                pe_max_lag = lag(pe_max), 
                v_chosen_change = v_chosen - lag(v_chosen),
                trial_neg_inv_sc = scale(-1/run_trial),
                rt_lag2_sc = lag(rt_csv_sc, 2),
                rt_lag3_sc = lag(rt_csv_sc, 3),
  ) %>% ungroup() %>% mutate(rt_vmax_lag_sc = scale(rt_vmax_lag),
                             abs_pe_sc = scale(abs_pe),
                             abs_pe_lag_sc = scale(abs_pe_lag)) %>%
  group_by(id) %>% mutate(total_earnings = sum(score_csv)) %>% ungroup() %>% mutate(id = as.integer(substr(id, 1, 5))) %>%
  select(id, run, trial, run_trial, score_csv,  rt_lag, rewFunc, rt_next,
         swing_above_median, first10,reward, reward_lag, rt_above_1s, rt_csv, entropy, entropy_lag, abs_pe_f, abs_pe, rt_vmax, rt_vmax_lag, rt_vmax_change,
         ev,rt_vmax_change_next, trial_neg_inv_sc, rt_csv_sc, rt_lag_sc, rt_vmax_lag_sc, rt_vmax_change, 
         v_entropy_wi, v_entropy_wi_lead, v_entropy_wi_change_lag, v_entropy_wi_change,
         v_max_wi, outcome) %>% mutate(rewom=if_else(score_csv > 0, "rew", "om"))
trial_df$rewFunc <- as.factor(trial_df$rewFunc)
levels(trial_df$rewFunc) <- c("DEV", "IEV", "CEV", "CEVR")


# make cluster ----
if (parallel) {
  f <- Sys.getenv('PBS_NODEFILE')
  library(parallel)
  ncores <- detectCores()
  nodelist <- if (nzchar(f)) readLines(f) else rep('localhost', ncores)
  
  cat("Node list allocated to this job\n")
  print(nodelist)
  
  cl <- makePSOCKcluster(nodelist, outfile='')
  print(cl) ##; print(unclass(cl))
  registerDoParallel(cl)
}
# loop over sensors ----

# define frequencies for the loop
test <- as_tibble(readRDS(paste0("MEG", all_sensors[1], "_tf.rds"))) %>% filter(Time>start_time) %>%
  rename(id = Subject, trial = Trial, run = Run, evt_time = Time, pow = Pow) %>%
  mutate(pow = scale2(pow)) # scale signal across subjects
freqs <- as.list(unique(test$Freq))
freqs <- freqs[1] # TEST ONLY
message("Decoding: analyzing censor data")
pb <- txtProgressBar(0, max = length(all_sensors)*length(freqs), style = 3)

# DECODING ANALYSES ----
s = 0

if(decode) {
  biglist <- list()
  for (sensor in all_sensors) {
    s = s + 1
    setwd(medusa_dir)
    rt <- as_tibble(readRDS(paste0("MEG", sensor, "_tf.rds"))) %>% filter(Time>start_time & Time < end_time) %>%
      rename(id = Subject, trial = Trial, run = Run, evt_time = Time, pow = Pow) %>%
      mutate(
        pow = log(pow + 10^-25)) %>% group_by(Freq) %>% mutate(pow = scale2(pow)) %>% ungroup()
    rt$sensor <- as.character(sensor) 
    # combine with behavior ----
    # message("Merging with behavior")
    rt_comb <- trial_df %>% 
      group_by(id, run) %>% ungroup() %>% inner_join(rt) %>% arrange(id, run, run_trial, evt_time)
    rt_comb$evt_time_f <- as.factor(rt_comb$evt_time)
    timepoints = as.character(unique(rt_comb$evt_time))
    if (s %% 1 == 0) {setTxtProgressBar(pb, s)}
    newlist <- list()
    f = 1
    for (freq in freqs) {
      f = f + 1
      rt_comb_f <- rt_comb %>% filter(Freq == freq)
      rt_comb_t <- split(rt_comb_f, rt_comb_f$evt_time)
      
      dd <- foreach(d = iter(rt_comb_t), j = icount(), .combine='rbind',.packages=c("lme4", "tidyverse", "broom.mixed", "car"), .noexport = c("rt_comb", "rt")) %dopar% {
        t <- timepoints[[j]]
        md <-  lmerTest::lmer(pow ~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + scale(rt_vmax_lag)  + scale(rt_vmax_change) + 
                                v_entropy_wi + v_entropy_wi_change  + 
                                v_max_wi  + scale(abs_pe) + outcome + 
                                (1|id), d, control=lmerControl(optimizer = "nloptwrap"))
        while (any(grepl("failed to converge", md@optinfo$conv$lme4$messages) )) {
          print(md@optinfo$conv$lme4$conv)
          ss <- getME(md,c("theta","fixef"))
          md <- update(md, start=ss, control=lmerControl(optimizer = "bobyqa"))}  
        
        dm <- tidy(md)
        dm$sensor <- sensor
        dm$t <- t
        dm$freq <- freq
        # print(dm)
        dm}
      newlist[[freq]] <- dd}
    df <-  do.call(rbind,newlist)
    biglist[[sensor]] <- df}
  
  ddf <-  do.call(rbind,biglist)
  
  # format statistics ----
  
  ddf <- ddf %>% mutate(stat_order = as.factor(case_when(abs(statistic) < 2 ~ '1', 
                                                         abs(statistic) > 2 & abs(statistic) < 3 ~ '2', 
                                                         abs(statistic) > 3 ~ '3')),
                        p_value = as.factor(case_when(`p.value` > .05 ~ '1',
                                                      `p.value` < .05 & `p.value` > .01 ~ '2',
                                                      `p.value` < .01 & `p.value` > .001 ~ '3',
                                                      `p.value` <.001 ~ '4')))
  ddf$t <- as.numeric(ddf$t)
  # ddf$label <- as.factor(sub("_[^_]+$", "", ddf$label))
  ddf$stat_order <- factor(ddf$stat_order, labels = c("NS", "|t| > 2", "|t| > 3"))
  ddf$p_value <- factor(ddf$p_value, labels = c("NS", "p < .05", "p < .01", "p < .001"))
  terms <- unique(ddf$term[ddf$effect=="fixed"])
  # FDR correction ----
  ddf <- ddf  %>% group_by(term, sensor) %>% mutate(p_fdr = p.adjust(p.value, method = 'fdr'),
                                                    p_level_fdr = as.factor(case_when(
                                                      # p_fdr > .1 ~ '0',
                                                      # p_fdr < .1 & p_fdr > .05 ~ '1',
                                                      p_fdr > .05 ~ '1',
                                                      p_fdr < .05 & p_fdr > .01 ~ '2',
                                                      p_fdr < .01 & p_fdr > .001 ~ '3',
                                                      p_fdr <.001 ~ '4'))
  ) %>% ungroup() #%>% mutate(side = substr(as.character(label), nchar(as.character(label)), nchar(as.character(label))),
  #          region = substr(as.character(label), 1, nchar(as.character(label))-2))
  ddf$p_level_fdr <- factor(ddf$p_level_fdr, levels = c('1', '2', '3', '4'), labels = c("NS","p < .05", "p < .01", "p < .001"))
  ddf$`p, FDR-corrected` = ddf$p_level_fdr
  # save model stats ----
  setwd(decode_plot_dir)
  save(file = "TESTmeg_freq_medusa_decode_output_nested.Rdata", ddf)
  
}

message("RT prediction: analyzing censor data")
pb <- txtProgressBar(0, max = length(all_sensors), style = 3)

# RT ANALYSES ----
s = 0
if(rt_predict) {
  biglist <- list()
  for (sensor in all_sensors) {
    s = s + 1
    setwd(medusa_dir)
    rt <- as_tibble(readRDS(paste0("MEG", sensor, "_tf.rds"))) %>% filter(Time>start_time & Time < end_time) %>%
      rename(id = Subject, trial = Trial, run = Run, evt_time = Time, pow = Pow) %>%
      mutate(
        pow = log(pow + 10^-25)) %>% group_by(Freq) %>% mutate(pow = scale2(pow)) %>% ungroup()
    # ggplot(rt, aes(pow)) + geom_histogram() + facet_wrap(~Freq)
    rt$sensor <- as.character(sensor) 
    ## Diagnostics on power distributions
    setwd("~/OneDrive/collected_letters/papers/meg/plots/diagnostics/")
    t <- rt %>% group_by(evt_time, Freq) %>% summarise(.groups = "keep", mean_pow = mean(pow, na.rm = T)) # %>%
    #   group_by(Freq) %>%  mutate(mean_pow = scale2(mean_pow)) %>% filter(evt_time < .95)
    pdf(paste0(sensor, "_log_power_heatmap_scaled.pdf"))
    ggplot(t, aes(evt_time, Freq,  fill = mean_pow)) + geom_tile() + scale_fill_viridis(option = "plasma")  + scale_color_grey()
    dev.off()
    # combine with behavior ----
    # message("Merging with behavior")
    rt_comb <- trial_df %>% 
      group_by(id, run) %>% ungroup() %>% inner_join(rt) %>% arrange(id, run, run_trial, evt_time)
    rt_comb$evt_time_f <- as.factor(rt_comb$evt_time)
    timepoints = as.character(unique(rt_comb$evt_time))
    if (s %% 1 == 0) {setTxtProgressBar(pb, s)}
    newlist <- list()
    
    for (freq in freqs) {
      rt_comb_f <- rt_comb %>% filter(Freq == freq)
      rt_comb_t <- split(rt_comb_f, rt_comb_f$evt_time)
      message(paste(s, freq))
      
      dd <- foreach(d = iter(rt_comb_l), j = icount(), .combine='rbind',.packages=c("lme4", "tidyverse", "broom.mixed", "car"), .noexport = c("rt_comb", "rt")) %dopar% {
        t <- timepoints[[j]]
        # load data ----
        # message("Loading")
        # wrangle into wide -----
        # message("Wranging")
        rt_wide <- d %>% filter(evt_time==t) %>%  group_by(id, run, run_trial) %>% 
          pivot_wider(names_from = c(evt_time_f), values_from = pow) 
        # analysis -----
        # timepoints <- timepoints[1]# TEST ONLY
        # message(paste("Analyzing timepoint", t,  sep = " "))
        rt_wide$h<-rt_wide[[t]]
        if (random) {
          md <-  lmerTest::lmer(scale(rt_next) ~ scale(pow) * rt_csv_sc * outcome  + scale(pow) * scale(rt_vmax)  +
                                  scale(pow) * rt_lag_sc + 
                                  (rt_csv_sc + rt_lag_sc + scale(rt_vmax)|id), d, control=lmerControl(optimizer = "nloptwrap"))
          while (any(grepl("failed to converge", md@optinfo$conv$lme4$messages) )) {
            print(md@optinfo$conv$lme4$conv)
            ss <- getME(md,c("theta","fixef"))
            md <- update(md, start=ss, control=lmerControl(optimizer = "bobyqa"))}  
        } else {
          md <-  lmerTest::lmer(scale(rt_next) ~ scale(pow) * rt_csv_sc * outcome  + scale(pow) * scale(rt_vmax)  +
                                  scale(pow) * rt_lag_sc + 
                                  (1|id), d, control=lmerControl(optimizer = "nloptwrap"))
          while (any(grepl("failed to converge", md@optinfo$conv$lme4$messages) )) {
            print(md@optinfo$conv$lme4$conv)
            ss <- getME(md,c("theta","fixef"))
            md <- update(md, start=ss, control=lmerControl(optimizer = "bobyqa"))}  
        }
        
        dm <- tidy(md)
        dm$sensor <- sensor
        dm$t <- t
        dm$freq <- freq
        # print(dm)
        dm}
      newlist[[freq]] <- dd}
    df <-  do.call(rbind,newlist)
    biglist[[sensor]] <- df}
  
  rdf <-  do.call(rbind,biglist)
  # format statistics ----
  
  rdf <- rdf %>% mutate(stat_order = as.factor(case_when(abs(statistic) < 2 ~ '1', 
                                                         abs(statistic) > 2 & abs(statistic) < 3 ~ '2', 
                                                         abs(statistic) > 3 ~ '3')),
                        p_value = as.factor(case_when(`p.value` > .05 ~ '1',
                                                      `p.value` < .05 & `p.value` > .01 ~ '2',
                                                      `p.value` < .01 & `p.value` > .001 ~ '3',
                                                      `p.value` <.001 ~ '4')))
  rdf$t <- as.numeric(rdf$t)
  # rdf$label <- as.factor(sub("_[^_]+$", "", rdf$label))
  rdf$stat_order <- factor(rdf$stat_order, labels = c("NS", "|t| > 2", "|t| > 3"))
  rdf$p_value <- factor(rdf$p_value, labels = c("NS", "p < .05", "p < .01", "p < .001"))
  terms <- unique(rdf$term[rdf$effect=="fixed"])
  terms <- terms[grepl("(h)",terms)]
  # FDR correction ----
  rdf <- rdf  %>% group_by(term, sensor) %>% mutate(p_fdr = p.adjust(p.value, method = 'fdr'),
                                                    p_level_fdr = as.factor(case_when(
                                                      # p_fdr > .1 ~ '0',
                                                      # p_fdr < .1 & p_fdr > .05 ~ '1',
                                                      p_fdr > .05 ~ '1',
                                                      p_fdr < .05 & p_fdr > .01 ~ '2',
                                                      p_fdr < .01 & p_fdr > .001 ~ '3',
                                                      p_fdr <.001 ~ '4'))
  ) %>% ungroup() #%>% mutate(side = substr(as.character(label), nchar(as.character(label)), nchar(as.character(label))),
  #          region = substr(as.character(label), 1, nchar(as.character(label))-2))
  rdf$p_level_fdr <- factor(rdf$p_level_fdr, levels = c('1', '2', '3', '4'), labels = c("NS","p < .05", "p < .01", "p < .001"))
  rdf$`p, FDR-corrected` = rdf$p_level_fdr
  setwd(rt_plot_dir)
  save(file = "meg_freq_medusa_rt_predict_output_nested.Rdata", rdf)
}
stopCluster(cl)
