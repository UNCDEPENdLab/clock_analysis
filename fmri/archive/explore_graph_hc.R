library(dplyr)
load('./fmri/keuka_brain_behavior_analyses/explore_clock_rev_vba_out.rdata')
df <- vba_output
df <- df %>% group_by(id) %>% mutate(total_earnings = sum(score_csv))
df <- df %>%
  group_by(id, run) %>%  dplyr::mutate(rt_swing = abs(c(NA, diff(rt_csv))),
                                       rt_swing_lr = abs(log(rt_csv/lag(rt_csv))),
                                       rt_cs = scale(rt_csv),
                                       rt_lag_cs = lag(rt_cs),
                                       rt_lag = lag(rt_csv) ,
                                       rt_swing_lag = lag(rt_swing),
                                       omission_lag = lag(score_csv==0),
                                       rt_vmax_lag = lag(rt_vmax),
                                       v_max_wi = scale(v_max),
                                       v_max_wi_lag = lag(v_max_wi),
                                       v_entropy_wi = scale(v_entropy),
                                       v_max_b = mean(na.omit(v_max)),
                                       v_entropy_b = mean(na.omit(v_entropy)),
                                       rt_change = rt_csv - rt_lag,
                                       pe_max_lag = lag(pe_max),
                                       abs_pe_max_lag = abs(pe_max_lag),
                                       rt_vmax_change = rt_vmax - rt_vmax_lag) %>% ungroup() %>%
  dplyr::mutate(block=ceiling(trial/40),
                rev_trial = trial - (block-1)*40, rev_trial_ctr = rev_trial - 20) %>% ungroup()
load('./fmri/keuka_brain_behavior_analyses/explore_demo.rdata')
s <- explore_demo
s$id <- s$registration_redcapid
s$group <- factor(s$registration_group, levels = c("ATT", "HC", "DEP", "IDE"))
s <- s %>% select(id,group)
dfs <- merge(s,df)
dfs$earnings_sp<-ifelse(dfs$total_earnings >= median(dfs$total_earnings[!duplicated(dfs$id)]),"Total Earning >= Median","Total Earning < Median")
save(dfs,file = 'explore_prelim_behav_data.Rdata')



stop("STOP LOADING...")

