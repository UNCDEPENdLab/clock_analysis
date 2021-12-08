library(dplyr)
library(brms)
load("feedback_hipp_wide_ts.Rdata")

setwd(file.path(getMainDir(), "clock_analysis", "fmri", "hippo_voxelwise", "mplus_var_hippo"))
str(fb_wide_ex)
fb_wide_ex <- fb_wide_ex %>% mutate(block=factor(paste0(id, run))) #%>% dplyr::select(-evt_time)

#spread evt_time wide
fb_r <- fb_wide_ex %>% dplyr::select(id, block, run, run_trial, evt_time, v_entropy_wi, ends_with("_r"))

names(fb_r) <- sub("hipp_(\\d+)_r", "rh\\1", names(fb_r), perl=TRUE)
#shorten hippocampus labels to "rh" for left hippocampus (to avoid Mplus 8-character complaints)
fb_r_lags <- fb_r %>% na.omit() %>% arrange(id, run, run_trial, evt_time) %>%  group_by(id, run, run_trial) %>%
  mutate_at(vars(starts_with("rh")), list(l=~lag(., order_by=evt_time))) %>% #need unquoted order_by argument
  group_by(id, run) %>% mutate(v_entropy_next=lead(v_entropy_wi, n=1, order_by=run_trial)) %>%
  ungroup() %>% filter(evt_time > 2)

#mlvar analysis using Stan/brms

#remove correlations of intercepts and slopes
b3 <- bf(rh3 ~ rh3_l*v_entropy_next + rh6_l*v_entropy_next + rh9_l*v_entropy_next + rh12_l*v_entropy_next + (1|int|id) + (0 + rh3_l + rh6_l + rh9_l + rh12_l + v_entropy_next |slo|id))
b6 <- bf(rh6 ~ rh6_l*v_entropy_next + rh3_l*v_entropy_next + rh9_l*v_entropy_next + rh12_l*v_entropy_next + (1|int|id) + (0 + rh6_l + rh3_l + rh9_l + rh12_l + v_entropy_next |slo|id))
b9 <- bf(rh9 ~ rh9_l*v_entropy_next + rh3_l*v_entropy_next + rh6_l*v_entropy_next + rh12_l*v_entropy_next + (1|int|id) + (0 + rh9_l + rh3_l + rh6_l + rh12_l + v_entropy_next |slo|id))
b12 <- bf(rh12 ~ rh12_l*v_entropy_next + rh3_l*v_entropy_next + rh6_l*v_entropy_next + rh9_l*v_entropy_next + (1|int|id) + (0 + rh12_l + rh3_l + rh6_l + rh9_l + v_entropy_next |slo|id))

m5_rh_4slc_entropy_next <- brm(b3 + b6 + b9 + b12 + set_rescor(TRUE), data = fb_r_lags, chains = 4, cores = 4, autocor = cor_arma(~run_trial|id, 1), iter=2500)

summary(m5_rh_4slc_entropy_next)

save(m5_rh_4slc_entropy_next, file="m5_brms_rh_4slc_entropy_next.RData")
