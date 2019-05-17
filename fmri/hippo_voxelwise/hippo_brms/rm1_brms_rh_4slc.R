library(dplyr)
library(brms)
load("feedback_hipp_wide_ts.Rdata")

setwd(file.path(getMainDir(), "clock_analysis", "fmri", "hippo_voxelwise", "hippo_brms"))
str(fb_wide)
fb_wide <- fb_wide %>% mutate(block=factor(paste0(id, run))) #%>% dplyr::select(-evt_time)

#spread evt_time wide
fb_r <- fb_wide %>% dplyr::select(id, block, run, run_trial, evt_time, ends_with("_r"))

names(fb_r) <- sub("hipp_(\\d+)_r", "rh\\1", names(fb_r), perl=TRUE)
#shorten hippocampus labels to "rh" for left hippocampus (to avoid Mplus 8-character complaints)
fb_r_lags <- fb_r %>% na.omit() %>% arrange(id, run, run_trial, evt_time) %>%  group_by(id, run, run_trial) %>%
  mutate_at(vars(starts_with("rh")), list(l=~lag(., order_by=evt_time))) %>% #need unquoted order_by argument
  ungroup()

#mlvar analysis using Stan/brms

#remove correlations of intercepts and slopes
b2 <- bf(rh2 ~ rh2_l + rh5_l + rh9_l + rh12_l + (1|int|id) + (0 + rh2_l + rh5_l + rh9_l + rh12_l |slo|id))
b5 <- bf(rh5 ~ rh5_l + rh2_l + rh9_l + rh12_l + (1|int|id) + (0 + rh5_l + rh2_l + rh9_l + rh12_l |slo|id))
b9 <- bf(rh9 ~ rh9_l + rh2_l + rh5_l + rh12_l + (1|int|id) + (0 + rh9_l + rh2_l + rh5_l + rh12_l |slo|id))
b12 <- bf(rh12 ~ rh12_l + rh2_l + rh5_l + rh9_l + (1|int|id) + (0 + rh12_l + rh2_l + rh5_l + rh9_l |slo|id))

rm1_rh_4slc <- brm(b2 + b5 + b9 + b12 + set_rescor(TRUE), data = fb_r_lags, chains = 4, cores = 4, autocor = cor_arma(~run_trial|id, 1), iter=2500)

summary(rm1_rh_4slc)

save(rm1_rh_4slc, file="output/rm1_brms_rh_4slc.RData")
