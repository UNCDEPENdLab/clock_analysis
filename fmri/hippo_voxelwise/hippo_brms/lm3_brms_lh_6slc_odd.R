library(dplyr)
library(brms)
load("feedback_hipp_wide_ts.Rdata")

setwd(file.path(getMainDir(), "clock_analysis", "fmri", "hippo_voxelwise", "hippo_brms"))
str(fb_wide)
fb_wide <- fb_wide %>% mutate(block=factor(paste0(id, run))) #%>% dplyr::select(-evt_time)

#spread evt_time wide
fb_l <- fb_wide %>% dplyr::select(id, block, run, run_trial, evt_time, ends_with("_l"))

names(fb_l) <- sub("hipp_(\\d+)_l", "lh\\1", names(fb_l), perl=TRUE)
#shorten hippocampus labels to "lh" for left hippocampus (to avoid Mplus 8-character complaints)
fb_l_lags <- fb_l %>% na.omit() %>% arrange(id, run, run_trial, evt_time) %>%  group_by(id, run, run_trial) %>%
  mutate_at(vars(starts_with("lh")), list(l=~lag(., order_by=evt_time))) %>% #need unquoted order_by argument
  ungroup()

#mlvar analysis using Stan/brms

#remove correlations of intercepts and slopes
b1 <-  bf(lh1 ~  lh1_l + lh3_l + lh5_l + lh7_l + lh9_l + lh11_l + (1|int|id) + (0 + lh1_l + lh3_l + lh5_l + lh7_l + lh9_l + lh11_l |slo|id))
b3 <-  bf(lh3 ~  lh1_l + lh3_l + lh5_l + lh7_l + lh9_l + lh11_l + (1|int|id) + (0 + lh1_l + lh3_l + lh5_l + lh7_l + lh9_l + lh11_l |slo|id))
b5 <-  bf(lh5 ~  lh1_l + lh3_l + lh5_l + lh7_l + lh9_l + lh11_l + (1|int|id) + (0 + lh1_l + lh3_l + lh5_l + lh7_l + lh9_l + lh11_l |slo|id))
b7 <-  bf(lh7 ~  lh1_l + lh3_l + lh5_l + lh7_l + lh9_l + lh11_l + (1|int|id) + (0 + lh1_l + lh3_l + lh5_l + lh7_l + lh9_l + lh11_l |slo|id))
b9 <-  bf(lh9 ~  lh1_l + lh3_l + lh5_l + lh7_l + lh9_l + lh11_l + (1|int|id) + (0 + lh1_l + lh3_l + lh5_l + lh7_l + lh9_l + lh11_l |slo|id))
b11 <- bf(lh11 ~ lh1_l + lh3_l + lh5_l + lh7_l + lh9_l + lh11_l + (1|int|id) + (0 + lh1_l + lh3_l + lh5_l + lh7_l + lh9_l + lh11_l |slo|id))

lm3_lh_6slc_odd <- brm(b1 + b3 + b5 + b7 + b9 + b11 + set_rescor(TRUE), data = fb_l_lags,
  chains = 4, cores = 4, autocor = cor_arma(~run_trial|id, 1), iter=2500, init_r=0.9)

summary(lm3_lh_6slc_odd)

save(lm3_lh_6slc_odd, file="output/lm3_brms_lh_6slc_odd.RData")
