library(dplyr)
library(brms)
load("feedback_hipp_wide_ts.Rdata")

setwd(file.path(getMainDir(), "clock_analysis", "fmri", "hippo_voxelwise", "mplus_var_hippo"))
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
b2 <- bf(lh2 ~ lh2_l + lh5_l + lh9_l + lh12_l + (1|int|id) + (0 + lh2_l + lh5_l + lh9_l + lh12_l |slo|id))
b5 <- bf(lh5 ~ lh5_l + lh2_l + lh9_l + lh12_l + (1|int|id) + (0 + lh5_l + lh2_l + lh9_l + lh12_l |slo|id))
b9 <- bf(lh9 ~ lh9_l + lh2_l + lh5_l + lh12_l + (1|int|id) + (0 + lh9_l + lh2_l + lh5_l + lh12_l |slo|id))
b12 <- bf(lh12 ~ lh12_l + lh2_l + lh5_l + lh9_l + (1|int|id) + (0 + lh12_l + lh2_l + lh5_l + lh9_l |slo|id))

lm1_lh_4slc <- brm(b2 + b5 + b9 + b12 + set_rescor(TRUE), data = fb_l_lags, chains = 4, cores = 4, autocor = cor_arma(~run_trial|id, 1), iter=2500)

summary(lm1_lh_4slc)

save(lm1_lh_4slc, file="lm1_brms_lh_4slc.RData")
