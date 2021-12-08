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

library(brms)
#mlvar analysis using Stan/brms
## b1 <- bf(lh1 ~ lh1_l + lh4_l + lh8_l + lh11_l + (1+lh1_l|p|id))
## b4 <- bf(lh4 ~ lh4_l + lh1_l + lh8_l + lh11_l + (1+lh4_l|p|id))
## b8 <- bf(lh8 ~ lh8_l + lh1_l + lh4_l + lh11_l + (1+lh8_l|p|id))
## b11 <- bf(lh11 ~ lh11_l + lh1_l + lh4_l + lh8_l + (1+lh11_l|p|id))

## fit4 <- brm(b1 + b4 + b8 + b11 + set_rescor(TRUE), data = fb_l_lags, chains = 4, cores = 4, autocor = cor_arma(~run_trial|id, 1))

#remove correlations of intercepts and slopes
b1 <- bf(lh1 ~ lh1_l + lh4_l + lh8_l + lh11_l + (1|int|id) + (0 + lh1_l + lh4_l + lh8_l + lh11_l |slo|id))
b4 <- bf(lh4 ~ lh4_l + lh1_l + lh8_l + lh11_l + (1|int|id) + (0 + lh4_l + lh1_l + lh8_l + lh11_l |slo|id))
b8 <- bf(lh8 ~ lh8_l + lh1_l + lh4_l + lh11_l + (1|int|id) + (0 + lh8_l + lh1_l + lh4_l + lh11_l |slo|id))
b11 <- bf(lh11 ~ lh11_l + lh1_l + lh4_l + lh8_l + (1|int|id) + (0 + lh11_l + lh1_l + lh4_l + lh8_l |slo|id))

fit_bigre <- brm(b1 + b4 + b8 + b11 + set_rescor(TRUE), data = fb_l_lags, chains = 4, cores = 4, autocor = cor_arma(~run_trial|id, 1), iter=2000)

summary(fit_bigre)

save(fit_bigre, file="brms_test_4slc_bigre.RData")
