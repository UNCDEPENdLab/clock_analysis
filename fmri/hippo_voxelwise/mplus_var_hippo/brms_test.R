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
b1 <- bf(lh1 ~ lh2_l + (1|p|id))
b2 <- bf(lh2 ~ lh1_l + (1|p|id))
fit2 <- brm(b1 + b2, data = fb_l_lags, chains = 4, cores = 4, autocor = cor_arma(~run_trial|id, 1))

summary(fit2)

save(fit2, file="brms_test.RData")
