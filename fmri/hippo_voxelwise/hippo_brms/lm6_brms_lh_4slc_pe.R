library(dplyr)
library(brms)
load("feedback_hipp_wide_ts.Rdata")

setwd(file.path(getMainDir(), "clock_analysis", "fmri", "hippo_voxelwise", "mplus_var_hippo"))
str(fb_wide_ex)
fb_wide_ex <- fb_wide_ex %>% mutate(block=factor(paste0(id, run))) #%>% dplyr::select(-evt_time)

#spread evt_time wide
fb_l <- fb_wide_ex %>% dplyr::select(id, block, run, run_trial, evt_time, pe_max, ends_with("_l"))

names(fb_l) <- sub("hipp_(\\d+)_l", "lh\\1", names(fb_l), perl=TRUE)
#shorten hippocampus labels to "lh" for left hippocampus (to avoid Mplus 8-character complaints)
fb_l_lags <- fb_l %>% na.omit() %>% arrange(id, run, run_trial, evt_time) %>%  group_by(id, run, run_trial) %>%
  mutate_at(vars(starts_with("lh")), list(l=~lag(., order_by=evt_time))) %>% #need unquoted order_by argument
  ungroup() %>% filter(evt_time >= 0 & evt_time <= 2)

#mlvar analysis using Stan/brms

#remove correlations of intercepts and slopes
b3 <- bf(lh3 ~ lh3_l*pe_max + lh6_l*pe_max + lh9_l*pe_max + lh12_l*pe_max + (1|int|id) + (0 + lh3_l + lh6_l + lh9_l + lh12_l + pe_max |slo|id))
b6 <- bf(lh6 ~ lh6_l*pe_max + lh3_l*pe_max + lh9_l*pe_max + lh12_l*pe_max + (1|int|id) + (0 + lh6_l + lh3_l + lh9_l + lh12_l + pe_max |slo|id))
b9 <- bf(lh9 ~ lh9_l*pe_max + lh3_l*pe_max + lh6_l*pe_max + lh12_l*pe_max + (1|int|id) + (0 + lh9_l + lh3_l + lh6_l + lh12_l + pe_max |slo|id))
b12 <- bf(lh12 ~ lh12_l*pe_max + lh3_l*pe_max + lh6_l*pe_max + lh9_l*pe_max + (1|int|id) + (0 + lh12_l + lh3_l + lh6_l + lh9_l + pe_max |slo|id))

m6_lh_4slc_pe <- brm(b3 + b6 + b9 + b12 + set_rescor(TRUE), data = fb_l_lags, chains = 4, cores = 4, autocor = cor_arma(~run_trial|id, 1), iter=2500)

summary(m6_lh_4slc_pe)

save(m6_lh_4slc_pe, file="m6_brms_lh_4slc_pe.RData")
