library(dplyr)
library(brms)
load("feedback_hipp_wide_ts.Rdata")

setwd(file.path(getMainDir(), "clock_analysis", "fmri", "hippo_voxelwise", "hippo_brms"))
str(fb_wide)
fb_wide <- fb_wide %>% mutate(block=factor(paste0(id, run))) #%>% dplyr::select(-evt_time)

#spread evt_time wide
fb_r <- fb_wide %>% dplyr::select(id, block, run, trial, run_trial, evt_time, ends_with("_r"))

names(fb_r) <- sub("hipp_(\\d+)_r", "rh\\1", names(fb_r), perl=TRUE)
#shorten hippocampus labels to "rh" for left hippocampus (to avoid Mplus 8-character complaints)
fb_r_lags <- fb_r %>% na.omit() %>% arrange(id, run, run_trial, evt_time) %>%  group_by(id, run, run_trial) %>%
  mutate_at(vars(starts_with("rh")), list(l=~lag(., order_by=evt_time))) %>% #need unquoted order_by argument
  ungroup()

#mlvar analysis using Stan/brms

#remove correlations of intercepts and slopes
b1 <-  bf(rh1 ~  rh1_l + rh3_l + rh5_l + rh7_l + rh9_l + rh11_l + (1|int|id) + (0 + rh1_l + rh3_l + rh5_l + rh7_l + rh9_l + rh11_l |slo|id))
b3 <-  bf(rh3 ~  rh1_l + rh3_l + rh5_l + rh7_l + rh9_l + rh11_l + (1|int|id) + (0 + rh1_l + rh3_l + rh5_l + rh7_l + rh9_l + rh11_l |slo|id))
b5 <-  bf(rh5 ~  rh1_l + rh3_l + rh5_l + rh7_l + rh9_l + rh11_l + (1|int|id) + (0 + rh1_l + rh3_l + rh5_l + rh7_l + rh9_l + rh11_l |slo|id))
b7 <-  bf(rh7 ~  rh1_l + rh3_l + rh5_l + rh7_l + rh9_l + rh11_l + (1|int|id) + (0 + rh1_l + rh3_l + rh5_l + rh7_l + rh9_l + rh11_l |slo|id))
b9 <-  bf(rh9 ~  rh1_l + rh3_l + rh5_l + rh7_l + rh9_l + rh11_l + (1|int|id) + (0 + rh1_l + rh3_l + rh5_l + rh7_l + rh9_l + rh11_l |slo|id))
b11 <- bf(rh11 ~ rh1_l + rh3_l + rh5_l + rh7_l + rh9_l + rh11_l + (1|int|id) + (0 + rh1_l + rh3_l + rh5_l + rh7_l + rh9_l + rh11_l |slo|id))


#under brms 2.9.0, this cor_arma throws an error: "Error: Time points within groups must be unique."
#this is not corrected using ~run_trial|id:run
rm3_rh_6slc_odd <- brm(b1 + b3 + b5 + b7 + b9 + b11 + set_rescor(TRUE), data = fb_r_lags,
  chains = 4, cores = 4, autocor = cor_arma(~run_trial|id, p=1), iter=2500, init_r=0.9)

#switch to overall trial, which crosses the run boundary, but that seems fine for our simple goal
#this still throws an error for cor_arma under brms 2.9. Need to investigate
#rm3_rh_6slc_odd <- brm(b1 + b3 + b5 + b7 + b9 + b11 + set_rescor(TRUE), data = fb_r_lags,
#  chains = 4, cores = 4, autocor = cor_arma(~trial|id, 1), iter=2500, init_r=0.9)

summary(rm3_rh_6slc_odd)

save(rm3_rh_6slc, file="output/rm3_brms_rh_6slc_odd.RData")
