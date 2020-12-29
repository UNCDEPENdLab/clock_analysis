# prepares clock data for mixed-effects Cox models 
# where response hazard is predicted by lagged deconvolved signal
library(readr)
library(ggplot2)
library(tidyverse)
library(lme4)
library(survival)
library(coxme)
library(survminer)
library(ggpubr)
# devtools::install_github('junkka/ehahelper') # requires gfortran, $ brew cask install gfortran
library(ehahelper)
library(car)

# basedir <- "~/Data_Analysis"
basedir <- "~/code"
coxme_dir <- file.path(basedir, "clock_analysis/coxme")
medusa_dir = "~/Box/SCEPTIC_fMRI/dan_medusa"
cache_dir = "~/Box/SCEPTIC_fMRI/dan_medusa/cache"
repo_directory <- file.path(basedir,"clock_analysis")
cwd <- getwd()
setwd(coxme_dir)

# load survival and wide MEDUSA ----

# load coxme survival object
# save(file = file.path(cache_dir, "fMRI_coxme_objects_for_trial_inv_MEDUSA_Dec15_2020.Rdata"), bb, fbb)
print("Loading survival objects")
load (file.path(cache_dir, "fMRI_coxme_objects_for_trial_inv_MEDUSA_Dec15_2020.Rdata"))

# load wide MEDUSAs
print("Loading clock decons")
load(file.path(cache_dir, 'clock_dan_wide_ts.Rdata'))
print("Loading RT decons")
load(file.path(cache_dir, 'rt_dan_wide_ts.Rdata'))

print("Adding lags")
# U and V within- and between-trial
bb <- bb %>% group_by(ID, run, run_trial) %>% mutate(value_wi_t = scale(value),
                                        uncertainty_wi_t = scale(uncertainty),
                                        value_b_t = mean(value),
                                        uncertainty_b_t = mean(uncertainty)) %>% ungroup()
fbb <- fbb %>% group_by(ID, run, run_trial) %>% mutate(value_wi_t = scale(value),
                                                     uncertainty_wi_t = scale(uncertainty),
                                                     value_b_t = mean(value),
                                                     uncertainty_b_t = mean(uncertainty)) %>% ungroup()

# get lagged brain signal ----
# start with RT-locked decons, get lags for RT prediction
labels <- names(rt_wide[grepl("_R|_r|_L|_l", names(rt_wide))])
rt_wide <- rt_wide %>% group_by(id, run) %>% arrange(id, run, run_trial)
clock_wide <- clock_wide %>% group_by(id, run) %>% arrange(id, run, run_trial)

print("Preparing for merge")
if (lagged_decon) {
for (label in labels) {
  varname = paste0(label, "_lag")
  rt_wide <- rt_wide %>% 
    mutate(!!varname := lag(!!as.name(label)))
}
for (label in labels) {
  varname = paste0(label, "_lag")
  clock_wide <- clock_wide %>% 
    mutate(!!varname := lag(!!as.name(label)))
}

# # spot-check
# test <- rt_wide %>% filter(id==10637) %>% select(c(run, run_trial, !!label, !!varname))
# View(test)
# merge ----
# merge survival objects with lagged MEDUSA data
rt_wide <- rt_wide %>% select(-labels) %>% mutate(ID = id) %>% ungroup()
clock_wide <- clock_wide %>% select(-labels) %>% mutate(ID = id) %>% ungroup()
print("Merging")
rt_lag_bb <- inner_join(bb, rt_wide)
rt_lag_fbb <- inner_join(fbb, rt_wide)
clock_lag_bb <- inner_join(bb, clock_wide)
clock_lag_fbb <- inner_join(fbb, clock_wide)
print("Saving")
save(file = file.path(cache_dir,"fMRI_coxme_objects_with_wtrial_inv_medusa_lagged_Dec29_2020.RData"), clock_lag_bb, clock_lag_fbb, rt_lag_bb, rt_lag_fbb)
} else {
  rt_wide <- rt_wide %>%  mutate(ID = id) %>% ungroup()
  clock_wide <- clock_wide %>%  mutate(ID = id) %>% ungroup()
  print("Merging")
  rt_bb <- inner_join(bb, rt_wide)
  rt_fbb <- inner_join(fbb, rt_wide)
  clock_bb <- inner_join(bb, clock_wide)
  clock_fbb <- inner_join(fbb, clock_wide)
  print("Saving")
  save(file = file.path(cache_dir,"fMRI_coxme_objects_with_wtrial_inv_medusa_UNlagged_Dec29_2020.RData"), clock_bb, clock_fbb, rt_bb, rt_fbb)
}

# save ----

# reset working directory
setwd(cwd)
