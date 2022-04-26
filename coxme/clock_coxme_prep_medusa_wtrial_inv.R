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
medusa_dir = "/Volumes/GoogleDrive/My Drive/SCEPTIC_fMRI/dan_medusa"
cache_dir = "/Volumes/GoogleDrive/My Drive/SCEPTIC_fMRI/dan_medusa/cache"
repo_directory <- file.path(basedir,"clock_analysis")
cwd <- getwd()
setwd(coxme_dir)

# these are set inside clock_coxme_wv_medusa.R
censor_clock_post_rt = F # censor post-RT epoch when analyzing clock data
lagged_decon = F # whether to use lagged (by 1s) or unlagged decon
# load survival and wide MEDUSA ----

# load coxme survival object
# save(file = file.path(cache_dir, "fMRI_coxme_objects_for_trial_inv_MEDUSA_Dec15_2020.Rdata"), bb, fbb)
cat("Loading survival objects\n") # bb - fMRI survival object, # fbb - same, first 1s and last 0.5s filtered out
load (file.path(cache_dir, "fMRI_coxme_objects_for_trial_inv_MEDUSA_Dec15_2020.Rdata"))

# load wide (different variables for each within-trial timepoint) MEDUSAs
cat("Loading clock decons\n")
load(file.path(cache_dir, 'clock_dan_wide_ts.Rdata'))
cat("Loading RT decons\n")
load(file.path(cache_dir, 'rt_dan_wide_ts.Rdata'))

cat("Scaling behavioral variables\n")
# decompose U and V into within- and between-trial
bb <- bb %>% group_by(ID, run, run_trial) %>% mutate(value_wi_t = scale(value),
                                                     uncertainty_wi_t = scale(uncertainty),
                                                     value_b_t = mean(value),
                                                     uncertainty_b_t = mean(uncertainty)) %>% ungroup()
fbb <- fbb %>% group_by(ID, run, run_trial) %>% mutate(value_wi_t = scale(value),
                                                       uncertainty_wi_t = scale(uncertainty),
                                                       value_b_t = mean(value),
                                                       uncertainty_b_t = mean(uncertainty)) %>% ungroup()

# swap in censored clock-aligned data
# if (censor_clock_post_rt) {clock_wide <- clock_wide_cens}

# get lagged brain signal ----
# lags are for next RT prediction
rt_labels <- names(rt_wide[grepl("_R|_r|_L|_l", names(rt_wide))])
clock_labels <- names(clock_wide[grepl("_R|_r|_L|_l", names(clock_wide))])

rt_wide <- rt_wide %>% group_by(id, run) %>% arrange(id, run, run_trial)
clock_wide <- clock_wide %>% group_by(id, run) %>% arrange(id, run, run_trial)

cat("Getting decon lags\n")
if (lagged_decon) {
  for (label in rt_labels) {
    varname = paste0(label, "_lag")
    rt_wide <- rt_wide %>% 
      mutate(!!varname := lag(!!as.name(label)))
  }
  for (label in clock_labels) {
    varname = paste0(label, "_lag")
    clock_wide <- clock_wide %>% 
      mutate(!!varname := lag(!!as.name(label)))
  }
  
  # # spot-check
  # test <- rt_wide %>% filter(id==10637) %>% select(c(run, run_trial, !!label, !!varname))
  # View(test)
  # merge ----
  # merge survival objects with lagged MEDUSA data
  cat("Removing unlagged decons\n")
  rt_wide <- rt_wide %>% select(-all_of(rt_labels)) %>% mutate(ID = id) %>% ungroup()
  clock_wide <- clock_wide %>% select(-all_of(clock_labels)) %>% mutate(ID = id) %>% ungroup()

  cat("Merging\n")
  rt_lag_bb <- inner_join(bb, rt_wide)
  rt_lag_fbb <- inner_join(fbb, rt_wide)
  clock_lag_bb <- inner_join(bb, clock_wide)
  clock_lag_fbb <- inner_join(fbb, clock_wide)
  cat("Saving\n")
  save(file = file.path(cache_dir,"fMRI_coxme_objects_with_wtrial_inv_medusa_lagged_Dec29_2020.RData"), clock_lag_bb, clock_lag_fbb, rt_lag_bb, rt_lag_fbb)
} else {
  rt_wide <- rt_wide %>%  mutate(ID = id) %>% ungroup()
  clock_wide <- clock_wide %>%  mutate(ID = id) %>% ungroup()
  cat("Merging\n")
  rt_bb <- inner_join(bb, rt_wide)
  rt_fbb <- inner_join(fbb, rt_wide)
  clock_bb <- inner_join(bb, clock_wide)
  clock_fbb <- inner_join(fbb, clock_wide)
  cat("Saving\n")
  save(file = file.path(cache_dir,"fMRI_coxme_objects_with_wtrial_inv_medusa_UNlagged_Dec29_2020.RData"), clock_bb, clock_fbb, rt_bb, rt_fbb)
}

# save ----

# reset working directory
setwd(cwd)

