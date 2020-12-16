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
load (file.path(cache_dir, "fMRI_coxme_objects_for_trial_inv_MEDUSA_Dec15_2020.Rdata"))

# load wide MEDUSAs
load(file.path(cache_dir, 'clock_dan_wide_ts.Rdata'))
load(file.path(cache_dir, 'rt_dan_wide_ts.Rdata'))

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

if (lagged_decon) {
rt_wide <- rt_wide %>% group_by(id, run) %>% arrange(id, run, run_trial)
for (label in labels) {
  varname = paste0(label, "_lag")
  rt_wide <- rt_wide %>% 
    mutate(!!varname := lag(!!as.name(label)))
  }
# # spot-check
# test <- rt_wide %>% filter(id==10637) %>% select(c(run, run_trial, !!label, !!varname))
# View(test)

# merge ----
# merge survival objects with lagged MEDUSA data
rt_wide <- rt_wide %>% select(-labels) %>% mutate(ID = id) %>% ungroup()
} else {
  rt_wide <- rt_wide %>%  mutate(ID = id) %>% ungroup()
}
medlag_bb <- inner_join(bb, rt_wide)
medlag_fbb <- inner_join(fbb, rt_wide)

# save ----
save(file = file.path(cache_dir,"fMRI_coxme_objects_with_wtrial_inv_medusa_Dec15_2020.RData"), medlag_bb, medlag_fbb)
# reset working directory
setwd(cwd)
