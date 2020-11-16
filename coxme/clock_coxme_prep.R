# runs mixed-effects Cox models on clock data
# when running the first time, first run compute_sceptic_fmri_statistics.R
# basedir <- "~/Data_Analysis"
basedir <- "~/code"
setwd(file.path(basedir, "clock_analysis/coxme"))
library(readr)
library(ggplot2)
library(tidyr)
library(tibble)
library(reshape2)
library(emmeans)
library(dplyr)
library(lme4)
library(survival)
library(coxme)
library(survminer)
# library(OIsurv)
library(ggpubr)
# devtools::install_github('junkka/ehahelper') # requires gfortran, $ brew cask install gfortran
library(ehahelper)
library(broom)
library(broom.mixed)
library(car)
# load(file="clock_for_coxme_value_only_070518.RData")
# msdf <- read_csv("~/code/clock_analysis/coxme/sceptic_signals/mmclock_meg_fixed_uv_ureset_fixedparams_meg_ffx_trial_statistics.csv.gz")


#################
## Prep fMRI data
#################

load(file="clock_for_coxme_value_only_070518.RData")
# msdf <- read_csv("~/code/clock_analysis/coxme/sceptic_signals/mmclock_meg_fixed_uv_ureset_fixedparams_meg_ffx_trial_statistics.csv.gz")

# sanity checks on within-trial matrices

# add DAN and other betas
setwd(file.path(basedir, 'clock_analysis/fmri/keuka_brain_behavior_analyses/'))
load('trial_df_and_vh_pe_clusters_u.Rdata')
hdf <- df %>% select (ID, DAN, dan_parietal, dan_l_sfg, dan_r_sfg, general_entropy,
                      pe_f1_cort_hipp, pe_f2_cerebell, pe_f3_str, pe_PH_r, pe_ips, med_par, fef, entropy_vlPFC) %>% unique()
sdf <- sdf %>% inner_join(hdf, by = "ID")

# mark events (responses)
sdf$response <- round(sdf$rt/1000, digits = 1)==sdf$t2

# number bins
sdf <- sdf %>% group_by(ID, run, trial) %>% dplyr::mutate(bin = 1:n(), time = bin/10) %>% ungroup() %>% 
  group_by(ID) %>% mutate(value_wi = scale(value),
                          uncertainty_wi = scale(uncertainty),
                          value_b = mean(value),
                          uncertainty_b = mean(uncertainty), 
                          trial_neg_inv_sc = scale(-1/trial)) %>% ungroup() 
# filter out no-go zones at the edges
# censoring the entire first second seems too harsh, remove first 500 ms
fdf <- sdf %>% filter(bin >5 & bin <35)
sdf <- sdf %>% select(ID, run, trial, bin, everything())

#use trial-wise data frame for getting lagged timesteps
trialwise <- sdf %>% filter(bin==1) %>% 
  group_by(ID, run) %>%
  mutate(timesteplag1=dplyr::lag(timestep, n=1, order_by=trial),
         timesteplag2=dplyr::lag(timestep, n=2, order_by=trial),
         timesteplag3=dplyr::lag(timestep, n=3, order_by=trial),
         timesteplag4=dplyr::lag(timestep, n=4, order_by=trial)) %>%
  ungroup() %>% select(ID, run, trial, timesteplag1, timesteplag2, timesteplag3, timesteplag4)

sdf <- sdf %>% select(-timesteplag) %>% left_join(trialwise, by=c("ID", "run", "trial"))


get_wv_smile <- function(microdf, nlags=3, nbefore=0, nafter=1, spec=FALSE) {
  smile <- rep(0, nrow(microdf))
  lcols <- ifelse(isTRUE(spec), paste0("timesteplag", nlags), paste0("timesteplag", 1:nlags))
  wvprefix <- ifelse(isTRUE(spec), "wvs", "wv")
  vv <- microdf[1,lcols] #just first row (bin) of this trial is needed
  if (nrow(vv) > 0L) {
    allpos <- lapply(vv, function(x) {
      if (is.na(x)) {
        return(NULL)
      } else {
        return(seq.int(x-nbefore, x+nafter))
      }
    })
    
    topopulate <- unname(unlist(allpos))
    smile[topopulate[topopulate <= nrow(microdf)]] <- 1
  }
  microdf[[paste0(wvprefix, nlags, "b", nbefore, "a", nafter)]] <- smile
  return(microdf)
}

#get_wv_smile assumes bin-ordered data
sdf <- sdf %>% arrange(ID, run, trial, bin)
sdf$splitbasis <- with(sdf, paste(ID, run, trial, sep="."))

splitdf <- split(sdf, sdf$splitbasis)
splitdf <- lapply(splitdf, function(microdf) {
  microdf %>% get_wv_smile(n=3, nbefore=0, nafter=1) %>%
    get_wv_smile(n=3, nbefore=1, nafter=1) %>%
    get_wv_smile(n=3, nbefore=1, nafter=2) %>%
    get_wv_smile(n=4, nbefore=0, nafter=1) %>%
    get_wv_smile(n=4, nbefore=1, nafter=1) %>%
    get_wv_smile(n=4, nbefore=1, nafter=2) %>%
    get_wv_smile(n=1, nbefore=1, nafter=1, spec=TRUE) %>%
    get_wv_smile(n=2, nbefore=1, nafter=1, spec=TRUE) %>%
    get_wv_smile(n=3, nbefore=1, nafter=1, spec=TRUE) %>%
    get_wv_smile(n=4, nbefore=1, nafter=1, spec=TRUE)
})

bb <- bind_rows(splitdf)
bb$run_trial <- bb$trial

fbb <- bb %>% filter(bin >10 & bin <35)
#spot check
bb %>% filter(trial > 5) %>% select(bin, timestep, timesteplag1, timesteplag2, wv3b0a1)
table(bb$wv3b0a1)
# table(bb$wv4b0a1)

### end wv smiles

###############
# Prep MEG data
###############

# start with MEG stats produced by parse_sceptic_outputs.R
setwd(file.path(basedir, "clock_analysis/coxme"))
meg <- read_rds("sceptic_signals/mmclock_meg_sceptic_selective_detailedstats_Nov2020.rds")
megu <- read_rds("sceptic_signals/mmclock_meg_sceptic_fixeduv_detailedstats_Nov2020.rds")
mvdf <- meg$v_df %>% pivot_longer(-c(id, run, trial, rewFunc, y_chosen), names_to = "bin_string", values_to = "v") %>%
  mutate(bin = as.numeric(substr(bin_string, 3, 5))) %>% filter(bin<=y_chosen)
# # sanity check
# ggplot(mvdf, aes(bin, v, color = rewFunc)) + geom_smooth()
mudf <- megu$u_df %>% pivot_longer(-c(id, run, trial, rewFunc, y_chosen), names_to = "bin_string", values_to = "u") %>%
  mutate(bin = as.numeric(substr(bin_string, 3, 5))) %>% filter(bin<=y_chosen) %>% select(id, run, trial, u, bin)
# # sanity check
# ggplot(mudf, aes(bin, u, color = rewFunc)) + geom_smooth()
mvudf <- merge(mvdf, mudf) %>% arrange(id, run, trial, bin) %>% mutate (
  ID = substr(id, 1,5)
) %>% select(-id)
# ggplot(mdf, aes(bin, u, color = rewFunc)) + geom_smooth()
# ggplot(mdf, aes(bin, v, color = rewFunc)) + geom_smooth()

mvudf$t2 <- mvudf$bin/10
mvudf$t1 <- mvudf$t2 - .1
mvudf$response <- mvudf$y_chosen==mvudf$bin
# test coxme
# summary(test <- coxme(Surv(t1,t2,response) ~ v + u +
#                                 (1|id), mdf))

# add trial-level more variables

load("../fmri/keuka_brain_behavior_analyses/trial_df_and_vh_pe_clusters_u.Rdata")
msdf <- merge(mdf, mvudf)
# ggplot(msdf, aes(bin, v, color = rewFunc)) + geom_smooth()

msdf$trial <- as.numeric(msdf$trial)
msdf$timestep <- msdf$y_chosen
# scale value and uncertainty
msdf <- msdf %>% group_by(ID) %>% mutate(value_wi = scale(v),
                          uncertainty_wi = scale(u),
                          value_b = mean(v),
                          uncertainty_b = mean(u), 
                          trial_neg_inv_sc = scale(-1/run_trial)) %>% ungroup() 
# filter out no-go zones at the edges
# censoring the entire first second seems too harsh, remove first 500 ms
mfdf <- msdf %>% filter(bin >5 & bin <35)

# diagnostics on uncertainty and value distributions

# ggplot(msdf %>% filter(bin >10 & bin <35), aes(time, value_wi, color = rewFunc)) + geom_smooth(method = 'gam', formula = y~splines::ns(x,2)) 
# ggplot(msdf %>% filter(bin >10 & bin <35), aes(time, uncertainty_wi, color = rewFunc)) + geom_smooth(method = 'gam', formula = y~splines::ns(x,2))
# corr.test(fdf$value_wi, fdf$uncertainty_wi)
# corr.test(fdf$value, fdf$uncertainty)
# 


##Build 0/1 WV smiles (selection history with different buffer sizes)
# msdf <- msdf %>% select(ID, run, trial, bin, everything())

#use trial-wise data frame for getting lagged timesteps
trialwise <- msdf %>% filter(bin==1) %>% 
  group_by(ID, run) %>%
  mutate(timesteplag1=dplyr::lag(timestep, n=1, order_by=run_trial),
         timesteplag2=dplyr::lag(timestep, n=2, order_by=run_trial),
         timesteplag3=dplyr::lag(timestep, n=3, order_by=run_trial),
         timesteplag4=dplyr::lag(timestep, n=4, order_by=run_trial)) %>%
  ungroup() %>% select(ID, run, trial, timesteplag1, timesteplag2, timesteplag3, timesteplag4)

# msdf <- msdf %>% select(-timesteplag) %>% left_join(trialwise, by=c("ID", "run", "trial"))
msdf <- msdf %>% left_join(trialwise, by=c("ID", "run", "trial"))


get_wv_smile <- function(microdf, nlags=3, nbefore=0, nafter=1, spec=FALSE) {
  smile <- rep(0, nrow(microdf))
  lcols <- ifelse(isTRUE(spec), paste0("timesteplag", nlags), paste0("timesteplag", 1:nlags))
  wvprefix <- ifelse(isTRUE(spec), "wvs", "wv")
  vv <- microdf[1,lcols] #just first row (bin) of this trial is needed
  if (nrow(vv) > 0L) {
    allpos <- lapply(vv, function(x) {
      if (is.na(x)) {
        return(NULL)
      } else {
        return(seq.int(x-nbefore, x+nafter))
      }
    })
    
    topopulate <- unname(unlist(allpos))
    smile[topopulate[topopulate <= nrow(microdf)]] <- 1
  }
  microdf[[paste0(wvprefix, nlags, "b", nbefore, "a", nafter)]] <- smile
  return(microdf)
}

#get_wv_smile assumes bin-ordered data
msdf <- msdf %>% arrange(ID, run, trial, bin)

msdf$splitbasis <- with(msdf, paste(ID, run, trial, sep="."))

splitdf <- split(msdf, msdf$splitbasis)
splitdf <- lapply(splitdf, function(microdf) {
  microdf %>% get_wv_smile(n=3, nbefore=0, nafter=1) %>%
    get_wv_smile(n=3, nbefore=1, nafter=1) %>%
    get_wv_smile(n=3, nbefore=1, nafter=2) %>%
    get_wv_smile(n=4, nbefore=0, nafter=1) %>%
    get_wv_smile(n=4, nbefore=1, nafter=1) %>%
    get_wv_smile(n=4, nbefore=1, nafter=2) %>%
    get_wv_smile(n=1, nbefore=1, nafter=1, spec=TRUE) %>%
    get_wv_smile(n=2, nbefore=1, nafter=1, spec=TRUE) %>%
    get_wv_smile(n=3, nbefore=1, nafter=1, spec=TRUE) %>%
    get_wv_smile(n=4, nbefore=1, nafter=1, spec=TRUE)
})

mbb <- bind_rows(splitdf)

mfbb <- mbb %>% filter(bin >10 & bin <35)

#spot check
# mbb %>% filter(trial > 5) %>% select(bin, timestep, timesteplag1, timesteplag2, wv3b0a1)
#  table(mbb$wv3b0a1)
#  table(mbb$wv4b0a1)

# spot-check obsessively
pdf('meg_check.pdf', height = 6, width = 6)
ggplot(mbb, aes(run_trial, value_wi, color = rewFunc)) + geom_smooth()
dev.off()

pdf('fmri_check.pdf', height = 6, width = 6)
ggplot(bb, aes(run_trial, value_wi, color = rewFunc)) + geom_smooth()
dev.off()


save(file = "fMRI_MEG_coxme_objects_Nov15_2020", bb, fbb, mbb, mfbb)
 