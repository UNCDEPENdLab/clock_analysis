library(modelr)
library(tidyverse)
library(lme4)
library(afex)
library(broom)
library(broom.mixed) #plays will with afex p-values in lmer wrapper
library(ggpubr)
library(car)
library(viridis)
library(psych)
library(corrplot)
library(foreach)
library(doParallel)
library(psych)
repo_directory <- "~/code/clock_analysis"
medusa_dir = "~/Box/SCEPTIC_fMRI/MEG_20Hz_n63/"
diag_dir =  "~/Box/SCEPTIC_fMRI/MEG_20Hz_n63/diags/"
behavioral_data_file = "~/code/clock_analysis/meg/MEG_n63_behavioral_data_preprocessed_trial_df.RDS"
source("~/code/fmri.pipeline/R/mixed_by.R")
stopifnot(dir.exists(medusa_dir))  

setwd(medusa_dir)

# options, files ----
plots = F
decode = F  # main analysis analogous to Fig. 4 E-G in NComm 2020
rt = T # predicts next response based on signal and behavioral variables
online = F # whether to analyze clock-aligned ("online") or RT-aligned ("offline") responses
exclude_first_run = F
reg_diagnostics = F
start_time = -3
domain = "tf" # "time"
label_sensors = F
test = F
scale_winsor = T
# # Kai’s guidance on sensors is: ‘So for FEF, I say focus on 612/613, 543/542, 1022/1023, 
# # For IPS, 1823, 1822, 2222,2223.’
# fef_sensors <- c("0612","0613", "0542", "0543","1022")
# ips_sensors <- c("1823", "1822", "2222","2223")
# dan_sensors <- c(fef_sensors,ips_sensors)

if (domain == "time") {
  files <-  gsub("//", "/", list.files(medusa_dir,pattern = "20Hz", full.names = T))
  files <- files[grepl("MEG", files)]
  if (test) {files <- files[1:4]}
} else if (domain == "tf") {
  files <-  gsub("//", "/", list.files(medusa_dir,pattern = "tf", full.names = T))
  files <- files[grepl("MEG", files)]
}

if (label_sensors | scale_winsor) {
  # make cluster ----
  library(parallel)
  ncores <- detectCores()
  if (ncores > 1L) {
    cl <- makeCluster(ncores)
    registerDoParallel(cl)
    on.exit(try(stopCluster(cl)))
  } else {
    registerDoSEQ()
  }
}

if (label_sensors) {
  for (this_file in files) {
    d <- readRDS(this_file)
    d$sensor <- str_extract(this_file, "[[:digit:]]{4}")
    saveRDS(d, file = this_file)
  }
}

scale2 <- function(x, na.rm = FALSE) (x - mean(x, na.rm = na.rm)) / sd(x, na.rm)
if (scale_winsor & domain == "time") {
  foreach(i = 1:length(files), .packages=c("tidyverse", "psych")) %dopar% {
    d <- readRDS(files[i])
    d$signal_scaled <- winsor(scale2(d$Signal), trim = .01)
    saveRDS(d, file = files[i])
  } 
  stopCluster(cl)} else if (scale_winsor & domain == "tf") {
  foreach(i = 1:length(files), .packages=c("tidyverse", "psych")) %dopar% {
    scale2 <- function(x, na.rm = FALSE) (x - mean(x, na.rm = na.rm)) / sd(x, na.rm)
    d <- readRDS(files[i])
    d$pow_scaled <- scale2(winsor(d$Pow, trim = .01))
    saveRDS(d, file = files[i])
  }
  stopCluster(cl)}
# files <- files[1] # TEST ONLY

# # take first few for testing
# all_sensors <- all_sensors[1:4]


# 
# sample_data <- readRDS(files[2])
# sample_data$pow_scaled <- winsor((sample_data$Pow), trim = .01)
# 
# ggplot(sample_data, aes(pow_scaled)) + geom_histogram() + facet_wrap(~Subject)
# ggplot(sample_data, aes(Pow)) + geom_histogram() + facet_wrap(~Subject)
# 
# check data - looks good!
# sensor <- all_sensors[[1]]
# rt <- as_tibble(readRDS(paste0("MEG", sensor, "_20Hz.rds"))) %>% filter(Time>start_time) %>%
#   rename(id = Subject, trial = Trial, run = Run, evt_time = Time, signal = Signal) %>%
#   mutate(signal = winsor(scale2(signal), trim = .01))  # scale signal across subjects
# cd(diag_dir)
# pdf("0111_signal_distributions.pdf", height = 30, width = 30)
# ggplot(rt, aes(signal)) + geom_histogram() + facet_wrap(~id)
# dev.off()

# check that merging variables are named identically in MEG and behavior files

# mixed_by call
trial_df <- readRDS(behavioral_data_file)
trial_df$Subject <- trial_df$id
trial_df$Run <- trial_df$run
trial_df$Trial <- trial_df$trial
decode_formula = formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + scale(rt_vmax_lag)  + scale(rt_vmax_change) + 
                           v_entropy_wi + v_entropy_wi_change + v_max_wi  + scale(abs_pe) + outcome + (1|Subject))
if (domain == "tf") {
  splits = c("Time", "sensor", "Freq")
  outcome = "pow_scaled"} else if (domain == "time") {
    splits = c("Time", "sensor")
    outcome = "signal_scaled"} 
ddf <- mixed_by(files, outcomes = outcome, rhs_model_formulae = decode_formula , split_on = splits, external_df = trial_df,
                padjust_by = "term", padjust_method = "fdr", ncores = 20, refit_on_nonconvergence = 3)

ddf <- as_tibble(ddf)

# save output
setwd("~/OneDrive/collected_letters/papers/meg/plots/rt_decode/")
saveRDS(ddf, file = "meg_mixed_by_tf_ddf.RDS")
