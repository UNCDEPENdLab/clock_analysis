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
encode = T  # main analysis analogous to Fig. 4 E-G in NComm 2020

rt = F # predicts next response based on signal and behavioral variables
online = F # whether to analyze clock-aligned ("online") or RT-aligned ("offline") responses
exclude_first_run = F
reg_diagnostics = F
start_time = -3
domain = "tf" # "time"
label_sensors = F
test = F
scale_winsor = F
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
    d$pow_scaled <- scale2(winsor(d$Pow, trim = .005))
    saveRDS(d, file = files[i])
  }
  stopCluster(cl)}
# files <- files[1] # TEST ONLY

# # take first few for testing
# all_sensors <- all_sensors[1:4]


# 
sample_data <- readRDS(files[2])
sample_data$pow_scaled <- winsor((sample_data$Pow), trim = .005)

ggplot(sample_data, aes(pow_scaled)) + geom_histogram() + facet_wrap(~Subject)
ggplot(sample_data, aes(Pow)) + geom_histogram() + facet_wrap(~Subject)
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
# trial_df$Subject <- trial_df$id
# trial_df$Run <- trial_df$run
# trial_df$Trial <- trial_df$trial
# trial_df$rt_next_sc <- scale(trial_df$rt_next)
# # save behavioral data with new variables
# saveRDS(trial_df, behavioral_data_file)

# encode_formula = formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + scale(rt_vmax_lag)  + scale(rt_vmax_change) + 
#                            v_entropy_wi + v_entropy_wi_change + v_max_wi  + scale(abs_pe) + outcome + (1|Subject))
encode_formula_e = formula(~ scale(rt_vmax_lag)*echange_f1_early + scale(rt_vmax_lag)*echange_f2_late + scale(rt_vmax_lag)*e_f1 +
                           scale(abs_pe)*echange_f1_early + scale(abs_pe)*echange_f2_late + scale(abs_pe)*e_f1 +
                           outcome*echange_f1_early + outcome*echange_f2_late + outcome*e_f1 +
                           rt_csv_sc*echange_f1_early + rt_csv_sc*echange_f2_late + rt_csv_sc*e_f1 +
                           trial_neg_inv_sc*echange_f1_early + trial_neg_inv_sc*echange_f2_late + trial_neg_inv_sc*e_f1 +
                           v_entropy_wi_change*echange_f1_early + v_entropy_wi_change*echange_f2_late + v_entropy_wi*e_f1 + rt_lag_sc*e_f1 + (1|Subject))
encode_formula_pe = formula(~ scale(rt_vmax_lag)*abs_pe_f2_early + scale(rt_vmax_lag)*abs_pabs_pe_f3_late_mid + scale(rt_vmax_lag)*abs_pe_f3_late +
                              scale(abs_pe)*abs_pe_f2_early + scale(abs_pe)*abs_pabs_pe_f3_late_mid + scale(abs_pe)*abs_pe_f3_late +
                              outcome*abs_pe_f2_early + outcome*abs_pabs_pe_f3_late_mid + outcome*abs_pe_f3_late +
                              rt_csv_sc*abs_pe_f2_early + rt_csv_sc*abs_pabs_pe_f3_late_mid + rt_csv_sc*abs_pe_f3_late +
                              trial_neg_inv_sc*abs_pe_f2_early + trial_neg_inv_sc*abs_pabs_pe_f3_late_mid + trial_neg_inv_sc*abs_pe_f3_late +
                              v_entropy_wi_change*abs_pe_f2_early + v_entropy_wi_change*abs_pabs_pe_f3_late_mid + v_entropy_wi*abs_pe_f3_late + rt_lag_sc*abs_pe_f3_late + (1|Subject))

rt_tf_formula = formula( ~ pow_scaled * rt_csv_sc * outcome  + pow_scaled * scale(rt_vmax)  +
                        pow_scaled * rt_lag_sc + 
                        (1|id))
rt_time_formula = formula( ~ signal_scaled * rt_csv_sc * outcome  + signal_scaled * scale(rt_vmax)  +
                             signal_scaled * rt_lag_sc + 
                           (1|id))
rt_outcome = "rt_next_sc"
if (domain == "tf") {
  splits = c("Time", "sensor", "Freq")
  signal_outcome = "pow_scaled"} else if (domain == "time") {
    splits = c("Time", "sensor")
    signal_outcome = "signal_scaled"} 
if (encode) {
  ddf <- as_tibble(mixed_by(files, outcomes = signal_outcome, rhs_model_formulae = encode_formula_pe, split_on = splits, external_df = trial_df,
                            padjust_by = "term", padjust_method = "fdr", ncores = 20, refit_on_nonconvergence = 5))
  # save output
  setwd("~/OneDrive/collected_letters/papers/meg/plots/rt_decode/")
  if (domain == "time") {
    saveRDS(ddf, file = "meg_mixed_by_time_ranefs_mult_interactions_ddf.RDS")
  } else if (domain == "tf") {
    saveRDS(ddf, file = "meg_mixed_by_tf_ranefs_mult_interactions_e_ddf.RDS")
  }
}

if (rt) {
  rdf <- as_tibble(mixed_by(files, outcomes = rt_outcome, rhs_model_formulae = rt_time_formula , split_on = splits, external_df = trial_df,
                            padjust_by = "term", padjust_method = "fdr", ncores = 20, refit_on_nonconvergence = 3))
  # save output
  setwd("~/OneDrive/collected_letters/papers/meg/plots/rt_rt/")
  if (domain == "time") {
    saveRDS(rdf, file = "meg_mixed_by_time_rdf.RDS")
  } else if (domain == "tf") {
    saveRDS(rdf, file = "meg_mixed_by_tf_rdf.RDS")
  }
}


