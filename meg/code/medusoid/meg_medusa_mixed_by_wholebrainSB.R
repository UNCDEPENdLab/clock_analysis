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
sensor_list <- read.table("~/code/clock_analysis/meg/code/meg_sensors_annotated.txt", header=TRUE, colClasses="character") %>%
  filter(dan != "no") %>% pull(sensor)
sensor_map <- read.table("~/code/clock_analysis/meg/code/meg_sensors_annotated.txt", header=TRUE, colClasses="character")
repo_directory <- "~/code/clock_analysis"
behavioral_data_file <- "~/code/clock_analysis/meg/MEG_n63_behavioral_data_preprocessed_trial_df.RDS"
source("~/code/fmri.pipeline/R/mixed_by.R")
# main analysis analogous to Fig. 4 E-G in NComm 2020
<<<<<<< HEAD:meg/code/medusoid/meg_medusa_mixed_by_wholebrain.R
alignment <- Sys.getenv("epoch")
regressor <- Sys.getenv("regressor")
message(paste0("Regressor: ", regressor))
debug = F
if (regressor=="entropy_change" | regressor=="abs_pe" | regressor=="reward" | regressor=="entropy_kld") {
=======

debug <- F
>>>>>>> a1dcac9d442691376262b6946cb7e273893b9b6f:meg/code/medusoid/meg_medusa_mixed_by_wholebrainSB.R
encode  <- T
rt_predict <- F
} else if (regressor=="rt") {
rt_predict <- T
encode <- F}
finish <- F
cat("Run encoding model: ", as.character(encode), "\n")
cat("Run rt prediction model: ", as.character(rt_predict), "\n")
domain = "tf"

if (debug) {alignment = "RT"}
# bin <- Sys.getenv("bin")
# 
stopifnot(alignment %in% c("RT", "clock", "feedback"))
if (whoami::username()=="ayd1") {
medusa_dir <- paste0("/bgfs/adombrovski/tfr_rds1/", alignment)
} else if (whoami::username()=="dnpl") {
  medusa_dir <- paste0("/proj/mnhallqlab/projects/Clock_MEG/atfr_rds/", alignment)
}


cat("Alignment: ", as.character(alignment), "\n")
cat("Domain: ", as.character(domain), "\n")

stopifnot(dir.exists(medusa_dir))  
setwd(medusa_dir)
message("Directory ", medusa_dir)
label_sensors = F
test = F
scale_winsor = F

ncores <- as.numeric(future::availableCores())
# core_prop <- as.numeric(Sys.getenv("core_prop"))
# if (whoami::username()=="dombax" | whoami::username() == "alexdombrovski") {
#   ncores <- floor(detectCores()*core_prop)
# }
# # Kai’s guidance on sensors is: ‘So for FEF, I say focus on 612/613, 543/542, 1022/1023, 
# # For IPS, 1823, 1822, 2222,2223.’
# fef_sensors <- c("0612","0613", "0542", "0543","1022")
# ips_sensors <- c("1823", "1822", "2222","2223")
# dan_sensors <- c(fef_sensors,ips_sensors)

sourcefilestart <- as.numeric(Sys.getenv("sourcefilestart"))
# incrementby <- as.numeric(Sys.getenv("incrementby"))
if (debug) {
sourcefilestart = 1
#  files <- "/bgfs/adombrovski/tfr_rds1/RT/temporal_l_group_freq_split_f_03.536.rds"
}
#if (finish) {
#  setwd("/bgfs/adombrovski/tfr_rds1/RT/results")
#  temp <- readRDS("meg_ddf_wholebrain_ec_rs.rds")
#  finished_files <- unique(temp$.filename)
#  setwd("/bgfs/adombrovski/tfr_rds1/RT")
#  all_files <- as.vector(list.files(pattern = "freq_t_all"))
#  remaining_files <- setdiff(all_files, finished_files)
#  files <- remaining_files[sourcefilestart:min((sourcefilestart + incrementby), length(remaining_files))]
#} else {
  setwd(medusa_dir)
  all_files <- list.files(pattern = "freq_t", full.names = T)
  files <- all_files[sourcefilestart]
#}


message(paste0("Processing files "))
cat(files)
# get, preprocess behavioral data
# for simplicity, change all matrix columns (wi-centering and scaling) to numeric
trial_df <- readRDS(behavioral_data_file) %>% as.data.frame(lapply(trial_df, function(x) {
  if (inherits(x, "matrix")) { x <- as.vector(x) }
  return(x)
})) %>% mutate(Subject = as.integer(id), Trial = trial, Run = run) %>% group_by(id, run) %>%
  mutate(abs_pe_lag = lag(abs_pe),
         v_entropy_wi_lag = lag(v_entropy_wi),
         rt_lag2_sc = lag(rt_lag_sc)
  ) %>% ungroup() 

get_kldsum <- function(v1, v2) {
  require(LaplacesDemon)
  stopifnot(length(v1) == length(v2))
  if (any(is.na(v1)) || any(is.na(v2))) { return(NA_real_) }
  kk <- KLD(v1, v2)
  return(kk$sum.KLD.px.py)
}

trial_df <- trial_df %>% group_by(id, run) %>% arrange(id, run, run_trial) %>% mutate(
  rt_lag2 = lag(rt_lag),
  rt_lag3 = lag(rt_lag2),
  rt_lag4 = lag(rt_lag3),
  rt_lag5 = lag(rt_lag4)) %>% ungroup() %>%
  rowwise() %>% mutate(
    kld4 = get_kldsum(c(rt_lag4, rt_lag3, rt_lag2, rt_lag), c(rt_lag5, rt_lag4, rt_lag3, rt_lag2)),
    kld3 = get_kldsum(c(rt_lag3, rt_lag2, rt_lag), c(rt_lag4, rt_lag3, rt_lag2))) %>%
  ungroup() %>% group_by(id, run) %>% mutate(kld3_lag = lag(kld3),
                                             kld4_lag = lag(kld4),
                                             kld3_cum2 = kld3 + kld3_lag,
                                             kld4_cum2 = kld4 + kld4_lag
  ) %>%
  ungroup() 

if (alignment=="RT" | alignment=="feedback") {
  # basic encoding with no random slopes
  encode_formula = formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + scale(rt_vmax_lag)  + scale(rt_vmax_change) + 
                             v_entropy_wi + v_entropy_wi_change + v_max_wi  + scale(abs_pe) + outcome + (1|Subject) + (1|Sensor))
  # random slopes of selected regressor
  if (regressor=="entropy") {
  encode_formula_rs = formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + 
                                  v_entropy_wi + scale(abs_pe) + outcome + (v_entropy_wi|Subject) + (v_entropy_wi|Sensor))
<<<<<<< HEAD:meg/code/medusoid/meg_medusa_mixed_by_wholebrain.R
  } else if (regressor=="entropy_change") {
  encode_formula_rs = formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + 
=======
# random slope of v_entropy
  encode_formula_rs_e_kld_sub = formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + kld3 +
                                  v_entropy_wi + scale(abs_pe) + outcome + (v_entropy_wi|Subject) + (1|Sensor))
  # random slope of v_entropy
  encode_formula_rs_ec = formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + 
>>>>>>> a1dcac9d442691376262b6946cb7e273893b9b6f:meg/code/medusoid/meg_medusa_mixed_by_wholebrainSB.R
                                   v_entropy_wi_change + scale(abs_pe) + outcome + (v_entropy_wi_change|Subject) + (v_entropy_wi_change|Sensor))
  } else if (regressor=="abs_pe") {
  encode_formula_rs = formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + 
                                   v_entropy_wi_change + scale(abs_pe) + outcome + (scale(abs_pe)|Subject) + (scale(abs_pe)|Sensor))
  } else if (regressor=="reward") {
  encode_formula_rs = formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + 
                                   v_entropy_wi_change + scale(abs_pe) + outcome + (outcome|Subject) + (outcome|Sensor))
  } else if (regressor=="entropy_kld") {
      encode_formula_rs_e = formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + kld3 +
                                  v_entropy_wi + scale(abs_pe) + outcome + (1|Subject) + (v_entropy_wi|Sensor))
  }
  rt_predict_formula = formula( ~ scale(Pow) * rt_csv_sc * outcome  + scale(Pow) * scale(rt_vmax)  +
                                  scale(Pow) * rt_lag_sc + (1|id) + (1|Sensor))
  # random slopes or RT_csv and RT_Vmax, note: no random slope of power
  rt_predict_formula_rs = formula( ~ scale(Pow) * rt_csv_sc * outcome  + scale(Pow) * scale(rt_vmax)  +
                                     scale(Pow) * rt_lag_sc + 
                                     (rt_csv_sc + scale(rt_vmax)|id))
  rt_outcome = "rt_next"
  c
} else if (alignment=="clock") {
  encode_formula = formula(~ rt_vmax + reward_lag + rt_csv_sc + rt_lag_sc + v_max_wi + trial_neg_inv_sc + 
                             v_entropy_wi + v_entropy_wi_change_lag + (1|Subject) + (1|Sensor))
<<<<<<< HEAD:meg/code/medusoid/meg_medusa_mixed_by_wholebrain.R
  if (regressor=="entropy") {
  encode_formula_rs = formula(~ reward_lag + rt_csv_sc + rt_lag_sc + trial_neg_inv_sc + 
                                  v_entropy_wi + (v_entropy_wi|Subject) + (v_entropy_wi|Sensor))
  } else if (regressor=="entropy_change") {
  encode_formula_rs = formula(~ reward_lag + rt_csv_sc + rt_lag_sc + trial_neg_inv_sc + 
=======
  # random slope of v_entropy
  encode_formula_rs_e =  formula(~ reward_lag + rt_csv_sc + rt_lag_sc + trial_neg_inv_sc +
                                  v_entropy_wi + (v_entropy_wi|Subject) + (v_entropy_wi|Sensor))
  encode_formula_rs_e_kld_sub =  formula(~ reward_lag + rt_csv_sc + rt_lag_sc + trial_neg_inv_sc + kld3 +
                                  v_entropy_wi + (v_entropy_wi|Subject) + (1|Sensor))

  # random slope of entropy_change
  encode_formula_rs_ec = formula(~ reward_lag + rt_csv_sc + rt_lag_sc + trial_neg_inv_sc + 
>>>>>>> a1dcac9d442691376262b6946cb7e273893b9b6f:meg/code/medusoid/meg_medusa_mixed_by_wholebrainSB.R
                                   v_entropy_wi_change_lag + (v_entropy_wi_change_lag|Subject) + (v_entropy_wi_change_lag|Sensor))
  } else if (regressor=="abs_pe") {
  encode_formula_rs = formula(~ reward_lag + rt_csv_sc + rt_lag_sc + trial_neg_inv_sc + scale(abs_pe_lag) +
                                   v_entropy_wi_change_lag + (scale(abs_pe_lag)|Subject) + (scale(abs_pe_lag)|Sensor))
  } else if (regressor=="reward") {
  encode_formula_rs = formula(~ reward_lag + rt_csv_sc + rt_lag_sc + trial_neg_inv_sc + scale(abs_pe_lag) +
                                   v_entropy_wi_change_lag + (reward_lag|Subject) + (reward_lag|Sensor))
  } else if (regressor=="entropy_kld") {
    encode_formula_rs =  formula(~ reward_lag + rt_csv_sc + rt_lag_sc + trial_neg_inv_sc + kld3 +
                                  v_entropy_wi + (1|Subject) + (v_entropy_wi|Sensor))
  }
  rt_predict_formula = formula( ~ scale(Pow) * rt_lag_sc * reward_lag  + scale(Pow) * scale(rt_vmax)  +
                                  (1|id) + (1|Sensor))
  # random slopes
  rt_predict_formula_rs = formula( ~ scale(Pow) * rt_lag_sc * reward_lag  + scale(Pow) * rt_vmax_lag_sc  + scale(Pow) * rt_lag2_sc +
                                     (rt_lag_sc + rt_vmax_lag_sc|id))
  rt_outcome = "rt_csv_sc"
}




#signal_outcome = "pow_scaled"
signal_outcome = "Pow"
#new approach: transform outcome variable at the time of computation
trans_func <- function(x) { DescTools::Winsorize(x, probs=c(.005, 1), na.rm=TRUE) }
#only drop bottom 0.5%

if (encode) {
  splits = c("Time", ".filename", "Freq")
  gc()
<<<<<<< HEAD:meg/code/medusoid/meg_medusa_mixed_by_wholebrain.R
  message(paste0("Using RHS formula: ", encode_formula_rs))
  ddf <- as_tibble(mixed_by(files, outcomes = signal_outcome, rhs_model_formulae = encode_formula_rs, split_on = splits,
                            external_df = trial_df, external_merge_by=c("Subject", "Run", "Trial"), padjust_by = "term", padjust_method = "BY", ncores = ncores,
                            refit_on_nonconvergence = 5, outcome_transform=trans_func, tidy_args=list(effects=c("fixed", "ran_vals", "ran_pars", "ran_coefs"), conf.int=TRUE)))
  saveRDS(ddf, file = paste0("meg_mixed_by_tf_ddf_wholebrain_", regressor, "_rs_single_", alignment, sourcefilestart))
=======
  ddf <- as_tibble(mixed_by(files, outcomes = signal_outcome, rhs_model_formulae = encode_formula_rs_e_kld, split_on = splits,
                            external_df = trial_df, external_merge_by=c("Subject", "Run", "Trial"), padjust_by = "term", padjust_method = "BY", ncores = ncores,
                            refit_on_nonconvergence = 5, outcome_transform=trans_func, tidy_args=list(effects=c("fixed", "ran_vals", "ran_pars", "ran_coefs"), conf.int=TRUE)))
  saveRDS(ddf, file = paste0("meg_mixed_by_tf_ddf_wholebrain_entropy_kld_sub_rs_single_", alignment, sourcefilestart))
>>>>>>> a1dcac9d442691376262b6946cb7e273893b9b6f:meg/code/medusoid/meg_medusa_mixed_by_wholebrainSB.R
}
# ddf <- as_tibble(mixed_by(files, outcomes = signal_outcome, rhs_model_formulae = encode_formula_rs_e, split_on = splits,
#                          external_df = trial_df, external_merge_by=c("Subject", "Run", "Trial"), padjust_by = "term", padjust_method = "BY", ncores = ncores,
#                          refit_on_nonconvergence = 5, outcome_transform=trans_func, tidy_args=list(effects=c("fixed", "ran_vals", "ran_pars", "ran_coefs"), conf.int=TRUE)))
# saveRDS(ddf, file = paste0("meg_mixed_by_tf_ddf_wholebrain_entropy_rs_single", alignment, sourcefilestart))
#}

if (rt_predict) {
  splits = c("Time", ".filename", "Freq", "Sensor")
  gc()
  message(paste0("Using RHS formula: ", rt_predict_formula_rs))
  rdf <- as_tibble(mixed_by(files, outcomes = rt_outcome, rhs_model_formulae = rt_predict_formula_rs , split_on = splits, external_df = trial_df,
                            padjust_by = "term", padjust_method = "BY", ncores = ncores, refit_on_nonconvergence = 5, outcome_transform=trans_func, 
                            tidy_args=list(effects=c("fixed", "ran_vals", "ran_pars", "ran_coefs"), conf.int=TRUE)))
  # rdf$sensor <- readr::parse_number(rdf$.filename)
  saveRDS(rdf, file = paste0("meg_tf_rdf_wholebrain_rt_rs_single_sensor_", alignment, sourcefilestart))
}
