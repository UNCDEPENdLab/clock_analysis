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
library(data.table)
sensor_list <- read.table("~/code/clock_analysis/meg/code/meg_sensors_annotated.txt", header=TRUE, colClasses="character") %>%
  filter(dan != "no") %>% pull(sensor)
sensor_map <- read.table("~/code/clock_analysis/meg/code/meg_sensors_annotated.txt", header=TRUE, colClasses="character")
repo_directory <- "~/code/clock_analysis"
behavioral_data_file <- "~/code/clock_analysis/meg/MEG_n63_behavioral_data_preprocessed_trial_df.RDS"
source("~/code/fmri.pipeline/R/mixed_by.R")

# main analysis analogous to Fig. 4 E-G in NComm 2020
debug = F #VERY CAREFUL, THIS MAKES IT RUN ON THE FIRST FILE ONLY
if (debug) {
  Sys.setenv(epoch = "RT")
  Sys.setenv(regressor = "abspe_by_rew")
}
alignment <- Sys.getenv("epoch")
regressor <- Sys.getenv("regressor")
message(paste0("Regressor: ", regressor))

<<<<<<< HEAD
if (regressor=="entropy_change" | regressor=="entropy_change_ri" | regressor == "entropy" | regressor=="abs_pe" | regressor == "entropy_change_full" | regressor == "entropy_change_sel" |
    regressor=="reward" | regressor=="reward_ri" | regressor=="entropy_kld" | regressor == "entropy_change_pos" | regressor == "entropy_change_neg" | regressor == "v_max" | regressor == "abspe_by_rew") {
=======
if (regressor=="entropy_change" | regressor == "entropy" | regressor=="abs_pe" | regressor == "entropy_change_full" | regressor == "entropy_change_sel" |
    regressor=="reward" | regressor=="entropy_kld" | regressor == "entropy_change_pos" | regressor == "entropy_change_neg" | regressor == "v_max" | 
    regressor == "abspe_by_rew" | regressor = "signed_pe") {
>>>>>>> e8da4b277597590eccabf5913765b757f3c426d1
  encode  <- T
  rt_predict <- F
} else if (regressor=="rt") {
  rt_predict <- T
  encode <- F}
finish <- F
cat("Run encoding model: ", as.character(encode), "\n")
cat("Run rt prediction model: ", as.character(rt_predict), "\n")
domain = "tf"

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
incrementby <- as.numeric(Sys.getenv("incrementby"))
if (debug) {
  sourcefilestart = 9
}
setwd(medusa_dir)
all_files <- list.files(pattern = "freq_t", full.names = T)
files <- all_files[sourcefilestart]
message(paste0("Starting at number ", sourcefilestart))
message(paste0("Processing file "))
cat(files)
# get behavioral data
## preprocessing has already been done, just historic here for documentation
## for simplicity, change all matrix columns (wi-centering and scaling) to numeric
# get_kldsum <- function(v1, v2) {
#   require(LaplacesDemon)
#   stopifnot(length(v1) == length(v2))
#   if (any(is.na(v1)) || any(is.na(v2))) { return(NA_real_) }
#   kk <- KLD(v1, v2)
#   return(kk$sum.KLD.px.py)
# }
trial_df <- readRDS(behavioral_data_file)
# back-calculate PE_max
# trial_df <- trial_df %>% group_by(id, run) %>% arrange(id, run, run_trial) %>% mutate(pe_max = abs_pe*reward_centered*2,
#                                                                                       pe_max_sc = scale(pe_max)) %>% ungroup()
# trial_df <- trial_df %>% 
#   as.data.table(lapply(trial_df, function(x) {
#   if (inherits(x, "matrix")) { x <- as.vector(x) }
#   return(x)
# })) %>%
#   mutate(Subject = as.integer(id), Trial = trial, Run = run) %>% group_by(id, run) %>%
#   mutate(abs_pe_lag = lag(abs_pe),
#          v_entropy_wi_lag = lag(v_entropy_wi),
#          rt_lag2_sc = lag(rt_lag_sc)
#   ) %>% ungroup() %>% group_by(id, run) %>% arrange(id, run, run_trial) %>% mutate(
#     rt_lag2 = lag(rt_lag),
#     rt_lag3 = lag(rt_lag2),
#     rt_lag4 = lag(rt_lag3),
#     rt_lag5 = lag(rt_lag4),
#     entropy_change_pos_lag = lag(entropy_change_pos_wi),
#     entropy_change_neg_lag = lag(entropy_change_neg_wi),
#     entropy_wi_change_lag_full = lag(v_entropy_wi_change_full)
#   ) %>% ungroup() %>%
#   rowwise() %>% mutate(
#     kld4 = get_kldsum(c(rt_lag4, rt_lag3, rt_lag2, rt_lag), c(rt_lag5, rt_lag4, rt_lag3, rt_lag2)),
#     kld3 = get_kldsum(c(rt_lag3, rt_lag2, rt_lag), c(rt_lag4, rt_lag3, rt_lag2))) %>%
#   ungroup()


if (alignment=="RT" | alignment=="feedback") {
  # basic encoding with no random slopes
  encode_formula = formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + scale(rt_vmax_lag)  + scale(rt_vmax_change) + 
                             v_entropy_wi + v_entropy_wi_change + v_max_wi  + scale(abs_pe) + outcome + (1|Subject) + (1|Sensor))
  # random slopes of selected regressor
  if (regressor=="entropy") {
    encode_formula_rs = formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + 
                                  v_entropy_wi + scale(abs_pe) + outcome + (v_entropy_wi|Subject) + (v_entropy_wi|Sensor))
  } else if (regressor=="entropy_change") {
    encode_formula_rs = formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + 
                                  v_entropy_wi_change + scale(abs_pe) + outcome + (v_entropy_wi_change|Subject) + (v_entropy_wi_change|Sensor))
  } else if (regressor=="entropy_change_ri") {
    encode_formula_ri = formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + 
                                  v_entropy_wi_change + scale(abs_pe) + outcome + (1|Subject) + (v_entropy_wi_change|Sensor))
  } else if (regressor=="abs_pe") {
    encode_formula_rs = formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + 
                                  v_entropy_wi_change + scale(abs_pe) + outcome + (scale(abs_pe)|Subject) + (scale(abs_pe)|Sensor))
  } else if (regressor=="abspe_by_rew") {
    encode_formula_rs = formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + 
                                  v_entropy_wi_change + abs_pe_sc * reward_centered + (1|Subject) + (reward_centered + abs_pe_sc|Sensor))
    emtrend_encode = "abs_pe_sc"
    emtrend_reward_centered = "reward_centered"
  } else if (regressor=="signed_pe") {
    encode_formula_rs = formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + 
                                  v_entropy_wi_change + pe_max_sc + (pe_max_sc|Subject) + (pe_max_sc|Sensor))
    emtrend_encode = "pe_max_sc"
  } else if (regressor=="reward") {
    encode_formula_rs = formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + 
                                  v_entropy_wi_change + scale(abs_pe) + outcome + (outcome|Subject) + (outcome|Sensor))
  } else if (regressor=="reward_ri") {
    encode_formula_rs = formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + 
                                  v_entropy_wi_change + scale(abs_pe) + outcome + (1|Subject) + (outcome|Sensor))
  } else if (regressor=="entropy_kld") {
    encode_formula_ri = formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + kld3 +
                                  v_entropy_wi + scale(abs_pe) + outcome + (1|Subject) + (v_entropy_wi|Sensor))
  } else if (regressor=="entropy_change_pos") {
    encode_formula_ri = formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + kld3 +
                                  entropy_change_pos_wi + scale(abs_pe) + outcome + (1|Subject) + (entropy_change_pos_wi|Sensor))
  } else if (regressor=="entropy_change_neg") {
    encode_formula_ri = formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + kld3 +
                                  entropy_change_neg_wi + scale(abs_pe) + outcome + (1|Subject) + (entropy_change_neg_wi|Sensor))
  } else if (regressor == "entropy_change_sel") { # version without subject random slope for speed
    encode_formula_rs = formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + 
                                  v_entropy_wi_change + scale(abs_pe) + outcome + (1|Subject) + (v_entropy_wi_change|Sensor))
  } else if (regressor == "entropy_change_full") { # version without subject random slope for speed
    encode_formula_rs = formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + 
                                  v_entropy_wi_change_full + scale(abs_pe) + outcome + (1|Subject) + (v_entropy_wi_change_full|Sensor))
  } else if (regressor == "v_max") {
    # run the strongest version of the model
    encode_formula_rs = formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + v_max_wi + reward_lag +
                                   v_entropy_wi_change + scale(abs_pe) + outcome + (v_max_wi|Subject) + (v_max_wi|Sensor))
  }
  rt_predict_formula = formula( ~ scale(Pow) * rt_csv_sc * outcome  + scale(Pow) * scale(rt_vmax)  +
                                  scale(Pow) * rt_lag_sc + (1|id) + (1|Sensor))
  # random slopes or RT_csv and RT_Vmax, note: no random slope of power
  rt_predict_formula_ri = formula( ~ Pow * rt_csv_sc * outcome  + Pow * rt_vmax  +
                                     Pow * rt_lag_sc + 
                                     (1|id))

  rt_predict_formula_rs = formula( ~ Pow * rt_csv_sc * outcome  + Pow * rt_vmax  +
                                     Pow * rt_lag_sc + 
                                     (rt_csv_sc + rt_vmax|id))

  rt_outcome = "rt_next"
  emtrend_rt = "rt_csv_sc"
  emtrend_rtvmax = "rt_vmax"
  emtrend_reward = "outcome"
} else if (alignment=="clock") {
  encode_formula = formula(~ rt_vmax + reward_lag + rt_csv_sc + rt_lag_sc + v_max_wi + trial_neg_inv_sc + 
                             v_entropy_wi + v_entropy_wi_change_lag + (1|Subject) + (1|Sensor))
  if (regressor=="entropy") {
    encode_formula_rs = formula(~ reward_lag + rt_csv_sc + rt_lag_sc + trial_neg_inv_sc + 
                                  v_entropy_wi + (v_entropy_wi|Subject) + (v_entropy_wi|Sensor))
  } else if (regressor=="entropy_change") {
    encode_formula_rs = formula(~ reward_lag + rt_csv_sc + rt_lag_sc + trial_neg_inv_sc + 
                                  v_entropy_wi_change_lag + (v_entropy_wi_change_lag|Subject) + (v_entropy_wi_change_lag|Sensor))
  } else if (regressor=="entropy_change_ri") {
    encode_formula_rs = formula(~ reward_lag + rt_csv_sc + rt_lag_sc + trial_neg_inv_sc + 
                                  v_entropy_wi_change_lag + (1|Subject) + (v_entropy_wi_change_lag|Sensor))
  } else if (regressor=="abs_pe") {
    encode_formula_rs = formula(~ reward_lag + rt_csv_sc + rt_lag_sc + trial_neg_inv_sc + scale(abs_pe_lag) +
                                  v_entropy_wi_change_lag + (scale(abs_pe_lag)|Subject) + (scale(abs_pe_lag)|Sensor))
  } else if (regressor=="reward") {
    encode_formula_rs = formula(~ reward_lag + rt_csv_sc + rt_lag_sc + trial_neg_inv_sc + scale(abs_pe_lag) +
                                  v_entropy_wi_change_lag + (reward_lag|Subject) + (reward_lag|Sensor))
  } else if (regressor=="entropy_kld") {
    encode_formula_rs =  formula(~ reward_lag + rt_csv_sc + rt_lag_sc + trial_neg_inv_sc + kld3 +
                                   v_entropy_wi + (1|Subject) + (v_entropy_wi|Sensor))
  } else if (regressor=="entropy_change_pos") {
    encode_formula_rs =  formula(~ reward_lag + rt_csv_sc + rt_lag_sc + trial_neg_inv_sc + kld3 +
                                   entropy_change_pos_lag + (1|Subject) + (entropy_change_pos_lag|Sensor))
  } else if (regressor=="entropy_change_neg") {
    encode_formula_rs =  formula(~ reward_lag + rt_csv_sc + rt_lag_sc + trial_neg_inv_sc + kld3 +
                                   entropy_change_neg_lag + (1|Subject) + (entropy_change_neg_lag|Sensor))
                                   
  } else if (regressor=="entropy_change_sel") {
    encode_formula_rs = formula(~ reward_lag + rt_csv_sc + rt_lag_sc + trial_neg_inv_sc + 
                                  v_entropy_wi_change_lag + (1|Subject) + (v_entropy_wi_change_lag|Sensor))
  } else if (regressor=="entropy_change_full") {
    encode_formula_rs = formula(~ reward_lag + rt_csv_sc + rt_lag_sc + trial_neg_inv_sc + 
                                  entropy_wi_change_lag_full + (1|Subject) + (entropy_wi_change_lag_full|Sensor))
  } else if (regressor == "v_max") {
    encode_formula_rs = formula(~ reward_lag + rt_csv_sc + rt_lag_sc + trial_neg_inv_sc + scale(abs_pe_lag) + v_max_wi + 
                                  v_entropy_wi_change_lag + (v_max_wi|Subject) + (v_max_wi|Sensor))
  }
  
  rt_predict_formula = formula( ~ Pow * rt_lag_sc * reward_lag  + Pow * rt_vmax  +
                                  (1|id) + (1|Sensor))
  # random slopes
  rt_predict_formula_ri = formula( ~ Pow * rt_lag_sc * reward_lag  + Pow * rt_vmax_lag_sc  + Pow * rt_lag2_sc +
                                     (1|id))
  rt_predict_formula_rs = formula( ~ Pow * rt_lag_sc * reward_lag  + Pow * rt_vmax_lag_sc  + Pow * rt_lag2_sc +
                                     (rt_lag_sc + rt_vmax_lag_sc|id))
  rt_outcome = "rt_csv_sc"
  emtrend_rt = "rt_lag_sc"
  emtrend_rtvmax = "rt_vmax_lag_sc"
  emtrend_reward = "reward_lag"
}



#signal_outcome = "pow_scaled"
signal_outcome = "Pow"
#new approach: transform outcome variable at the time of computation
trans_func <- function(x) { DescTools::Winsorize(x, probs=c(.005, 1), na.rm=TRUE) }
#only drop bottom 0.5%

if (encode) {
  splits = c("Time", ".filename", "Freq")
  gc()
  if (str_detect(regressor, "_ri")) {
      message(paste0("Using RHS formula: ", encode_formula_ri))
  } else {
  message(paste0("Using RHS formula: ", encode_formula_rs))}
  ddf <- mixed_by(files, outcomes = signal_outcome, rhs_model_formulae = list(ri = encode_formula_ri), split_on = splits,
                            external_df = trial_df, external_merge_by=c("Subject", "Run", "Trial"), padjust_by = "term", padjust_method = "BY", ncores = ncores,
                            refit_on_nonconvergence = 5, outcome_transform=trans_func, tidy_args=list(effects=c("fixed", "ran_vals", "ran_pars", "ran_coefs"), conf.int=TRUE,
                                                                                                      calculate =c("parameter_estimates_reml","fit_statistics")) #,
                            #emtrends_spec = list(
                            #  list(outcome=signal_outcome, model_name="ri", var=emtrend_encode, specs=c(emtrend_reward_centered), at = list(reward_centered = c(-0.5, 0.5))))
                              )
    if (str_detect(regressor, "_ri")) {
  saveRDS(ddf, file = paste0("meg_mixed_by_tf_ddf_wholebrain_", regressor, "_single_", alignment, sourcefilestart))} else {
  saveRDS(ddf, file = paste0("meg_mixed_by_tf_ddf_wholebrain_", regressor, "_rs_single_", alignment, sourcefilestart))}


if (rt_predict) {
  splits = c("Time", ".filename", "Freq", "Sensor")
  gc()
  message(paste0("Using RHS formula: ", rt_predict_formula_rs))
  rsdf <- mixed_by(files, outcomes = rt_outcome, rhs_model_formulae = list(rt_single_sensor_rs=rt_predict_formula_rs), split_on = splits, external_df = trial_df,
                            padjust_by = "term", padjust_method = "BY", ncores = ncores, refit_on_nonconvergence = 5, outcome_transform=trans_func, 
                            tidy_args=list(effects=c("fixed", "ran_vals", "ran_pars", "ran_coefs"), conf.int=TRUE, calculate=c("parameter_estimates_reml")), 
                            emtrends_spec = list(
                              list(outcome=rt_outcome, model_name="rt_single_sensor_rs", var=emtrend_rt, specs=c("Pow", emtrend_reward), at = list(Pow = c(-2,2))), 
                           list(outcome=rt_outcome, model_name="rt_single_sensor_rs", var=emtrend_rtvmax, specs=c("Pow", emtrend_reward), at = list(Pow = c(-2,2)))
                          ), scale_predictors = "Pow")
  saveRDS(rsdf, file = paste0("meg_tf_rdf_wholebrain_rt_rs_single_sensor_", alignment, sourcefilestart))
  message(paste0("Using RHS formula: ", rt_predict_formula_ri))
  ridf <- mixed_by(files, outcomes = rt_outcome, rhs_model_formulae = list(rt_single_sensor_ri=rt_predict_formula_ri), split_on = splits, external_df = trial_df,
                  padjust_by = "term", padjust_method = "BY", ncores = ncores, refit_on_nonconvergence = 5, outcome_transform=trans_func, 
                  tidy_args=list(effects=c("fixed", "ran_vals", "ran_pars", "ran_coefs"), conf.int=TRUE, calculate=c("parameter_estimates_reml")), 
                  emtrends_spec = list(
                    list(outcome=rt_outcome, model_name="rt_single_sensor_ri", var=emtrend_rt, specs=c("Pow", emtrend_reward), at = list(Pow = c(-2,2))), 
                    list(outcome=rt_outcome, model_name="rt_single_sensor_ri", var=emtrend_rtvmax, specs=c("Pow", emtrend_reward), at = list(Pow = c(-2,2)))
                  ), scale_predictors = "Pow")
  saveRDS(ridf, file = paste0("meg_tf_rdf_wholebrain_rt_ri_single_sensor_", alignment, sourcefilestart))
  
}
