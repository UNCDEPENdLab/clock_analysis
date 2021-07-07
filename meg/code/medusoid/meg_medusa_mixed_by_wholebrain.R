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

debug = F
encode  <- T
rt_predict <- F
finish <- T
cat("Run encoding model: ", as.character(encode), "\n")
cat("Run rt prediction model: ", as.character(rt_predict), "\n")
domain = "tf"
alignment <- Sys.getenv("epoch")
if (debug) {alignment = "RT"}
# bin <- Sys.getenv("bin")
# 
stopifnot(alignment %in% c("RT", "clock", "feedback"))

medusa_dir <- paste0("/bgfs/adombrovski/tfr_rds1/", alignment)


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

if (finish) {
  setwd("/bgfs/adombrovski/tfr_rds1/RT/results")
  temp <- readRDS("meg_ddf_wholebrain_ec_rs.rds")
  finished_files <- unique(temp$.filename)
  setwd("/bgfs/adombrovski/tfr_rds1/RT")
  all_files <- as.vector(list.files(pattern = "freq_t_all"))
  remaining_files <- setdiff(all_files, finished_files)
  files <- remaining_files[sourcefilestart:min((sourcefilestart + incrementby), length(remaining_files))]
} else {
  all_files <- list.files(pattern = "freq_t", full.names = T)
   files <- all_files[sourcefilestart:min((sourcefilestart + incrementby), length(all_files))]
}

if (debug) {
  files <- "/bgfs/adombrovski/tfr_rds1/RT/temporal_l_group_freq_split_f_03.536.rds"
}
message(paste0("Processing files "))
cat(files)
# get behavioral data
trial_df <- readRDS(behavioral_data_file)

#for simplicity, change all matrix columns (wi-centering and scaling) to numeric
trial_df <- as.data.frame(lapply(trial_df, function(x) {
  if (inherits(x, "matrix")) { x <- as.vector(x) }
  return(x)
}))

# write random subject-level vectors for sanity checks
# preproc trial-df
trial_df <- trial_df %>% mutate(Subject = as.integer(id), Trial = trial, Run = run) %>% group_by(id, run) %>%
  mutate(abs_pe_lag = lag(abs_pe),
         v_entropy_wi_lag = lag(v_entropy_wi),
         rt_lag2_sc = lag(rt_lag_sc)
  ) %>% ungroup() %>% group_by(id) %>%
  mutate(rand1 = rnorm(1),
         rand2 = rnorm(1),
         rand3 = rnorm(1),
         rand4 = rnorm(1))
#message(str(trial_df))
# trial_df$Subject <- trial_df$id
# trial_df$Run <- trial_df$run
# trial_df$Trial <- trial_df$trial
# trial_df$rt_next_sc <- scale(trial_df$rt_next)
# # # save behavioral data with new variables
# saveRDS(trial_df, behavioral_data_file)


if (alignment=="RT" | alignment=="feedback") {
  encode_formula = formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + scale(rt_vmax_lag)  + scale(rt_vmax_change) + 
                             v_entropy_wi + v_entropy_wi_change + v_max_wi  + scale(abs_pe) + outcome + (1|Subject) + (1|Sensor))
  # random slope of v_entropy
  encode_formula_rs_e = formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + 
                                  v_entropy_wi + scale(abs_pe) + outcome + (v_entropy_wi|Subject) + (v_entropy_wi|Sensor))
  # random slope of v_entropy
  encode_formula_rs_ec = formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + 
                                   v_entropy_wi_change + scale(abs_pe) + outcome + (v_entropy_wi_change|Subject) + (v_entropy_wi_change|Sensor))
  rt_predict_formula = formula( ~ scale(Pow) * rt_csv_sc * outcome  + scale(Pow) * scale(rt_vmax)  +
                                  scale(Pow) * rt_lag_sc + (1|id) + (1|Sensor))
  # random slopes or RT_csv and RT_Vmax
  rt_predict_formula_rs = formula( ~ scale(Pow) * rt_csv_sc * outcome  + scale(Pow) * scale(rt_vmax)  +
                                     scale(Pow) * rt_lag_sc + 
                                     (rt_csv_sc + scale(rt_vmax)|id) + (scale(Pow) * rt_csv_sc + scale(Pow) * scale(rt_vmax)|Sensor))
  rt_outcome = "rt_next"
  
} else if (alignment=="clock") {
  encode_formula = formula(~ rt_vmax + reward_lag + rt_csv_sc + rt_lag_sc + v_max_wi + trial_neg_inv_sc + 
                             v_entropy_wi + v_entropy_wi_change_lag + (1|Subject) + (1|Sensor))
  # random slope of v_entropy
  encode_formula_rs_e = formula(~ reward_lag + rt_csv_sc + rt_lag_sc + trial_neg_inv_sc + 
                                  v_entropy_wi + (v_entropy_wi|Subject) + (v_entropy_wi|Sensor))
  # random slope of entropy_change
  encode_formula_rs_ec = formula(~ reward_lag + rt_csv_sc + rt_lag_sc + trial_neg_inv_sc + 
                                   v_entropy_wi_change_lag + (v_entropy_wi_change_lag|Subject) + (v_entropy_wi_change_lag|Sensor))
  rt_predict_formula = formula( ~ scale(Pow) * rt_lag_sc * reward_lag  + scale(Pow) * scale(rt_vmax)  +
                                  (1|id) + (1|Sensor))
  # random slopes
  rt_predict_formula_rs = formula( ~ scale(Pow) * rt_lag_sc * reward_lag  + scale(Pow) * rt_vmax_lag_sc  + scale(Pow) * rt_lag2_sc +
                                     (rt_lag_sc + rt_vmax_lag_sc|id) + (scale(Pow) * rt_lag_sc + scale(Pow) * rt_vmax_lag_sc|Sensor))
  rt_outcome = "rt_csv_sc"
}




splits = c("Time", ".filename", "Freq")
#signal_outcome = "pow_scaled"
signal_outcome = "Pow"
#new approach: transform outcome variable at the time of computation
trans_func <- function(x) { DescTools::Winsorize(x, probs=c(.005, 1), na.rm=TRUE) }
#only drop bottom 0.5%

if (encode) {
  gc()
  ddf <- as_tibble(mixed_by(files, outcomes = signal_outcome, rhs_model_formulae = encode_formula_rs_ec, split_on = splits,
                            external_df = trial_df, external_merge_by=c("Subject", "Run", "Trial"), padjust_by = "term", padjust_method = "BY", ncores = ncores,
                            refit_on_nonconvergence = 5, outcome_transform=trans_func, tidy_args=list(effects=c("fixed", "ran_vals", "ran_pars", "ran_coefs"), conf.int=TRUE)))
  saveRDS(ddf, file = paste0("meg_mixed_by_tf_ddf_wholebrain_entropy_change_rs_", alignment, sourcefilestart))
}
# ddf <- as_tibble(mixed_by(files, outcomes = signal_outcome, rhs_model_formulae = encode_formula_rs_e, split_on = splits,
#                          external_df = trial_df, external_merge_by=c("Subject", "Run", "Trial"), padjust_by = "term", padjust_method = "BY", ncores = ncores,
#                          refit_on_nonconvergence = 5, outcome_transform=trans_func, tidy_args=list(effects=c("fixed", "ran_vals", "ran_pars", "ran_coefs"), conf.int=TRUE)))
# saveRDS(ddf, file = paste0("meg_mixed_by_tf_ddf_combined_entropy_rs_", alignment, basename(files)))


if (rt_predict) {
  gc()
  rdf <- as_tibble(mixed_by(files, outcomes = rt_outcome, rhs_model_formulae = rt_predict_formula_rs , split_on = splits, external_df = trial_df,
                            padjust_by = "term", padjust_method = "BY", ncores = ncores, refit_on_nonconvergence = 5, outcome_transform=trans_func, 
                            tidy_args=list(effects=c("fixed", "ran_vals", "ran_pars", "ran_coefs"), conf.int=TRUE)))
  # rdf$sensor <- readr::parse_number(rdf$.filename)
  saveRDS(rdf, file = paste0("meg_tf_rdf_wholebrain_rt_rs_int", alignment, sub(".rds", "", sub("_group_freq_split", "", sourcefilestart))))
}
