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
encode  <- TRUE
rt_predict <- TRUE
cat("Run encoding model: ", as.character(encode), "\n")
cat("Run rt prediction model: ", as.character(rt_predict), "\n")
domain = "tf"
alignment <- Sys.getenv("epoch")
bin <- Sys.getenv("bin")

 stopifnot(alignment %in% c("RT", "clock", "feedback"))
      
if (whoami::username() == "ayd1") {medusa_dir <- paste0("/bgfs/adombrovski/tfr_rds1/", alignment, "/grouped_tf")}
      

cat("Alignment: ", as.character(alignment), "\n")
cat("Domain: ", as.character(domain), "\n")

stopifnot(dir.exists(medusa_dir))  
setwd(medusa_dir)

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

files <- paste0(Sys.getenv("sourcefile"))
message(paste0("Processing file ",files))

# get behavioral data
trial_df <- readRDS(behavioral_data_file)

#for simplicity, change all matrix columns (wi-centering and scaling) to numeric
trial_df <- as.data.frame(lapply(trial_df, function(x) {
  if (inherits(x, "matrix")) { x <- as.vector(x) }
  return(x)
}))

# write random subject-level vectors for sanity checks
trial_df <- trial_df %>% group_by(id) %>%
  mutate(rand1 = rnorm(1),
         rand2 = rnorm(1),
         rand3 = rnorm(1),
         rand4 = rnorm(1))

# trial_df$Subject <- trial_df$id
# trial_df$Run <- trial_df$run
# trial_df$Trial <- trial_df$trial
# trial_df$rt_next_sc <- scale(trial_df$rt_next)
# # # save behavioral data with new variables
# saveRDS(trial_df, behavioral_data_file)


if (alignment=="RT" | alignment=="feedback") {
  encode_formula = formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + scale(rt_vmax_lag)  + scale(rt_vmax_change) + 
                             v_entropy_wi + v_entropy_wi_change + v_max_wi  + scale(abs_pe) + outcome + (1|Subject) + (1|sensor))
  # random slope of v_entropy
  encode_formula_rs_e = formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + 
                             v_entropy_wi + scale(abs_pe) + outcome + (v_entropy_wi|Subject) + (v_entropy_wi|sensor))
  
} else if (alignment=="clock") {
  encode_formula = formula(~ rt_vmax + reward_lag + rt_csv_sc + rt_lag_sc + v_max_wi + trial_neg_inv_sc + 
                             v_entropy_wi + v_entropy_wi_change_lag + (1|Subject) + (1|sensor))
  # random slope of v_entropy
  encode_formula_rs_e = formula(~ reward_lag + rt_csv_sc + rt_lag_sc + trial_neg_inv_sc + 
                             v_entropy_wi + (v_entropy_wi|Subject) + (v_entropy_wi|sensor))
  
}


# preproc trial-df
trial_df <- trial_df %>% mutate(Subject = as.integer(id), Trial = trial, Run = run) %>% group_by(id, run) %>%
  mutate(abs_pe_lag = lag(abs_pe),
         v_entropy_wi_lag = lag(v_entropy_wi)) %>% ungroup()
message(str(trial_df))
if (alignment == "RT" | alignment == "feedback") {
  rt_outcome = "rt_next"
  encode_formula_e = formula(~ scale(rt_vmax_lag)*echange_f1_early + scale(rt_vmax_lag)*echange_f2_late + scale(rt_vmax_lag)*e_f1 + scale(rt_vmax_lag)*abs_pe_f2_early + scale(rt_vmax_lag)*abs_pe_f1_mid + scale(rt_vmax_lag)*abs_pe_f3_late +
                               scale(abs_pe)*echange_f1_early + scale(abs_pe)*echange_f2_late + scale(abs_pe)*e_f1 +
                               outcome*echange_f1_early + outcome*echange_f2_late + outcome*e_f1 + outcome*abs_pe_f2_early + outcome*abs_pe_f1_mid + outcome*abs_pe_f3_late +
                               rt_csv_sc*echange_f1_early + rt_csv_sc*echange_f2_late + rt_csv_sc*e_f1 + rt_csv_sc*abs_pe_f2_early + rt_csv_sc*abs_pe_f1_mid + rt_csv_sc*abs_pe_f3_late +
                               trial_neg_inv_sc*echange_f1_early + trial_neg_inv_sc*echange_f2_late + trial_neg_inv_sc*e_f1 + trial_neg_inv_sc*abs_pe_f2_early + trial_neg_inv_sc*abs_pe_f1_mid + trial_neg_inv_sc*abs_pe_f3_late +
                               v_entropy_wi_change*echange_f1_early + v_entropy_wi_change*echange_f2_late + v_entropy_wi*e_f1 + 
                               scale(abs_pe)*abs_pe_f2_early + scale(abs_pe)*abs_pe_f1_mid + scale(abs_pe)*abs_pe_f3_late +(1|Subject) + (1|sensor))
  encode_formula_pe = formula(~ scale(rt_vmax_lag)*abs_pe_f2_early + scale(rt_vmax_lag)*abs_pe_f1_mid + scale(rt_vmax_lag)*abs_pe_f3_late +
                                scale(abs_pe)*abs_pe_f2_early + scale(abs_pe)*abs_pe_f1_mid + scale(abs_pe)*abs_pe_f3_late +
                                outcome*abs_pe_f2_early + outcome*abs_pe_f1_mid + outcome*abs_pe_f3_late +   
                                rt_csv_sc*abs_pe_f2_early + rt_csv_sc*abs_pe_f1_mid + rt_csv_sc*abs_pe_f3_late +
                                trial_neg_inv_sc*abs_pe_f2_early + trial_neg_inv_sc*abs_pe_f1_mid + trial_neg_inv_sc*abs_pe_f3_late +
                                v_entropy_wi_change*abs_pe_f2_early + v_entropy_wi_change*abs_pe_f1_mid + v_entropy_wi_change*abs_pe_f3_late + 
                                v_entropy_wi*abs_pe_f2_early + v_entropy_wi*abs_pe_f1_mid + v_entropy_wi*abs_pe_f3_late + 
                                rt_lag_sc*abs_pe_f2_early + rt_lag_sc*abs_pe_f1_mid + rt_lag_sc*abs_pe_f3_late + (1|Subject) + (1|sensor))
  encode_formula_r = formula(~ scale(rt_vmax_lag)*rand1 + scale(rt_vmax_lag)*rand2 + scale(rt_vmax_lag)*rand3 +
                                scale(abs_pe)*rand1 + scale(abs_pe)*rand2 + scale(abs_pe)*rand3 +
                                outcome*rand1 + outcome*rand2 + outcome*rand3 +   
                                rt_csv_sc*rand1 + rt_csv_sc*rand2 + rt_csv_sc*rand3 +
                                trial_neg_inv_sc*rand1 + trial_neg_inv_sc*rand2 + trial_neg_inv_sc*rand3 +
                                v_entropy_wi_change*rand1 + v_entropy_wi_change*rand2 + v_entropy_wi_change*rand3 + 
                                v_entropy_wi*rand1 + v_entropy_wi*rand2 + v_entropy_wi*rand3 + 
                                rt_lag_sc*rand1 + rt_lag_sc*rand2 + rt_lag_sc*rand3 + (1|Subject) + (1|sensor))
  
} else if (alignment == "clock") {
  rt_outcome = "rt_csv_sc"
  encode_formula_e = formula(~ scale(rt_vmax_lag)*echange_f1_early + scale(rt_vmax_lag)*echange_f2_late + scale(rt_vmax_lag)*e_f1 + scale(rt_vmax_lag)*abs_pe_f2_early + scale(rt_vmax_lag)*abs_pe_f1_mid + scale(rt_vmax_lag)*abs_pe_f3_late +
                               scale(abs_pe_lag)*echange_f1_early + scale(abs_pe_lag)*echange_f2_late + scale(abs_pe_lag)*e_f1 +
                               reward_lag*echange_f1_early + reward_lag*echange_f2_late + reward_lag*e_f1 + reward_lag*abs_pe_f2_early + reward_lag*abs_pe_f1_mid + reward_lag*abs_pe_f3_late +
                               rt_lag_sc*echange_f1_early + rt_lag_sc*echange_f2_late + rt_lag_sc*e_f1 + rt_lag_sc*abs_pe_f2_early + rt_lag_sc*abs_pe_f1_mid + rt_lag_sc*abs_pe_f3_late +
                               trial_neg_inv_sc*echange_f1_early + trial_neg_inv_sc*echange_f2_late + trial_neg_inv_sc*e_f1 + trial_neg_inv_sc*abs_pe_f2_early + trial_neg_inv_sc*abs_pe_f1_mid + trial_neg_inv_sc*abs_pe_f3_late +
                               v_entropy_wi_change*echange_f1_early + v_entropy_wi_change*echange_f2_late + v_entropy_wi*e_f1 + 
                               scale(abs_pe_lag)*abs_pe_f2_early + scale(abs_pe_lag)*abs_pe_f1_mid + scale(abs_pe_lag)*abs_pe_f3_late +(1|Subject) + (1|sensor))
  
encode_formula_pe = formula(~ scale(rt_vmax_lag)*abs_pe_f2_early + scale(rt_vmax_lag)*abs_pe_f1_mid + scale(rt_vmax_lag)*abs_pe_f3_late +
                                scale(abs_pe)*abs_pe_f2_early + scale(abs_pe)*abs_pe_f1_mid + scale(abs_pe)*abs_pe_f3_late +
                                outcome*abs_pe_f2_early + outcome*abs_pe_f1_mid + outcome*abs_pe_f3_late +   
                                rt_csv_sc*abs_pe_f2_early + rt_csv_sc*abs_pe_f1_mid + rt_csv_sc*abs_pe_f3_late +
                                trial_neg_inv_sc*abs_pe_f2_early + trial_neg_inv_sc*abs_pe_f1_mid + trial_neg_inv_sc*abs_pe_f3_late +
                                v_entropy_wi_change*abs_pe_f2_early + v_entropy_wi_change*abs_pe_f1_mid + v_entropy_wi_change*abs_pe_f3_late + 
                                v_entropy_wi*abs_pe_f2_early + v_entropy_wi*abs_pe_f1_mid + v_entropy_wi*abs_pe_f3_late + 
                                rt_lag_sc*abs_pe_f2_early + rt_lag_sc*abs_pe_f1_mid + rt_lag_sc*abs_pe_f3_late + (1|Subject) + (1|sensor))
encode_formula_r = formula(~ scale(rt_vmax_lag)*rand1 + scale(rt_vmax_lag)*rand2 + scale(rt_vmax_lag)*rand3 +
                             scale(abs_pe)*rand1 + scale(abs_pe)*rand2 + scale(abs_pe)*rand3 +
                             outcome*rand1 + outcome*rand2 + outcome*rand3 +   
                             rt_csv_sc*rand1 + rt_csv_sc*rand2 + rt_csv_sc*rand3 +
                             trial_neg_inv_sc*rand1 + trial_neg_inv_sc*rand2 + trial_neg_inv_sc*rand3 +
                             v_entropy_wi_change*rand1 + v_entropy_wi_change*rand2 + v_entropy_wi_change*rand3 + 
                             v_entropy_wi*rand1 + v_entropy_wi*rand2 + v_entropy_wi*rand3 + 
                             rt_lag_sc*rand1 + rt_lag_sc*rand2 + rt_lag_sc*rand3 + (1|Subject) + (1|sensor))

}
if (domain == "tf") {
  splits = c("Time", ".filename", "Freq")
  #signal_outcome = "pow_scaled"
  signal_outcome = "Pow"
  #new approach: transform outcome variable at the time of computation
  trans_func <- function(x) { DescTools::Winsorize(x, probs=c(.005, 1), na.rm=TRUE) }
  #only drop bottom 0.5%
  if (alignment == "RT" | alignment == "feedback") {
    rt_predict_formula = formula( ~ scale(Pow) * rt_csv_sc * outcome  + scale(Pow) * scale(rt_vmax)  +
                                    scale(Pow) * rt_lag_sc + (1|id) + (1|sensor))
  } else {
    rt_predict_formula = formula( ~ scale(Pow) * rt_lag_sc * outcome  + scale(Pow) * scale(rt_vmax)  +
                                    (1|id) + (1|sensor))}
} 
  

gc()
if (encode) {
  ddf <- as_tibble(mixed_by(files, outcomes = signal_outcome, rhs_model_formulae = encode_formula_rs_e, split_on = splits,
                            external_df = trial_df, external_merge_by=c("Subject", "Run", "Trial"), padjust_by = "term", padjust_method = "BY", ncores = ncores,
                            refit_on_nonconvergence = 5, outcome_transform=trans_func, tidy_args=list(effects=c("fixed", "random"), conf.int=TRUE)))
  # refit_on_nonconvergence = 5))
  # ddf$sensor <- readr::parse_number(ddf$.filename)
  # save output
  #setwd("~/OneDrive/collected_letters/papers/meg/plots/clock_encode/")
    # saveRDS(ddf, file = "meg_mixed_by_tf_ranefs_mult_interactions_e_ddf.RDS")
    # saveRDS(ddf, file = "meg_mixed_by_tf_clock_ddf.RDS")
    saveRDS(ddf, file = paste0("meg_mixed_by_tf_ddf_combined_entropy_rs", alignment,basename(files)))
  
}  
gc()
if (rt_predict) {
  rdf <- as_tibble(mixed_by(files, outcomes = rt_outcome, rhs_model_formulae = rt_predict_formula , split_on = splits, external_df = trial_df,
                            padjust_by = "term", padjust_method = "BY", ncores = ncores, refit_on_nonconvergence = 5, outcome_transform=trans_func))
  # rdf$sensor <- readr::parse_number(rdf$.filename)
  

    saveRDS(rdf, file = paste0("meg_mixed_by_tf_ddf_combined_rt_rs", alignment, basename(files)))
  }




### Deprecated stuff
## if (label_sensors || scale_winsor) {
##   # make cluster ----
##   library(parallel)
##   ncores <- detectCores()
##   if (ncores > 1L) {
##     cl <- makeCluster(ncores)
##     registerDoParallel(cl)
##     on.exit(try(stopCluster(cl)))
##   } else {
##     registerDoSEQ()
##   }
## }

## if (label_sensors) {
##   for (this_file in files) {
##     d <- readRDS(this_file)
##     d$sensor <- str_extract(this_file, "[[:digit:]]{4}")
##     saveRDS(d, file = this_file)
##     rm(d)
##   }
## }

## scale2 <- function(x, na.rm = FALSE) (x - mean(x, na.rm = na.rm)) / sd(x, na.rm)
## if (scale_winsor) {
##   if (domain == "time") {
##     foreach(i = 1:length(files), .packages=c("tidyverse", "psych")) %dopar% {
##       d <- readRDS(files[i])
##       d$signal_scaled <- winsor(scale2(d$Signal), trim = .01)
##       saveRDS(d, file = files[i])
##       return(NULL)
##     }    
##   } else if (domain == "tf") {
##     foreach(i = 1:length(files), .packages=c("tidyverse", "psych")) %dopar% {
##       scale2 <- function(x, na.rm = FALSE) (x - mean(x, na.rm = na.rm)) / sd(x, na.rm)
##       d <- readRDS(files[i])
##       d$pow_scaled <- scale2(winsor(d$Pow, trim = .005))
##       saveRDS(d, file = files[i])
##       return(NULL)
##     }
##   }
##   try(stopCluster(cl))
## }


# files <- files[1] # TEST ONLY

# # take first few for testing
# all_sensors <- all_sensors[1:4]


# 
# sample_data <- readRDS(files[2])
# sample_data$pow_scaled <- winsor((sample_data$Pow), trim = .005)
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
