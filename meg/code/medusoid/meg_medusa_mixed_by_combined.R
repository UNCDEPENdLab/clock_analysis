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

if (whoami::username() == "dombax" || whoami::username() == "alexdombrovski") {
  sensor_list <- read.table("~/code/clock_analysis/meg/code/meg_sensors_annotated.txt", header=TRUE, colClasses="character") %>%
    filter(dan != "no") %>% pull(sensor)
  sensor_map <- read.table("~/code/clock_analysis/meg/code/meg_sensors_annotated.txt", header=TRUE, colClasses="character")
  repo_directory <- "~/code/clock_analysis"
  behavioral_data_file <- "~/code/clock_analysis/meg/MEG_n63_behavioral_data_preprocessed_trial_df.RDS"
  source("~/code/fmri.pipeline/R/mixed_by.R")
} else {
  sensor_list <- read.table("/proj/mnhallqlab/projects/clock_analysis/meg/code/meg_sensors_annotated.txt", header=TRUE, colClasses="character") %>%
    filter(dan != "no") %>% pull(sensor)
  sensor_map <- read.table("/proj/mnhallqlab/projects/clock_analysis/meg/code/meg_sensors_annotated.txt", header=TRUE, colClasses="character") 
  repo_directory <- "/proj/mnhallqlab/projects/clock_analysis"
  behavioral_data_file <- "/proj/mnhallqlab/projects/clock_analysis/meg/code/medusoid/MEG_n63_behavioral_data_preprocessed_trial_df.RDS"
  
  source("/proj/mnhallqlab/users/michael/fmri.pipeline/R/mixed_by.R")
}

#repo_directory <- "~/code/clock_analysis"
# <- _dir = "~/Box/SCEPTIC_fMRI/MEG_20Hz_n63_clockalign/"
#diag_dir =  "~/Box/SCEPTIC_fMRI/MEG_20Hz_n63_clockalign/diags/"
#behavioral_data_file = "~/Box/SCEPTIC_fMRI/sceptic_model_fits/MEG_n63_behavioral_data_preprocessed_trial_df.RDS"
#source("~/code/fmri.pipeline/R/mixed_by.R")

# options, files ----
plots = F

# main analysis analogous to Fig. 4 E-G in NComm 2020
encode <- Sys.getenv("encode")
if (encode == "") {
  encode <- T
} else {
  encode <- as.logical(encode) #support T/F or TRUE/FALSE
}
cat("Run encoding model: ", as.character(encode), "\n")


# predicts next response based on signal and behavioral variables
rt_predict <- Sys.getenv("rt_predict")
if (rt_predict == "") {
  rt_predict <- T
} else {
  rt_predict <- as.logical(rt_predict) #support T/F or TRUE/FALSE
}
cat("Run rt prediction model: ", as.character(rt_predict), "\n")
domain = Sys.getenv("domain") # "time", "tf"
alignment <- Sys.getenv("alignment")
if (encode == "") {
  alignment <- "RT"
} else {
  stopifnot(alignment %in% c("RT", "clock", "feedback"))
  medusa_dir <- paste0("/proj/mnhallqlab/projects/Clock_MEG/tfr_rds/", alignment, "/time_freq")
  if (whoami::username() == "dombax") {
    if (domain == "tf") {medusa_dir <- paste0("/Users/Shared/tfr_rds/", alignment, "/grouped_tf")
    orig_dir <- paste0("/Users/Shared/tfr_rds/", alignment, "/time_freq")
    } else if (domain == "time") {
      medusa_dir <- paste0("/Users/Shared/tfr_rds/", alignment, "/grouped_time")
      orig_dir <- paste0("/Users/Shared/tfr_rds/", alignment, "/downsamp_20Hz_mean")
    } } else if (whoami::username() == "alexdombrovski" & domain == "time") {
    if (alignment == "clock") {medusa_dir <- "~/Box/SCEPTIC_fMRI/MEG_20Hz_n63_clockalign/"} else if (alignment == "RT") {
      medusa_dir <- "~/Box/SCEPTIC_fMRI/MEG_20Hz_n63/"
    }
    } else {stop()
    geterrmessage("Cannot find the data")}
}
cat("Alignment: ", as.character(alignment), "\n")
cat("Domain: ", as.character(domain), "\n")

stopifnot(dir.exists(medusa_dir))  
setwd(medusa_dir)

label_sensors = F
test = F
scale_winsor = F

ncores <- 36
if (whoami::username()=="dombax" | whoami::username() == "alexdombrovski") {
  ncores <- foor(detectCores()*.7)
}
# # Kai’s guidance on sensors is: ‘So for FEF, I say focus on 612/613, 543/542, 1022/1023, 
# # For IPS, 1823, 1822, 2222,2223.’
# fef_sensors <- c("0612","0613", "0542", "0543","1022")
# ips_sensors <- c("1823", "1822", "2222","2223")
# dan_sensors <- c(fef_sensors,ips_sensors)

if (domain == "time") {
  files <-  gsub("//", "/", list.files(medusa_dir,pattern = "group", full.names = T))
  # files <- files[grepl("MEG", files)]
  # files <- files[!endsWith(suffix = "1_20Hz.rds",x = files)]
  if (test) {files <- files[1]}
} else if (domain == "tf") {
  files <-  gsub("//", "/", list.files(medusa_dir,pattern = "group", full.names = T))
  #files <- files[grepl("MEG", files)]
  #build list from preferred sensors
  # files <- file.path(medusa_dir, paste0("MEG", sensor_list, "_tf.rds"))
}
filenum <- as.numeric(Sys.getenv("filenum"))
files <- files[filenum]
group_sensors = Sys.getenv("group_sensors")
# group sensors into DAN nodes
if (group_sensors) {
  sensor_map <- sensor_map %>% filter(sensor %in% sensor_list)
  if (domain=="time") {sensor_map$filename <- file.path(orig_dir, paste0("MEG", sensor_list, "_20Hz.rds"))} else if 
  (domain=="tf") {sensor_map$filename <- file.path(orig_dir, paste0("MEG", sensor_list, "_tf.rds"))}
  
  sensor_map <- sensor_map %>% mutate(node = paste(lobe, hemi, sep = "_")) 
  sensor_map <- split(sensor_map, sensor_map$node)
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  on.exit(try(stopCluster(cl)))
  foreach(d = iter(sensor_map), .packages=c("data.table", "tidyverse")) %dopar% {
    message(paste(unique(d$node)))
    l <- lapply(d$filename, readRDS)
    names(l) <- d$sensor
    group_data <- data.table::rbindlist(l, idcol = "sensor")
    setwd(medusa_dir)
    saveRDS(group_data, file = paste0(unique(d$node), "_group.RDS"))
    return(NULL)}
  files <- paste0(names(sensor_map), "_group.RDS")
}
# check that merging variables are named identically in MEG and behavior files

# mixed_by call
trial_df <- readRDS(behavioral_data_file)

#for simplicity, change all matrix columns (wi-centering and scaling) to numeric
trial_df <- as.data.frame(lapply(trial_df, function(x) {
  if (inherits(x, "matrix")) { x <- as.vector(x) }
  return(x)
}))


# trial_df$Subject <- trial_df$id
# trial_df$Run <- trial_df$run
# trial_df$Trial <- trial_df$trial
# trial_df$rt_next_sc <- scale(trial_df$rt_next)
# # # save behavioral data with new variables
# saveRDS(trial_df, behavioral_data_file)


if (alignment=="RT" | alignment=="feedback") {
  encode_formula = formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + scale(rt_vmax_lag)  + scale(rt_vmax_change) + 
                             v_entropy_wi + v_entropy_wi_change + v_max_wi  + scale(abs_pe) + outcome + (1|Subject) + (1|sensor))
} else if (alignment=="clock") {
  encode_formula = formula(~ rt_vmax + reward_lag + rt_csv_sc + rt_lag_sc + v_max_wi + trial_neg_inv_sc + 
                             v_entropy_wi + v_entropy_wi_change_lag + (1|Subject) + (1|sensor))
}

encode_formula_e = formula(~ scale(rt_vmax_lag)*echange_f1_early + scale(rt_vmax_lag)*echange_f2_late + scale(rt_vmax_lag)*e_f1 +
                             scale(abs_pe)*echange_f1_early + scale(abs_pe)*echange_f2_late + scale(abs_pe)*e_f1 +
                             outcome*echange_f1_early + outcome*echange_f2_late + outcome*e_f1 +
                             rt_csv_sc*echange_f1_early + rt_csv_sc*echange_f2_late + rt_csv_sc*e_f1 +
                             trial_neg_inv_sc*echange_f1_early + trial_neg_inv_sc*echange_f2_late + trial_neg_inv_sc*e_f1 +
                             v_entropy_wi_change*echange_f1_early + v_entropy_wi_change*echange_f2_late + v_entropy_wi*e_f1 + rt_lag_sc*e_f1 + (1|Subject) + (1|sensor))
encode_formula_pe = formula(~ scale(rt_vmax_lag)*abs_pe_f2_early + scale(rt_vmax_lag)*abs_pe_f1_mid + scale(rt_vmax_lag)*abs_pe_f3_late +
                              scale(abs_pe)*abs_pe_f2_early + scale(abs_pe)*abs_pe_f1_mid + scale(abs_pe)*abs_pe_f3_late +
                              outcome*abs_pe_f2_early + outcome*abs_pe_f1_mid + outcome*abs_pe_f3_late +   
                              rt_csv_sc*abs_pe_f2_early + rt_csv_sc*abs_pe_f1_mid + rt_csv_sc*abs_pe_f3_late +
                              trial_neg_inv_sc*abs_pe_f2_early + trial_neg_inv_sc*abs_pe_f1_mid + trial_neg_inv_sc*abs_pe_f3_late +
                              v_entropy_wi_change*abs_pe_f2_early + v_entropy_wi_change*abs_pe_f1_mid + v_entropy_wi_change*abs_pe_f3_late + 
                              v_entropy_wi*abs_pe_f2_early + v_entropy_wi*abs_pe_f1_mid + v_entropy_wi*abs_pe_f3_late + 
                              rt_lag_sc*abs_pe_f2_early + rt_lag_sc*abs_pe_f1_mid + rt_lag_sc*abs_pe_f3_late + (1|Subject) + (1|sensor))


# preproc trial-df
trial_df <- trial_df %>% mutate(Subject = as.integer(id), Trial = trial, Run = run)
if (alignment == "RT" | alignment == "feedback") {
  rt_outcome = "rt_next"} else if (alignment == "clock") {
    rt_outcome = "rt_csv_sc"
  }
if (domain == "tf") {
  splits = c("Time", ".filename", "Freq")
  #signal_outcome = "pow_scaled"
  signal_outcome = "Pow"
  #new approach: transform outcome variable at the time of computation
  trans_func <- function(x) { DescTools::Winsorize(x, probs=c(.005, 1), na.rm=TRUE) } #only drop bottom 0.5%
  rt_predict_formula = formula( ~ scale(Pow) * rt_csv_sc * outcome  + scale(Pow) * scale(rt_vmax)  +
                                  scale(Pow) * rt_lag_sc + (1|id) + (1|sensor))
  
} else if (domain == "time") {
  splits = c("Time", ".filename")
  #signal_outcome = "signal_scaled"
  signal_outcome = "signal_scaled"
  trans_func <- function(x) { DescTools::Winsorize(x, probs=c(.01, .99), na.rm=TRUE) } #drop top and bottom 1%
  rt_predict_formula = formula( ~ signal_scaled * rt_csv_sc * outcome  + signal_scaled * scale(rt_vmax)  +
                                  signal_scaled * rt_lag_sc + (1|id) + (1|sensor))
  
}
gc()
if (encode) {
  ddf <- as_tibble(mixed_by(files, outcomes = signal_outcome, rhs_model_formulae = encode_formula_e, split_on = splits,
                            external_df = trial_df, external_merge_by=c("Subject", "Run", "Trial"), padjust_by = "term", padjust_method = "BY", ncores = ncores,
                            refit_on_nonconvergence = 5, outcome_transform=trans_func))
  # refit_on_nonconvergence = 5))
  # ddf$sensor <- readr::parse_number(ddf$.filename)
  # save output
  #setwd("~/OneDrive/collected_letters/papers/meg/plots/clock_encode/")
  if (domain == "time") {
    # saveRDS(ddf, file = "meg_mixed_by_time_clock_ddf.RDS")
    # saveRDS(ddf, file = "meg_mixed_by_time_ranefs_mult_interactions_pe_ddf.RDS")
    saveRDS(ddf, file = paste0("meg_mixed_by_time_ddf_combined_ranefs", alignment, ".RDS"))    
  } else if (domain == "tf") {
    # saveRDS(ddf, file = "meg_mixed_by_tf_ranefs_mult_interactions_e_ddf.RDS")
    # saveRDS(ddf, file = "meg_mixed_by_tf_clock_ddf.RDS")
    saveRDS(ddf, file = paste0("meg_mixed_by_tf_ddf_combined_ranefs", alignment,files, ".RDS"))    
  }
}  
gc()
if (rt_predict) {
  rdf <- as_tibble(mixed_by(files, outcomes = rt_outcome, rhs_model_formulae = rt_predict_formula , split_on = splits, external_df = trial_df,
                            padjust_by = "term", padjust_method = "BY", ncores = ncores, refit_on_nonconvergence = 5, outcome_transform=trans_func))
  # rdf$sensor <- readr::parse_number(rdf$.filename)
  
  # save output
  #setwd("~/OneDrive/collected_letters/papers/meg/plots/rt_rt/")
  if (domain == "time") {
    saveRDS(rdf, file = paste0("meg_mixed_by_time_rdf_combined_", alignment, ".RDS"))
  } else if (domain == "tf") {
    saveRDS(rdf, file = paste0("meg_mixed_by_tf_rdf_combined_", alignment, filenum, ".RDS"))
  }
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
