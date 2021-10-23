# test merge of within-trial UV with MEG data
library(foreach)
library(doParallel)
library(tidyverse)
medusa_dir <- "/bgfs/adombrovski/tfr_rds1/clock"
code_dir <- "~/code/clock_analysis/meg/data/"
merge_timepoints = F

uv_censored_data_file <- "~/code/clock_analysis/meg/data/meg_n63_uv_behav_ends_censored_Jul2021.Rds"
uv_uncensored_data_file <- "~/code/clock_analysis/meg/data/meg_n63_uv_behav_uncensored_Jul2021.Rds"
uv <- readRDS(uv_uncensored_data_file) 
uv <- uv %>% as.data.table(lapply(uv, function(x) {
  if (inherits(x, "matrix")) { x <- as.vector(x) }
  return(x)
}))

# sample file for testing
# meg <- readRDS(file.path(medusa_dir, "freq_t_all_sensors02.102_-0.006.rds"))

# merge files across all timepoints
all_files <- grep(list.files(path=medusa_dir, pattern = "freq_t_all_sensors", full.names = F), pattern='ddf|-', invert=TRUE, value=TRUE)
# freqs <- as.character(unique(parse_number(all_files)))
freqs <- unique(substr(all_files, 19, 24))
# the idea is to loop over all frequencies, read files into a list and write them
if (merge_timepoints) {
  # ncores <- detectCores()/8
  # ncores <- 2
  # cl <- makeCluster(ncores)
  # registerDoParallel(cl)
  # on.exit(try(stopCluster(cl)))
  setwd(medusa_dir)
  # 10 Hz file is double the size, need to recalculate
  # foreach(f = c(10,18), .packages=c("data.table", "tidyverse")) %dopar% {
  for (f in 1:length(freqs)) {
    gc()
    freq <- freqs[f]
    message(freq)
    freq_files <- grep(all_files, pattern = freq, value = T)
    message("Merging files")
    message(freq_files)
    # freq_files <- freq_files[1:8]
    fl <- lapply(freq_files, readRDS)
    message("Merged, binding")
    fdf <- data.table::rbindlist(fl, use.names = T)
    timestep_breaks <- c(0, unique(uv$timestep)) # these are really ends of time bins
    labels <- timestep_breaks[1:length(timestep_breaks)-1]
    fdf$timestep <- as.numeric(as.character(cut(fdf$Time, breaks = timestep_breaks, labels = labels)))
    message("Saving")
    saveRDS(fdf, file = paste0("freq_split_", freq, ".Rds"))
    # return(NULL)}
  }
}

# save subset of 40Hz data as a prototype
sdf <- fdf %>% filter(Sensor == "1822" | Sensor == "1823")
sdf <- sdf %>% mutate(timestep =  as.numeric(as.character(timestep)))
# get U and V leads
uv <- uv %>% arrange(id, trial, timestep) %>% group_by(id, trial) %>%
  mutate(
    v_lead1_wi = lead(value_wi, 1),
    u_lead1_wi = lead (uncertainty_wi, 1),
    v_lead2_wi = lead(value_wi, 2),
    u_lead2_wi = lead (uncertainty_wi, 2),
    v_lead4_wi = lead(value_wi, 4),
    u_lead4_wi = lead (uncertainty_wi, 4),
    v_lead1_wi_t = lead(value_wi_t, 1),
    u_lead1_wi_t = lead (uncertainty_wi_t, 1),
    v_lead2_wi_t = lead(value_wi_t, 2),
    u_lead2_wi_t = lead (uncertainty_wi_t, 2),
  ) %>% ungroup()
setwd(code_dir)
saveRDS(uv, file = "meg_n63_uv_behav_uncensored_leads_Jul2021.Rds")
uvc <- uv %>% filter(timestep > 0.4 & timestep < 3.6)
saveRDS(uv, file = "meg_n63_uv_behav_censored_leads_Jul2021.Rds")

# prototype a mixed_by on just 2 sensors

signal_outcome = "Pow"
trans_func <- function(x) { DescTools::Winsorize(x, probs=c(.005, 1), na.rm=TRUE) }
source("~/code/fmri.pipeline/R/mixed_by.R")
library(splines)

encode_formula_uv =  formula(~ bs(timestep, df = 4) + value_wi_t + uncertainty_wi_t + value_b_t + uncertainty_b_t + (1|Subject) + (value_wi_t + uncertainty_wi_t|Sensor))
for (f in 1:length(freqs)) {
  gc()
  freq <- freqs[f]
  message(freq)
  all_fsplit_files <- grep("group", list.files(pattern = "freq_split"), value = T, invert = T)
  freq_files <- grep(all_fsplit_files, pattern = freq, value = T)
  message("Processing file")
  message(freq_files)
  fdf <- readRDS(freq_files)
  ddf <- mixed_by(fdf, outcomes = signal_outcome, rhs_model_formulae = encode_formula_uv,
                  external_df = uvc, external_merge_by=c("Subject", "Trial", "timestep"), 
                  padjust_by = "term", padjust_method = "BY", ncores = 1, calculate=c("parameter_estimates_reml"),
                  refit_on_nonconvergence = 5, outcome_transform=trans_func, tidy_args=list(effects=c("fixed", "ran_vals", "ran_pars", "ran_coefs"), conf.int=TRUE))
  ddf$coef_df_reml[ddf$coef_df_reml$effect=="fixed"]
  
  saveRDS(ddf, file = paste0("meg_mixed_by_tf_ddf_wholebrain_uv_basic_", freq ,"_.Rds"))
}
encode_formula_uv =  formula(~ bs(timestep, df = 4) + value_wi_t + uncertainty_wi_t + value_b_t + uncertainty_b_t + (value_wi_t + uncertainty_wi_t + value_b_t + uncertainty_b_t |Subject) + (value_wi_t + uncertainty_wi_t + value_b_t + uncertainty_b_t|Sensor))
for (f in 1:length(freqs)) {
  gc()
  freq <- freqs[f]
  message(freq)
  all_fsplit_files <- grep("group", list.files(pattern = "freq_split"), value = T, invert = T)
  freq_files <- grep(all_fsplit_files, pattern = freq, value = T)
  message("Processing file")
  message(freq_files)
  fdf <- readRDS(freq_files)
  ddf <- mixed_by(fdf, outcomes = signal_outcome, rhs_model_formulae = encode_formula_uv,
                  external_df = uvc, external_merge_by=c("Subject", "Trial", "timestep"), 
                  padjust_by = "term", padjust_method = "BY", ncores = 1, calculate=c("parameter_estimates_reml"),
                  refit_on_nonconvergence = 5, outcome_transform=trans_func, tidy_args=list(effects=c("fixed", "ran_vals", "ran_pars", "ran_coefs"), conf.int=TRUE))
  ddf$coef_df_reml[ddf$coef_df_reml$effect=="fixed"]
  
  saveRDS(ddf, file = paste0("meg_mixed_by_tf_ddf_wholebrain_uv_rs_", freq ,"_.Rds"))
}