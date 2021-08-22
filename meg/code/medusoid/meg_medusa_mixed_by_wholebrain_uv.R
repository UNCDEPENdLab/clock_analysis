library(tidyverse)

uv_censored_data_file <- "~/code/clock_analysis/meg/data/meg_n63_uv_behav_censored_leads_Jul2021.Rds"
uv_uncensored_data_file <- "~/code/clock_analysis/meg/data/meg_n63_uv_behav_uncensored_leads_Jul2021.Rds"
source("~/code/fmri.pipeline/R/mixed_by.R")
window = censored
debug = F
alignment = "clock"
if (whoami::username()=="ayd1") {
  medusa_dir <- paste0("/bgfs/adombrovski/tfr_rds1/", alignment)
} else if (whoami::username()=="dnpl") {
  medusa_dir <- paste0("/proj/mnhallqlab/projects/Clock_MEG/atfr_rds/", alignment)
}

stopifnot(dir.exists(medusa_dir))  
setwd(medusa_dir)
message("Directory ", medusa_dir)

ncores <- as.numeric(future::availableCores())
sourcefilestart <- as.numeric(Sys.getenv("sourcefilestart"))
if (debug) {
  sourcefilestart = 1
}
setwd(medusa_dir)
all_files <- list.files(pattern = "freq_s", full.names = T)
files <- all_files[sourcefilestart]
message(paste0("Processing files "))
cat(files)
if (window == "censored") {uv_df <- readRDS(uv_censored_data_file)} else if (window == "uncensored") {uv_df <- readRDS(uv_uncensored_data_file)}
encode_formula_uv =  formula(~ bs(timestep, df = 4) + value_wi_t + uncertainty_wi_t + value_b_t + uncertainty_b_t + (value_wi_t + uncertainty_wi_t|Subject) + (value_wi_t + uncertainty_wi_t|Sensor))
signal_outcome = "Pow"
trans_func <- function(x) { DescTools::Winsorize(x, probs=c(.005, 1), na.rm=TRUE) }
gc()

ddf <- as_tibble(mixed_by(files, outcomes = signal_outcome, rhs_model_formulae = encode_formula_uv,
                          external_df = trial_df, external_merge_by=c("Subject", "Run", "Trial"), padjust_by = "term", padjust_method = "BY", ncores = ncores,
                          refit_on_nonconvergence = 5, outcome_transform=trans_func, tidy_args=list(effects=c("fixed", "ran_vals", "ran_pars", "ran_coefs"), conf.int=TRUE)))
saveRDS(ddf, file = paste0("meg_mixed_by_uv_", window, sourcefilestart))
