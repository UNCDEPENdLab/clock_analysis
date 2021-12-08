# filter MEG TF files by time point
# remove < -.4s and 
library(tidyverse)
library(stringr)
# setwd("/proj/mnhallqlab/projects/Clock_MEG/atfr_rds/RT")
# setwd("~/OneDrive/collected_letters/papers/meg/plots/wholebrain/output/test")
files <- list.files(pattern = "freq_t")
list <- strex::str_extract_numbers(files, negs = T, decimals = T)
times <- do.call(rbind,list)[,2]
freqs <- do.call(rbind,list)[,1]
bad_files <- files[times < -0.4 | freqs > 40 | times > 1.4]
target_bad_files <- paste0("bad/", bad_files)
file.rename(from = bad_files, to = target_bad_files)
