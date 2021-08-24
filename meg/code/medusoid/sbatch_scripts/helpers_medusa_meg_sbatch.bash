# count jobs


#!/bin/bash
squeue -u ayd1  -h -t pending -r | wc -l
squeue -u ayd1  -h -t running -r | wc -l

squeue --start -u ayd1

# split MEG TF files by timepoint
library(tidyverse)
# loop through files sequentially
for (alignment in c("RT", "clock")) {
  medusa_dir <- paste0("/bgfs/adombrovski/tfr_rds1/", alignment)
  setwd(medusa_dir)
  files <- list.files(pattern = "freq_split")
  for (filename in files) {
    l <- readRDS(filename)
    saveRDS(file = filename, data.table::rbindlist(l))
    }
  }
  

