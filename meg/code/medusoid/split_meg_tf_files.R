# split MEG TF files by timepoint
library(tidyverse)
library(parallel)
library(foreach)
library(doParallel)

# loop through files sequentially
repo_directory <- "~/code/clock_analysis"
for (alignment in c("RT", "clock")) {
  medusa_dir <- paste0("/Users/Shared/tfr_rds/", alignment, "/grouped_tf")
  setwd(medusa_dir)
  files <- list.files(pattern = "group.RDS")
  for (file in files) {
    df <- readRDS(file)
    df <- df %>% group_by(Freq)
    freqs <- unique(df$Freq)
    l <- group_split(df)
    names(l) <- freqs
    for (freq in freqs) {
      saveRDS(file = paste0(str_remove(file, ".RDS"), "_freq_split_", freq), data.table::rbindlist(l[freq]))
    }
  }
  
}
# load file
# split MEG TF files by timepoint
library(tidyverse)
# loop through files sequentially
for (alignment in c("RT", "clock")) {
  medusa_dir <- paste0("/Users/Shared/tfr_rds/", alignment, "/grouped_tf")
  setwd(medusa_dir)
  files <- list.files(pattern = "freq_split")
  ncores <- detectCores()/2
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  on.exit(try(stopCluster(cl)))
  foreach(filename %in% files, .packages=c("data.table", "tidyverse")) %dopar% {
    message(paste(filename))
    l <- readRDS(filename)
    saveRDS(file = filename, data.table::rbindlist(l))
    return(NULL)}
  # for (filename in files) {
  #   l <- readRDS(filename)
  #   saveRDS(file = filename, data.table::rbindlist(l))
  # }
}


