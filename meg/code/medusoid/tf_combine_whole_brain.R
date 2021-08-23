library(tidyverse)
library(foreach)
library(doParallel)
library(data.table)
library(psych)
library(stringr)

alignment = "RT"
orig_dir <- paste0("/bgfs/adombrovski/tfr_rds/", alignment, "/time_freq")
medusa_dir <- paste0("/bgfs/adombrovski/tfr_rds1/", alignment)
setwd(orig_dir)
files <- list.files(pattern = "MEG")
sensors <- unique(substr(files,4,7))
sensors <- sensors[105:length(sensors)]
for (sensor in sensors) {
  message(sensor)
  setwd(orig_dir)
  sensorfile <- files[grepl(sensor, files)]
  sdf <- readRDS(sensorfile)
  sdf$Sensor <- sensor
  setwd(medusa_dir)
  message(paste0("Read, saving ", sensor))
  sdf[, write.table(.SD, paste("freqTsplit", Freq, "_", Time, "_", sensor,  ".rds", sep = "")), by = c("Freq", "Time")]
}
setwd(medusa_dir)
lobefiles <- list.files(pattern = "freq_split")
freqs <- unique(str_remove_all(str_extract(lobefiles, "_[^_]+$"), pattern = c(".rds|_")))
times <- unique(sdf$Time)
files <- list.files(pattern = "freqTsplit")
for (freq in freqs) {
  message(freq)
  for (t in times) {
    message(t)
    t_files <- files[grepl(t, files)]
    tf_files <- t_files[grepl(freq, t_files)]
    tf <- do.call(rbind, lapply(tf_files, read.table, fill = T))  
    tf$Freq <- freq
    tf$Time <- t
    saveRDS(tf, paste("freq_t_all_sensors", freq, "_", t, ".rds", sep = ""))
  }
}
