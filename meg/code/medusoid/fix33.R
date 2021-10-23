# the idea is to go through freqTsplit by sensor files, merge by frequency
library(tidyverse)
for (epoch in c("clock", "RT")) {
  orig_dir <- paste0("/bgfs/adombrovski/tfr_rds/", epoch, "/time_freq")  
  target_dir <- paste0("/bgfs/adombrovski/tfr_rds1/", epoch, "/temp33")
  
  setwd(orig_dir)
  all_files <- list.files(pattern = "MEG")
  #freqs <- as.character(unique(parse_number(all_files)))
  #list33 <- list.files(pattern = '33_')
  #sensors <- unique(as.character(gsub(".*_", "", list33)))
  #sensors <- gsub(".rds", "", sensors)
  for (file in all_files) {
    setwd(orig_dir)
    df <- readRDS(file)
    df <- df %>% filter(Time == unique(df$Time)[49]) # I don't understand why it can't be 0.33
    freqs <- as.character(unique(df$Freq))
    sensor <- substr(file, 4, 7)
    message(sensor)
    for (f in 1:length(freqs)) {
      message(freqs[f])
      fdf <- df %>% filter(Freq == unique(df$Freq)[f])
      fdf$Sensor <- sensor
      setwd(target_dir)
      saveRDS(fdf, file = paste0("freqTsplit_", as.character(unique(fdf$Freq)), "_", as.character(unique(fdf$Time)), "_", sensor, ".Rds"))
    }
  }
}

for (epoch in c("clock", "RT")) {
  target_dir <- paste0("/bgfs/adombrovski/tfr_rds1/", epoch, "/temp33")
  setwd(target_dir)
  for (freq in freqs) {
    nfreq <- substr(freq, 3, 8)
    message(paste0(epoch, "_", nfreq))
            freq_files <- list.files(pattern = nfreq)
            fls <- lapply(freq_files, readRDS)
            ftdf <- data.table::rbindlist(fls)
            saveRDS(ftdf, file = paste0("freq_t_all_sensors", nfreq, "_", "0.33.rds" ))
  }
}