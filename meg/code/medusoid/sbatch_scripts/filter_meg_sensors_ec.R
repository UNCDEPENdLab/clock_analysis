# filter only sensors responsive to entropy change for the purpose of fMRI ->  MEG analyses

library(tidyverse)
library(foreach)
#library(doParallel)
 library(doSNOW)
library(parallel)
library(parallelly)

ec_sensors <- c(1642, 1842, 1912, 1913, 1923, 1922,
                2013, 2012, 2022, 2032, # from Kai
                2033, 2042, 2043, 2113, 2112, 2313, 2343)
ncores <- detectCores()/2
cl <- makeCluster(ncores)
registerDoSNOW(cl)
# registerDoParallel(cl)
on.exit(try(stopCluster(cl)))

# grab all the files
basedir <- "/bgfs/adombrovski/tfr_rds1/RT"
flist <- list.files(pattern = "^freq_t_all_sensors", path = basedir)
flist <- flist[1:3]
# flist <- list(flist)
setwd(basedir)
temp <- foreach(i = 1:(length(flist)), .combine='+', .packages=c("tidyverse")) %dopar% {
  # .GlobalEnv$f <- f
  # for (f in flist) {
  # f <- flist[i]
  # setwd(basedir)
  # df <- readRDS(f) %>% filter(Sensor %in% ec_sensors)
  # saveRDS(f, file=file.path(paste0("freq_t_ec_sensors_", unique(df$Freq), unique(df$Time), ".rds")))
  # return(NULL)
  }
registerDoSEQ()

# unregister_dopar <- function() {
#   env <- foreach:::.foreachGlobals
#   rm(list=ls(name=env), pos=env)
# }
# unregister_dopar()
