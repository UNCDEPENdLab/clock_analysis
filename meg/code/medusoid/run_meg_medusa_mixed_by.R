# not functional for TF yet

# run to avoid feeding extra data to workers
for (i in 1:6) {
  Sys.setenv(filenum = i)
  rm(list=ls())
  gc()
  # call, options
  Sys.setenv(encode = T, rt_predict = T, alignment = c("RT"), domain = "tf", group_sensors = F)
  source("~/code/clock_analysis/meg/code/medusoid/meg_medusa_mixed_by_combined.R")
  
  gc()
  Sys.setenv(encode = T, rt_predict = T, alignment = c("clock"), domain = "tf", group_sensors = Fi)
  source("~/code/clock_analysis/meg/code/medusoid/meg_medusa_mixed_by_combined.R")
}
