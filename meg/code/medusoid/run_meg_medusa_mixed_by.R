# not functional for TF yet

# run to avoid feeding extra data to workers
rm(list=ls())
# crashes on 4, Right Parietal
# rerun starting after it
# for (i in 1:6) {
for (i in 4) {
  Sys.setenv(filenum = i)
  gc()
  # call, options
  Sys.setenv(encode = F, rt_predict = T, alignment = c("RT"), domain = "tf", group_sensors = F)
  source("~/code/clock_analysis/meg/code/medusoid/meg_medusa_mixed_by_combined.R")
  
  gc()
  Sys.setenv(encode = F, rt_predict = F, alignment = c("clock"), domain = "tf", group_sensors = F)
  source("~/code/clock_analysis/meg/code/medusoid/meg_medusa_mixed_by_combined.R")
}
