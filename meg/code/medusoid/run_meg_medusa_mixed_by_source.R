  gc()
  Sys.setenv(encode = F, rt_predict = F, alignment = c("RT"), domain = "time", group_rois = F, bin = T)
  source("~/code/clock_analysis/meg/code/medusoid/meg_medusa_mixed_by_source.R")
  
  gc()
  Sys.setenv(encode = T, rt_predict = T, alignment = c("clock"), domain = "time", group_rois = T, bin = T)
  source("~/code/clock_analysis/meg/code/medusoid/meg_medusa_mixed_by_source.R")

