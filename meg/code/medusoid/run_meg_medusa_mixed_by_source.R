rm(list=ls())
for (i in 3:4) {
  Sys.setenv(filenum = i)
  gc()
  Sys.setenv(encode = T, rt_predict = T, alignment = c("RT"), domain = "tf", group_rois = F, bin = F)
  source("~/code/clock_analysis/meg/code/medusoid/meg_medusa_mixed_by_source.R")
  gc()
  Sys.setenv(encode = T, rt_predict = T, alignment = c("clock"), domain = "tf", group_rois = F, bin = F)
  source("~/code/clock_analysis/meg/code/medusoid/meg_medusa_mixed_by_source.R")
}
