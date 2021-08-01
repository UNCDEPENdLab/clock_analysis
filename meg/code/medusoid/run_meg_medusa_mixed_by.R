library(varhandle)
rm(list=ls())

# more cores can be used on smaller frontal (1,2) and temporal (5,6) files. 
for (i in c(1,2,5,6)) {
  rm.all.but(keep = "i")
  Sys.setenv(filenum = i)
  gc()
  # call, options
  Sys.setenv(encode = T, rt_predict = T, alignment = c("RT"), domain = "tf", group_sensors = F, core_prop = 0.8)
  source("~/code/clock_analysis/meg/code/medusoid/meg_medusa_mixed_by_combined.R")
  rm.all.but(keep = "i")
  gc()
  Sys.setenv(encode = T, rt_predict = T, alignment = c("clock"), domain = "tf", group_sensors = F, core_prop = 0.8)
  source("~/code/clock_analysis/meg/code/medusoid/meg_medusa_mixed_by_combined.R")
}

for (i in c(3:4)) {
  rm.all.but(keep = "i")
  Sys.setenv(filenum = i)
  gc()
  # call, options
  Sys.setenv(encode = T, rt_predict = T, alignment = c("RT"), domain = "tf", group_sensors = F, core_prop = 0.4)
  source("~/code/clock_analysis/meg/code/medusoid/meg_medusa_mixed_by_combined.R")
  rm.all.but(keep = "i")
  gc()
  Sys.setenv(encode = T, rt_predict = T, alignment = c("clock"), domain = "tf", group_sensors = F, core_prop = 0.4)
  source("~/code/clock_analysis/meg/code/medusoid/meg_medusa_mixed_by_combined.R")
}
