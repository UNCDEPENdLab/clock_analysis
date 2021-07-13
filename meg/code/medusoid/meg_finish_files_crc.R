library(tidyverse)

setwd("/bgfs/adombrovski/tfr_rds1/RT/results")
# ddf <- readRDS("meg_mixed_by_tf_ddf_wholebrain_entropy_change_rs_RTfreq_t_all_sensors80.000_1.992.rds")
# ecdf <- ddf %>% filter(term == "v_entropy_wi_change" & group == "Sensor")
# ggplot(ecdf, aes(level, estimate)) + geom_point() + geom_errorbar(aes(ymin = conf.low, ymax = conf.high))


# finish up whole-brain analyses that timed out

temp <- readRDS("meg_ddf_wholebrain_ec_rs.rds")
finished_files <- unique(temp$.filename)
setwd("/bgfs/adombrovski/tfr_rds1/RT")
all_files <- as.vector(list.files(pattern = "freq_t_all"))
remaining_files <- setdiff(all_files, finished_files)
