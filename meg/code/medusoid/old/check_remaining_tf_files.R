# inspect mixed_by results from whole-brain analyses of MEG
library(tidyverse)

setwd("~/OneDrive/collected_letters/papers/meg/plots/wholebrain")
ddf <- readRDS("meg_mixed_by_tf_ddf_wholebrain_entropy_change_rs_RTfreq_t_all_sensors80.000_1.992.rds")
ecdf <- ddf %>% filter(term == "v_entropy_wi_change" & group == "Sensor" & effect == "ran_coefs")
ggplot(ecdf %>% filter(effect == "ran_vals"), aes(level, estimate)) + geom_point() + geom_errorbar(aes(ymin = conf.low, ymax = conf.high))
ggplot(ecdf %>% filter(effect == "ran_coefs"), aes(level, estimate)) + geom_point() + geom_errorbar(aes(ymin = conf.low, ymax = conf.high))
fdf <- ddf %>% filter(effect=="fixed")


# deal with unfinished TF points
# run meg_tf_dual_wholebrain.R first
setwd(data_dir)
all_files <- readRDS("all_tf_files_rt.rds")
finished_files <- unique(rddf$.filename)
remaining_files <- setdiff(all_files, finished_files)

# what the hell is in those "finish" files?

file_pattern <- "meg_mixed_by_tf_ddf_wholebrain_entropy_change_rs_finishRT"
files <-  gsub("//", "/", list.files(data_dir, pattern = file_pattern, full.names = F))
l <- lapply(files, readRDS)
df <- data.table::rbindlist(l)
last_files <- unique(df$.filename)

ddf1 <- readRDS("meg_ddf_wholebrain_ec_rs_1.rds")
# ddf2 <- readRDS("meg_ddf_wholebrain_ec_rs_2.rds")
ddf3 <- readRDS("meg_ddf_wholebrain_ec_rs_3.rds")
# ddf4 <- readRDS("meg_ddf_wholebrain_ec_rs_4.rds")
# ddf5 <- readRDS("meg_ddf_wholebrain_ec_rs_5.rds")

ddf1 <- ddf1 %>% filter(alignment == "rt")
# ddf2 <- ddf2 %>% filter(alignment == "rt")
ddf3 <- ddf3 %>% filter(alignment == "rt")
# ddf4 <- ddf4 %>% filter(alignment == "rt")
# ddf5 <- ddf5 %>% filter(alignment == "rt")

files1 <- unique(ddf1$.filename)
# files2 <- unique(ddf2$.filename)
files3 <- unique(ddf3$.filename)
# files4 <- unique(ddf4$.filename)
# files5 <- unique(ddf5$.filename)

# merge two unique output dataframes
mddf <- unique(rbind(ddf1, ddf3))
# write almost-complete file
saveRDS(mddf, "meg_rddf_wholebrain_ec_rs_almost_complete.rds")

edf <- mddf %>% filter(effect=="fixed" & term == "(Intercept)")
termstr <- "intercept"
message(termstr)
fname = paste("meg_tf_combined_uncorrected_", termstr, ".pdf", sep = "")
pdf(fname, width = 10, height = 7.5)
print(ggplot(edf, aes(Time - 5, Freq)) + geom_tile(aes(fill = estimate), size = .01) +
        geom_vline(xintercept = 0, lty = "dashed", color = "black", size = 2) +
        geom_vline(xintercept = -5, lty = "dashed", color = "white", size = 2) +
        geom_vline(xintercept = -5.3, lty = "dashed", color = "white", size = 1) +
        geom_vline(xintercept = -2.5, lty = "dotted", color = "grey", size = 1) +
        scale_fill_viridis(option = "plasma") +  xlab(rt_epoch_label) + ylab("Frequency") +
        # facet_wrap( ~ node, ncol = 2) +
        geom_text(data = edf, x = -5.5, y = 5,aes(label = "Response(t)"), size = 2.5, color = "white", angle = 90) +
        geom_text(data = edf, x = -4.5, y = 5,aes(label = "Outcome(t)"), size = 2.5, color = "white", angle = 90) +
        geom_text(data = edf, x = 0.5, y = 6 ,aes(label = "Clock onset (t+1)"), size = 2.5, color = "black", angle = 90) +
        labs(alpha = expression(italic(p)[uncorrected])) + ggtitle(paste(termstr)) + theme_dark())
dev.off()
