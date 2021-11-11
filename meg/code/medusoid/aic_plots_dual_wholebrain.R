# AIC comparison between full and selective maintenance models
# uses MEG TF data

library(tidyverse)
library(lme4)
library(ggpubr)
library(car)
library(viridis)
library(ggnewscale)
library(RColorBrewer)
source("~/code/Rhelpers/theme_black.R")

repo_directory <- "~/code/clock_analysis"
data_dir <- "~/OneDrive/collected_letters/papers/meg/plots/wholebrain/output"
plot_dir <- "~/OneDrive/collected_letters/papers/meg/plots/wholebrain/"

rt_epoch_label = "Time relative to event, seconds"
encode = T
rt_predict = F
# regressors = c("reward")
p_adjust_method = "fdr"
regressor = c("entropy_change")
# regressors = c("entropy", "kld", "entropy_change", "entropy_change_neg", "entropy_change_pos", "reward")
print_filenames = T
fixed_only = F
reprocess = T
plots = T
aic_threshold = 100
diags = F
noclock = T
freq_threshold = 18 # set to 40 Hz (18th band) for full-spectrum output
setwd(data_dir)
# plots ----

message(paste0("Processing ", regressor))
epoch_label = "Time relative to feedback, seconds"
if (reprocess) {
  if (!noclock) {
    epoch_label = "Time relative to clock onset, seconds"
    # get clock-aligned data
    if (regressor=="entropy_change") {
      sel_file_pattern <- "*_change_sel_.*clock"} 
    files <-  gsub("//", "/", list.files(data_dir, pattern = sel_file_pattern, full.names = F))
    message(paste0("Found ", length(files), " files."))
    csl <- lapply(files, function(x) {
      if (print_filenames) { print(x) }
      df <- readRDS(x) 
      # if (class(df)=="list") {
      df <- df$fit_df
      # } else if (ncol(df)==3) {
      # df <- df$fit_df
      # }
      return(df)
    })
    csdf <- data.table::rbindlist(csl) %>% unique() %>% distinct(Time, Freq, .keep_all = TRUE)
    csdf <- csdf %>% mutate(t  = as.numeric(Time), alignment = "clock", model = "selective",
    )
    message("Processed clock-aligned selective. \n")
    
    if (regressor=="entropy_change") {
      full_file_pattern <- "*_change_full_.*clock"} 
    files <-  gsub("//", "/", list.files(data_dir, pattern = full_file_pattern, full.names = F))
    message(paste0("Found ", length(files), " files."))
    cfl <- lapply(files, function(x) {
      if (print_filenames) { print(x) }
      df <- readRDS(x) 
      # if (ncol(df)==3) {
      df <- df$fit_df
      # }
      return(df)
    })
    cfdf <- data.table::rbindlist(cfl) %>% unique() %>% distinct(Time, Freq, .keep_all = TRUE)
    cfdf <- cfdf %>% mutate(t  = as.numeric(Time), alignment = "clock", model = "full",
    )
    message("Processed clock-aligned full. Merging. \n")
    cddf <- rbind(csdf, cfdf)
  }
  # get RT-aligned
  if (regressor=="entropy_change") {
    sel_file_pattern <- "*_change_sel_.*RT"} 
  files <-  gsub("//", "/", list.files(data_dir, pattern = sel_file_pattern, full.names = F))
  message(paste0("Found ", length(files), " files."))
  rsl <- lapply(files, function(x) {
    if (print_filenames) { print(x) }
    df <- readRDS(x) 
    # if (ncol(df)==3) {
    df <- df$fit_df
    # }
    #        df <- df %>% filter(effect=="fixed")
    return(df)
  })
  rsdf <- data.table::rbindlist(rsl) %>% unique() %>% distinct(Time, Freq, .keep_all = TRUE)
  rsdf <- rsdf %>% mutate(t  = as.numeric(Time), alignment = "rt", model = "selective",
  )
  message("Processed clock-aligned selective. \n")
  if (regressor=="entropy_change") {
    full_file_pattern <- "*_change_full_.*RT"} 
  files <-  gsub("//", "/", list.files(data_dir, pattern = full_file_pattern, full.names = F))
  message(paste0("Found ", length(files), " files."))
  rfl <- lapply(files, function(x) {
    if (print_filenames) { print(x) }
    df <- readRDS(x) 
    # if (ncol(df)==3) {
    df <- df$fit_df
    # }
    #        df <- df %>% filter(effect=="fixed")
    return(df)
  })
  rfdf <- data.table::rbindlist(rfl) %>% unique() %>% distinct(Time, Freq, .keep_all = TRUE)
  rfdf <- rfdf %>% mutate(t  = as.numeric(Time), alignment = "rt", model = "full",
  )
  message("Processed RT-aligned full. \n")
  rddf <- rbind(rfdf, rsdf)
  
  rddf <- rddf %>% mutate(t  = Time - offset
  )
  if (!noclock) {offset = 4.3
  ddf <- rbind(rddf, cddf)
  } else {offset = 0.3}
  ddf <- rddf
  ddf$Freq <- gsub("f_", "", ddf$Freq)
  ddf$Freq <- ordered(as.numeric(substr(as.character(ddf$Freq), 1,4)))    
  ddf <- ddf %>% arrange(Time, Freq, alignment, model) %>% mutate(row = row_number()) 
  if (!noclock) {ddf <- rbind(cddf, rddf)} else {
    ddf <- rddf  
    wddf <- pivot_wider(ddf %>% select(Time, t, Freq, AIC, BIC, model, alignment), names_from = model, names_sep = "_", values_from = c(AIC, BIC)) %>%
      mutate(AIC_diff = AIC_full - AIC_selective,
             BIC_diff = BIC_full - BIC_selective)
  }
  # convert to wide
  wddf <- pivot_wider(ddf %>% select(Time, t, Freq, AIC, BIC, model, alignment), names_from = model, names_sep = "_", values_from = c(AIC, BIC)) %>%
    mutate(AIC_diff = AIC_full - AIC_selective,
           BIC_diff = BIC_full - BIC_selective)
  saveRDS(ddf, file = paste0("meg_ddf_wholebrain_compression_", regressor, ".rds"))
  saveRDS(wddf, file = paste0("meg_wddf_wholebrain_compression_", regressor, ".rds"))
} else {
  setwd(paste0(plot_dir, "/", "compression"))
  wddf <- readRDS(paste0("meg_wddf_wholebrain_compression_", regressor, ".rds"))
}

if (plots) {
  message("Processed output, plotting.  \n")
  # p1 <- ggplot(wddf %>% filter(alignment=="clock" & abs(AIC_diff) > aic_threshold), aes(Time, Freq)) + geom_tile(aes(fill = AIC_diff)) + 
  #   scale_fill_viridis(option = "plasma") + theme_dark()
  # p2 <- ggplot(wddf %>% filter(alignment=="rt" & abs(AIC_diff) > aic_threshold), aes(Time, Freq)) + geom_tile(aes(fill = AIC_diff)) + 
  #   scale_fill_viridis(option = "plasma") + theme_dark()
  # ggarrange(p2, p1)
  fname = paste0("meg_tf_compression_", regressor, ".pdf")
  wddf <- wddf %>% mutate(AIC_diff_thresholded  = case_when(
    abs(AIC_diff) > aic_threshold ~ AIC_diff,
    TRUE ~ NA_real_
  )) %>% filter(as.numeric(Freq)<freq_threshold)
  if (!noclock){
  pdf(fname, width = 10, height = 5)
  ggplot(wddf , aes(t, Freq)) + geom_tile(aes(fill = AIC_diff_thresholded)) + 
    geom_vline(xintercept = 0, lty = "dashed", color = "black", size = 2) +
    geom_vline(xintercept = -offset + 0.3, lty = "dashed", color = "white", size = 2) +
    geom_vline(xintercept = -offset, lty = "dashed", color = "white", size = 1) +
    geom_vline(xintercept = -2, lty = "dotted", color = "grey", size = 1) +
    scale_fill_viridis(option = "plasma") +  xlab(rt_epoch_label) + ylab("Frequency") +
    geom_text(data = wddf, x = -offset-.1, y = 5,aes(label = "Response(t)"), size = 3, color = "white", angle = 90) +
    geom_text(data = wddf, x = -offset+.4, y = 5,aes(label = "Outcome(t)"), size = 3, color = "white", angle = 90) +
    geom_text(data = wddf, x = 0.15, y = 6 ,aes(label = "Clock onset (t+1)"), size = 3, color = "black", angle = 90) +
    scale_x_continuous(limits = c(-5.3,1.7), breaks = c(0-offset+0.3, 0.20-offset+0.3, 0.4-offset+0.3, 0.60-offset+0.3, 1-offset+0.3, 1.5-offset+0.3, 0, 1), labels = c("0", "0.2", "0.4", "0.6", "1", "1.5",  "0", "1")) +
    labs(fill = expression("AIC"["full"] - AIC["selective"])) + ggtitle("Evidence of value compression") + theme_dark()
  dev.off()
  } else { # RT only
    pdf(fname, width = 7, height = 4.5)
    ggplot(wddf , aes(t, Freq)) + geom_tile(aes(fill = AIC_diff_thresholded)) + 
      geom_vline(xintercept = 0, lty = "dashed", color = "white", size = 2) +
      geom_vline(xintercept = -0.3, lty = "dashed", color = "white", size = 1) +
      scale_fill_viridis(option = "plasma") +  xlab(rt_epoch_label) + ylab("Frequency") +
      # facet_wrap( ~ node, ncol = 2) + 
      geom_text(data = wddf, x = -0.4, y = 5,aes(label = "Response(t)"), size = 2.5, color = "white", angle = 90) +
      geom_text(data = wddf, x = 0.1, y = 5,aes(label = "Outcome(t)"), size = 2.5, color = "white", angle = 90) +
      scale_x_continuous(limits = c(-1,1.3), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1, 1.2)) +
      labs(fill = expression("AIC"["full"] - AIC["selective"])) + ggtitle("Evidence of value compression") + theme_dark()
    dev.off()
  }
}
