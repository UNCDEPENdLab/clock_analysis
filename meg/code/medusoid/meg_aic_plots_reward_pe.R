# AIC comparison between alternative PE models
# rew + abs(PE) will be the reference, based on intuitive coefficient maps
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
# data_dir <- "~/OneDrive/collected_letters/papers/meg/plots/wholebrain/output"
data_dir <- "/Volumes/GoogleDrive/.shortcut-targets-by-id/1ukjK6kTlaR-LXIqX6nylYOPWu1j3XGyF/SCEPTIC_fMRI/MEG/output"
plot_dir <- "~/OneDrive/collected_letters/papers/meg/plots/wholebrain/"

rt_epoch_label = "Time relative to event, seconds"
encode = T
rt_predict = F
# regressors = c("reward")
p_adjust_method = "fdr"

reference = "reward"
comparator = "abs_pe_rs"
# regressors = c("entropy", "kld", "entropy_change", "entropy_change_neg", "entropy_change_pos", "reward")
print_filenames = T
fixed_only = F
reprocess = F
plots = T
aic_threshold = 1000
diags = F
noclock = T
freq_threshold = 18 # set to 40 Hz (18th band) for full-spectrum output
setwd(data_dir)
# plots ----

message(paste0("Processing ", reference))
epoch_label = "Time relative to feedback, seconds"
if (reprocess) {
  if (reference=="reward") {
    ref_file_pattern <- "*_reward_rs_.*RT"} 
  files <-  gsub("//", "/", list.files(data_dir, pattern = ref_file_pattern, full.names = F))
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
  rsdf <- rsdf %>% mutate(t  = as.numeric(Time), alignment = "rt", model = "reference",
  )
  message("Processed reference. \n")
  if (comparator=="signed_pe_rs") {
    comp_file_pattern <- "*signed_pe_rs*"} else if (comparator=="abs_pe_rs") {
      comp_file_pattern <- "*abs_pe_rs*"
    }
  files <-  gsub("//", "/", list.files(data_dir, pattern = comp_file_pattern, full.names = F))
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
  rfdf <- rfdf %>% mutate(t  = as.numeric(Time), alignment = "rt", model = "comparator",
  )
  message("Processed RT-aligned full. \n")
  ddf <- rbind(rfdf, rsdf)
  
  
  offset = 0.3
  ddf <- ddf %>% mutate(t  = Time - offset)
  ddf$Freq <- gsub("f_", "", ddf$Freq)
  ddf$Freq <- ordered(as.numeric(substr(as.character(ddf$Freq), 1,4)))    
  ddf <- ddf %>% arrange(Time, Freq, alignment, model) %>% mutate(row = row_number()) 
  wddf <- pivot_wider(ddf %>% select(Time, t, Freq, AIC, BIC, model, alignment), names_from = model, names_sep = "_", values_from = c(AIC, BIC)) %>%
      mutate(AIC_diff = AIC_reference - AIC_comparator,
             BIC_diff = BIC_reference - BIC_comparator)
  
  # convert to wide
  wddf <- pivot_wider(ddf %>% select(Time, t, Freq, AIC, BIC, model, alignment), names_from = model, names_sep = "_", values_from = c(AIC, BIC)) %>%
    mutate(AIC_diff = AIC_reference - AIC_comparator,
           BIC_diff = BIC_reference - BIC_comparator)
  setwd(paste0(plot_dir, "/", "aic_comparisons"))
  saveRDS(ddf, file = paste0("meg_ddf_wholebrain_aic_comparison_", comparator, "VS", reference, ".rds"))
  saveRDS(wddf, file = paste0("meg_wddf_wholebrain_aic_comparison_", comparator, "VS", reference, ".rds"))
} else {
  setwd(paste0(plot_dir, "/", "aic_comparisons"))
  wddf <- readRDS(paste0("meg_wddf_wholebrain_aic_comparison_", comparator, "VS", reference, ".rds"))
}

if (plots) {
  message("Processed output, plotting.  \n")
  # p1 <- ggplot(wddf %>% filter(alignment=="clock" & abs(AIC_diff) > aic_threshold), aes(Time, Freq)) + geom_tile(aes(fill = AIC_diff)) + 
  #   scale_fill_viridis(option = "plasma") + theme_dark()
  # p2 <- ggplot(wddf %>% filter(alignment=="rt" & abs(AIC_diff) > aic_threshold), aes(Time, Freq)) + geom_tile(aes(fill = AIC_diff)) + 
  #   scale_fill_viridis(option = "plasma") + theme_dark()
  # ggarrange(p2, p1)
  fname = paste0("meg_tf_aic_comparison_", comparator, "VS", reference, ".pdf")
  wddf <- wddf %>% mutate(AIC_diff_thresholded  = case_when(
    abs(AIC_diff) > aic_threshold ~ AIC_diff,
    TRUE ~ NA_real_
  )) %>% filter(as.numeric(Freq)<freq_threshold)
  
  pdf(fname, width = 7, height = 4.5)
  print(
    ggplot(wddf , aes(t, Freq)) + geom_tile(aes(fill = AIC_diff_thresholded)) + 
    geom_vline(xintercept = 0, lty = "dashed", color = "white", size = 2) +
    geom_vline(xintercept = -0.3, lty = "dashed", color = "white", size = 1) +
    scale_fill_viridis(option = "plasma") +  xlab(rt_epoch_label) + ylab("Frequency") +
    # facet_wrap( ~ node, ncol = 2) + 
    geom_text(data = wddf, x = -0.4, y = 5,aes(label = "Response(t)"), size = 2.5, color = "white", angle = 90) +
    geom_text(data = wddf, x = 0.1, y = 5,aes(label = "Outcome(t)"), size = 2.5, color = "white", angle = 90) +
    scale_x_continuous(limits = c(-1,1.3), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1, 1.2)) +
    labs(fill = expression("AIC"[reference] - "AIC"[comparator])) + ggtitle("Model comparison") + theme_dark()
  )
  dev.off()
  
}
