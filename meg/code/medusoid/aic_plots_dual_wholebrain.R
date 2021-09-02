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

clock_epoch_label = "Time relative to clock onset, seconds"
rt_epoch_label = "Time relative to outcome, seconds"
encode = T
rt_predict = F
# regressors = c("reward")
p_adjust_method = "fdr"
regressors = c("entropy_change")
# regressors = c("entropy", "kld", "entropy_change", "entropy_change_neg", "entropy_change_pos", "reward")
print_filenames = T
fixed_only = F
reprocess = F
plots = T
diags = F
average = F
noclock = F
setwd(data_dir)
# plots ----
for (regressor in regressors) {
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
        if (ncol(df)==3) {
          df <- df$fit_df
        }
        #        df <- df %>% filter(effect=="fixed")
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
        if (ncol(df)==3) {
          df <- df$fit_df
        }
        #        df <- df %>% filter(effect=="fixed")
        return(df)
      })
      cfdf <- data.table::rbindlist(cfl) %>% unique() %>% distinct(Time, Freq, .keep_all = TRUE)
      cfdf <- cfdf %>% mutate(t  = as.numeric(Time), alignment = "clock", model = "full",
      )
      message("Processed clock-aligned full. \n")
      
    }
    # get RT-aligned
    if (regressor=="entropy_change") {
      file_pattern <- "*_change_rs_single_.*RT"} else if (regressor=="entropy") {
        file_pattern <- "_entropy_rs.*RT"} else if (regressor=="kld") {
          file_pattern <- ".*kld.*RT"} else if (regressor=="entropy_change_pos") {
            file_pattern <- ".*entropy_change_pos_rs.*RT"} else if (regressor=="entropy_change_neg") {
              file_pattern <- ".*entropy_change_neg_rs.*RT"}  else if (regressor=="reward") {
                file_pattern <- ".*reward_rs.*RT"} 
    # file_pattern <- "ddf_combined_entropy_rsRT|ddf_combined_entropy_change_rs_RT"
    # file_pattern <- "meg_mixed_by_tf_ddf_wholebrain_entropy_change_rs_RT|meg_mixed_by_tf_ddf_wholebrain_entropy_change_rs_finishRT"
    # file_pattern <- "entropy_rs_singleRT"
    files <-  gsub("//", "/", list.files(data_dir, pattern = file_pattern, full.names = F))
    message(paste0("Found ", length(files), " files."))
    rl <- lapply(files, function(x) {
      if (print_filenames) { print(x) }
      df <- readRDS(x) 
      if (ncol(df)==3) {
        df <- df$coef_df_reml
      }
      #      df <- df %>% filter(effect=="fixed")
      return(df)
    })
    rddf <- data.table::rbindlist(rl)  %>% unique()  %>% distinct(Time, Freq, term, effect, group, level, .keep_all = TRUE)
    # rddf$node <- sub("_group.*", "", rddf$.filename)
    rddf$alignment <- "RT"
    if (!noclock) {offset = 4.3} else {offset = 0.3}
    rddf <- rddf %>% mutate(t  = Time - offset, 
                            alignment = "rt",
                            term = str_replace(term, "rt_csv_sc", "RT_t"),
                            term = str_replace(term, "outcomeReward", "reward_t"),
                            term = str_replace(term, "v_entropy_wi_change", "entropy_change_t"),
                            term = str_replace(term, "entropy_change_neg_wi", "entropy_change_neg_t"),
                            term = str_replace(term, "entropy_change_pos_wi", "entropy_change_pos_t")
    )
    # saveRDS(rddf, file = "meg_ddf_wholebrain_ec_rs_rt.rds")
    message("Processed RT-aligned, merging.  \n")
    if (!noclock) {ddf <- rbind(cddf, rddf)} else {
      ddf <- rddf  
    }
    ddf <- ddf %>% mutate(p_value = as.factor(case_when(`p.value` > .05 ~ '1',
                                                        `p.value` < .05 & `p.value` > .01 ~ '2',
                                                        `p.value` < .01 & `p.value` > .001 ~ '3',
                                                        `p.value` <.001 ~ '4'))                        # sensor = as.character(sensor)
    )
    ddf$p_value <- factor(ddf$p_value, labels = c("NS", "p < .05", "p < .01", "p < .001"))
    ddf <- ddf  %>% ungroup() %>% group_by(term, alignment) %>% mutate(p_fdr = p.adjust(p.value, method = p_adjust_method),
                                                                       p_level_fdr = as.factor(case_when(
                                                                         # p_fdr > .1 ~ '0',
                                                                         # p_fdr < .1 & p_fdr > .05 ~ '1',
                                                                         p_fdr > .05 ~ '1',
                                                                         p_fdr < .05 & p_fdr > .01 ~ '2',
                                                                         p_fdr < .01 & p_fdr > .001 ~ '3',
                                                                         p_fdr <.001 ~ '4'))) %>% ungroup()
    
    ddf$p_level_fdr <- factor(ddf$p_level_fdr, levels = c('1', '2', '3', '4'), labels = c("NS","p < .05", "p < .01", "p < .001"))
    ddf$`p, FDR-corrected` = ddf$p_level_fdr
    ddf$Freq <- gsub("f_", "", ddf$Freq)
    ddf$Freq <- ordered(as.numeric(substr(as.character(ddf$Freq), 1,4)))    
    setwd(paste0(plot_dir, "/", regressor))
    # saveRDS(ddf, file = "meg_ddf_e_ec_rs.rds")
    saveRDS(ddf, file = paste0("meg_ddf_wholebrain_", regressor, ".rds"))} else {
      setwd(paste0(plot_dir, "/", regressor))
      ddf <- readRDS(paste0("meg_ddf_wholebrain_", regressor, ".rds"))
      ddf <- ddf  %>% ungroup() %>% 
        filter((Time < 1.5 & Time > -1 & alignment=="rt") | (Time < 1.5 & Time > -2 & alignment=="clock")) %>% 
        group_by(term, alignment) %>% mutate(p_fdr = p.adjust(p.value, method = p_adjust_method),
                                             p_level_fdr = as.factor(case_when(
                                               # p_fdr > .1 ~ '0',
                                               # p_fdr < .1 & p_fdr > .05 ~ '1',
                                               p_fdr > .05 ~ '1',
                                               p_fdr < .05 & p_fdr > .01 ~ '2',
                                               p_fdr < .01 & p_fdr > .001 ~ '3',
                                               p_fdr <.001 ~ '4'))) %>% ungroup()
      ddf$p_level_fdr <- factor(ddf$p_level_fdr, levels = c('1', '2', '3', '4'), labels = c("NS","p < .05", "p < .01", "p < .001"))
      ddf$`p, FDR-corrected` = ddf$p_level_fdr
    }
  terms <- unique(ddf$term[ddf$effect=="fixed"])
  if (plots) {
    message("Plotting.  \n")
    if (!noclock) {offset = 4.3} else {offset = 0.3}
    for (fe in terms) {
      edf <- ddf %>% filter(term == paste(fe) & effect=="fixed")
      termstr <- str_replace_all(fe, "[^[:alnum:]]", "_")
      message(termstr)
      fname = paste("meg_tf_combined_uncorrected_", termstr, ".pdf", sep = "")
      if (!noclock) { # plots both
        pdf(fname, width = 10, height = 5)
        print(ggplot(edf, aes(t, Freq)) + geom_tile(aes(fill = estimate, alpha = p_value), size = .01) +
                geom_vline(xintercept = 0, lty = "dashed", color = "black", size = 2) +
                geom_vline(xintercept = -offset + 0.3, lty = "dashed", color = "white", size = 2) +
                geom_vline(xintercept = -offset, lty = "dashed", color = "white", size = 1) +
                geom_vline(xintercept = -2, lty = "dotted", color = "grey", size = 1) +
                scale_fill_viridis(option = "plasma") +  xlab(rt_epoch_label) + ylab("Frequency") +
                geom_text(data = edf, x = -offset-.1, y = 5,aes(label = "Response(t)"), size = 2.5, color = "white", angle = 90) +
                geom_text(data = edf, x = -offset+.4, y = 5,aes(label = "Outcome(t)"), size = 2.5, color = "white", angle = 90) +
                geom_text(data = edf, x = 0.5, y = 6 ,aes(label = "Clock onset (t+1)"), size = 2.5, color = "black", angle = 90) +
                scale_x_continuous(limits = c(-5.3,1.7), breaks = c(0-offset+0.3, 0.20-offset+0.3, 0.4-offset+0.3, 0.60-offset+0.3, 1-offset+0.3, 1.5-offset+0.3, 0, 1), labels = c("0", "0.2", "0.4", "0.6", "1", "1.5",  "0", "1")) +
                labs(alpha = expression(italic(p)[uncorrected])) + ggtitle(paste(termstr)) + theme_dark())
        dev.off()
        fname = paste("meg_tf_rt_all_dan_FDR_", termstr, ".pdf", sep = "")
        pdf(fname, width = 10, height = 5)
        print(ggplot(edf, aes(t, Freq)) + geom_tile(aes(fill = estimate, alpha = p_level_fdr), size = .01) +
                geom_vline(xintercept = 0, lty = "dashed", color = "black", size = 2) +
                geom_vline(xintercept = -offset + 0.3, lty = "dashed", color = "white", size = 2) +
                geom_vline(xintercept = -offset, lty = "dashed", color = "white", size = 1) +
                geom_vline(xintercept = -2, lty = "dotted", color = "grey", size = 1) +
                scale_fill_viridis(option = "plasma") +  xlab(rt_epoch_label) + ylab("Frequency") +
                geom_text(data = edf, x = -offset-.1, y = 5,aes(label = "Response(t)"), size = 2.5, color = "white", angle = 90) +
                geom_text(data = edf, x = -offset+.4, y = 5,aes(label = "Outcome(t)"), size = 2.5, color = "white", angle = 90) +
                geom_text(data = edf, x = 0.5, y = 6 ,aes(label = "Clock onset (t+1)"), size = 2.5, color = "black", angle = 90) +
                scale_x_continuous(limits = c(-5.3,1.7), breaks = c(0-offset+0.3, 0.20-offset+0.3, 0.4-offset+0.3, 0.60-offset+0.3, 1-offset+0.3, 1.5-offset+0.3, 0, 1), labels = c("0", "0.2", "0.4", "0.6", "1", "1.5",  "0", "1")) +
                labs(alpha = expression(italic(p)[corrected])) + ggtitle(paste(termstr)) + theme_dark())    # 
        dev.off() }
      else { # plots only feedback-aligned
        fname = paste("meg_tf_combined_uncorrected_", termstr, ".pdf", sep = "")
        pdf(fname, width = 7, height = 4.5)
        print(ggplot(edf, aes(t, Freq)) + geom_tile(aes(fill = estimate, alpha = p_value), size = .01) +
                geom_vline(xintercept = 0, lty = "dashed", color = "white", size = 2) +
                geom_vline(xintercept = -0.3, lty = "dashed", color = "white", size = 1) +
                scale_fill_viridis(option = "plasma") +  xlab(rt_epoch_label) + ylab("Frequency") +
                geom_text(data = edf, x = -0.4, y = 5,aes(label = "Response(t)"), size = 2.5, color = "white", angle = 90) +
                geom_text(data = edf, x = 0.1, y = 5,aes(label = "Outcome(t)"), size = 2.5, color = "white", angle = 90) +
                scale_x_continuous(limits = c(-1,1.7), breaks = c(0, 0.2, 0.4, 0.6, 0.6, 1, 1.5)) +
                labs(alpha = expression(italic(p)[uncorrected])) + ggtitle(paste(termstr)) + theme_dark())
        dev.off() 
        fname = paste("meg_tf_rt_all_dan_FDR_", termstr, ".pdf", sep = "")
        pdf(fname, width = 7, height = 4.5)
        print(ggplot(edf, aes(t, Freq)) + geom_tile(aes(fill = estimate, alpha = p_level_fdr), size = .01) +
                geom_vline(xintercept = 0, lty = "dashed", color = "white", size = 2) +
                geom_vline(xintercept = -0.3, lty = "dashed", color = "white", size = 1) +
                scale_fill_viridis(option = "plasma") +  xlab(rt_epoch_label) + ylab("Frequency") +
                # facet_wrap( ~ node, ncol = 2) + 
                geom_text(data = edf, x = -0.4, y = 5,aes(label = "Response(t)"), size = 2.5, color = "white", angle = 90) +
                geom_text(data = edf, x = 0.1, y = 5,aes(label = "Outcome(t)"), size = 2.5, color = "white", angle = 90) +
                scale_x_continuous(limits = c(-1,1.7), breaks = c(0, 0.2, 0.4, 0.6, 0.6, 1, 1.5)) +
                labs(alpha = expression(italic(p)[corrected])) + ggtitle(paste(termstr)) + theme_dark())    # 
        dev.off() }
    }
  }
}

