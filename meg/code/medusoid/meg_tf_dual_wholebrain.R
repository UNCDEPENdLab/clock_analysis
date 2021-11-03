# plots MLM coefficients by sensor across timepoints and frequencies
# first run meg_time_freq_load_analyze_medusa.R

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
p_adjust_method = "fdr"
regressors = c("abspe_by_rew")
# regressors = c("entropy", "kld", "entropy_change", "entropy_change_neg", "entropy_change_pos", "reward")
print_filenames = T
fixed_only = F
reprocess = T
plots = T
diags = F
average = F
noclock = T
freq_threshold = 19 # set to 40 Hz (19th band) for full-spectrum output
setwd(data_dir)
# plots ----
if (encode) {  
  message("Processing encoding results")
  for (regressor in regressors) {
    if (regressor=="reward") {noclock = T}
    message(paste0("Processing ", regressor))
    epoch_label = "Time relative to feedback, seconds"
    if (reprocess) {
      if (!noclock) {
        epoch_label = "Time relative to clock onset, seconds"
        # get clock-aligned data
        if (regressor=="entropy_change") {
          file_pattern <- "*_change_rs_single_.*clock"} else if (regressor=="entropy") {
            file_pattern <- "_entropy_rs.*clock"} else if (regressor=="kld") {
              file_pattern <- ".*kld.*clock"} else if (regressor=="entropy_change_pos") {
                file_pattern <- ".*entropy_change_pos_rs.*clock"} else if (regressor=="entropy_change_neg") {
                  file_pattern <- ".*entropy_change_neg_rs.*clock"} else if (regressor=="reward") {
                    file_pattern <- ".*reward_rs.*clock"} else if (regressor == "v_max") {
                      file_pattern <- ".*v_max_rs.*clock"} 
        files <-  gsub("//", "/", list.files(data_dir, pattern = file_pattern, full.names = F))
        message(paste0("Found ", length(files), " files."))
        cl <- lapply(files, function(x) {
          if (print_filenames) { print(x) }
          df <- readRDS(x) 
          # if (ncol(df)==3) {
            df <- df$coef_df_reml
          # }
          #        df <- df %>% filter(effect=="fixed")
          return(df)
        })
        cddf <- data.table::rbindlist(cl) %>% unique() %>% distinct(Time, Freq, term, effect, group, level, .keep_all = TRUE)
        cddf$alignment <- "clock"
        cddf <- cddf %>% mutate(t  = as.numeric(Time), alignment = "clock",
                                term = str_replace(term, "rt_lag_sc", "RT_t"),
                                term = str_replace(term, "reward_lagReward", "reward_t"),
                                term = str_replace(term, "v_entropy_wi_change_lag", "entropy_change_t"),
                                term = str_replace(term, "entropy_change_neg_lag", "entropy_change_neg_t"),
                                term = str_replace(term, "entropy_change_pos_lag", "entropy_change_pos_t")
        )
        message("Processed clock-aligned. \n")}
      # get RT-aligned
      if (regressor=="entropy_change") {
        file_pattern <- "*_change_rs_single_.*RT"} else if (regressor=="entropy") {
          file_pattern <- ".*_entropy_rs.*RT"} else if (regressor=="kld") {
            file_pattern <- ".*kld.*RT"} else if (regressor=="entropy_change_pos") {
              file_pattern <- ".*entropy_change_pos_rs.*RT"} else if (regressor=="entropy_change_neg") {
                file_pattern <- ".*entropy_change_neg_rs.*RT"}  else if (regressor=="reward") {
                  file_pattern <- ".*reward_rs.*RT"} else if (regressor=="v_max"){
                    file_pattern <- ".*v_max_rs.*RT"} else if (regressor=="abs_pe") {
                      file_pattern <- ".*abs_pe.*RT"} else if(regressor =="signed_pe") {
                        file_pattern <- ".*signed_pe.*"} else if(regressor =="abspe_by_rew") {
                          file_pattern <- ".*abspe_by_rew.*"
                      }
      # file_pattern <- "ddf_combined_entropy_rsRT|ddf_combined_entropy_change_rs_RT"
      # file_pattern <- "meg_mixed_by_tf_ddf_wholebrain_entropy_change_rs_RT|meg_mixed_by_tf_ddf_wholebrain_entropy_change_rs_finishRT"
      # file_pattern <- "entropy_rs_singleRT"
      files <-  gsub("//", "/", list.files(data_dir, pattern = file_pattern, full.names = F))
      message(paste0("Found ", length(files), " files."))
      rl <- lapply(files, function(x) {
        if (print_filenames) { print(x) }
        df <- readRDS(x) 
        if(class(df) == "list") {df <- df$coef_df_reml}
        # if (ncol(df)<4) {
        #   df <- df$coef_df_reml
        # }
        #      df <- df %>% filter(effect=="fixed")
        return(df)
      })
      rddf <- data.table::rbindlist(rl)  %>% unique()  %>% distinct(Time, Freq, term, effect, group, level, rhs, .keep_all = TRUE)
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
        if (freq_threshold>0) {
        edf <- ddf %>% filter(term == paste(fe) & effect=="fixed" & as.numeric(Freq) < freq_threshold)} else {
          edf <- ddf %>% filter(term == paste(fe) & effect=="fixed")
        }
        termstr <- str_replace_all(fe, "[^[:alnum:]]", "_")
        message(termstr)
        fname = paste("meg_tf_combined_uncorrected_", termstr, ".pdf", sep = "")
        if (!noclock) { # plots both
          pdf(fname, width = 6, height = 3)
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
          pdf(fname, width = 6, height = 3)
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
                  scale_x_continuous(limits = c(-1,1.3), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1, 1.2)) +
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
                  scale_x_continuous(limits = c(-1,1.3), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1, 1.2)) +
                  labs(alpha = expression(italic(p)[corrected])) + ggtitle(paste(termstr)) + theme_dark())    # 
          dev.off() }
      }
    }
  }
}

# system("for i in *FDR*.pdf; do sips -s format png $i --out $i.png; done")

if (rt_predict) {
  setwd(data_dir)
  # plots ----
  epoch_label = "Time relative to clock onset, seconds"
  # get clock-aligned data
  file_pattern <- ".*rdf.*clock.*"
  message("Processing ")
  files <-  gsub("//", "/", list.files(data_dir, pattern = file_pattern, full.names = F))
  cl <- lapply(files, function(x) {
    print(x)
    df <- readRDS(x) %>% dplyr::filter(effect=="fixed")
    return(df)
  })
  # what's the deal with 0.33?  They have the "f_" prefix in Freq
  # df33 <- data.table::rbindlist(cl) %>% filter(grepl('Pow', term) & (Time > 0.329 & Time < .331))
  # need to remove the prefix
  crdf <- data.table::rbindlist(cl) %>% filter(grepl('Pow', term))
  # crdf <- data.table::rbindlist(cl) %>% filter(grepl('Pow', term))& (Time < 0.329 | Time > .331))
  # crdf$node <- sub("_group.*", "", crdf$.filename)
  crdf <- crdf %>% mutate(t  = as.numeric(Time), alignment = "clock",
                          term = str_replace_all(term, "[^[:alnum:]]", ""),
                          term = str_replace(term, "scalePow", "Power"),
                          term = str_replace(term, "rtlagsc", "*RT_t"),
                          term = str_replace(term, "rtvmaxlagsc", "*RT_Vmax_t"),
                          term = str_replace(term, "rewardlagReward", "*Reward_t"),
                          term = str_replace(term, "rtlag2sc", "*RT_tMINUS1"), 
                          Freq = gsub("f_", "", Freq) # removes prefix for 0.33 s
  )
  # get RT-aligned
  message("Processed clock-aligned")
  file_pattern <-  ".*rdf.*RT.*"
  files <-  gsub("//", "/", list.files(data_dir, pattern = file_pattern, full.names = F))
  rl <- lapply(files, function(x) {
    print(x)
    df <- readRDS(x) %>% dplyr::filter(effect=="fixed")
    return(df)
  })
  rrdf <- data.table::rbindlist(rl) %>% filter(grepl('Pow', term))
  rrdf$node <- sub("_group.*", "", rrdf$.filename)
  rrdf <- rrdf %>% mutate(t  = as.numeric(Time) - 5, alignment = "RT",
                          term = str_replace_all(term, "[^[:alnum:]]", ""),
                          term = str_replace(term, "scalePow", "Power"),
                          term = str_replace(term, "rtcsvsc", "*RT_t"),
                          term = str_replace(term, pattern = "scalertvmax", replacement = "*RT_Vmax_t"),
                          term = str_replace(term, "rtlagsc", "*RT_tMINUS1"),
                          term = str_replace(term, "outcomeReward", "*Reward_t"),
                          Freq = gsub("f_", "", Freq) # removes prefix for 0.33 s
  )
  message("Processed RT-aligned, merging")
  rdf <- rbind(droplevels(crdf), droplevels(rrdf))
  terms <- unique(rdf$term[rdf$effect=="fixed"]) 
  terms <- terms[grepl("Power", terms)]
  # rdf$sensor <- readr::parse_number(rdf$.filename, trim_ws = F)
  # rdf$sensor <- stringr::str_pad(rdf$sensor, 4, "0",side = "left")
  rdf <- rdf %>% mutate(p_value = as.factor(case_when(`p.value` > .05 ~ '1',
                                                      `p.value` < .05 & `p.value` > .01 ~ '2',
                                                      `p.value` < .01 & `p.value` > .001 ~ '3',
                                                      `p.value` <.001 ~ '4'))                        # sensor = as.character(sensor)
  )
  rdf$p_value <- factor(rdf$p_value, labels = c("NS", "p < .05", "p < .01", "p < .001"))
  
  # FDR labeling ----
  rdf <- rdf  %>% mutate(
    # p_fdr = padj_fdr_term,
    p_level_fdr = as.factor(case_when(
      # p_fdr > .1 ~ '0',
      # p_fdr < .1 & p_fdr > .05 ~ '1',
      padj_BY_term > .05 ~ '1',
      padj_BY_term < .05 & padj_BY_term > .01 ~ '2',
      padj_BY_term < .01 & padj_BY_term > .001 ~ '3',
      padj_BY_term <.001 ~ '4'))
  ) %>% ungroup() #%>% mutate(side = substr(as.character(label), nchar(as.character(label)), nchar(as.character(label))),
  # region = substr(as.character(label), 1, nchar(as.character(label))-2))
  rdf$p_level_fdr <- factor(rdf$p_level_fdr, levels = c('1', '2', '3', '4'), labels = c("NS","p < .05", "p < .01", "p < .001"))
  rdf$`p, FDR-corrected` = rdf$p_level_fdr
  # levels(rdf$Freq) <- signif(as.numeric(substr(levels(rdf$Freq), 1,3)), 3)
  rdf$Freq <- ordered(as.numeric(substr(as.character(rdf$Freq), 1,4)))
  setwd(data_dir)
  saveRDS(rdf, "meg_rt_predict_rdf_Aug17.Rds")
  # rdf$Freq <- fct_rev(rdf$Freq)
  # rdf$numFreq <- as.numeric(paste(rdf$Freq))
  # drop magnetometers (if any)
  # rdf <- rdf %>% filter(!grepl("1$", sensor))
  # add sensor labels
  
  fef_sensors <- c("612","613", "542", "543","1022")
  ips_sensors <- c("1823", "1822", "2222","2223")
  
  dan_sensors <- c(fef_sensors, ips_sensors)
  dan_rdf <- rdf %>% filter(Sensor %in% dan_sensors)
  
  setwd(plot_dir)
  
  
  for (fe in terms) {
    edf <- dan_rdf %>% filter(term == paste(fe) & effect=="fixed") 
    termstr <- str_replace_all(fe, "[^[:alnum:]]", "")
    message(termstr)
    fname = paste("RT_predict_uncorrected_", termstr, ".pdf", sep = "")
    if (termstr=="Power") {lolim = -.1; hilim = .075; hilabel = "Suppression"; lolabel = "Synchronization"
    } else if (termstr=="PowerRTt") {lolim = -.05; hilim = .2; hilabel = "Suppression"; lolabel = "Synchronization"
    } else if (termstr=="PowerRTVmaxt") {lolim = -.01; hilim = .023; hilabel = "Synchronization"; lolabel = "Suppression"
    } else if (termstr=="PowerRTtMINUS1") {lolim = -.01; hilim = .023; hilabel = "Synchronization"; lolabel = "Suppression"
    } else if (termstr=="PowerRewardt") {lolim = -.01; hilim = .023; hilabel = "Synchronization"; lolabel = "Suppression"
    } else if (termstr=="PowerRTtRewardt") {lolim = -.02; hilim = .04; hilabel = "Synchronization"; lolabel = "Suppression"
    } else  {lolim = -.03; hilim = .02}
    pdf(fname, width = 10, height = 7.5)
    print(ggplot(edf %>% filter(estimate < 0), aes(t, Freq)) + geom_tile(aes(fill = estimate, alpha = p_value), size = .01) +
            scale_fill_distiller(palette = "Oranges", direction = 1, name = lolabel, limits = c(lolim, 0)) + scale_x_continuous(breaks = pretty(edf$t, n = 5)) + labs(fill = lolabel) +
            new_scale_fill() +
            geom_tile(data = edf %>% filter(estimate > 0), aes(t, Freq, fill = estimate, alpha = p_value), size = .01) +
            scale_y_discrete(limits = levels(edf$Freq)) +
            scale_fill_distiller(palette = "YlGnBu", direction = -1, name = hilabel, limits = c(0, hilim)) +
            scale_color_grey() + xlab(epoch_label) + ylab("Frequency") +
            # facet_wrap( ~ node, ncol = 2) +
            facet_wrap( ~ Sensor) +
            geom_vline(xintercept = 0, lty = "dashed", color = "white", size = 1) + theme_black() +
            geom_vline(xintercept = -5, lty = "dashed", color = "white", size = 1) +
            geom_vline(xintercept = -5.3, lty = "dashed", color = "white", size = .5) +
            geom_vline(xintercept = -2.5, lty = "dotted", color = "grey", size = .5) +
            geom_text(data = edf, x = -6, y = 5,aes(label = "Response(t)"), size = 3, color = "grey", angle = 90) +
            geom_text(data = edf, x = -4.5, y = 5,aes(label = "Outcome(t)"), size = 3, color = "grey", angle = 90) +
            geom_text(data = edf, x = -0.75, y = 7.5 ,aes(label = "Clock onset (t+1)"), size = 3, color = "grey", angle = 90) +
            labs(alpha = expression(italic(p)[uncorrected]), fill = hilabel) + ggtitle(paste(termstr)))
    
    dev.off()
    # fname = paste("RT_predict_FDR_", termstr, ".pdf", sep = "")
    # if (termstr=="Power") {lolim = -.1; hilim = .075; hilabel = "Suppression"; lolabel = "Synchronization"
    # } else if (termstr=="PowerRTt") {lolim = -.03; hilim = .1; hilabel = "Suppression"; lolabel = "Synchronization"
    # } else if (termstr=="PowerRTVmaxt") {lolim = -.01; hilim = .023; hilabel = "Synchronization"; lolabel = "Suppression"
    # } else if (termstr=="PowerRTtMINUS1") {lolim = -.01; hilim = .023; hilabel = "Synchronization"; lolabel = "Suppression"
    # } else if (termstr=="PowerRewardt") {lolim = -.01; hilim = .023; hilabel = "Synchronization"; lolabel = "Suppression"
    # } else if (termstr=="PowerRTtRewardt") {lolim = -.02; hilim = .04; hilabel = "Synchronization"; lolabel = "Suppression"
    # } else  {lolim = -.03; hilim = .02}
    # pdf(fname, width = 10, height = 9)
    # print(ggplot(edf %>% filter(estimate < 0), aes(t, Freq)) + geom_tile(aes(fill = estimate, alpha = p_level_fdr), size = .01) + 
    #         scale_fill_distiller(palette = "Oranges", direction = 1, name = lolabel, limits = c(lolim, 0)) + scale_x_continuous(breaks = pretty(edf$t, n = 5)) + labs(fill = lolabel) +
    #         new_scale_fill() +
    #         geom_tile(data = edf %>% filter(estimate > 0), aes(t, Freq, fill = estimate, alpha = p_level_fdr), size = .01) +
    #         scale_y_discrete(limits = levels(edf$Freq)) +
    #         scale_fill_distiller(palette = "YlGnBu", direction = -1, name = hilabel, limits = c(0, hilim)) + 
    #         scale_color_grey() + xlab(epoch_label) + ylab("Frequency") +
    #         # facet_wrap( ~ node, ncol = 2) + 
    #         facet_wrap( ~ Sensor) +
    #         geom_vline(xintercept = 0, lty = "dashed", color = "white", size = 1) + theme_black() + 
    #         geom_vline(xintercept = -5, lty = "dashed", color = "white", size = 1) +
    #         geom_vline(xintercept = -5.3, lty = "dashed", color = "white", size = .5) +
    #         geom_vline(xintercept = -2.5, lty = "dotted", color = "grey", size = .5) +
    #         geom_text(data = edf, x = -6, y = 5,aes(label = "Response(t)"), size = 3, color = "grey", angle = 90) +
    #         geom_text(data = edf, x = -4.5, y = 5,aes(label = "Outcome(t)"), size = 3, color = "grey", angle = 90) +
    #         geom_text(data = edf, x = -0.75, y = 7.5 ,aes(label = "Clock onset (t+1)"), size = 3, color = "grey", angle = 90) +
    #         labs(alpha = expression(italic(p)[FDR-corrected]), fill = hilabel) + ggtitle(paste(termstr)))
    # dev.off()
  }
  # legend:
  # Power           - main effect of power predicting sooner/later responses
  # PowerRTt        - power predicting RT swings
  # PowerRewardt    - power predicting post-reward RT shortening
  # PowerRTVmaxt    - power predicting convergence on RT_Vmax
  # PowerRTtMINUS1  - power predicting reselection/swing away from RT(t-1)
  # PowerRTtRewardt - power predicting win-stay/lose-shift
  
  # convert to PNG for Word
  # system("for i in *meg_tf*.pdf; do sips -s format png $i --out $i.png; done")
  # save rdf
  saveRDS(rdf, "meg_tf_rs_RT_prediction_rdf.rds")
}

# diagnostics on random slopes
if (diags) {
  # encoding
  ddfe_sensor <- ddf %>% filter(effect=="ran_coefs" & term=="v_entropy_wi" & group=="sensor") 
  # dimensions are sensor*time*frequency*alignment {RT, clock}, 72*84*22*2
  ggplot(ddfe_sensor, aes(estimate)) + geom_histogram() + facet_wrap(~node, ncol = 2)
  
  ddfe_subject <- ddf %>% filter(effect=="ran_coefs" & term=="v_entropy_wi" & group=="Subject") 
  # dimensions are sensor*time*frequency*alignment {RT, clock}, 72*84*22*2
  ggplot(ddfe_subject, aes(estimate)) + geom_histogram() + facet_wrap(~node, ncol = 2)
  
  ddfec_sensor <- ddf %>% filter(effect=="ran_coefs" & term=="entropy_change_t" & group=="sensor") 
  # dimensions are sensor*time*frequency*alignment {RT, clock}, 72*84*22*2
  ggplot(ddfec_sensor, aes(estimate)) + geom_histogram() + facet_wrap(~node, ncol = 2)
  
  ddfec_subject <- ddf %>% filter(effect=="ran_coefs" & term=="entropy_change_t" & group=="Subject") 
  # dimensions are sensor*time*frequency*alignment {RT, clock}, 72*84*22*2
  ggplot(ddfec_subject, aes(estimate)) + geom_histogram() + facet_wrap(~node, ncol = 2)
  
  # RT prediction
  rdfp_sensor <- rdf %>% filter(effect=="ran_coefs" & term=="Power" & group=="sensor") 
  # dimensions are sensor*time*frequency*alignment {RT, clock}, 72*84*22*2
  ggplot(rdfp_sensor, aes(estimate)) + geom_histogram() + facet_wrap(~node, ncol = 2)
  
  rdfp_subject <- rdf %>% filter(effect=="ran_coefs" & term=="Power" & group=="id") 
  # dimensions are sensor*time*frequency*alignment {RT, clock}, 72*84*22*2
  ggplot(rdfp_subject, aes(estimate)) + geom_histogram() + facet_wrap(~node, ncol = 2)
  
  rdfx_sensor <- rdf %>% filter(effect=="ran_coefs" & term=="RT_t * Power" & group=="sensor") 
  # dimensions are sensor*time*frequency*alignment {RT, clock}, 72*84*22*2
  ggplot(rdfx_sensor, aes(estimate)) + geom_histogram() + facet_wrap(~node, ncol = 2)
  
  rdfx_subject <- rdf %>% filter(effect=="ran_coefs" & term=="RT_t * Power" & group=="id") 
  # dimensions are sensor*time*frequency*alignment {RT, clock}, 72*84*22*2
  ggplot(rdfx_subject, aes(estimate)) + geom_histogram() + facet_wrap(~node, ncol = 2)
  
  rdfv_sensor <- rdf %>% filter(effect=="ran_coefs" & term=="RT_Vmax_t * Power" & group=="sensor") 
  # dimensions are sensor*time*frequency*alignment {RT, clock}, 72*84*22*2
  ggplot(rdfx_sensor, aes(estimate)) + geom_histogram() + facet_wrap(~node, ncol = 2)
  
  rdfv_subject <- rdf %>% filter(effect=="ran_coefs" & term=="RT_Vmax_t * Power" & group=="id") 
  # dimensions are sensor*time*frequency*alignment {RT, clock}, 72*84*22*2
  ggplot(rdfv_subject, aes(estimate)) + geom_histogram() + facet_wrap(~node, ncol = 2)
  
}

if (average) {
  # average RT*power effect across rewards and punishments
  # take the fixed effects of interest
  rtdf <- rdf %>% filter(effect=="fixed" & (term == "Power*RT_t" |  term == "Power*RT_trewardlagReward" ))
  # make it wide
  
  wdf <- rtdf %>% pivot_wider(names_from = "term", values_from = c(estimate, conf.low, conf.high, p_level_fdr), id_cols = c(Freq, t, alignment, node))
  wdf <- wdf %>% rowwise() %>% mutate(
    coef_power_rt_mean = as.matrix(mean(`estimate_Power*RT_t`, `estimate_Power*RT_trewardlagReward`, trim = 0, na.rm = TRUE)),
    p_level_fdr = min(ordered(`p_level_fdr_Power*RT_t`), ordered(`p_level_fdr_Power*RT_trewardlagReward`), na.rm = T)
  )
  # crude plot
  lolim = -.025; hilim = .04; hilabel = "Suppression"; lolabel = "Synchronization"
  ggplot(wdf %>% filter(coef_power_rt_mean < 0) , aes(t, Freq)) + geom_tile(aes(fill = coef_power_rt_mean, alpha = p_level_fdr), size = .01) + 
    scale_fill_distiller(palette = "Oranges", direction = 1, name = lolabel, limits = c(lolim, 0)) + scale_x_continuous(breaks = pretty(edf$t, n = 5)) + labs(fill = lolabel) +
    new_scale_fill() +
    geom_tile(data = wdf %>% filter(coef_power_rt_mean > 0), aes(t, Freq, fill = coef_power_rt_mean, alpha = p_level_fdr), size = .01) + theme_dark() +
    geom_vline(xintercept = 0, lty = "dashed", color = "white", size = 1) + theme_black() + 
    scale_fill_distiller(palette = "YlGnBu", direction = -1, name = hilabel, limits = c(0, hilim))+ scale_color_grey() + xlab(epoch_label) + ylab("Frequency") + 
    geom_vline(xintercept = -5, lty = "dashed", color = "white", size = 1) +
    geom_vline(xintercept = -5.3, lty = "dashed", color = "white", size = .5) +
    geom_vline(xintercept = -2.5, lty = "dotted", color = "grey", size = .5) +
    facet_wrap( ~ node, ncol = 2) + 
    geom_text(data = edf, x = -6, y = 5,aes(label = "Response(t)"), size = 3, color = "grey", angle = 90) +
    geom_text(data = edf, x = -4.5, y = 5,aes(label = "Outcome(t)"), size = 3, color = "grey", angle = 90) +
    geom_text(data = edf, x = -0.75, y = 7.5 ,aes(label = "Clock onset (t+1)"), size = 3, color = "grey", angle = 90) + 
    labs(alpha = expression(italic(p)[FDR-corrected]))
  
  
  
}

