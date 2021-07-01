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

# library(psych)
repo_directory <- "~/code/clock_analysis"
data_dir <- "~/OneDrive/collected_letters/papers/meg/plots/tf_combined_rs/output"
plot_dir <- "~/OneDrive/collected_letters/papers/meg/plots/tf_combined_rs/"
# rt_encode_plot_dir = "~/OneDrive/collected_letters/papers/meg/plots/rt_encode/"  
# clock_encode_plot_dir = "~/OneDrive/collected_letters/papers/meg/plots/clock_encode/"  
# dual_encode_plot_dir = "~/OneDrive/collected_letters/papers/meg/plots/dual_encode/"  
# 
# rt_rt_plot_dir = "~/OneDrive/collected_letters/papers/meg/plots/rt_rt/"
# clock_rt_plot_dir = "~/OneDrive/collected_letters/papers/meg/plots/clock_rt/"
# dual_rt_plot_dir = "~/OneDrive/collected_letters/papers/meg/plots/dual_rt/"

clock_epoch_label = "Time relative to clock onset, seconds"
rt_epoch_label = "Time relative to outcome, seconds"
encode = T
rt_predict = F
plots = T
diags = F
average = F
noclock = T
setwd(data_dir)
# plots ----
if (encode) {  
  message("Processing encoding results")
  epoch_label = "Time relative to feedback, seconds"
  
  # encode_formula_e = formula(~ scale(rt_vmax_lag)*echange_f1_early + scale(rt_vmax_lag)*echange_f2_late + scale(rt_vmax_lag)*e_f1 +
  #                              scale(abs_pe)*echange_f1_early + scale(abs_pe)*echange_f2_late + scale(abs_pe)*e_f1 +
  #                              outcome*echange_f1_early + outcome*echange_f2_late + outcome*e_f1 +
  #                              rt_csv_sc*echange_f1_early + rt_csv_sc*echange_f2_late + rt_csv_sc*e_f1 +
  #                              trial_neg_inv_sc*echange_f1_early + trial_neg_inv_sc*echange_f2_late + trial_neg_inv_sc*e_f1 +
  #                              v_entropy_wi_change*echange_f1_early + v_entropy_wi_change*echange_f2_late + v_entropy_wi*e_f1 + rt_lag_sc*e_f1 + (1|Subject) + (1|sensor))
  if (!noclock) {
    epoch_label = "Time relative to clock onset, seconds"
    # get clock-aligned data
    # file_pattern <- "ddf_combined_entropy_rs_clock|ddf_combined_entropy_change_rs_clock"
    file_pattern <- "*wholebrain*clock*"
    files <-  gsub("//", "/", list.files(data_dir, pattern = file_pattern, full.names = F))
    cl <- lapply(files, readRDS)
    # names(cl) <- 
    # names(cl) <- node_list[parse_number(files)]
    # # names(l) <- node_list[1:length(files)]
    # cddf <- data.table::rbindlist(cl, idcol = "node")
    cddf <- data.table::rbindlist(cl)
    cddf$node <- sub("_group.*", "", cddf$.filename)
    cddf$alignment <- "clock"
    cddf <- cddf %>% mutate(t  = as.numeric(Time), alignment = "clock",
                            term = str_replace(term, "rt_lag_sc", "RT_t"),
                            term = str_replace(term, "reward_lagReward", "reward_t"),
                            term = str_replace(term, "v_entropy_wi_change_lag", "entropy_change_t")
    )
    message("Processed clock-aligned")}
  # get RT-aligned
  # file_pattern <- "ddf_combined_entropy_rsRT|ddf_combined_entropy_change_rs_RT"
  file_pattern <- "meg_mixed_by_tf_ddf_wholebrain_entropy_change_rs_RT"
  files <-  gsub("//", "/", list.files(data_dir, pattern = file_pattern, full.names = F))
  rl <- lapply(files, readRDS)
  rddf <- data.table::rbindlist(rl)
  rddf$node <- sub("_group.*", "", rddf$.filename)
  rddf$alignment <- "RT"
  if (!noclock) {offset = 5} else {offset = 0}
  rddf <- rddf %>% mutate(t  = Time - offset, 
                          alignment = "rt",
                          term = str_replace(term, "rt_csv_sc", "RT_t"),
                          term = str_replace(term, "outcomeReward", "reward_t"),
                          term = str_replace(term, "v_entropy_wi_change", "entropy_change_t")
  )
  message("Processed RT-aligned, merging")
  if (!noclock) {ddf <- rbind(cddf, rddf)} else {
    ddf <- rddf  
  }
  
  terms <- unique(ddf$term[ddf$effect=="fixed"])
  # ddf$sensor <- readr::parse_number(ddf$.filename, trim_ws = F)
  # ddf$sensor <- stringr::str_pad(ddf$sensor, 4, "0",side = "left")
  ddf <- ddf %>% mutate(p_value = as.factor(case_when(`p.value` > .05 ~ '1',
                                                      `p.value` < .05 & `p.value` > .01 ~ '2',
                                                      `p.value` < .01 & `p.value` > .001 ~ '3',
                                                      `p.value` <.001 ~ '4'))                        # sensor = as.character(sensor)
  )
  ddf$p_value <- factor(ddf$p_value, labels = c("NS", "p < .05", "p < .01", "p < .001"))
  
  # FDR labeling ----
  ddf <- ddf  %>% mutate(
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
  ddf$p_level_fdr <- factor(ddf$p_level_fdr, levels = c('1', '2', '3', '4'), labels = c("NS","p < .05", "p < .01", "p < .001"))
  ddf$`p, FDR-corrected` = ddf$p_level_fdr
  
  levels(ddf$Freq) <- signif(as.numeric(substr(levels(ddf$Freq), 3,6)), 2)
  ddf$Freq <- fct_rev(ddf$Freq)
  # drop magnetometers (if any)
  # ddf <- ddf %>% filter(!grepl("1$", sensor))
  # add sensor labels
  setwd(plot_dir)
  # saveRDS(ddf, file = "meg_ddf_e_ec_rs.rds")
  saveRDS(ddf, file = "meg_ddf_wholebrain_ec_rs.rds")
  if (plots) {
    message("Plotting encoding results")
    setwd(plot_dir)
    for (fe in terms) {
      edf <- ddf %>% filter(term == paste(fe) & effect=="fixed")
      termstr <- str_replace_all(fe, "[^[:alnum:]]", "_")
      message(termstr)
      fname = paste("meg_tf_combined_uncorrected_", termstr, ".pdf", sep = "")
      pdf(fname, width = 10, height = 7.5)
      print(ggplot(edf, aes(t, Freq)) + geom_tile(aes(fill = estimate, alpha = p_value), size = .01) +
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
      
      fname = paste("meg_tf_rt_all_dan_FDR_", termstr, ".pdf", sep = "")
      pdf(fname, width = 10, height = 7.5)
      print(ggplot(edf, aes(t, Freq)) + geom_tile(aes(fill = estimate, alpha = p_level_fdr), size = .01) +
              geom_vline(xintercept = 0, lty = "dashed", color = "black", size = 2) +
              geom_vline(xintercept = -5, lty = "dashed", color = "white", size = 2) +
              geom_vline(xintercept = -5.3, lty = "dashed", color = "white", size = 1) +
              geom_vline(xintercept = -2.5, lty = "dotted", color = "grey", size = 1) +
              scale_fill_viridis(option = "plasma") +  xlab(rt_epoch_label) + ylab("Frequency") +
              # facet_wrap( ~ node, ncol = 2) + 
              geom_text(data = edf, x = -5.5, y = 5,aes(label = "Response(t)"), size = 2.5, color = "white", angle = 90) +
              geom_text(data = edf, x = -4.5, y = 5,aes(label = "Outcome(t)"), size = 2.5, color = "white", angle = 90) +
              geom_text(data = edf, x = 0.5, y = 6 ,aes(label = "Clock onset (t+1)"), size = 2.5, color = "black", angle = 90) +
              labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + theme_dark())    # 
      dev.off()
    }}
} 
# system("for i in *FDR*.pdf; do sips -s format png $i --out $i.png; done")

if (rt_predict) {
  setwd(data_dir)
  # plots ----
  epoch_label = "Time relative to clock onset, seconds"
  # get clock-aligned data
  file_pattern <- "meg_tf_rdf_combined_rt_rs_intclock"
  files <-  gsub("//", "/", list.files(data_dir, pattern = file_pattern, full.names = F))
  cl <- lapply(files, readRDS)
  crdf <- data.table::rbindlist(cl) %>% filter(grepl('Pow', term))
  crdf$node <- sub("_group.*", "", crdf$.filename)
  crdf <- crdf %>% mutate(t  = as.numeric(Time), alignment = "clock",
                          
                          term = str_replace_all(term, "[^[:alnum:]]", ""),
                          term = str_replace(term, "scalePow", "Power"),
                          term = str_replace(term, "rtlagsc", "*RT_t"),
                          term = str_replace(term, "rtvmaxlagsc", "*RT_Vmax_t"),
                          term = str_replace(term, "rewardlagReward", "*Reward_t"),
                          term = str_replace(term, "rtlag2sc", "*RT_tMINUS1")
  )
  # get RT-aligned
  message("Processed clock-aligned")
  file_pattern <- "meg_tf_rdf_combined_rt_rs_intRT"
  files <-  gsub("//", "/", list.files(data_dir, pattern = file_pattern, full.names = F))
  rl <- lapply(files, readRDS)
  rrdf <- data.table::rbindlist(rl) %>% filter(grepl('Pow', term))
  rrdf$node <- sub("_group.*", "", rrdf$.filename)
  rrdf <- rrdf %>% mutate(t  = as.numeric(Time) - 5, alignment = "RT",
                          term = str_replace_all(term, "[^[:alnum:]]", ""),
                          term = str_replace(term, "scalePow", "Power"),
                          term = str_replace(term, "rtcsvsc", "*RT_t"),
                          term = str_replace(term, pattern = "scalertvmax", replacement = "*RT_Vmax_t"),
                          term = str_replace(term, "rtlagsc", "*RT_tMINUS1"),
                          term = str_replace(term, "outcomeReward", "*Reward_t")
  )
  message("Processed RT-aligned, merging")
  rdf <- rbind(crdf, rrdf)
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
  levels(rdf$Freq) <- signif(as.numeric(substr(levels(rdf$Freq), 3,6)), 2)
  rdf$Freq <- fct_rev(rdf$Freq)
  rdf$numFreq <- as.numeric(paste(rdf$Freq))
  # drop magnetometers (if any)
  # rdf <- rdf %>% filter(!grepl("1$", sensor))
  # add sensor labels
  setwd(plot_dir)
  
  
  for (fe in terms) {
    edf <- rdf %>% filter(term == paste(fe) & effect=="fixed") 
    termstr <- str_replace_all(fe, "[^[:alnum:]]", "")
    message(termstr)
    # fname = paste("RT_predict_uncorrected_", termstr, ".pdf", sep = "")      
    # pdf(fname, width = 10, height = 7.5)
    # print(ggplot(edf %>% filter(estimate < 0), aes(t, Freq)) + geom_tile(aes(fill = estimate, alpha = p_value), size = .01) + 
    #         scale_fill_distiller(palette = "OrRd", direction = 1, name = "Exploration", limits = c(-.03, 0)) + scale_x_continuous(breaks = pretty(edf$t, n = 20)) + labs(fill = "Exploration") +
    #         labs(alpha = expression(italic(p)[uncorrected])) + ggtitle(paste(termstr)) +
    #         new_scale_fill() +
    #         geom_tile(data = edf %>% filter(estimate > 0), aes(t, Freq, fill = estimate, alpha = p_value), size = .01) + theme_dark() +
    #         geom_vline(xintercept = 0, lty = "dashed", color = "black", size = 2) + theme_bw() + 
    #         scale_fill_distiller(palette = "GnBu", direction = -1, name = "Exploitation", limits = c(0, .05))+ scale_color_grey() + xlab(epoch_label) + ylab("Frequency") + 
    #         geom_vline(xintercept = 0, lty = "dashed", color = "black", size = 2) +
    #         geom_vline(xintercept = -5, lty = "dashed", color = "black", size = 2) +
    #         geom_vline(xintercept = -5.3, lty = "dashed", color = "black", size = 1) +
    #         geom_vline(xintercept = -2.5, lty = "dotted", color = "grey", size = 1) +
    #         facet_wrap( ~ node, ncol = 2) + 
    #         geom_text(data = edf, x = -5.5, y = 5,aes(label = "Response(t)"), size = 2.5, color = "black", angle = 90) +
    #         geom_text(data = edf, x = -4.5, y = 5,aes(label = "Outcome(t)"), size = 2.5, color = "black", angle = 90) +
    #         geom_text(data = edf, x = 0.5, y = 6 ,aes(label = "Clock onset (t+1)"), size = 2.5, color = "black", angle = 90) +
    #         labs(alpha = expression(italic(p)[uncorrected])) + ggtitle(paste(termstr))) + scale_y_discrete(limits = rev(levels(rdf$p_level_fdr))) 
    # 
    # dev.off()
    fname = paste("RT_predict_FDR_", termstr, ".pdf", sep = "")
    if (termstr=="Power") {lolim = -.1; hilim = .075; hilabel = "Suppression"; lolabel = "Synchronization"
    } else if (termstr=="PowerRTt") {lolim = -.03; hilim = .1; hilabel = "Suppression"; lolabel = "Synchronization"
    } else if (termstr=="PowerRTVmaxt") {lolim = -.01; hilim = .023; hilabel = "Synchronization"; lolabel = "Suppression"
    } else if (termstr=="PowerRTtMINUS1") {lolim = -.01; hilim = .023; hilabel = "Synchronization"; lolabel = "Suppression"
    } else if (termstr=="PowerRewardt") {lolim = -.01; hilim = .023; hilabel = "Synchronization"; lolabel = "Suppression"
    } else if (termstr=="PowerRTtRewardt") {lolim = -.02; hilim = .04; hilabel = "Synchronization"; lolabel = "Suppression"
    } else  {lolim = -.03; hilim = .02}
    pdf(fname, width = 10, height = 9)
    print(ggplot(edf %>% filter(estimate < 0), aes(t, Freq)) + geom_tile(aes(fill = estimate, alpha = p_level_fdr), size = .01) + 
            scale_fill_distiller(palette = "Oranges", direction = 1, name = lolabel, limits = c(lolim, 0)) + scale_x_continuous(breaks = pretty(edf$t, n = 5)) + labs(fill = lolabel) +
            new_scale_fill() +
            geom_tile(data = edf %>% filter(estimate > 0), aes(t, Freq, fill = estimate, alpha = p_level_fdr), size = .01) +
            scale_y_discrete(limits = levels(edf$Freq)) +
            scale_fill_distiller(palette = "YlGnBu", direction = -1, name = hilabel, limits = c(0, hilim)) + 
            scale_color_grey() + xlab(epoch_label) + ylab("Frequency") +
            facet_wrap( ~ node, ncol = 2) + 
            geom_vline(xintercept = 0, lty = "dashed", color = "white", size = 1) + theme_black() + 
            geom_vline(xintercept = -5, lty = "dashed", color = "white", size = 1) +
            geom_vline(xintercept = -5.3, lty = "dashed", color = "white", size = .5) +
            geom_vline(xintercept = -2.5, lty = "dotted", color = "grey", size = .5) +
            geom_text(data = edf, x = -6, y = 5,aes(label = "Response(t)"), size = 3, color = "grey", angle = 90) +
            geom_text(data = edf, x = -4.5, y = 5,aes(label = "Outcome(t)"), size = 3, color = "grey", angle = 90) +
            geom_text(data = edf, x = -0.75, y = 7.5 ,aes(label = "Clock onset (t+1)"), size = 3, color = "grey", angle = 90) +
            labs(alpha = expression(italic(p)[FDR-corrected]), fill = hilabel) + ggtitle(paste(termstr)))
    dev.off()
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

