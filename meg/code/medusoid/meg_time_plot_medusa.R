# plots MLM coefficients by sensor across timepoints and frequencies
# first run meg_time_freq_load_analyze_medusa.R

library(tidyverse)
library(lme4)
library(ggpubr)
library(car)
library(viridis)
# library(psych)

# what to run
# you can run all options at once
online = F
encode = F  # main analysis analogous to Fig. 4 E-G in NComm 2020
rt_predict = T # predicts next response based on signal and behavioral variables
# random = T # whether to use data from analyses where behavioral variables have both fixed and random effects
# uncorrected = F # whether to plot uncorrected data (FDR-corrected always plotted)

repo_directory <- "~/code/clock_analysis"
if (online) {encode_plot_dir = "~/OneDrive/collected_letters/papers/meg/plots/clock_encode/"} else {
  encode_plot_dir = "~/OneDrive/collected_letters/papers/meg/plots/rt_encode/"  
}

# encode_plot_dir = "~/OneDrive/collected_letters/papers/meg/plots/clock_encode/"
# rt_plot_dir = "~/OneDrive/collected_letters/papers/meg/plots/rt_rt/"
if (online) {rt_plot_dir = "~/OneDrive/collected_letters/papers/meg/plots/clock_rt/"} else {
rt_plot_dir = "~/OneDrive/collected_letters/papers/meg/plots/rt_rt/"}
if (online) {
  epoch_label = "Time relative to clock onset, seconds"} else {epoch_label = "Time relative to outcome, seconds"}


sensor_map <- read.table("~/OneDrive/collected_letters/papers/meg/plots/meg_sensors_annotated.txt", header = T, colClasses = "character") %>%
  mutate(lobe=ordered(lobe, levels=c("frontal", "temporal", "parietal", "occipital")),
         dan=factor(dan, levels=c("yes", "maybe", "no"), labels=c("1yes", "2maybe", "3no"))) %>%
  #mutate(dan_grad=ordered(dan, levels=c("no", "maybe", "yes")))
  group_by(hemi, dan, lobe) %>% mutate(newlab=paste(substr(dan[1], 1,2), substr(lobe[1], 1,1), sprintf("%02d", 1:n()), sensor, sep="_")) %>% ungroup() #hemi[1], 
  #group_by(dan, hemi) %>% mutate(newlab=paste(dan[1], 1:n(), sep="_")) %>% ungroup()
# plots ----
if (encode) {  
  message("Plotting decoding results")
  setwd(encode_plot_dir)
  epoch_label = "Time relative to clock onset, seconds"
# basic models on n = 63
  ddf <- readRDS("meg_mixed_by_time_clock_ddf.RDS")
    # ddf <- readRDS("meg_mixed_by_time_ranefs_mult_interactions_pe_ddf.RDS")
  terms <- unique(ddf$term[ddf$effect=="fixed"])
  
  ddf <- ddf %>% mutate(p_value = as.factor(case_when(`p.value` > .05 ~ '1',
                                                      `p.value` < .05 & `p.value` > .01 ~ '2',
                                                      `p.value` < .01 & `p.value` > .001 ~ '3',
                                                      `p.value` <.001 ~ '4')),
                        t  = as.numeric(Time),
                        sensor = as.character(sensor))
  ddf$p_value <- factor(ddf$p_value, labels = c("NS", "p < .05", "p < .01", "p < .001"))
  # terms <- unique(ddf$term[ddf$effect=="fixed"])
  # terms <- terms[grepl(!"_f", terms)]
  # # FDR labeling ----
  ddf <- ddf  %>% mutate(p_fdr = padj_fdr_term,
                                            p_level_fdr = as.factor(case_when(
                                              # p_fdr > .1 ~ '0',
                                              # p_fdr < .1 & p_fdr > .05 ~ '1',
                                              p_fdr > .05 ~ '1',
                                              p_fdr < .05 & p_fdr > .01 ~ '2',
                                              p_fdr < .01 & p_fdr > .001 ~ '3',
                                              p_fdr <.001 ~ '4'))
  ) %>% ungroup() #%>% mutate(side = substr(as.character(label), nchar(as.character(label)), nchar(as.character(label))),
  #          region = substr(as.character(label), 1, nchar(as.character(label))-2))
  ddf$p_level_fdr <- factor(ddf$p_level_fdr, levels = c('1', '2', '3', '4'), labels = c("NS","p < .05", "p < .01", "p < .001"))
  ddf$`p, FDR-corrected` = ddf$p_level_fdr
  
  # drop the *1 sensors
  ddf <- ddf %>% filter(!grepl("1$", sensor))
  # add sensor labels
  ddf <- ddf %>% merge(sensor_map)
  # ggplot(ddf %>% filter(term==terms[1]), aes(t, newlab)) + geom_tile(aes(fill = abs(estimate), alpha = p_value)) + 
    facet_grid(lobe  ~ hemi, scales = "free")
  for (fe in terms) {
    edf <- ddf %>% filter(term == paste(fe)) 
    termstr <- str_replace_all(fe, "[^[:alnum:]]", "_")
    fname = paste("meg_time_clock_uncorrected_n63_test", termstr, ".pdf", sep = "")
    message(fname)
    pdf(fname, width = 30, height = 16)
    print(ggplot(edf, aes(t, newlab)) + geom_tile(aes(fill = abs(estimate), alpha = p_value)) + 
            facet_wrap(~lobe * hemi, nrow = 2, scales = "free") +
            geom_vline(xintercept = 0, lty = "dashed", color = "black", size = 2) +
            scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) + ylab("Sensor") +
            labs(alpha = expression(italic(p)[uncorrected])) + ggtitle(paste(termstr)))
    dev.off()
    fname = paste("meg_time_clock_FDR_n63_test", termstr, ".pdf", sep = "")
    pdf(fname, width = 30, height = 16)
    print(ggplot(edf, aes(t, newlab)) + geom_tile(aes(fill = abs(estimate), alpha = p_level_fdr)) + 
            facet_wrap(~lobe * hemi, nrow = 2, scales = "free") +
            geom_vline(xintercept = 0, lty = "dashed", color = "black", size = 2) +
            scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) + ylab("Sensor") +
            labs(alpha = expression(italic(p)[uncorrected])) + ggtitle(paste(termstr)))
    dev.off()
  }
} 
# system("for i in *_time_clock_*.pdf; do sips -s format png $i --out $i.png; done")

if(rt_predict) {
  # plots ----
  setwd(rt_plot_dir)
  # rdf <- readRDS("meg_mixed_by_time_rdf.RDS")
  rdf <- as_tibble(readRDS("meg_mixed_by_time_rdf_combined_feedback.RDS"))
  terms <- unique(rdf$term[rdf$effect=="fixed"])
  terms <- terms[grepl("signal", terms)]
  rdf <- rdf %>% mutate(p_value = as.factor(case_when(`p.value` > .05 ~ '1',
                                                      `p.value` < .05 & `p.value` > .01 ~ '2',
                                                      `p.value` < .01 & `p.value` > .001 ~ '3',
                                                      `p.value` <.001 ~ '4')),
                        t  = as.numeric(Time),
                        lobe = stringr::str_extract(.filename, "[^_]+"),
                        hemi = substr(stringr::str_extract(.filename, "_l_|_r_"), 2,2))
                        # sensor = as.character(sensor))
  rdf$p_value <- factor(rdf$p_value, labels = c("NS", "p < .05", "p < .01", "p < .001"))
  terms <- unique(rdf$term[rdf$effect=="fixed"])
  terms <- terms[grepl("signal", terms)]
  # FDR labeling ----
  rdf <- rdf  %>% filter(t>-4) %>% mutate(p_fdr = padj_BY_term,
                         p_level_fdr = as.factor(case_when(
                           # p_fdr > .1 ~ '0',
                           # p_fdr < .1 & p_fdr > .05 ~ '1',
                           p_fdr > .05 ~ '1',
                           p_fdr < .05 & p_fdr > .01 ~ '2',
                           p_fdr < .01 & p_fdr > .001 ~ '3',
                           p_fdr <.001 ~ '4'))
  ) %>% ungroup() #%>% mutate(side = substr(as.character(label), nchar(as.character(label)), nchar(as.character(label))),
  #          region = substr(as.character(label), 1, nchar(as.character(label))-2))
  rdf$p_level_fdr <- factor(rdf$p_level_fdr, levels = c('1', '2', '3', '4'), labels = c("NS","p < .05", "p < .01", "p < .001"))
  rdf$`p, FDR-corrected` = rdf$p_level_fdr
  
  # drop the *1 sensors
  # rdf <- rdf %>% filter(!grepl("1$", sensor))
  # add sensor labels
  # rdf <- rdf %>% merge(sensor_map)
  # ggplot(rdf %>% filter(term==terms[1]), aes(t, newlab)) + geom_tile(aes(fill = abs(estimate), alpha = p_value)) + 
  # facet_grid(lobe  ~ hemi, scales = "free")
  for (fe in terms) {
    edf <- rdf %>% filter(term == paste(fe)) 
    termstr <- str_replace_all(fe, "[^[:alnum:]]", "_")
    fname = paste("meg_time_combined_feedback_uncorrected_n63", termstr, ".pdf", sep = "")
    pdf(fname, width = 10, height = 4)
    print(ggplot(edf, aes(t, lobe)) + geom_tile(aes(fill = estimate, alpha = p_value)) + 
            facet_wrap(~hemi, nrow = 1, scales = "free") +
            geom_vline(xintercept = 0, lty = "dashed", color = "black", size = 2) +
            scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) + ylab("Lobe") +
            labs(alpha = expression(italic(p)[uncorrected])) + ggtitle(paste(termstr)))
    dev.off()
    fname = paste("meg_time_combined_feedback_FDR_n63", termstr, ".pdf", sep = "")
    pdf(fname, width = 10, height = 4)
    print(ggplot(edf, aes(t, lobe)) + geom_tile(aes(fill = estimate, alpha = p_level_fdr)) + 
            facet_wrap(~hemi, nrow = 1, scales = "free") +
            geom_vline(xintercept = 0, lty = "dashed", color = "black", size = 2) +
            scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) + ylab("Lobe") +
            labs(alpha = expression(italic(p)[uncorrected])) + ggtitle(paste(termstr)))
    dev.off()
  }
  
  # convert to PNG for Word
  # system("for i in *63*.pdf; do sips -s format png $i --out $i.png; done")
}
