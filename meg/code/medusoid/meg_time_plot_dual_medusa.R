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
online = T
encode = T  # main analysis analogous to Fig. 4 E-G in NComm 2020
rt_predict = F # predicts next response based on signal and behavioral variables
# random = T # whether to use data from analyses where behavioral variables have both fixed and random effects
# uncorrected = F # whether to plot uncorrected data (FDR-corrected always plotted)

repo_directory <- "~/code/clock_analysis"
rt_encode_plot_dir = "~/OneDrive/collected_letters/papers/meg/plots/rt_encode/"  
clock_encode_plot_dir = "~/OneDrive/collected_letters/papers/meg/plots/clock_encode/"  
dual_encode_plot_dir = "~/OneDrive/collected_letters/papers/meg/plots/dual_encode/"  

rt_rt_plot_dir = "~/OneDrive/collected_letters/papers/meg/plots/rt_rt/"
clock_rt_plot_dir = "~/OneDrive/collected_letters/papers/meg/plots/clock_rt/"
dual_rt_plot_dir = "~/OneDrive/collected_letters/papers/meg/plots/dual_rt/"

clock_epoch_label = "Time relative to clock onset, seconds"
rt_epoch_label = "Time relative to outcome, seconds"


sensor_map <- read.table("~/OneDrive/collected_letters/papers/meg/plots/meg_sensors_annotated.txt", header = T, colClasses = "character") %>%
  mutate(lobe=ordered(lobe, levels=c("frontal", "temporal", "parietal", "occipital")),
         dan=factor(dan, levels=c("yes", "maybe", "no"), labels=c("1yes", "2maybe", "3no"))) %>%
  #mutate(dan_grad=ordered(dan, levels=c("no", "maybe", "yes")))
  group_by(hemi, dan, lobe) %>% mutate(newlab=paste(substr(dan[1], 1,2), substr(lobe[1], 1,1), sprintf("%02d", 1:n()), sensor, sep="_")) %>% ungroup() #hemi[1], 
  #group_by(dan, hemi) %>% mutate(newlab=paste(dan[1], 1:n(), sep="_")) %>% ungroup()
# plots ----
if (encode) {  
  message("Plotting decoding results")
  setwd(rt_encode_plot_dir)
  rddf <- readRDS("meg_mixed_by_time_ddf.RDS") %>% filter(Time > -2 & Time < 2) %>%
    mutate(t  = as.numeric(Time) - 5, alignment = "rt")
  setwd(clock_encode_plot_dir)
  cddf <- readRDS("meg_mixed_by_time_clock_ddf.RDS") %>% filter(Time > -2 & Time < 2) %>%
    mutate(t  = as.numeric(Time), alignment = "clock")
  ddf <- rbind(rddf, cddf)
    # ddf <- readRDS("meg_mixed_by_time_ranefs_mult_interactions_pe_ddf.RDS")
  terms <- unique(ddf$term[ddf$effect=="fixed"])
  
  ddf <- ddf %>% mutate(p_value = as.factor(case_when(`p.value` > .05 ~ '1',
                                                      `p.value` < .05 & `p.value` > .01 ~ '2',
                                                      `p.value` < .01 & `p.value` > .001 ~ '3',
                                                      `p.value` <.001 ~ '4')),
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
  # ggplot(ddf %>% filter(term==terms[1]), aes(t, newlab)) + geom_tile(aes(fill = estimate, alpha = p_value)) + 
    # facet_grid(lobe  ~ hemi, scales = "free")
  setwd(dual_encode_plot_dir)
  for (fe in terms) {
    edf <- ddf %>% filter(term == paste(fe)) 
    termstr <- str_replace_all(fe, "[^[:alnum:]]", "_")
    fname = paste("meg_time_dual_uncorrected_n63", termstr, ".pdf", sep = "")
    message(fname)
    pdf(fname, width = 30, height = 16)
    print(ggplot(edf, aes(t, newlab)) + geom_tile(aes(fill = estimate, alpha = p_value)) + 
            facet_wrap(~lobe * hemi, nrow = 2, scales = "free") +
            geom_vline(xintercept = 0, lty = "dashed", color = "black", size = 2) +
            geom_vline(xintercept = -5, lty = "dashed", color = "red", size = 2) +
            geom_vline(xintercept = -5.3, lty = "dashed", color = "red", size = 1) +
            scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(clock_epoch_label) + ylab("Sensor") +
            labs(alpha = expression(italic(p)[uncorrected])) + ggtitle(paste(termstr)))
    dev.off()
    fname = paste("meg_time_dual_FDR_n63", termstr, ".pdf", sep = "")
    pdf(fname, width = 30, height = 16)
    print(ggplot(edf, aes(t, newlab)) + geom_tile(aes(fill = estimate, alpha = p_level_fdr)) + 
            facet_wrap(~lobe * hemi, nrow = 2, scales = "free") +
            geom_vline(xintercept = 0, lty = "dashed", color = "black", size = 2) +
            geom_vline(xintercept = -5, lty = "dashed", color = "red", size = 2) +
            geom_vline(xintercept = -5.3, lty = "dashed", color = "red", size = 1) +
            scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(clock_epoch_label) + ylab("Sensor") +
            labs(alpha = expression(italic(p)[uncorrected])) + ggtitle(paste(termstr)))
    dev.off()
  }
} 
# system("for i in *.pdf; do sips -s format jpeg $i --out $i.jpeg; done")

if(rt_predict) {
  message("Plotting decoding results")
  setwd(rt_rt_plot_dir)
  rrdf <- readRDS("meg_mixed_by_time_rdf.RDS") %>% filter(Time > -2 & Time < 2) %>%
    mutate(t  = as.numeric(Time) - 5, alignment = "rt")
  setwd(clock_rt_plot_dir)
  crdf <- readRDS("meg_mixed_by_time_rdf.RDS") %>% filter(Time > -2 & Time < 2) %>%
    mutate(t  = as.numeric(Time), alignment = "clock")
  rdf <- rbind(rrdf, crdf)
  # ddf <- readRDS("meg_mixed_by_time_ranefs_mult_interactions_pe_ddf.RDS")
  terms <- unique(ddf$term[rdf$effect=="fixed"])
  terms <- terms[grepl("signal", terms)]
  rdf <- rdf %>% mutate(p_value = as.factor(case_when(`p.value` > .05 ~ '1',
                                                      `p.value` < .05 & `p.value` > .01 ~ '2',
                                                      `p.value` < .01 & `p.value` > .001 ~ '3',
                                                      `p.value` <.001 ~ '4')),
                        t  = t,
                        sensor = as.character(sensor))
  rdf$p_value <- factor(rdf$p_value, labels = c("NS", "p < .05", "p < .01", "p < .001"))
  terms <- unique(rdf$term[rdf$effect=="fixed"])
  terms <- terms[grepl("signal", terms)]
  # FDR labeling ----
  rdf <- rdf  %>% filter(t>-1) %>% mutate(p_fdr = padj_fdr_term,
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
  rdf <- rdf %>% filter(!grepl("1$", sensor))
  # add sensor labels
  rdf <- rdf %>% merge(sensor_map)
  # ggplot(rdf %>% filter(term==terms[1]), aes(t, newlab)) + geom_tile(aes(fill = estimate, alpha = p_value)) + 
  facet_grid(lobe  ~ hemi, scales = "free")
  for (fe in terms) {
    edf <- rdf %>% filter(term == paste(fe)) 
    termstr <- str_replace_all(fe, "[^[:alnum:]]", "_")
    fname = paste("meg_time_clock_uncorrected_n63", termstr, ".pdf", sep = "")
    pdf(fname, width = 20, height = 16)
    print(ggplot(edf, aes(t, newlab)) + geom_tile(aes(fill = estimate, alpha = p_value)) + 
            facet_wrap(~lobe * hemi, nrow = 2, scales = "free") +
            geom_vline(xintercept = 0, lty = "dashed", color = "black", size = 2) +
            scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) + ylab("Sensor") +
            labs(alpha = expression(italic(p)[uncorrected])) + ggtitle(paste(termstr)))
    dev.off()
    fname = paste("meg_time_clock_FDR_n63", termstr, ".pdf", sep = "")
    pdf(fname, width = 20, height = 16)
    print(ggplot(edf, aes(t, newlab)) + geom_tile(aes(fill = estimate, alpha = p_level_fdr)) + 
            facet_wrap(~lobe * hemi, nrow = 2, scales = "free") +
            geom_vline(xintercept = 0, lty = "dashed", color = "black", size = 2) +
            scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) + ylab("Sensor") +
            labs(alpha = expression(italic(p)[uncorrected])) + ggtitle(paste(termstr)))
    dev.off()
  }
  
  # convert to PNG for Word
  # system("for i in *63*.pdf; do sips -s format png $i --out $i.png; done")
}
