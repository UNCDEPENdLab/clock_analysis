# plots MLM coefficients by sensor across timepoints and frequencies
# first run meg_time_freq_load_analyze_medusa.R

library(tidyverse)
library(lme4)
library(ggpubr)
library(car)
library(viridis)
# library(psych)
repo_directory <- "~/code/clock_analysis"
data_dir <- "~/Box/SCEPTIC_fMRI/MEG_TimeFreq/"
rt_encode_plot_dir = "~/OneDrive/collected_letters/papers/meg/plots/rt_encode/"  
clock_encode_plot_dir = "~/OneDrive/collected_letters/papers/meg/plots/clock_encode/"  
dual_encode_plot_dir = "~/OneDrive/collected_letters/papers/meg/plots/dual_encode/"  

rt_rt_plot_dir = "~/OneDrive/collected_letters/papers/meg/plots/rt_rt/"
clock_rt_plot_dir = "~/OneDrive/collected_letters/papers/meg/plots/clock_rt/"
dual_rt_plot_dir = "~/OneDrive/collected_letters/papers/meg/plots/dual_rt/"

clock_epoch_label = "Time relative to clock onset, seconds"
rt_epoch_label = "Time relative to outcome, seconds"
uncorrected_input = F
online = F
encode = T

sensor_map <- read.table("~/OneDrive/collected_letters/papers/meg/plots/meg_sensors_annotated.txt", header = T, colClasses = "character") %>%
  mutate(lobe=ordered(lobe, levels=c("frontal", "temporal", "parietal", "occipital")),
         dan=factor(dan, levels=c("yes", "maybe", "no"), labels=c("1yes", "2maybe", "3no"))) %>%
  #mutate(dan_grad=ordered(dan, levels=c("no", "maybe", "yes")))
  group_by(hemi, dan, lobe) %>% mutate(sensor_lab=paste(substr(dan[1], 1,2), substr(lobe[1], 1,1), sprintf("%02d", 1:n()), sensor, sep="_"),
                                       short_lab=paste(substr(dan[1], 1,2), substr(lobe[1], 1,1), sprintf("%02d", 1:n()), sep="_")) %>% ungroup() #hemi[1], 
#group_by(dan, hemi) %>% mutate(newlab=paste(dan[1], 1:n(), sep="_")) %>% ungroup()

setwd(data_dir)
# plots ----
if (encode) {  
  message("Plotting encoding results")
  # epoch_label = "Time relative to outcome, seconds"
  
  if (online) {
    epoch_label = "Time relative to clock onset, seconds"
    # add clock-aligned file
  } else {
    epoch_label = "Time relative to feedback, seconds"
    ddf <- as_tibble(readRDS("meg_mixed_by_tf_allDAN_rtencode_ddf.RDS"))}
  
  terms <- unique(ddf$term[ddf$effect=="fixed"])
  ddf$sensor <- readr::parse_number(ddf$.filename, trim_ws = F)
  ddf$sensor <- stringr::str_pad(ddf$sensor, 4, "0",side = "left")
  ddf <- ddf %>% mutate(p_value = as.factor(case_when(`p.value` > .05 ~ '1',
                                                      `p.value` < .05 & `p.value` > .01 ~ '2',
                                                      `p.value` < .01 & `p.value` > .001 ~ '3',
                                                      `p.value` <.001 ~ '4')),
                        t  = as.numeric(Time) - .3,
                        sensor = as.character(sensor))
  ddf$p_value <- factor(ddf$p_value, labels = c("NS", "p < .05", "p < .01", "p < .001"))
  terms <- unique(ddf$term[ddf$effect=="fixed"])
  
  if (uncorrected_input) {
    # FDR correction ----
    ddf <- ddf  %>% filter(t>-4) %>% group_by(term) %>% mutate(p_fdr = p.adjust(p.value, method = 'fdr'),
                                                               p_level_fdr = as.factor(case_when(
                                                                 # p_fdr > .1 ~ '0',
                                                                 # p_fdr < .1 & p_fdr > .05 ~ '1',
                                                                 p_fdr > .05 ~ '1',
                                                                 p_fdr < .05 & p_fdr > .01 ~ '2',
                                                                 p_fdr < .01 & p_fdr > .001 ~ '3',
                                                                 p_fdr <.001 ~ '4'))
    ) %>% ungroup() #%>% mutate(side = substr(as.character(label), nchar(as.character(label)), nchar(as.character(label))),
    
  } 
  # FDR labeling ----
  ddf <- ddf  %>% mutate(
    # p_fdr = padj_fdr_term,
                         p_level_fdr = as.factor(case_when(
                           # p_fdr > .1 ~ '0',
                           # p_fdr < .1 & p_fdr > .05 ~ '1',
                           padj_fdr_term > .05 ~ '1',
                           padj_fdr_term < .05 & padj_fdr_term > .01 ~ '2',
                           padj_fdr_term < .01 & padj_fdr_term > .001 ~ '3',
                           padj_fdr_term <.001 ~ '4'))
  ) %>% ungroup() #%>% mutate(side = substr(as.character(label), nchar(as.character(label)), nchar(as.character(label))),
           # region = substr(as.character(label), 1, nchar(as.character(label))-2))
  ddf$p_level_fdr <- factor(ddf$p_level_fdr, levels = c('1', '2', '3', '4'), labels = c("NS","p < .05", "p < .01", "p < .001"))
  ddf$`p, FDR-corrected` = ddf$p_level_fdr
  
  levels(ddf$Freq) <- substr(unique(ddf$Freq), 3,6)
  ddf$freq <- fct_rev(ddf$Freq)
  # drop magnetometers (if any)
  ddf <- ddf %>% filter(!grepl("1$", sensor))
  # add sensor labels
  ddf <- ddf %>% merge(sensor_map)
  
  setwd(rt_encode_plot_dir)
  for (fe in terms) {
    edf <- ddf %>% filter(term == paste(fe) & (ddf$lobe == "frontal" | ddf$lobe == "parietal") & dan == "1yes") 
    termstr <- str_replace_all(fe, "[^[:alnum:]]", "_")
    fname = paste("meg_tf_rt_all_dan_uncorrected_", termstr, ".pdf", sep = "")
    pdf(fname, width = 45, height = 45)
    print(ggplot(edf, aes(t, freq)) + geom_tile(aes(fill = estimate, alpha = p_value), size = .01) +
            geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + #scale_y_discrete(limits = rev) +
            geom_vline(xintercept = -.3, lty = "dashed", color = "#FF0000", size = 1) + #scale_y_discrete(limits = rev) +
            scale_fill_viridis(option = "plasma") +  xlab(epoch_label) + ylab("Frequency") +
            facet_grid(lobe*short_lab ~ hemi) + geom_text(data = edf, x = 1.5,y = 2,aes(label = sensor), size = 8, color = "white") +
            labs(alpha = expression(italic(p)[uncorrected])) + ggtitle(paste(termstr)) + theme_dark())
    dev.off()

    fname = paste("meg_tf_rt_all_dan_FDR_", termstr, ".pdf", sep = "")
    pdf(fname, width = 45, height = 45)
    print(ggplot(edf, aes(t, freq)) + geom_tile(aes(fill = estimate, alpha = p_level_fdr), size = .01) + 
            geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + #scale_y_discrete(limits = rev) +
            geom_vline(xintercept = -.3, lty = "dashed", color = "#FF0000", size = 1) + #scale_y_discrete(limits = rev) +
            scale_fill_viridis(option = "plasma") +  xlab(epoch_label) + ylab("Frequency") + 
            facet_grid(lobe*short_lab ~ hemi) + geom_text(data = edf, x = 1.5,y = 2,aes(label = sensor), size = 8, color = "white") +
            labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + theme_dark())
    dev.off()
  }
} 
system("for i in *scaled*.pdf; do sips -s format png $i --out $i.png; done")

if(rt_predict) {
  # plots ----
  setwd('~/OneDrive/collected_letters/papers/meg/plots/rt_rt')
  epoch_label = "Time relative to outcome, seconds"
  if (random) {rt_results_fname = "meg_freq_medusa_rt_predict_output_random.Rdata"} else {
    rt_results_fname = "meg_mixed_by_tf_rdf.RDS"  
  }
  rdf <- readRDS(rt_results_fname)
  terms <- unique(rdf$term[rdf$effect=="fixed"])
  terms <- terms[grepl("(pow)",terms)]
  rdf <- rdf %>% mutate(p_value = as.factor(case_when(`p.value` > .05 ~ '1',
                                                      `p.value` < .05 & `p.value` > .01 ~ '2',
                                                      `p.value` < .01 & `p.value` > .001 ~ '3',
                                                      `p.value` <.001 ~ '4')),
                        t  = as.numeric(Time),
                        sensor = as.character(sensor))
  rdf$p_value <- factor(rdf$p_value, labels = c("NS", "p < .05", "p < .01", "p < .001"))
  terms <- unique(rdf$term[rdf$effect=="fixed"])
  # FDR labeling ----
  rdf <- rdf  %>% filter(t>-2) %>% mutate(p_fdr = padj_fdr_term,
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
  
  levels(rdf$Freq) <- substr(unique(rdf$Freq), 3,6)
  rdf$freq <- fct_rev(rdf$Freq)
  library(ggnewscale)

  for (fe in terms) {
    edf <- rdf %>% filter(term == paste(fe)) 
    termstr <- str_replace_all(fe, "[^[:alnum:]]", "_")

    if (uncorrected) {
    if (random) {fname = paste("meg_tf_dan_uncorrected_random_", termstr, ".pdf", sep = "")  
    } else {fname = paste("meg_tf_dan_uncorrected_", termstr, ".pdf", sep = "")}
    
    pdf(fname, width = 20, height = 8)
    print(ggplot(edf %>% filter(estimate < 0), aes(t, freq)) + geom_tile(aes(fill = estimate, alpha = p_value), size = .01) + 
            scale_fill_distiller(palette = "OrRd", direction = 1, limits = c(-.1, 0)) + scale_color_grey(limits = rev(levels(rdf$p_level_fdr))) + xlab(epoch_label) + ylab("Frequency") + facet_wrap(~sensor) +
            labs(alpha = expression(italic(p)[uncorrected])) + ggtitle(paste(termstr)) +
            new_scale_fill() +
            geom_tile(data = edf %>% filter(estimate > 0), aes(t, freq, fill = estimate, alpha = p_value), size = .01) + theme_dark() +
            geom_vline(xintercept = 0, lty = "dashed", color = "black", size = 2) + theme_bw() + 
            scale_fill_distiller(palette = "GnBu", direction = -1, limits = c(0, .15))+ scale_color_grey() + xlab(epoch_label) + ylab("Frequency") + facet_wrap(~sensor, ncol = 3) +
            labs(alpha = expression(italic(p)[uncorrected])) + ggtitle(paste(termstr))) + scale_y_discrete(limits = rev(levels(rdf$p_level_fdr)))
    dev.off()
    }
    if (random) {fname = paste("meg_tf_dan_FDR_random_", termstr, ".pdf", sep = "")  
    } else {fname = paste("meg_tf_dan_FDR_", termstr, ".pdf", sep = "")}
    pdf(fname, width = 20, height = 8)
    print(ggplot(edf %>% filter(estimate < 0), aes(t, freq)) + geom_tile(aes(fill = estimate, alpha = p_level_fdr), size = .01) + 
            scale_fill_distiller(palette = "OrRd", direction = 1, name = "Exploration", limits = c(-.1, 0)) + scale_x_continuous(breaks = pretty(edf$t, n = 20)) + labs(fill = "Exploration") +
            new_scale_fill() +
            geom_tile(data = edf %>% filter(estimate > 0), aes(t, freq, fill = estimate, alpha = p_level_fdr), size = .01) +
            geom_vline(xintercept = 0, lty = "dashed", color = "black", size = 2) + theme_bw() + 
            scale_fill_distiller(palette = "GnBu", direction = -1, name = "Exploitation", limits = c(0, .15))+ scale_color_grey() + xlab(epoch_label) + ylab("Frequency") + facet_wrap(~sensor, ncol = 3) +
            labs(alpha = expression(italic(p)[FDR-corrected])) + ggtitle(paste(termstr))) + labs(fill = "Exploitation")
    dev.off()
  }
  
  # convert to PNG for Word
  system("for i in *meg_tf*.pdf; do sips -s format png $i --out $i.png; done")
  
  # TESTING THE "UNIFIED" VERSION -- currently not functional
  unified = F
  if (unified) {
    for (fe in terms) {
      edf <- rdf %>% filter(term == paste(fe)) 
      termstr <- str_replace_all(fe, "[^[:alnum:]]", "_")
      if (uncorrected) {
        if (random) {fname = paste("meg_tf_dan_uncorrected_unified_random_", termstr, ".pdf", sep = "")  
        } else {fname = paste("meg_tf_dan_uncorrected_", termstr, ".pdf", sep = "")}
        
        pdf(fname, width = 10, height = 4)
        print(ggplot(edf %>% filter(estimate < 0), aes(t, freq)) + geom_tile(aes(fill = estimate, alpha = p_value), size = .01) + 
                scale_fill_distiller(palette = "OrRd", direction = 1, limits = c(-.2, 0)) + scale_color_grey(limits = rev(levels(rdf$p_level_fdr))) + xlab(epoch_label) + ylab("Frequency") + 
                labs(alpha = expression(italic(p)[uncorrected])) + ggtitle(paste(termstr)) +
                new_scale_fill() +
                geom_tile(data = edf %>% filter(estimate > 0), aes(t, freq, fill = estimate, alpha = p_value), size = .01) + theme_dark() +
                geom_vline(xintercept = 0, lty = "dashed", color = "black", size = 2) + theme_bw() + 
                scale_fill_distiller(palette = "GnBu", direction = -1, limits = c(0, .2))+ scale_color_grey() + xlab(epoch_label) + ylab("Frequency") + 
                labs(alpha = expression(italic(p)[uncorrected])) + ggtitle(paste(termstr))) + scale_y_discrete(limits = rev(levels(rdf$p_level_fdr)))
        dev.off()
      }
      if (random) {fname = paste("meg_tf_dan_FDR_unified_random_", termstr, ".pdf", sep = "")  
      } else {fname = paste("meg_tf_dan_FDR_", termstr, ".pdf", sep = "")}
      pdf(fname, width = 10, height = 4)
      print(ggplot(edf %>% filter(estimate < 0), aes(t, freq)) + geom_tile(aes(fill = estimate, alpha = p_level_fdr), size = .01) + 
              scale_fill_distiller(palette = "OrRd", direction = 1, name = "Exploration", limits = c(-.15, 0)) + scale_x_continuous(breaks = pretty(edf$t, n = 20)) + labs(fill = "Exploration") +
              new_scale_fill() +
              geom_tile(data = edf %>% filter(estimate > 0), aes(t, freq, fill = estimate, alpha = p_level_fdr), size = .01) +
              geom_vline(xintercept = 0, lty = "dashed", color = "black", size = 2) + theme_bw() + 
              scale_fill_distiller(palette = "GnBu", direction = -1, name = "Exploitation")+ scale_color_grey() + xlab(epoch_label) + ylab("Frequency") + 
              labs(alpha = expression(italic(p)[FDR-corrected])) + ggtitle(paste(termstr))) + labs(fill = "Exploration")
      
      dev.off()
    }
  }
}
