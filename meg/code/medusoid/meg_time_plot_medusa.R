# plots MLM coefficients by sensor across timepoints and frequencies
# first run meg_time_freq_load_analyze_medusa.R

library(tidyverse)
library(lme4)
library(ggpubr)
library(car)
library(viridis)
# library(psych)
repo_directory <- "~/code/clock_analysis"
encode_plot_dir = "~/OneDrive/collected_letters/papers/meg/plots/rt_decode/"
rt_plot_dir = "~/OneDrive/collected_letters/papers/meg/plots/rt_rt//"


# what to run
# you can run all options at once
encode = T  # main analysis analogous to Fig. 4 E-G in NComm 2020
rt_predict = F # predicts next response based on signal and behavioral variables
random = T # whether to use data from analyses where behavioral variables have both fixed and random effects
uncorrected = F # whether to plot uncorrected data (FDR-corrected always plotted)

# plots ----
if (encode) {  
  message("Plotting decoding results")
  setwd(encode_plot_dir)
  epoch_label = "Time relative to outcome, seconds"
  ddf <- readRDS("meg_mixed_by_time_ranefs_mult_interactions_ddf.RDS")
  terms <- unique(ddf$term[ddf$effect=="fixed"])
  
  ddf <- ddf %>% mutate(p_value = as.factor(case_when(`p.value` > .05 ~ '1',
                                                      `p.value` < .05 & `p.value` > .01 ~ '2',
                                                      `p.value` < .01 & `p.value` > .001 ~ '3',
                                                      `p.value` <.001 ~ '4')),
                        t  = as.numeric(Time),
                        sensor = as.character(sensor))
  ddf$p_value <- factor(ddf$p_value, labels = c("NS", "p < .05", "p < .01", "p < .001"))
  terms <- unique(ddf$term[ddf$effect=="fixed"])
  # FDR labeling ----
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
  
  for (fe in terms) {
    edf <- ddf %>% filter(term == paste(fe)) 
    termstr <- str_replace_all(fe, "[^[:alnum:]]", "_")
    fname = paste("meg_time_uncorrected_n63", termstr, ".pdf", sep = "")
    pdf(fname, width = 16, height = 30)
    print(ggplot(edf, aes(t, sensor)) + geom_tile(aes(fill = abs(estimate), alpha = p_value)) + 
            geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
            scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) + ylab("Sensor") +
            labs(alpha = expression(italic(p)[uncorrected])) + ggtitle(paste(termstr)))
    dev.off()
    fname = paste("meg_time_FDR_n63", termstr, ".pdf", sep = "")
    pdf(fname, width = 16, height = 30)
    print(ggplot(edf, aes(t, sensor)) + geom_tile(aes(fill = abs(estimate), alpha = p_level_fdr)) + 
            geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
            scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) + ylab("Sensor") +
            labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)))
    dev.off()
  }
} 
system("for i in *n63*.pdf; do sips -s format png $i --out $i.png; done")

if(rt_predict) {
  # plots ----
  setwd('~/OneDrive/collected_letters/papers/meg/plots/rt_rt')
  epoch_label = "Time relative to outcome, seconds"
  rdf <- readRDS("meg_mixed_by_time_rdf.RDS")
  terms <- unique(rdf$term[rdf$effect=="fixed"])
  terms <- terms[grepl("signal", terms)]
  rdf <- rdf %>% mutate(p_value = as.factor(case_when(`p.value` > .05 ~ '1',
                                                      `p.value` < .05 & `p.value` > .01 ~ '2',
                                                      `p.value` < .01 & `p.value` > .001 ~ '3',
                                                      `p.value` <.001 ~ '4')),
                        t  = as.numeric(Time),
                        sensor = as.character(sensor))
  rdf$p_value <- factor(rdf$p_value, labels = c("NS", "p < .05", "p < .01", "p < .001"))
  terms <- unique(rdf$term[rdf$effect=="fixed"])
  # FDR labeling ----
  rdf <- rdf  %>% mutate(p_fdr = padj_fdr_term,
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
  
  for (fe in terms) {
    edf <- rdf %>% filter(term == paste(fe)) 
    termstr <- str_replace_all(fe, "[^[:alnum:]]", "_")
    fname = paste("meg_time_uncorrected_n63", termstr, ".pdf", sep = "")
    pdf(fname, width = 16, height = 30)
    print(ggplot(edf, aes(t, sensor)) + geom_tile(aes(fill = abs(estimate), alpha = p_value)) + 
            geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
            scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) + ylab("Sensor") +
            labs(alpha = expression(italic(p)[uncorrected])) + ggtitle(paste(termstr)))
    dev.off()
    fname = paste("meg_time_FDR_n63", termstr, ".pdf", sep = "")
    pdf(fname, width = 16, height = 30)
    print(ggplot(edf, aes(t, sensor)) + geom_tile(aes(fill = abs(estimate), alpha = p_level_fdr)) + 
            geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
            scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) + ylab("Sensor") +
            labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)))
    dev.off()
  }
  
  # convert to PNG for Word
  system("for i in *63*.pdf; do sips -s format png $i --out $i.png; done")
}
