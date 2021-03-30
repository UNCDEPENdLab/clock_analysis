# plots MLM coefficients by sensor across timepoints and frequencies
# first run meg_time_freq_load_analyze_medusa.R

library(tidyverse)
library(lme4)
library(ggpubr)
library(car)
library(viridis)
# library(psych)
repo_directory <- "~/code/clock_analysis"

# what to run
# you can run all options at once
decode = T  # main analysis analogous to Fig. 4 E-G in NComm 2020
rt_predict = T # predicts next response based on signal and behavioral variables

# plots ----
if (decode) {  
  message("Plotting decoding results")
  setwd('~/OneDrive/collected_letters/papers/meg/plots/rt_decode')
  epoch_label = "Time relative to outcome, seconds"
  decode_results_fname = "meg_freq_medusa_rt_predict_output_all.Rdata"
  load(decode_results_fname)
  terms <- unique(ddf$term[ddf$effect=="fixed"])
  for (fe in terms) {
    edf <- ddf %>% filter(term == paste(fe)) 
    termstr <- str_replace_all(fe, "[^[:alnum:]]", "_")
    fname = paste("meg_tf_dan_uncorrected_", termstr, ".pdf", sep = "")
    pdf(fname, width = 30, height = 30)
    print(ggplot(edf, aes(t, freq)) + geom_tile(aes(fill = estimate, alpha = p_value)) + 
            geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
            scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) + ylab("Frequency") + facet_wrap(~sensor) +
            labs(alpha = expression(italic(p)[uncorrected])) + ggtitle(paste(termstr)))
    dev.off()
    fname = paste("meg_tf_dan_FDR_", termstr, ".pdf", sep = "")
    pdf(fname, width = 30, height = 30)
    print(ggplot(edf, aes(t, freq)) + geom_tile(aes(fill = estimate, alpha = p_level_fdr)) + 
            geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
            scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) + ylab("Frequency") + facet_wrap(~sensor) +
            labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)))
    dev.off()
  }
} 

if(rt) {
  # plots ----
  setwd('~/OneDrive/collected_letters/papers/meg/plots/rt_rt')
  epoch_label = "Time relative to outcome, seconds"
  rt_results_fname = "meg_freq_medusa_rt_predict_output_all.Rdata"
  decode_results_fname = "meg_freq_medusa_rt_predict_output_all.Rdata"
  load(decode_results_fname)
  terms <- unique(rdf$term[rdf$effect=="fixed"])
  terms <- terms[grepl("(h)",terms)]
  for (fe in terms) {
    edf <- rdf %>% filter(term == paste(fe)) 
    termstr <- str_replace_all(fe, "[^[:alnum:]]", "_")
    fname = paste("meg_tf_dan_uncorrected_", termstr, ".pdf", sep = "")
    pdf(fname, width = 30, height = 30)
    print(ggplot(edf, aes(t, freq)) + geom_tile(aes(fill = estimate, alpha = p_value)) + 
            geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
            scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) + ylab("Frequency") + facet_wrap(~sensor) +
            labs(alpha = expression(italic(p)[uncorrected])) + ggtitle(paste(termstr)))
    dev.off()
    fname = paste("meg_tf_dan_FDR_", termstr, ".pdf", sep = "")
    pdf(fname, width = 30, height = 30)
    print(ggplot(edf, aes(t, freq)) + geom_tile(aes(fill = estimate, alpha = p_level_fdr)) + 
            geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
            scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) + ylab("Frequency") + facet_wrap(~sensor) +
            labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)))
    dev.off()
  }
}