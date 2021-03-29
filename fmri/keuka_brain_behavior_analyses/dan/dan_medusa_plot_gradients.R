# plots MLM coefficients by stream (dorso-dorsal, oculomotor, dorso-ventral) and visuomotor gradient
# first run dan_decode_rt_prediction.R

# library(modelr)
library(tidyverse)
library(lme4)
# library(afex)
# library(broom)
# library(broom.mixed) #plays will with afex p-values in lmer wrapper
library(ggpubr)
library(car)
library(viridis)
library(psych)
# library(corrplot)
# library(foreach)
# library(doParallel)
library(readxl)
repo_directory <- "~/code/clock_analysis"

# what to run
decode = T  # main analysis analogous to Fig. 4 E-G in NComm 2020
rt_predict = F # predicts next response based on signal and behavioral variables
online = T # whether to analyze clock-aligned ("online") or RT-aligned ("offline") responses

# directories
setwd('~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/')
cache_dir <- "~/Box/SCEPTIC_fMRI/dan_medusa/cache/"
repo_dir <- "~/code/clock_analysis"


if(decode) {
  message("\nPlotting decoding data")
  if (online) {
    setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/dan/plots/clock_decode')
    epoch_label = "Time relative to clock onset, seconds"
    decode_results_fname = "clock_decode_output.Rdata"
  } else {setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/dan/plots/rt_decode')
    epoch_label = "Time relative to outcome, seconds"
    decode_results_fname = "rt_decode_output.Rdata"}
  load(decode_results_fname)
  terms <- unique(ddf$term[ddf$effect=="fixed"])
  for (fe in terms) {
    # fe <- terms[1] # test only
    edf <- ddf %>% filter(term == paste(fe) & t < 8) 
    termstr <- str_replace_all(fe, "[^[:alnum:]]", "_")
    # plot visuomotor gradients
    fname = paste("visuomotor_", termstr, ".pdf", sep = "")
    pdf(fname, width = 9, height = 3.5)
    print(ggplot(edf, aes(t, as.factor(visuomotor_grad))) + geom_tile(aes(fill = estimate, alpha = `p, FDR-corrected`), size = 1) +  
            geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + facet_wrap(~side) +
            scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) + 
            labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab("") +
      scale_y_discrete(labels=c("1" = "MT+, control", "2" = "Parieto-occipital", "3" = "Post. parietal", "4" = "Frontal")))
    dev.off()
    # plot stream gradients
    fname = paste("streams_", termstr, ".pdf", sep = "")
    pdf(fname, width = 9, height = 3.5)
    print(ggplot(edf, aes(t, stream)) + geom_tile(aes(fill = estimate, alpha = `p, FDR-corrected`), size = 1) +  
            geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + facet_wrap(~side) +
            scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) + 
            labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab("") +
      scale_y_discrete(labels=c("visual-motion" = "MT+,\ncontrol", "ventro-dorsal" = "Ventro-dorsal\nstream",
                                "oculomotor" = "Oculomotor\nstream", "dorso-dorsal" = "Dorso-dorsal\nstream")))
    dev.off()
    
  }
}
# add labels
## RT prediction ----

if(rt_predict) {
  # plots ----
  message("\nPlotting RT prediction")
  if (online) {
    setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/dan/plots/clock_rt')
    epoch_label = "Time relative to clock onset, seconds"  
    rt_results_fname = "clock_rt_predict_output.Rdata"
  } else {
    setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/dan/plots/rt_rt')
    epoch_label = "Time relative to outcome, seconds"
    rt_results_fname = "rt_rt_predict_output.Rdata"
  }
  load(rt_results_fname)
  terms <- unique(rdf$term[rdf$effect=="fixed"])
  terms <- terms[grepl("(h)",terms)]
  for (fe in terms) {
    edf <- rdf %>% filter(term == paste(fe) & t < 8) 
    termstr <- str_replace_all(fe, "[^[:alnum:]]", "_")
    # plot visuomotor gradients
    fname = paste("visuomotor_", termstr, ".pdf", sep = "")
    pdf(fname, width = 9, height = 3.5)
    print(ggplot(edf, aes(t, as.factor(visuomotor_grad))) + geom_tile(aes(fill = estimate, alpha = `p, FDR-corrected`), size = 1) +  
            geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + facet_wrap(~side) +
            scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) + 
            labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab("") +
            scale_y_discrete(labels=c("1" = "MT+, control", "2" = "Parieto-occipital", "3" = "Post. parietal", "4" = "Frontal")))
    dev.off()
    # plot stream gradients
    fname = paste("streams_", termstr, ".pdf", sep = "")
    pdf(fname, width = 9, height = 3.5)
    print(ggplot(edf, aes(t, stream)) + geom_tile(aes(fill = estimate, alpha = `p, FDR-corrected`), size = 1) +  
            geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + facet_wrap(~side) +
            scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) + 
            labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab("") +
            scale_y_discrete(labels=c("visual-motion" = "MT+,\ncontrol", "ventro-dorsal" = "Ventro-dorsal\nstream",
                                      "oculomotor" = "Oculomotor\nstream", "dorso-dorsal" = "Dorso-dorsal\nstream")))
    dev.off()
  }
}
  



