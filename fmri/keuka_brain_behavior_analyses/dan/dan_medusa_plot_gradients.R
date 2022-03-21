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
#repo_directory <- "~/code/clock_analysis"
repo_directory <- "~/Data_Analysis/clock_analysis"

# what to run
# you can run all options at once
decode = T  # main analysis analogous to Fig. 4 E-G in NComm 2020
rt_predict = T # predicts next response based on signal and behavioral variables
alignments = c(T,F) # whether to analyze clock-aligned ("online") or RT-aligned ("offline") responses

# directories
setwd(file.path(repo_directory, 'fmri/keuka_brain_behavior_analyses/'))
cache_dir <- "~/Box/SCEPTIC_fMRI/dan_medusa/cache/"

for (online in alignments)
{
  if(decode) {
    message("\nPlotting decoding data")
    message("\nPlotting visuomotor decoding")
    if (online) {
      setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/dan/plots/clock_decode')
      epoch_label = "Time relative to clock onset, seconds"
      decode_results_fname = "clock_decode_output_visuomotor.Rdata"
    } else {setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/dan/plots/rt_decode')
      epoch_label = "Time relative to outcome, seconds"
      decode_results_fname = "rt_decode_output_visuomotor.Rdata"}
    load(decode_results_fname)
    terms <- unique(ddf$term[ddf$effect=="fixed"])
    for (fe in terms) {
      # fe <- terms[1] # test only
      edf <- ddf %>% filter(term == paste(fe) & t < 8) 
      termstr <- str_replace_all(fe, "[^[:alnum:]]", "_")
      # plot visuomotor gradients
      fname = paste("visuomotor_", termstr, ".pdf", sep = "")
      pdf(fname, width = 9, height = 3.5)
      print(ggplot(edf, aes(t, as.factor(zone))) + geom_tile(aes(fill = estimate, alpha = `p, FDR-corrected`), size = 1) +
              geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + facet_wrap(~side) +
              scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) +
              labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab("") +
              scale_y_discrete(labels=c("1" = "MT+, control", "2" = "Parieto-occipital", "3" = "Post. parietal", "4" = "Frontal")))
      dev.off()
    }
    message("\nPlotting streams decoding")
    if (online) {
      setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/dan/plots/clock_decode')
      epoch_label = "Time relative to clock onset, seconds"
      decode_results_fname = "clock_decode_output_streams.Rdata"
    } else {setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/dan/plots/rt_decode')
      epoch_label = "Time relative to outcome, seconds"
      decode_results_fname = "rt_decode_output_streams.Rdata"}
    load(decode_results_fname)
    terms <- unique(ddf$term[ddf$effect=="fixed"])
    for (fe in terms) {
      # fe <- terms[1] # test only
      edf <- ddf %>% filter(term == paste(fe) & t < 8) 
      termstr <- str_replace_all(fe, "[^[:alnum:]]", "_")
      # plot stream gradients
      fname = paste("streams_", termstr, ".pdf", sep = "")
      pdf(fname, width = 9, height = 3.5)
      print(ggplot(edf, aes(t, zone)) + geom_tile(aes(fill = estimate, alpha = `p, FDR-corrected`), size = 1) +
              geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + facet_wrap(~side) +
              scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) +
              labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab("") +
              scale_y_discrete(labels=c("visual-motion" = "MT+,\ncontrol", "ventro-dorsal" = "Ventro-dorsal\nstream",
                                        "oculomotor" = "Oculomotor\nstream", "dorso-dorsal" = "Dorso-dorsal\nstream")))
      dev.off()

      fname = paste("streams_line_", termstr, ".pdf", sep = "")
      pdf(fname, width = 9, height = 3.5)
      gg <- ggplot(edf, aes(x=t, y=estimate, ymin=estimate-std.error, ymax=estimate+std.error, color=zone, size=`p, FDR-corrected`)) + 
              geom_line(size = 1) + geom_point() +
              geom_errorbar() +
              geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + facet_wrap(~side) +
              scale_color_brewer(palette="Set1") + xlab(epoch_label) + 
              labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab("")
      print(gg)
      dev.off()

    }
  }
  # add labels
  ## RT prediction ----
  
  if(rt_predict) {
    # plots ----
    message("\nPlotting RT prediction")
    message("\nPlotting RT visuomotor")
    if (online) {
      setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/dan/plots/clock_rt')
      epoch_label = "Time relative to clock onset, seconds"  
      rt_results_fname = "clock_rt_output_visuomotor.Rdata"
    } else {
      setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/dan/plots/rt_rt')
      epoch_label = "Time relative to outcome, seconds"
      rt_results_fname = "rt_rt_output_visuomotor.Rdata"
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
      print(ggplot(edf, aes(t, as.factor(region))) + geom_tile(aes(fill = estimate, alpha = `p, FDR-corrected`), size = 1) +  
              geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + facet_wrap(~side) +
              scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) + 
              labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab("") +
              scale_y_discrete(labels=c("1" = "MT+, control", "2" = "Parieto-occipital", "3" = "Post. parietal", "4" = "Frontal")))
      dev.off()
    }
    message("\nPlotting RT streams")
    if (online) {
      setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/dan/plots/clock_rt')
      epoch_label = "Time relative to clock onset, seconds"  
      rt_results_fname = "clock_rt_output_streams.Rdata"
    } else {
      setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/dan/plots/rt_rt')
      epoch_label = "Time relative to outcome, seconds"
      rt_results_fname = "rt_rt_output_streams.Rdata"
    }
    load(rt_results_fname)
    terms <- unique(rdf$term[rdf$effect=="fixed"])
    terms <- terms[grepl("(h)",terms)]
    for (fe in terms) {
      edf <- rdf %>% filter(term == paste(fe) & t < 8) 
      termstr <- str_replace_all(fe, "[^[:alnum:]]", "_")
      # plot stream gradients
      fname = paste("streams_", termstr, ".pdf", sep = "")
      pdf(fname, width = 9, height = 3.5)
      print(ggplot(edf, aes(t, region)) + geom_tile(aes(fill = estimate, alpha = `p, FDR-corrected`), size = 1) +  
              geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + facet_wrap(~side) +
              scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) + 
              labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab("") +
              scale_y_discrete(labels=c("visual-motion" = "MT+,\ncontrol", "ventro-dorsal" = "Ventro-dorsal\nstream",
                                        "oculomotor" = "Oculomotor\nstream", "dorso-dorsal" = "Dorso-dorsal\nstream")))
      dev.off()
      
      fname = paste("streams_line_", termstr, ".pdf", sep = "")
      pdf(fname, width = 9, height = 3.5)
      gg <- ggplot(edf, aes(x=t, y=estimate, ymin=estimate-std.error, ymax=estimate+std.error, color=region, size=`p, FDR-corrected`)) + 
        geom_line(size = 1) + geom_point() +
        geom_errorbar() +
        geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + facet_wrap(~side) +
        scale_color_brewer(palette="Set1") + xlab(epoch_label) + 
        labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)) + ylab("")
      print(gg)
      dev.off()
      
    }
  }
}




