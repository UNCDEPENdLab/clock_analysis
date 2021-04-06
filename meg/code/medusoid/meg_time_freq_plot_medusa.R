# plots MLM coefficients by sensor across timepoints and frequencies
# first run meg_time_freq_load_analyze_medusa.R

library(tidyverse)
library(lme4)
library(ggpubr)
library(car)
library(viridis)
# library(psych)
repo_directory <- "~/code/clock_analysis"
decode_plot_dir = "~/OneDrive/collected_letters/papers/meg/plots/rt_decode/"
rt_plot_dir = "~/OneDrive/collected_letters/papers/meg/plots/rt_rt//"


# what to run
# you can run all options at once
decode = F  # main analysis analogous to Fig. 4 E-G in NComm 2020
rt_predict = T # predicts next response based on signal and behavioral variables
random = T # whether to use data from analyses where behavioral variables have both fixed and random effects
uncorrected = F # whether to plot uncorrected data (FDR-corrected always plotted)

# plots ----
if (decode) {  
  message("Plotting decoding results")
  setwd(decode_plot_dir)
  epoch_label = "Time relative to outcome, seconds"
  load("meg_freq_medusa_decode_output_all.Rdata")
  terms <- unique(ddf$term[ddf$effect=="fixed"])
  
  # # within-sensor FDR correction
  # ddf <- ddf  %>% group_by(term, sensor) %>% mutate(p_fdr = p.adjust(p.value, method = 'fdr'),
  #                                                   p_level_fdr = as.factor(case_when(
  #                                                     # p_fdr > .1 ~ '0',
  #                                                     # p_fdr < .1 & p_fdr > .05 ~ '1',
  #                                                     p_fdr > .05 ~ '1',
  #                                                     p_fdr < .05 & p_fdr > .01 ~ '2',
  #                                                     p_fdr < .01 & p_fdr > .001 ~ '3',
  #                                                     p_fdr <.001 ~ '4'))
  # ) %>% ungroup() #%>% mutate(side = substr(as.character(label), nchar(as.character(label)), nchar(as.character(label))),
  # #          region = substr(as.character(label), 1, nchar(as.character(label))-2))
  # ddf$p_level_fdr <- factor(ddf$p_level_fdr, levels = c('1', '2', '3', '4'), labels = c("NS","p < .05", "p < .01", "p < .001"))
  # ddf$`p, FDR-corrected` = ddf$p_level_fdr
  
  
  for (fe in terms) {
    edf <- ddf %>% filter(term == paste(fe)) 
    termstr <- str_replace_all(fe, "[^[:alnum:]]", "_")
    fname = paste("meg_tf_dan_uncorrected_", termstr, ".pdf", sep = "")
    pdf(fname, width = 18, height = 8)
    print(ggplot(edf, aes(t, freq)) + geom_tile(aes(fill = estimate, alpha = p_value), size = .01) + 
            geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
            scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) + ylab("Frequency") + facet_wrap(~sensor) +
            labs(alpha = expression(italic(p)[uncorrected])) + ggtitle(paste(termstr)))
    dev.off()
    fname = paste("meg_tf_dan_FDR_", termstr, ".pdf", sep = "")
    pdf(fname, width = 18, height = 8)
    print(ggplot(edf, aes(t, freq)) + geom_tile(aes(fill = estimate, alpha = p_level_fdr), size = .01) + 
            geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
            scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) + ylab("Frequency") + facet_wrap(~sensor) +
            labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)))
    dev.off()
  }
} 

if(rt_predict) {
  # plots ----
  setwd('~/OneDrive/collected_letters/papers/meg/plots/rt_rt')
  epoch_label = "Time relative to outcome, seconds"
  if (random) {rt_results_fname = "meg_freq_medusa_rt_predict_output_nested.Rdata"} else {
    rt_results_fname = "meg_freq_medusa_rt_predict_output_nested.Rdata"  
  }
  load(rt_results_fname)
  terms <- unique(rdf$term[rdf$effect=="fixed"])
  terms <- terms[grepl("(pow)",terms)]
  levels(rdf$freq) <- substr(unique(rdf$freq), 3,6)
  rdf$freq <- fct_rev(rdf$freq)
  # # within-sensor FDR correction
  # rdf <- rdf  %>% group_by(term, sensor) %>% mutate(p_fdr = p.adjust(p.value, method = 'fdr'),
  #                                                   p_level_fdr = as.factor(case_when(
  #                                                     # p_fdr > .1 ~ '0',
  #                                                     # p_fdr < .1 & p_fdr > .05 ~ '1',
  #                                                     p_fdr > .05 ~ '1',
  #                                                     p_fdr < .05 & p_fdr > .01 ~ '2',
  #                                                     p_fdr < .01 & p_fdr > .001 ~ '3',
  #                                                     p_fdr <.001 ~ '4'))
  # ) %>% ungroup() #%>% mutate(side = substr(as.character(label), nchar(as.character(label)), nchar(as.character(label))),
  # #          region = substr(as.character(label), 1, nchar(as.character(label))-2))
  # rdf$p_level_fdr <- factor(rdf$p_level_fdr, levels = c('1', '2', '3', '4'), labels = c("NS","p < .05", "p < .01", "p < .001"))
  # rdf$`p, FDR-corrected` = rdf$p_level_fdr
  library(ggnewscale)

  for (fe in terms) {
    edf <- rdf %>% filter(term == paste(fe)) 
    termstr <- str_replace_all(fe, "[^[:alnum:]]", "_")
<<<<<<< HEAD
    if (uncorrected) {
    if (random) {fname = paste("meg_tf_dan_uncorrected_random_", termstr, ".pdf", sep = "")  
    } else {fname = paste("meg_tf_dan_uncorrected_", termstr, ".pdf", sep = "")}
    
    pdf(fname, width = 20, height = 8)
=======
    fname = paste("meg_tf_dan_uncorrected_random_", termstr, ".pdf", sep = "")
    pdf(fname, width = 18, height = 8)
>>>>>>> 720825a3639b026074ddbbf98398e8dd3d71ec23
    print(ggplot(edf %>% filter(estimate < 0), aes(t, freq)) + geom_tile(aes(fill = estimate, alpha = p_value), size = .01) + 
            scale_fill_distiller(palette = "OrRd", direction = 1, limits = c(-.2, 0)) + scale_color_grey(limits = rev(levels(rdf$p_level_fdr))) + xlab(epoch_label) + ylab("Frequency") + facet_wrap(~sensor) +
            labs(alpha = expression(italic(p)[uncorrected])) + ggtitle(paste(termstr)) +
            new_scale_fill() +
            geom_tile(data = edf %>% filter(estimate > 0), aes(t, freq, fill = estimate, alpha = p_value), size = .01) + theme_dark() +
            geom_vline(xintercept = 0, lty = "dashed", color = "black", size = 2) + theme_bw() + 
            scale_fill_distiller(palette = "GnBu", direction = -1, limits = c(0, .2))+ scale_color_grey() + xlab(epoch_label) + ylab("Frequency") + facet_wrap(~sensor, ncol = 3) +
            labs(alpha = expression(italic(p)[uncorrected])) + ggtitle(paste(termstr))) + scale_y_discrete(limits = rev(levels(rdf$p_level_fdr)))
    dev.off()
<<<<<<< HEAD
    }
    if (random) {fname = paste("meg_tf_dan_FDR_random_", termstr, ".pdf", sep = "")  
    } else {fname = paste("meg_tf_dan_FDR_", termstr, ".pdf", sep = "")}
    pdf(fname, width = 20, height = 8)
=======
    
    fname = paste("meg_tf_dan_FDR_random_", termstr, ".pdf", sep = "")
    pdf(fname, width = 18, height = 8)
>>>>>>> 720825a3639b026074ddbbf98398e8dd3d71ec23
    print(ggplot(edf %>% filter(estimate < 0), aes(t, freq)) + geom_tile(aes(fill = estimate, alpha = p_level_fdr), size = .01) + 
            scale_fill_distiller(palette = "OrRd", direction = 1, name = "Exploration", limits = c(-.3, 0)) + scale_x_continuous(breaks = pretty(edf$t, n = 20)) + labs(fill = "Exploration") +
            new_scale_fill() +
<<<<<<< HEAD
            geom_tile(data = edf %>% filter(estimate > 0), aes(t, freq, fill = estimate, alpha = p_level_fdr), size = .01) +
            geom_vline(xintercept = 0, lty = "dashed", color = "black", size = 2) + theme_bw() + 
            scale_fill_distiller(palette = "GnBu", direction = -1, name = "Exploitation", limits = c(0, .3))+ scale_color_grey() + xlab(epoch_label) + ylab("Frequency") + facet_wrap(~sensor, ncol = 3) +
            labs(alpha = expression(italic(p)[FDR-corrected])) + ggtitle(paste(termstr))) + labs(fill = "Exploitation")

=======
            geom_tile(data = edf %>% filter(estimate > 0), aes(t, freq, fill = estimate, alpha = p_value), size = .01) +
            geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + theme_dark() + 
            scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) + ylab("Frequency") + facet_wrap(~sensor) +
            labs(alpha = expression(italic(p)[FDR-corrected])) + ggtitle(paste(termstr))) 
>>>>>>> 720825a3639b026074ddbbf98398e8dd3d71ec23
    dev.off()
  }
  
  # TESTING THE "UNIFIED" VERSION
  unified = T
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
              scale_fill_distiller(palette = "OrRd", direction = 1, name = "Exploration") + scale_x_continuous(breaks = pretty(edf$t, n = 20)) + labs(fill = "Exploration") +
              new_scale_fill() +
              geom_tile(data = edf %>% filter(estimate > 0), aes(t, freq, fill = estimate, alpha = p_level_fdr), size = .01) +
              geom_vline(xintercept = 0, lty = "dashed", color = "black", size = 2) + theme_bw() + 
              scale_fill_distiller(palette = "GnBu", direction = -1, name = "Exploitation")+ scale_color_grey() + xlab(epoch_label) + ylab("Frequency") + 
              labs(alpha = expression(italic(p)[FDR-corrected])) + ggtitle(paste(termstr))) + labs(fill = "Exploration")
      
      dev.off()
    }
  }
}
