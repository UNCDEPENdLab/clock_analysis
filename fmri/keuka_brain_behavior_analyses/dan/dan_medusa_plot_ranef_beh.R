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
# repo_directory <- "~/Data_Analysis/clock_analysis"

# what to run
# you can run all options at once
decode = T  # main analysis analogous to Fig. 4 E-G in NComm 2020
rt_predict = F # predicts next response based on signal and behavioral variables
alignments = c(F) # whether to analyze clock-aligned ("online", T) or RT-aligned ("offline", F) responses
meg = "_MEG"
# directories
# setwd(file.path(repo_directory, 'fmri/keuka_brain_behavior_analyses/'))
plot_dir <- "~/OneDrive/collected_letters/papers/sceptic_fmri/dan/plots/ranef_beh_rt/"
for (online in alignments)
{
  # plots ----
  message("\nPlotting RT prediction")
  message("\nPlotting RT visuomotor")
  if (online) {
    setwd(plot_dir)
    epoch_label = "Time relative to clock onset, seconds"  
    clock_results_fname = "clock_rt_output_visuomotor.Rdata"
  } else {
    setwd(plot_dir)
    epoch_label = "Time relative to outcome, seconds"
    rt_results_fname = paste0("mixed_by_entropy_pe_ranef_effects", meg, ".rds")
  }
  rdf <- readRDS(rt_results_fname)
  terms <- unique(rdf$term[grepl("blup", rdf$term)])
  rdf <- rdf %>% mutate(p_BY = padj_BY_term,
                                          p_level_BY = as.factor(case_when(
                                            # p_BY > .1 ~ '0',
                                            # p_BY < .1 & p_BY > .05 ~ '1',
                                            p_BY > .05 ~ '1',
                                            p_BY < .05 & p_BY > .01 ~ '2',
                                            p_BY < .01 & p_BY > .001 ~ '3',
                                            p_BY <.001 ~ '4')),
                        t = evt_time,
                        ranef = case_when(
                          brain_ranef=="scale.abs_pe." ~ "abs(PE)",
                          brain_ranef=="v_entropy_wi" ~ "entropy"
                        )
  ) %>% ungroup() #%>% mutate(side = substr(as.character(label), nchar(as.character(label)), nchar(as.character(label))),
  #          region = substr(as.character(label), 1, nchar(as.character(label))-2))
  rdf$p_level_BY <- factor(rdf$p_level_BY, levels = c('1', '2', '3', '4'), labels = c("NS","p < .05", "p < .01", "p < .001"))
  rdf$`p, BY-corrected` = rdf$p_level_BY
  
  for (ranef in unique(rdf$ranef)) {
    for (fe in terms) {
      edf <- rdf %>% filter(term == paste(fe) & ranef==!!ranef) 
      termstr <- str_replace_all(fe, "[^[:alnum:]]", "_")
      # plot stream gradients
      fname = paste(meg, ranef, "streams", termstr, ".pdf", sep = "_")
      pdf(fname, width = 9, height = 3.5)
      print(ggplot(edf, aes(t, stream)) + geom_tile(aes(fill = estimate, alpha = `p, BY-corrected`), size = 1) +  
              geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + facet_wrap(~side) +
              scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) + 
              labs(alpha = expression(italic(p)[BY])) + ggtitle(paste(termstr)) + ylab("") +
              scale_y_discrete(labels=c("visual-motion" = "MT+,\ncontrol", "ventro-dorsal" = "Ventro-dorsal\nstream",
                                        "oculomotor" = "Oculomotor\nstream", "dorso-dorsal" = "Dorso-dorsal\nstream")))
      dev.off()
    }
  } 
}




