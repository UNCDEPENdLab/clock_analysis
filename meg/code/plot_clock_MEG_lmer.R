setwd("/gpfs/group/mnh5174/default/Michael/Clock_MEG")
library(cowplot)
library(ggplot2)
library(data.table)
library(dplyr)
library(viridis)

plot_lmers <- function(allres_combined, clustername) {
  models <- sort(unique(allres_combined$model))
  frequencies <- sort(unique(allres_combined$Freq))
  freq_labels <- substr(frequencies, 1, 4)

  for (m in models) {
    mdf <- allres_combined %>% filter(model==m)

    #residuals check
    residcheck <- mdf %>% filter(group=="fixed" & term=="(Intercept)" & Time >= 0)
    
    #plot overall effects
    allanova <- mdf %>% filter(group=="anova")
    
    toplot <- sort(unique(allanova$term))
    pdf(paste0("figures/", clustername, "_", m, "_anova_figures.pdf"), width=12, height=10)
    for (n in toplot) {
      thiseff <- allanova %>% filter(term==n & p.value < .1 & Time >= 0)
      thiseff$Freq <- ordered(thiseff$Freq, levels=frequencies, labels=freq_labels)
      g <- ggplot(thiseff, aes(x=Time, y=Freq, fill=p.value, color=p.value)) + 
        geom_tile() + scale_color_viridis("pvalue", direction=-1) + 
        scale_fill_viridis("pvalue", direction=-1) + ggtitle(paste("Model: ", m, ", Effect: ", n)) #coord_fixed() +
      plot(g)
    }
    dev.off()
    
    allbeta <- mdf %>% filter(group=="fixed")
    toplot <- sort(unique(allbeta$term))
    pdf(paste0("figures/", clustername, "_", m, "_t_figures.pdf"), width=12, height=10)
    for (n in toplot) {
      thiseff <- allbeta %>% filter(term==n & Time >= 0) #& abs(statistic) > 0.5
      thiseff$Freq <- ordered(thiseff$Freq, levels=frequencies, labels=freq_labels)
      g <- ggplot(thiseff, aes(x=Time, y=Freq, color=statistic, fill=statistic)) + geom_tile() + 
        #scale_fill_viridis("t") + scale_color_viridis("t") + ggtitle(paste("Model: ", m, ", Effect: ", n)) #coord_fixed() +
        scale_fill_distiller("t", palette="RdBu") + scale_color_distiller("t", palette="RdBu") + ggtitle(paste("Model: ", m, ", Effect: ", n)) #coord_fixed() +
      plot(g)
    }
    dev.off()

  }

  #aic plot. we have one lmer per time x freq combination
  #we currently have 10 models: m0--m9
  allaic <- allres_combined %>% filter(group=="fixed" & term=="(Intercept)" & Time >= 0) %>% group_by(Time, Freq) %>%
    do({
      df <- .
      df <- df[which.min(df$AIC),]
      df
    }) %>% ungroup() %>%
    mutate(Freq = ordered(Freq, levels=frequencies, labels=freq_labels))

  g <- ggplot(allaic, aes(x=Time, y=Freq, color=model, fill=model)) + geom_tile() +
    scale_color_brewer("Best model", palette="Set3") + scale_fill_brewer("Best model", palette="Set3") +  
    ggtitle("Best model by time and frequency")

  pdf(paste0("figures/", clustername, "_model_aic_comparison.pdf"))
  plot(g)
  dev.off()
}


output_files <- list.files("cluster\\d+_.*_lmerfits.RData", path="output", full.names=TRUE)
#load("lmer_results_m1-m9.RData")

clusters <- unique(sub("output/(cluster\\d+_MEG\\d+)_.*\\.RData", "\\1", output_files, perl=TRUE))
frequencies <- unique(sub("output/cluster\\d+_MEG\\d+_([0-9\\.]+)_data_.*\\.RData", "\\1", output_files, perl=TRUE))

for (cl in clusters) {
  allres <- list()
  for (fr in frequencies) {
    load(paste0("output/", cl, "_", fr, "_data_lmerfits.RData"))
    allres[[make.names(fr)]] <- res    
  }
  allres_combined <- rbindlist(lapply(allres, rbindlist)) #stack into one huge data frame
  plot_lmers(allres_combined, cl)
}



#allres_combined <- rbindlist(lapply(allres, rbindlist))

#forgot to put $Freq on the first time...
# allres_combined <- rbindlist(
#   lapply(1:length(allres), function(fbin) {
#     withfreq <- lapply(allres[[fbin]], function(df) {
#       df$Freq <- round(as.numeric(names(allres[fbin])), 3)
#       return(df)
#     })
#     rbindlist(withfreq)
#   })
# )



