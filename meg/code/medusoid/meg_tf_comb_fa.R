# plots MLM coefficients by sensor across timepoints and frequencies
# first run meg_time_freq_load_analyze_medusa.R

library(tidyverse)
library(lme4)
library(ggpubr)
library(car)
library(viridis)
library(ggnewscale)

# library(psych)
repo_directory <- "~/code/clock_analysis"
data_dir <- "~/OneDrive/collected_letters/papers/meg/plots/tf_combined/"  
# rt_encode_plot_dir = "~/OneDrive/collected_letters/papers/meg/plots/rt_encode/"  
# clock_encode_plot_dir = "~/OneDrive/collected_letters/papers/meg/plots/clock_encode/"  
# dual_encode_plot_dir = "~/OneDrive/collected_letters/papers/meg/plots/dual_encode/"  
# 
# rt_rt_plot_dir = "~/OneDrive/collected_letters/papers/meg/plots/rt_rt/"
# clock_rt_plot_dir = "~/OneDrive/collected_letters/papers/meg/plots/clock_rt/"
# dual_rt_plot_dir = "~/OneDrive/collected_letters/papers/meg/plots/dual_rt/"

encode = T
rt_predict = F

sensor_map <- read.table("~/OneDrive/collected_letters/papers/meg/plots/meg_sensors_annotated.txt", header = T, colClasses = "character") %>%
  mutate(lobe=ordered(lobe, levels=c("frontal", "temporal", "parietal", "occipital")),
         dan=factor(dan, levels=c("yes", "maybe", "no"), labels=c("1yes", "2maybe", "3no"))) %>%
  #mutate(dan_grad=ordered(dan, levels=c("no", "maybe", "yes")))
  group_by(hemi, dan, lobe) %>% mutate(sensor_lab=paste(substr(dan[1], 1,2), substr(lobe[1], 1,1), sprintf("%02d", 1:n()), sensor, sep="_"),
                                       short_lab=paste(substr(dan[1], 1,2), substr(lobe[1], 1,1), sprintf("%02d", 1:n()), sep="_")) %>% ungroup() %>% 
  mutate(node = paste(lobe, hemi, sep = "_")) 

#group_by(dan, hemi) %>% mutate(newlab=paste(dan[1], 1:n(), sep="_")) %>% ungroup()
node_list <- unique(sensor_map$node[!grepl(pattern = "occip", x = sensor_map$node)]) %>% sort()
setwd(data_dir)
# plots ----
  alignment = "clock"
  # get clock-aligned data
  file_pattern <- "ddf_combined_RT"
  files <-  gsub("//", "/", list.files(data_dir, pattern = file_pattern, full.names = F))
  cl <- lapply(files, readRDS)
  names(cl) <- node_list[parse_number(files)]
  # names(l) <- node_list[1:length(files)]
  cddf <- data.table::rbindlist(cl, idcol = "node")
  cddf <- cddf %>% mutate(t  = as.numeric(Time), alignment = "clock",
                          term = str_replace(term, "rt_lag_sc", "RT_t"),
                          term = str_replace(term, "reward_lagReward", "reward_t"),
                          term = str_replace(term, "scale\\(rt_vmax_lag\\)", "RT_Vmax_t*"),
                          term = str_replace(term, "scale\\(abs_pe_lag\\)", "absPE_t")
  )
  # get RT-aligned
  file_pattern <- "ddf_combined_clock"
  files <-  gsub("//", "/", list.files(data_dir, pattern = file_pattern, full.names = F))
  rl <- lapply(files, readRDS)
  names(rl) <- node_list[parse_number(files)]
  rddf <- data.table::rbindlist(rl, idcol = "node")
  rddf <- rddf %>% mutate(t  = Time - 5, 
                          alignment = "rt",
                          term = str_replace(term, "rt_csv_sc", "RT_t"),
                          term = str_replace(term, "outcomeReward", "reward_t"),
                          term = str_replace(term, "scale\\(rt_vmax_lag\\)", "RT_Vmax_t*"),
                          term = str_replace(term, "scale\\(abs_pe\\)", "absPE_t")
  )
  
  ddf <- rbind(cddf, rddf)
  terms <- unique(ddf$term[ddf$effect=="fixed"])
  # ddf$sensor <- readr::parse_number(ddf$.filename, trim_ws = F)
  # ddf$sensor <- stringr::str_pad(ddf$sensor, 4, "0",side = "left")
  ddf <- ddf %>% mutate(p_value = as.factor(case_when(`p.value` > .05 ~ '1',
                                                      `p.value` < .05 & `p.value` > .01 ~ '2',
                                                      `p.value` < .01 & `p.value` > .001 ~ '3',
                                                      `p.value` <.001 ~ '4'))                        # sensor = as.character(sensor)
  )
  ddf$p_value <- factor(ddf$p_value, labels = c("NS", "p < .05", "p < .01", "p < .001"))
  
  # FDR labeling ----
  ddf <- ddf  %>% mutate(
    # p_fdr = padj_fdr_term,
    p_level_fdr = as.factor(case_when(
      # p_fdr > .1 ~ '0',
      # p_fdr < .1 & p_fdr > .05 ~ '1',
      padj_BY_term > .05 ~ '1',
      padj_BY_term < .05 & padj_BY_term > .01 ~ '2',
      padj_BY_term < .01 & padj_BY_term > .001 ~ '3',
      padj_BY_term <.001 ~ '4'))
  ) %>% ungroup() #%>% mutate(side = substr(as.character(label), nchar(as.character(label)), nchar(as.character(label))),
  # region = substr(as.character(label), 1, nchar(as.character(label))-2))
  ddf$p_level_fdr <- factor(ddf$p_level_fdr, levels = c('1', '2', '3', '4'), labels = c("NS","p < .05", "p < .01", "p < .001"))
  ddf$`p, FDR-corrected` = ddf$p_level_fdr
  
  levels(ddf$Freq) <- substr(unique(ddf$Freq), 3,6)
  ddf$freq <- fct_rev(ddf$Freq)
  # drop magnetometers (if any)
  # ddf <- ddf %>% filter(!grepl("1$", sensor))
  # add sensor labels
  setwd(data_dir)
  
  # arrange estimates into variables defined by time and frequency
  wdf <- ddf %>% select(node, t, freq, estimate) %>% pivot_wider(names_from = c(node,t,freq), values_from = estimate)
  