# plots MLM coefficients by sensor across timepoints and frequencies
# first run meg_time_freq_load_analyze_medusa.R

library(tidyverse)
library(lme4)
library(ggpubr)
library(car)
library(viridis)
library(ggnewscale)

# library(psych)

# what to run
# you can run all options at once
online = T
encode = T  # main analysis analogous to Fig. 4 E-G in NComm 2020
rt_predict = F # predicts next response based on signal and behavioral variables
# random = T # whether to use data from analyses where behavioral variables have both fixed and random effects
# uncorrected = F # whether to plot uncorrected data (FDR-corrected always plotted)

repo_directory <- "~/code/clock_analysis"
rt_data_directory <- "/Volumes/GoogleDrive/My Drive/SCEPTIC_fMRI/MEG/dan_source/RT_time/"
clock_data_directory <- "/Volumes/GoogleDrive/My Drive/SCEPTIC_fMRI/MEG/dan_source/clock_time/"

dual_plot_dir = "~/OneDrive/collected_letters/papers/meg/plots/dual_source/"  

clock_epoch_label = "Time relative to clock onset, seconds"
rt_epoch_label = "Time relative to outcome, seconds"

labels <- as_tibble(readxl::read_excel("~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/MNH Dan Labels.xlsx")) %>%
  select(c("roinum", "plot_label", "Stream", "Visuomotor_Gradient", "Stream_Gradient", "lobe")) %>% filter(!is.na(lobe))

names(labels) <- c("roinum","label_short", "stream", "visuomotor_grad", "stream_grad", "lobe")
labels$stream_grad <- as.numeric(labels$stream_grad)
labels <- labels %>% arrange(visuomotor_grad, stream_grad) %>% mutate(
  side  = case_when(
    grepl("L_", label_short) ~ "L",
    grepl("R_", label_short) ~ "R"),
  label_short = substr(label_short, 3, length(label_short)),
  label = paste(visuomotor_grad, stream_grad, label_short, side, sep = "_"),
  stream_side = paste0(stream, "_", side),
  visuomotor_side = paste0(visuomotor_grad, "_", side)) %>% 
  select(c(label, label_short, side, roinum, stream, visuomotor_grad, stream_grad, stream_side, visuomotor_side, lobe))
labels <- labels %>% mutate(node = paste(lobe, side, sep = "_")) 
short_labels <- labels %>% mutate(roinum = as.factor(roinum)) %>% select(roinum, node)

## plots ----
if (encode) {  
  message("Plotting decoding results")
  setwd(rt_data_directory)
  rddf <- readRDS("meg_mixed_by_time_ddf_source_RT.RDS") %>% 
    mutate(Time = parse_number(gsub(",.*$", "", time_bin)),
           t  = as.numeric(Time) - 4.5, alignment = "rt",
           # get shared variables
           term = case_when(
             term=="rt_csv_sc" ~ "RT_t",
             term=="outcomeReward" ~ "reward_t",
             term=="scale(rt_vmax_lag)" ~ "RT_Vmax_t",
             term=="scale(rt_vmax_change)" ~ "RT_Vmax_change_t",
             term=="v_entropy_wi_change" ~ "V_entropy_change_t",
             TRUE ~ term
           )
    )
  rddf$node <- sub(".rds", "", sub("RT_|clock_", "", rddf$.filename))
  setwd(clock_data_directory)
  cddf <- readRDS("meg_mixed_by_time_ddf_source_clock.RDS") %>% 
    mutate(Time = parse_number(gsub(",.*$", "", time_bin)),
           t  = as.numeric(Time), alignment = "clock",
           # get shared variables
           term = case_when(
             term=="rt_lag_sc" ~ "RT_t",
             term=="reward_lagReward" ~ "reward_t",
             term=="rt_vmax" ~ "RT_Vmax_t",
             term=="scale(rt_vmax_change)" ~ "RT_Vmax_change_t",
             term=="v_entropy_wi_change_lag" ~ "V_entropy_change_t",
             TRUE ~ term
           )
    ) %>% filter(t < 3.5)
  cddf$node <- sub(".rds", "", sub("RT_|clock_", "", cddf$.filename))
  ddf <- rbind(rddf, cddf)
  # ddf <- rddf
  
  # ddf <- readRDS("meg_mixed_by_time_ranefs_mult_interactions_pe_ddf.RDS")
  terms <- unique(ddf$term[ddf$effect=="fixed"])
  ddf <- ddf %>% mutate(p_value = as.factor(case_when(`p.value` > .05 ~ '1',
                                                      `p.value` < .05 & `p.value` > .01 ~ '2',
                                                      `p.value` < .01 & `p.value` > .001 ~ '3',
                                                      `p.value` <.001 ~ '4')))
  ddf$p_value <- factor(ddf$p_value, labels = c("NS", "p < .05", "p < .01", "p < .001"))
  # terms <- unique(ddf$term[ddf$effect=="fixed"])
  # terms <- terms[grepl(!"_f", terms)]
  # # FDR labeling ----
  ddf <- ddf  %>% mutate(p_fdr = padj_BY_term,
                         p_level_fdr = as.factor(case_when(
                           # p_fdr > .1 ~ '0',
                           # p_fdr < .1 & p_fdr > .05 ~ '1',
                           p_fdr > .05 ~ '1',
                           p_fdr < .05 & p_fdr > .01 ~ '2',
                           p_fdr < .01 & p_fdr > .001 ~ '3',
                           p_fdr <.001 ~ '4')),
                         lobe = gsub("[_L|)R]", "", node),
                         side = gsub("frontal_|parietal_|temporal_", "", node))
  ddf$p_level_fdr <- factor(ddf$p_level_fdr, levels = c('1', '2', '3', '4'), labels = c("NS","p < .05", "p < .01", "p < .001"))
  ddf$`p, FDR-corrected` = ddf$p_level_fdr
  
  setwd(dual_plot_dir)
  for (fe in terms) {
    edf <- ddf %>% filter(term == paste(fe)) 
    termstr <- str_replace_all(fe, "[^[:alnum:]]", "_")
    fname = paste("meg_time_dual_uncorrected_abs", termstr, ".pdf", sep = "")
    message(fname)
    pdf(fname, width = 10, height = 3)
    print(ggplot(edf, aes(t, 1)) + geom_tile(aes(fill = abs(estimate), alpha = p_value)) + 
            facet_grid(lobe~side) + 
            geom_vline(xintercept = 0, lty = "dashed", color = "black", size = 2) +
            geom_vline(xintercept = -4.5, lty = "dashed", color = "red", size = 2) +
            geom_vline(xintercept = -4.8, lty = "dashed", color = "red", size = 1) +
            geom_vline(xintercept = -2.25, lty = "dotted", color = "grey", size = 1) +
            scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(clock_epoch_label) + 
            labs(alpha = expression(italic(p)[uncorrected])) + ggtitle(paste(termstr))) + theme(
              axis.text.y = element_blank())
    dev.off()
    fname = paste("meg_time_dual_FDR_abs", termstr, ".pdf", sep = "")
    pdf(fname, width = 10, height = 3)
    print(ggplot(edf, aes(t, 1)) + geom_tile(aes(fill = abs(estimate), alpha = p_level_fdr)) + 
            facet_grid(lobe~side) + 
            geom_vline(xintercept = 0, lty = "dashed", color = "black", size = 2) +
            geom_vline(xintercept = -4.5, lty = "dashed", color = "red", size = 2) +
            geom_vline(xintercept = -4.8, lty = "dashed", color = "red", size = 1) +
            
            geom_vline(xintercept = -2.25, lty = "dotted", color = "grey", size = 1) +
            scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(clock_epoch_label) +
            labs(alpha = expression(italic(p)[uncorrected])) + ggtitle(paste(termstr))) + theme(
              axis.text.y = element_blank())
    dev.off()
  }
} 
# system("for i in *FDR*short*.pdf; do sips -s format jpeg $i --out $i.jpeg; done")
# system("for i in *PE*.pdf; do sips -s format jpeg $i --out $i.jpeg; done")

if(rt_predict) {
  message("Plotting decoding results")
  # read in RT-aligned data
  setwd(rt_data_directory)
  rrdf <- readRDS("meg_mixed_by_time_rdf_source_RT.RDS") %>% 
    mutate(Time = parse_number(gsub(",.*$", "", time_bin)),
           t  = as.numeric(Time) - 4.5, alignment = "rt",
           term = case_when(
             term=="source_est:rt_csv_sc" ~ "source_est:RT_t",
             term=="source_est:rt_csv_sc:v_entropy_wi" ~ "source_est:RT_t:v_entropy_wi",
             TRUE ~ term))
  rrdf$node <- sub(".rds", "", sub("RT_|clock_", "", rrdf$.filename))
  # read in clock-aligned data
  setwd(clock_data_directory)
  crdf <- readRDS("meg_mixed_by_time_rdf_source_clock.RDS") %>% 
    mutate(Time = parse_number(gsub(",.*$", "", time_bin)),
           t  = as.numeric(Time), 
           alignment = "clock",
           term = case_when(
             term=="source_est:rt_lag_sc" ~ "source_est:RT_t",
             term=="source_est:rt_lag_sc:v_entropy_wi" ~ "source_est:RT_t:v_entropy_wi",
             TRUE ~ term)
    ) 
  crdf$node <- sub(".rds", "", sub("RT_|clock_", "", crdf$.filename))
  rdf <- rbind(rrdf, crdf)
  # ddf <- readRDS("meg_mixed_by_time_ranefs_mult_interactions_pe_ddf.RDS")
  terms <- unique(rdf$term[rdf$effect=="fixed"])
  terms <- terms[grepl("source", terms)]
  rdf <- rdf %>%  mutate(p_fdr = padj_BY_term,
                         p_level_fdr = as.factor(case_when(
                           # p_fdr > .1 ~ '0',
                           # p_fdr < .1 & p_fdr > .05 ~ '1',
                           p_fdr > .05 ~ '1',
                           p_fdr < .05 & p_fdr > .01 ~ '2',
                           p_fdr < .01 & p_fdr > .001 ~ '3',
                           p_fdr <.001 ~ '4')),
                         lobe = gsub("[_L|)R]", "", node),
                         side = gsub("frontal_|parietal_|temporal_", "", node),
                         p_value = as.factor(case_when(`p.value` > .05 ~ '1',
                                                      `p.value` < .05 & `p.value` > .01 ~ '2',
                                                      `p.value` < .01 & `p.value` > .001 ~ '3',
                                                      `p.value` <.001 ~ '4')))
  
  rdf$p_value <- factor(rdf$p_value, labels = c("NS", "p < .05", "p < .01", "p < .001"))
  terms <- unique(rdf$term[rdf$effect=="fixed"])
  terms <- terms[grepl("source", terms)]
    rdf$p_level_fdr <- factor(rdf$p_level_fdr, levels = c('1', '2', '3', '4'), labels = c("NS","p < .05", "p < .01", "p < .001"))
  rdf$`p, FDR-corrected` = rdf$p_level_fdr
  
  setwd(dual_plot_dir)
  for (fe in terms) {
    edf <- rdf %>% filter(term == paste(fe)) 
    termstr <- str_replace_all(fe, "[^[:alnum:]]", "_")
    fname = paste("meg_time_dual_uncorrected_source_", termstr, ".pdf", sep = "")
    message(fname)
    pdf(fname, width = 15, height = 6)
    print(ggplot(edf, aes(t, 1)) + geom_tile(aes(fill = estimate, alpha = p_value)) + 
            facet_grid(lobe~side) + 
            scale_fill_distiller(palette = "OrRd", direction = 1, name = "Exploration", limits = c(-.3, 0)) + scale_x_continuous(breaks = pretty(edf$t, n = 20)) + labs(fill = "Exploration") +
            labs(alpha = expression(italic(p)[uncorrected])) + ggtitle(paste(termstr)) +
            new_scale_fill() +
            geom_tile(data = edf %>% filter(estimate > 0), aes(t, 1, fill = estimate, alpha = p_value), size = .01) + theme_dark() +
            geom_vline(xintercept = 0, lty = "dashed", color = "black", size = 2) + theme_bw() + 
            scale_fill_distiller(palette = "GnBu", direction = -1, name = "Exploitation", limits = c(0, .3))+ scale_color_grey() + xlab(clock_epoch_label) + 
            geom_vline(xintercept = 0, lty = "dashed", color = "black", size = 2) +
            geom_vline(xintercept = -4.5, lty = "dashed", color = "black", size = 2) +
            geom_vline(xintercept = -4.8, lty = "dashed", color = "black", size = 1) +
            geom_vline(xintercept = -2.25, lty = "dotted", color = "grey", size = 1) +
            # scale_fill_distiller(palette = "RdBu", limits = c(-.05,.05)) + scale_color_grey() + xlab(clock_epoch_label) + 
            labs(alpha = expression(italic(p)[uncorrected])) + ggtitle(paste(termstr)))
    dev.off()
    fname = paste("meg_time_dual_FDR_source_", termstr, ".pdf", sep = "")
    pdf(fname, width = 15, height = 6)
    print(ggplot(edf, aes(t, 1)) + geom_tile(aes(fill = estimate, alpha = p_level_fdr)) + 
            facet_grid(lobe~side) + 
            geom_vline(xintercept = 0, lty = "dashed", color = "black", size = 2) +
            geom_vline(xintercept = -4.5, lty = "dashed", color = "black", size = 2) +
            geom_vline(xintercept = -4.8, lty = "dashed", color = "black", size = 1) +
            geom_vline(xintercept = -2.25, lty = "dotted", color = "grey", size = 1) +
            scale_fill_distiller(palette = "OrRd", direction = 1, name = "Exploration", limits = c(-.3, 0)) + scale_x_continuous(breaks = pretty(edf$t, n = 20)) + labs(fill = "Exploration") +
            labs(alpha = expression(italic(p)[uncorrected])) + ggtitle(paste(termstr)) +
            new_scale_fill() +
            geom_tile(data = edf %>% filter(estimate > 0), aes(t, 1, fill = estimate, alpha = p_value), size = .01) + theme_dark() +
            geom_vline(xintercept = 0, lty = "dashed", color = "black", size = 2) + theme_bw() + 
            scale_fill_distiller(palette = "GnBu", direction = -1, name = "Exploitation", limits = c(0, .3))+ scale_color_grey() + xlab(clock_epoch_label) + 
            
            labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)))
    dev.off()
  }
  
  # convert to PNG for Word
  # system("for i in *source_est*.pdf; do sips -s format png $i --out $i.png; done")
}
