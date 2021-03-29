# with 'decode = T' makes MEDUSA decoding plots for Fig. 4 E-G.
# loops over decoding and RT prediction multi-level models for various regions and time points
# first run medusa_event_locked_lmer.R

library(modelr)
library(tidyverse)
library(lme4)
library(afex)
library(broom)
library(broom.mixed) #plays will with afex p-values in lmer wrapper
library(ggpubr)
library(car)
library(viridis)
library(psych)
library(corrplot)
library(foreach)
library(doParallel)
library(readxl)
repo_directory <- "~/code/clock_analysis"

# data & options ----

# data loading options
reprocess = T # otherwise load data from cache
if (!reprocess) {
  wide_only = T # only load wide data (parcels and timepoints as variables)
}
replicate_compression = F
if(replicate_compression) {reprocess = T}
# load MEDUSA deconvolved data
source(file.path(repo_directory, "fmri/keuka_brain_behavior_analyses/dan/load_medusa_data_dan.R"))

# what to run
plots = T
decode = T  # main analysis analogous to Fig. 4 E-G in NComm 2020
rt_predict = T # predicts next response based on signal and behavioral variables
online = F # whether to analyze clock-aligned ("online") or RT-aligned ("offline") responses
# online_alignment <- c(T, F)
exclude_first_run = T
reg_diagnostics = F


setwd('~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/')

# read in behavioral data
cache_dir <- "~/Box/SCEPTIC_fMRI/dan_medusa/cache/"
repo_dir <- "~/code/clock_analysis"
load(file.path(repo_dir, '/fmri/keuka_brain_behavior_analyses/trial_df_and_vh_pe_clusters_u.Rdata'))

# select relevant columns for compactness
df <- df %>% select(id, run, run_trial, rewFunc,emotion, last_outcome, rt_csv, score_csv, rt_next, pe_max, rt_vmax, rt_vmax_lag,
                    rt_vmax_change, v_max_wi, v_entropy_wi, v_entropy_b, v_entropy, v_max_b, u_chosen_quantile, u_chosen_quantile_lag, u_chosen_quantile_change, 
                    rt_vmax_lag_sc, rt_lag_sc,rt_lag2_sc, rt_csv_sc, trial_neg_inv_sc, Age, Female, kld3, kld4)

if (online) {
  d <- merge(df, clock_wide, by = c("id", "run", "run_trial"))
} else { d <- merge(df, rt_wide, by = c("id", "run", "run_trial"))}

d <- d %>% group_by(id, run) %>% arrange(id, run, run_trial) %>% 
  mutate(rt_next = lead(rt_csv_sc),
         rt_change = rt_next - rt_csv_sc,
         rt_vmax_lead = lead(rt_vmax),
         rt_vmax_change_next = rt_vmax_lead - rt_vmax,
         v_entropy_wi_lead = lead(v_entropy_wi),
         v_entropy_wi_change = v_entropy_wi_lead-v_entropy_wi,
         v_entropy_wi_change_lag = lag(v_entropy_wi_change),
         u_chosen_quantile_next = lead(u_chosen_quantile),
         u_chosen_quantile_change_next = lead(u_chosen_quantile_change),
         kld3_lead = lead(kld3),
         kld3_lag = lag(kld3),
         outcome = case_when(
           score_csv>0 ~ 'Reward',
           score_csv==0 ~ "Omission"),
         abs_pe = abs(pe_max),
         abs_pe_lag = lag(abs_pe)
  ) %>% ungroup()
# remove first run
if (exclude_first_run) {
  d <- d %>% filter(run>1)
}

# diagnose regressor multicollinearity
if (reg_diagnostics) {
  regs <- d %>% select(rt_csv, rt_lag_sc, rt_vmax, rt_vmax_lag, rt_vmax_change, v_entropy_wi, v_entropy_wi_change, v_max_wi,
                       kld3, kld3_lag, abs_pe, score_csv, trial_neg_inv_sc)
  cormat <- corr.test(regs)
  corrplot(cormat$r, cl.lim=c(-1,1),
           method = "circle", tl.cex = 1.5, type = "upper", tl.col = 'black',
           order = "hclust", diag = FALSE,
           addCoef.col="black", addCoefasPercent = FALSE,
           p.mat = cormat$p, sig.level=0.05, insig = "blank")
}

scale2 <- function(x, na.rm = FALSE) (x - mean(x, na.rm = na.rm)) / sd(x, na.rm)
# scale decon across subjects as a predictor
# choice uncertainty prediction analyses run on scaled 'ds' dataframe instead of 'd'
# ds <- d %>% mutate_at(vars(starts_with("dan")), scale2, na.rm = TRUE) %>% ungroup()

## "Decoding" ----
# combined right and left hippocampus with side as a predictor
# if model does not converge, update with new starting values (not needed here)
labels <- names(d[grepl("_R_|_r_|_L_|_l_", names(d))])

# make cluster ----
f <- Sys.getenv('PBS_NODEFILE')
library(parallel)
ncores <- detectCores()
nodelist <- if (nzchar(f)) readLines(f) else rep('localhost', ncores)

cat("Node list allocated to this job\n")
print(nodelist)

cl <- makePSOCKcluster(nodelist, outfile='')
print(cl) ##; print(unclass(cl))
registerDoParallel(cl)
# loop over sensors ----
pb <- txtProgressBar(0, max = length(labels), style = 3)

# test
# labels <- labels[1:2]

if(decode) {
  message("\nDecoding: analyzing parcel data")
  ddf <- foreach(i = 1:length(labels), .packages=c("lme4", "tidyverse", "broom.mixed", "car"),
                 .combine='rbind', .noexport = c("clock_wide", "clock_wide_cens", "rt_wide")) %dopar% {
                   # message(paste("Analyzing timepoint", t,  sep = " "))
                   if (i %% 10 == 0) {setTxtProgressBar(pb, i)}
                   label <- as.character(labels[[i]])
                   d$h <- as.numeric(d[[label]])
                   if (online) {
                     md <-  lmerTest::lmer(h ~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + scale(rt_vmax_lag) + scale(rt_vmax_change) + 
                                             v_entropy_wi + v_entropy_wi_change_lag + v_entropy_wi_change  +
                                             kld3_lag  + v_max_wi  + scale(abs_pe_lag) + last_outcome + 
                                             (1|id), d, control=lmerControl(optimizer = "nloptwrap"))
                   } else {
                     md <-  lmerTest::lmer(h ~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + scale(rt_vmax_lag)  + scale(rt_vmax_change) + 
                                             v_entropy_wi + v_entropy_wi_change_lag + v_entropy_wi_change  + 
                                             kld3_lag  + v_max_wi  + scale(abs_pe) + outcome + 
                                             (1|id), d, control=lmerControl(optimizer = "nloptwrap")) }
                   while (any(grepl("failed to converge", md@optinfo$conv$lme4$messages) )) {
                     print(md@optinfo$conv$lme4$conv)
                     ss <- getME(md,c("theta","fixef"))
                     md <- update(md, start=ss)}
                   
                   dm <- tidy(md)
                   dm$label <- label
                   dm$t <- gsub(".*_", "\\1", label)
                   dm}
  # FDR correction ----
  message("\nFDR correction")
  ddf <- ddf %>% mutate(stat_order = as.factor(case_when(abs(statistic) < 2 ~ '1', 
                                                         abs(statistic) > 2 & abs(statistic) < 3 ~ '2', 
                                                         abs(statistic) > 3 ~ '3')),
                        p_value = as.factor(case_when(`p.value` > .05 ~ '1',
                                                      `p.value` < .05 & `p.value` > .01 ~ '2',
                                                      `p.value` < .01 & `p.value` > .001 ~ '3',
                                                      `p.value` <.001 ~ '4')))
  ddf$t <- as.numeric(ddf$t)
  ddf$label <- as.factor(sub("_[^_]+$", "", ddf$label))
  ddf$stat_order <- factor(ddf$stat_order, labels = c("NS", "|t| > 2", "|t| > 3"))
  ddf$p_value <- factor(ddf$p_value, labels = c("NS", "p < .05", "p < .01", "p < .001"))
  terms <- unique(ddf$term[ddf$effect=="fixed"])
  ddf <- ddf  %>% group_by(term) %>% mutate(p_fdr = p.adjust(p.value, method = 'fdr'),
                                            p_level_fdr = as.factor(case_when(
                                              # p_fdr > .1 ~ '0',
                                              # p_fdr < .1 & p_fdr > .05 ~ '1',
                                              p_fdr > .05 ~ '1',
                                              p_fdr < .05 & p_fdr > .01 ~ '2',
                                              p_fdr < .01 & p_fdr > .001 ~ '3',
                                              p_fdr <.001 ~ '4'))
  ) %>% ungroup() %>% mutate(side = substr(as.character(label), nchar(as.character(label)), nchar(as.character(label))),
                             region = substr(as.character(label), 1, nchar(as.character(label))-2))
  ddf$p_level_fdr <- factor(ddf$p_level_fdr, levels = c('1', '2', '3', '4'), labels = c("NS","p < .05", "p < .01", "p < .001"))
  ddf$`p, FDR-corrected` = ddf$p_level_fdr
  
  # plots ----
  message("\nPlotting")
  if (online) {
    setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/dan/plots/clock_decode')
    epoch_label = "Time relative to clock onset, seconds"
    decode_results_fname = "clock_decode_output.Rdata"
  } else {setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/dan/plots/rt_decode')
    epoch_label = "Time relative to outcome, seconds"
    decode_results_fname = "rt_decode_output.Rdata"}
  for (fe in terms) {
    edf <- ddf %>% filter(term == paste(fe) & t < 8) 
    termstr <- str_replace_all(fe, "[^[:alnum:]]", "_")
    fname = paste(termstr, ".pdf", sep = "")
    pdf(fname, width = 11, height = 6)
    print(ggplot(edf, aes(t, region)) + geom_tile(aes(fill = estimate, alpha = `p, FDR-corrected`), size = 1) +  
            geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + facet_wrap(~side) +
            scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) + ylab("Parcel") + 
            labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)))
    dev.off()
    # save model stats ----
  }
  # add labels
  message("\nLabelling results")
  all_labels <- as_tibble(read_excel("~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/MNH Dan Labels.xlsx")) %>%
    select(c("roinum", "plot_label", "Stream", "Visuomotor_Gradient", "Stream_Gradient"))
  names(all_labels) <- c("atlas_value","label_short", "stream", "visuomotor_grad", "stream_grad")
  all_labels$stream_grad <- as.numeric(all_labels$stream_grad)
  all_labels <- all_labels %>% arrange(visuomotor_grad, stream_grad) %>% mutate(
    side  = case_when(
      grepl("L_", label_short) ~ "L",
      grepl("R_", label_short) ~ "R"),
    label_short = substr(label_short, 3, length(label_short)),
    label = paste(visuomotor_grad, stream_grad, label_short, side, sep = "_")) %>% 
    select(c(label, label_short, side, atlas_value, stream, visuomotor_grad, stream_grad))
  ddf <- merge(ddf, all_labels)
  message("\nSaving results")
  save(file = decode_results_fname, ddf)
  gc()
}



## RT prediction ----

if(rt_predict) {
  message("\nRT prediction: analyzing parcel data")
  rdf <- foreach(i = 1:length(labels), .packages=c("lme4", "tidyverse", "broom.mixed", "car"),
                 .combine='rbind', .noexport = c("clock_wide", "clock_wide_cens", "rt_wide")) %dopar% {
                   # message(paste("Analyzing timepoint", t,  sep = " "))
                   if (i %% 10 == 0) {setTxtProgressBar(pb, i)}
                   label <- as.character(labels[[i]])
                   d$h <- as.numeric(d[[label]])
                   if (online) {
                     md <-  lmerTest::lmer(scale(rt_next) ~ scale(h) * scale(rt_vmax)  + 
                                   scale(h) * rt_csv_sc * last_outcome + scale(h) * rt_lag_sc + 
                                   (1|id), d, control=lmerControl(optimizer = "nloptwrap"))
                   } else {
                     md <-  lmerTest::lmer(scale(rt_next) ~ scale(h) * rt_csv_sc * outcome  + scale(h) * scale(rt_vmax)  +
                                   scale(h) * rt_lag_sc + 
                                   (1|id), d, control=lmerControl(optimizer = "nloptwrap"))}                   
                   while (any(grepl("failed to converge", md@optinfo$conv$lme4$messages) )) {
                     print(md@optinfo$conv$lme4$conv)
                     ss <- getME(md,c("theta","fixef"))
                     md <- update(md, start=ss)}
                   
                   dm <- tidy(md)
                   dm$label <- label
                   dm$t <- gsub(".*_", "\\1", label)
                   dm}
  # FDR correction ----
  message("\nFDR correction")
  rdf <- rdf %>% mutate(stat_order = as.factor(case_when(abs(statistic) < 2 ~ '1', 
                                                         abs(statistic) > 2 & abs(statistic) < 3 ~ '2', 
                                                         abs(statistic) > 3 ~ '3')),
                        p_value = as.factor(case_when(`p.value` > .05 ~ '1',
                                                      `p.value` < .05 & `p.value` > .01 ~ '2',
                                                      `p.value` < .01 & `p.value` > .001 ~ '3',
                                                      `p.value` <.001 ~ '4')))
  rdf$t <- as.numeric(rdf$t)
  rdf$label <- as.factor(sub("_[^_]+$", "", rdf$label))
  rdf$stat_order <- factor(rdf$stat_order, labels = c("NS", "|t| > 2", "|t| > 3"))
  rdf$p_value <- factor(rdf$p_value, labels = c("NS", "p < .05", "p < .01", "p < .001"))
  terms <- unique(rdf$term[rdf$effect=="fixed"])
  terms <- terms[grepl("(h)",terms)]
  rdf <- rdf  %>% group_by(term) %>% mutate(p_fdr = p.adjust(p.value, method = 'fdr'),
                                            p_level_fdr = as.factor(case_when(
                                              # p_fdr > .1 ~ '0',
                                              # p_fdr < .1 & p_fdr > .05 ~ '1',
                                              p_fdr > .05 ~ '1',
                                              p_fdr < .05 & p_fdr > .01 ~ '2',
                                              p_fdr < .01 & p_fdr > .001 ~ '3',
                                              p_fdr <.001 ~ '4'))
  ) %>% ungroup() %>% mutate(side = substr(as.character(label), nchar(as.character(label)), nchar(as.character(label))),
                             region = substr(as.character(label), 1, nchar(as.character(label))-2))
  rdf$p_level_fdr <- factor(rdf$p_level_fdr, levels = c('1', '2', '3', '4'), labels = c("NS","p < .05", "p < .01", "p < .001"))
  rdf$`p, FDR-corrected` = rdf$p_level_fdr
  
  # plots ----
  message("\nPlotting")
  if (online) {
    setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/dan/plots/clock_rt')
    epoch_label = "Time relative to clock onset, seconds"  
    rt_results_fname = "clock_rt_predict_output.Rdata"
  } else {
    setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/dan/plots/rt_rt')
    epoch_label = "Time relative to outcome, seconds"
    rt_results_fname = "rt_rt_predict_output.Rdata"
  }
  for (fe in terms) {
    edf <- rdf %>% filter(term == paste(fe) & t < 8) 
    termstr <- str_replace_all(fe, "[^[:alnum:]]", "_")
    fname = paste(termstr, ".pdf", sep = "")
    pdf(fname, width = 11, height = 6)
    print(ggplot(edf, aes(t, region)) + geom_tile(aes(fill = estimate, alpha = `p, FDR-corrected`), size = 1) +  
            geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + facet_wrap(~side) +
            scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) + ylab("Parcel") + 
            labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)))
    dev.off()
    # save output for inspection
    
  }
  # add labels
  message("\nLabelling results")
  all_labels <- as_tibble(read_excel("~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/MNH Dan Labels.xlsx")) %>%
    select(c("roinum", "plot_label", "Stream", "Visuomotor_Gradient", "Stream_Gradient"))
  names(all_labels) <- c("atlas_value","label_short", "stream", "visuomotor_grad", "stream_grad")
  all_labels$stream_grad <- as.numeric(all_labels$stream_grad)
  all_labels <- all_labels %>% arrange(visuomotor_grad, stream_grad) %>% mutate(
    side  = case_when(
      grepl("L_", label_short) ~ "L",
      grepl("R_", label_short) ~ "R"),
    label_short = substr(label_short, 3, length(label_short)),
    label = paste(visuomotor_grad, stream_grad, label_short, side, sep = "_")) %>% 
    select(c(label, label_short, side, atlas_value, stream, visuomotor_grad, stream_grad))
  rdf <- merge(rdf, all_labels)
  message("\nSaving results")
  save(file = rt_results_fname, rdf)
}
stopCluster(cl)
gc()



