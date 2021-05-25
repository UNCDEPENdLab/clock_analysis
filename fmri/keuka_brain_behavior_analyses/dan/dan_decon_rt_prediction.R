# with 'decode = T' makes MEDUSA decoding plots for Fig. 4 E-G.
# loops over decoding and uncertainty prediction multi-level models for various hippocampal slices and post-feedback time points
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
repo_directory <- "~/code/clock_analysis"

# data & options ----

# data loading options
reprocess = F # otherwise load data from cache
if (!reprocess) {
  wide_only = T # only load wide data (parcels and timepoints as variables)
}
# load MEDUSA deconvolved data
source(file.path(repo_directory, "fmri/keuka_brain_behavior_analyses/dan/load_medusa_data_dan.R"))

# what to run
plots = T
decode = F  # main analysis analogous to Fig. 4 E-G in NComm 2020
rt_predict = T # predicts next response based on signal and behavioral variables
online = T # whether to analyze clock-aligned ("online") or RT-aligned ("offline") responses
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
if (decode) {
  newlist <- list()
  for (label in labels) {print(paste("Processing parcel", label,  sep = " "))
    d$h<-d[[label]]
    if (online) {
      md <-  lmer(h ~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + scale(rt_vmax_lag) + 
                    v_entropy_wi + v_entropy_wi_change  + v_max_wi  + 
                    kld3_lag  + v_max_wi  + scale(abs_pe_lag) + last_outcome + 
                    (1|id), d, control=lmerControl(optimizer = "nloptwrap"))
    } else {
    md <-  lmer(h ~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + scale(rt_vmax_lag)  + scale(rt_vmax_change) + 
                  v_entropy_wi + v_entropy_wi_change  + v_max_wi  + 
                  kld3_lag  + v_max_wi  + scale(abs_pe) + outcome + 
                  (1|id), d, control=lmerControl(optimizer = "nloptwrap")) }
    while (any(grepl("failed to converge", md@optinfo$conv$lme4$messages) )) {
      print(md@optinfo$conv$lme4$conv)
      ss <- getME(md,c("theta","fixef"))
      md <- update(md, start=ss)}
    
    dm <- tidy(md)
    dm$label <- label
    dm$t <- gsub(".*_", "\\1", label)
    dm <- dm %>% mutate(stat_order = as.factor(case_when(abs(statistic) < 2 ~ '1', 
                                                         abs(statistic) > 2 & abs(statistic) < 3 ~ '2', 
                                                         abs(statistic) > 3 ~ '3')),
                        p_value = as.factor(case_when(p.value > .05 ~ '1',
                                                      p.value < .05 & p.value > .01 ~ '2',
                                                      p.value < .01 & p.value > .001 ~ '3',
                                                      p.value <.001 ~ '4'))
    )
    newlist[[label]]<-dm
  }
  
  ddf <- do.call(rbind,newlist)
  ddf$t <- as.numeric(ddf$t)
  ddf$label <- as.factor(sub("_[^_]+$", "", ddf$label))
  ddf$stat_order <- factor(ddf$stat_order, labels = c("NS", "|t| > 2", "|t| > 3"))
  ddf$p_value <- factor(ddf$p_value, labels = c("NS", "p < .05", "p < .01", "p < .001"))
  
  terms <- names(fixef(md))
  # FDR correction ----
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

    if (online) {
    setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/dan/plots/clock_decode')
          epoch_label = "Time relative to clock onset, seconds"
  } else {setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/dan/plots/rt_decode')
  epoch_label = "Time relative to outcome, seconds"}
  
  for (fe in terms) {
    edf <- ddf %>% filter(term == paste(fe) & t < 8) 
    termstr <- str_replace_all(fe, "[^[:alnum:]]", "_")
    pdf(paste(termstr, ".pdf", sep = ""), width = 11, height = 6)
    print(ggplot(edf, aes(t, region)) + geom_tile(aes(fill = estimate, alpha = `p, FDR-corrected`), size = 1) +  
            geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + facet_wrap(~side) +
            scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) + ylab("Parcel") + 
            labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)))
    dev.off()
    # save model stats ----
    save(file = "medusa_decode_output.Rdata", ddf)
  }
}

## RT prediction ----

if (rt_predict) {
  labels <- names(d[grepl("_R_|_r_|_L_|_l_", names(d))])
  newlist <- list()
  for (label in labels) {print(paste("Processing parcel", label,  sep = " "))
    d$h<-d[[label]]
    if (online) {
      md <-  lmer(scale(rt_next) ~ scale(h) * scale(rt_vmax)  + 
                    scale(h) * rt_csv_sc * last_outcome + scale(h) * rt_lag_sc + 
                    (1|id), d, control=lmerControl(optimizer = "nloptwrap"))
    } else {
    md <-  lmer(scale(rt_next) ~ scale(h) * rt_csv_sc * outcome  + scale(h) * scale(rt_vmax)  +
                  scale(h) * rt_lag_sc + 
                  (1|id), d, control=lmerControl(optimizer = "nloptwrap"))}
    
    # while (any(grepl("failed to converge", md@optinfo$conv$lme4$messages) )) {
    #   print(md@optinfo$conv$lme4$messages)
    #   ss <- getME(md,c("theta","fixef"))
    #   md <- update(md, start=ss)}
    
    dm <- tidy(md)
    dm$label <- label
    dm$t <- gsub(".*_", "\\1", label)
    dm <- dm %>% mutate(stat_order = as.factor(case_when(abs(statistic) < 2 ~ '1', 
                                                         abs(statistic) > 2 & abs(statistic) < 3 ~ '2', 
                                                         abs(statistic) > 3 ~ '3')),
                        p_value = as.factor(case_when(p.value > .05 ~ '1',
                                                      p.value < .05 & p.value > .01 ~ '2',
                                                      p.value < .01 & p.value > .001 ~ '3',
                                                      p.value <.001 ~ '4'))
    )
    newlist[[label]]<-dm
    # }
  }
  ddf <- do.call(rbind,newlist)
  ddf$t <- as.numeric(ddf$t)
  ddf$label <-  as.factor(sub("_[^_]+$", "", ddf$label))
  ddf$stat_order <- factor(ddf$stat_order, labels = c("NS", "|t| > 2", "|t| > 3"))
  ddf$p_value <- factor(ddf$p_value, labels = c("NS", "p < .05", "p < .01", "p < .001"))
  
  terms <- names(fixef(md))
  terms <- terms[grepl("(h)",terms)]
  
  # FDR correction within term, across regions, timepoints and both hemispheres ----
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
  #,
  # side_long = case_when(side=='l' ~ 'Left',
  #                       side=='r' ~ 'Right')
  
  ddf$p_level_fdr <- factor(ddf$p_level_fdr, levels = c('1', '2', '3', '4'), labels = c("NS","p < .05", "p < .01", "p < .001"))
  
  ddf$`p, FDR-corrected` = ddf$p_level_fdr
  
  # plots ----
  if (online) {
    setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/dan/plots/clock_rt')
    epoch_label = "Time relative to clock onset, seconds"  
  } else {
    setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/dan/plots/rt_rt')
    epoch_label = "Time relative to outcome, seconds"
  }
  for (fe in terms) {
    edf <- ddf %>% filter(term == paste(fe) & t < 8) 
    termstr <- str_replace_all(fe, "[^[:alnum:]]", "_")
    pdf(paste(termstr, ".pdf", sep = ""), width = 11, height = 6)
    print(ggplot(edf, aes(t, region)) + geom_tile(aes(fill = estimate, alpha = `p, FDR-corrected`), size = 1) +  
            geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + facet_wrap(~side) +
            scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) + ylab("Parcel") + 
            labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)))
    dev.off()
    # save output for inspection
    save(file = "medusa_rt_predict_output.Rdata", ddf)
  }
}


