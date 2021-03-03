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
library(foreach)
library(doParallel)
repo_directory <- "~/code/clock_analysis"

# data & options ----

# data loading options
reprocess = F # otherwise load data from cache
if (!reprocess) {
  wide_only = T # only load wide data (parcels and timepoints as variables)
}
replicate_compression = F
if(replicate_compression) {reprocess = T}
# load MEDUSA deconvolved data
source(file.path(repo_directory, "meg/code/medusoid/meg_load_medusa_data.R"))

# what to run
plots = T
decode = T  # main analysis analogous to Fig. 4 E-G in NComm 2020
rt_predict = F # predicts next response based on signal and behavioral variables
online = F # whether to analyze clock-aligned ("online") or RT-aligned ("offline") responses
exclude_first_run = T
reg_diagnostics = F


setwd('~/code/clock_analysis/meg/code/medusoid')

cache_dir <- "~/Box/SCEPTIC_fMRI/MEG_20Hz/cache/"
repo_dir <- "~/code/clock_analysis"
# load(file.path(repo_dir, '/fmri/keuka_brain_behavior_analyses/trial_df_and_vh_pe_clusters_u.Rdata'))

# select relevant columns for compactness
df <- trial_df %>% ungroup() %>%
  mutate(rt_lag_sc = scale(rt_lag),
         rt_csv_sc = scale(rt_csv),
         rt_vmax_lag_sc = scale(rt_vmax_lag)) %>%
  group_by(id, run) %>% arrange(id, run, run_trial) %>% 
  mutate(rt_next = lead(rt_csv_sc),
         rt_change = rt_next - rt_csv_sc,
         rt_vmax_lead = lead(rt_vmax),
         rt_vmax_change_next = rt_vmax_lead - rt_vmax,
         v_entropy_wi_lead = lead(v_entropy_wi),
         v_entropy_wi_change = v_entropy_wi_lead-v_entropy_wi,
         v_entropy_wi_change_lag = lag(v_entropy_wi_change),
         # kld3_lead = lead(kld3),
         # kld3_lag = lag(kld3),
         outcome = case_when(
           score_csv>0 ~ 'Reward',
           score_csv==0 ~ "Omission"),
         abs_pe = abs(pe_max),
         abs_pe_lag = lag(abs_pe),
         v_max_wi = scale(v_max),
         v_max_wi_lag = lag(v_max_wi),
         v_entropy_wi = scale(v_entropy),
         v_max_b = mean(na.omit(v_max)),
         v_entropy_b = mean(na.omit(v_entropy)),
         rt_change = rt_csv - rt_lag,
         pe_max_lag = lag(pe_max), 
         abs_pe_max_lag = abs(pe_max_lag), 
         rt_vmax_change = rt_vmax - rt_vmax_lag,
         v_chosen_change = v_chosen - lag(v_chosen),
         trial_neg_inv_sc = scale(-1/run_trial),
         rt_lag2_sc = lag(rt_csv_sc, 2),
         rt_lag3_sc = lag(rt_csv_sc, 3),) %>% ungroup() 

if (online) {
  d <- merge(df, clock_wide, by = c("id", "run", "run_trial"))
} else { d <- merge(df, rt_wide, by = c("id", "run", "run_trial"))}

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
# scale signal within subjects as a predictor
d <- d %>% mutate_at(vars(contains(".")), scale2, na.rm = TRUE) %>% ungroup()

# make cluster
f <- Sys.getenv('PBS_NODEFILE')
library(parallel)
ncores <- detectCores()
nodelist <- if (nzchar(f)) readLines(f) else rep('localhost', ncores)

cat("Node list allocated to this job\n")
print(nodelist)

cl <- makePSOCKcluster(nodelist, outfile='')
print(cl) ##; print(unclass(cl))
registerDoParallel(cl)


## "Decoding" ----
# combined right and left hippocampus with side as a predictor
# if model does not converge, update with new starting values (not needed here)

labels <- names(d[grepl("\\.", names(d))])
# labels <- labels[1:10]
if (decode) {
  # newlist <- list()
  # for (label in labels) {print(paste("Processing parcel", label,  sep = " "))
  ddf <- foreach(i = 1:length(labels), .packages=c("lme4", "tidyverse", "broom.mixed", "car"), 
                 .combine='rbind') %dopar% {
                   label <- labels[[i]]
                   d$h<-d[[label]]
                   if (online) {
                     md <-  lmer(h ~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + scale(rt_vmax_lag) + scale(rt_vmax_change) + 
                                   v_entropy_wi + v_entropy_wi_lead +  v_entropy_wi_change_lag + #v_entropy_wi_change  +
                                   v_max_wi  + scale(abs_pe_lag) + last_outcome + 
                                   (1|id), d, control=lmerControl(optimizer = "nloptwrap"))
                   } else {
                     md <-  lmer(h ~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + scale(rt_vmax_lag)  + scale(rt_vmax_change) + 
                                   v_entropy_wi + v_entropy_wi_lead + v_entropy_wi_change_lag + #v_entropy_wi_change  + 
                                   v_max_wi  + scale(abs_pe) + outcome + 
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
                   # newlist[[label]]<-dm
                   dm
                 }
  
  ddf <- do.call(rbind,newlist)
  ddf$t <- as.numeric(ddf$t)
  ddf$label <- as.factor(sub("_[^_]+$", "", ddf$label))
  ddf$stat_order <- factor(ddf$stat_order, labels = c("NS", "|t| > 2", "|t| > 3"))
  ddf$p_value <- factor(ddf$p_value, labels = c("NS", "p < .05", "p < .01", "p < .001"))
  
  terms <- names(fixef(md))
  # FDR correction ----
  ddf <- ddf  %>% filter(t>-2) %>% group_by(term) %>% mutate(p_fdr = p.adjust(p.value, method = 'fdr'),
                                                             p_level_fdr = as.factor(case_when(
                                                               # p_fdr > .1 ~ '0',
                                                               # p_fdr < .1 & p_fdr > .05 ~ '1',
                                                               p_fdr > .05 ~ '1',
                                                               p_fdr < .05 & p_fdr > .01 ~ '2',
                                                               p_fdr < .01 & p_fdr > .001 ~ '3',
                                                               p_fdr <.001 ~ '4'))
  ) %>% ungroup() #%>% mutate(side = substr(as.character(label), nchar(as.character(label)), nchar(as.character(label))),
  #          region = substr(as.character(label), 1, nchar(as.character(label))-2))
  ddf$p_level_fdr <- factor(ddf$p_level_fdr, levels = c('1', '2', '3', '4'), labels = c("NS","p < .05", "p < .01", "p < .001"))
  ddf$`p, FDR-corrected` = ddf$p_level_fdr
  
  # plots ----
  
  if (online) {
    setwd('~/OneDrive/collected_letters/papers/meg/plots/clock_decode')
    epoch_label = "Time relative to clock onset, seconds"
  } else {setwd('~/OneDrive/collected_letters/papers/meg/plots/rt_decode')
    epoch_label = "Time relative to outcome, seconds"}
  
  for (fe in terms) {
    edf <- ddf %>% filter(term == paste(fe)) 
    termstr <- str_replace_all(fe, "[^[:alnum:]]", "_")
    fname = paste("meg_", termstr, ".pdf", sep = "")
    if (replicate_compression){fname = paste(termstr,"_replicate_compression", ".pdf", sep = "")}
    pdf(fname, width = 16, height = 3)
    print(ggplot(edf, aes(t, label)) + geom_tile(aes(fill = estimate, alpha = p_value)) + 
            geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
            scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) + ylab("Sensor") +
            labs(alpha = expression(italic(p)[uncorrected])) + ggtitle(paste(termstr)))
    dev.off()
  }
  # save model stats ----
  save(file = "medusa_decode_output.Rdata", ddf)
  
}
# print(ggplot(edf, aes(t, label)) + geom_tile(aes(fill = estimate, alpha = p_value, size = 1)) +  
#         # print(ggplot(edf, aes(t, label)) + geom_tile(aes(fill = estimate, alpha = `p, FDR-corrected`), size = 1) +  
#         geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + 
#         scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) + ylab("Sensor") #+ labs(alpha = expression(italic(p)[uncorrected])) + ggtitle(paste(termstr))
#       # labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)))
#       dev.off()

## RT prediction ----

if (rt_predict) {
  labels <- names(d[grepl("\\.", names(d))])
  # newlist <- list()
  # for (label in labels) {print(paste("Processing parcel", label,  sep = " "))
  ddf <- foreach(i = 1:length(labels), .packages=c("lme4", "tidyverse", "broom.mixed", "car"), 
                 .combine='rbind') %dopar% {
                   label <- labels[[i]]
                   
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
                   # newlist[[label]]<-dm
                   # }
                 dm}
  # ddf <- do.call(rbind,newlist)
  ddf$t <- as.numeric(ddf$t)
  ddf$label <-  as.factor(sub("_[^_]+$", "", ddf$label))
  ddf$stat_order <- factor(ddf$stat_order, labels = c("NS", "|t| > 2", "|t| > 3"))
  ddf$p_value <- factor(ddf$p_value, labels = c("NS", "p < .05", "p < .01", "p < .001"))
  
  terms <- names(fixef(md))
  terms <- terms[grepl("(h)",terms)]
  
  # FDR correction ----
  ddf <- ddf  %>% filter(t>-2) %>% group_by(term) %>% mutate(p_fdr = p.adjust(p.value, method = 'fdr'),
                                                             p_level_fdr = as.factor(case_when(
                                                               # p_fdr > .1 ~ '0',
                                                               # p_fdr < .1 & p_fdr > .05 ~ '1',
                                                               p_fdr > .05 ~ '1',
                                                               p_fdr < .05 & p_fdr > .01 ~ '2',
                                                               p_fdr < .01 & p_fdr > .001 ~ '3',
                                                               p_fdr <.001 ~ '4'))
  ) %>% ungroup() #%>% mutate(side = substr(as.character(label), nchar(as.character(label)), nchar(as.character(label))),
  #          region = substr(as.character(label), 1, nchar(as.character(label))-2))
  ddf$p_level_fdr <- factor(ddf$p_level_fdr, levels = c('1', '2', '3', '4'), labels = c("NS","p < .05", "p < .01", "p < .001"))
  ddf$`p, FDR-corrected` = ddf$p_level_fdr
  
  # plots ----
  
  if (online) {
    setwd('~/OneDrive/collected_letters/papers/meg/plots/clock_rt')
    epoch_label = "Time relative to clock onset, seconds"
  } else {setwd('~/OneDrive/collected_letters/papers/meg/plots/rt_rt')
    epoch_label = "Time relative to outcome, seconds"}
  
  for (fe in terms) {
    edf <- ddf %>% filter(term == paste(fe)) 
    termstr <- str_replace_all(fe, "[^[:alnum:]]", "_")
    fname = paste("meg_", termstr, ".pdf", sep = "")
    if (replicate_compression){fname = paste(termstr,"_replicate_compression", ".pdf", sep = "")}
    pdf(fname, width = 16, height = 3)
    print(ggplot(edf, aes(t, label)) + geom_tile(aes(fill = estimate, alpha = p_value)) + 
            geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) +
            scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) + ylab("Sensor") +
            labs(alpha = expression(italic(p)[uncorrected])) + ggtitle(paste(termstr)))
    dev.off()
    # save output for inspection
    save(file = "medusa_rt_predict_output.Rdata", ddf)
  }
}
stopCluster(cl)


