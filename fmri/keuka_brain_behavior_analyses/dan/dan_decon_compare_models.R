# with 'decode = T' makes MEDUSA decoding plots for Fig. 4 E-G.
# loops over decoding and uncertainty prediction multi-level models for various hippocampal slices and post-feedback time points
# first run medusa_event_locked_lmer.R

library(modelr)
library(tidyverse)
library(lme4)
# library(afex)
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
replicate_compression = F
if(replicate_compression) {reprocess = T}
# load MEDUSA deconvolved data
source(file.path(repo_directory, "fmri/keuka_brain_behavior_analyses/dan/load_medusa_data_dan.R"))

# what to run
plots = T
decode = T  # main analysis analogous to Fig. 4 E-G in NComm 2020
rt_predict = T # predicts next response based on signal and behavioral variables
online = F # whether to analyze clock-aligned ("online") or RT-aligned ("offline") responses
exclude_first_run = T
reg_diagnostics = F


setwd('~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/')

# read in behavioral data
cache_dir <- "~/Box/SCEPTIC_fMRI/dan_medusa/cache/"
repo_dir <- "~/code/clock_analysis"
# load(file.path(repo_dir, '/fmri/keuka_brain_behavior_analyses/trial_df_and_vh_pe_clusters_u.Rdata'))
load(file.path(cache_dir, 'sceptic_trial_df_for_medusa.RData')) 
df <- trial_df
# select relevant columns for compactness
df <- df %>% select(id, run, run_trial, rewFunc,emotion, rt_csv, score_csv, rt_next, pe_max, rt_vmax, rt_vmax_lag,
                    rt_vmax_change, v_entropy_wi, v_entropy_b, v_entropy, v_max_b,
                    rt_vmax_lag_sc, rt_lag_sc,rt_lag2_sc, rt_csv_sc, trial_neg_inv_sc, Age, Female, v_entropy_wi,v_entropy_wi_change, v_entropy_wi_full,
                    v_entropy_wi_change_full, v_entropy_wi_full, rt_vmax_full, rt_vmax_change_full, pe_max_full)

if (online) {
  d <- merge(df, clock_wide, by = c("id", "run", "run_trial"))
} else { d <- merge(df, rt_wide, by = c("id", "run", "run_trial"))}

d <- d %>% group_by(id, run) %>% arrange(id, run, run_trial) %>% 
  mutate(pe_max_full_lag = lag(pe_max_full), 
         outcome = case_when(
           score_csv>0 ~ 'Reward',
           score_csv==0 ~ "Omission"),
         last_outcome = lag(outcome), 
         trial_neg_inv_sc  = scale(1/run_trial),
         rt_csv_sc = scale(rt_csv),
         rt_lag_sc = lag(rt_csv_sc),
         pe_max_lag_full = lag(pe_max_full)) %>% ungroup()
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
# labels <- labels[1:5] # to test
if (decode) {
  newlist <- list()
  for (label in labels) {print(paste("Processing parcel", label,  sep = " "))
    d$h<-d[[label]]
    if (online) {
      md <-  lme4::lmer(h ~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + scale(rt_vmax_lag) + scale(rt_vmax_change) + 
                          v_entropy_wi + v_entropy_wi_lead +  v_entropy_wi_change_lag + #v_entropy_wi_change  +
                          v_max_wi  + scale(abs_pe_lag) + last_outcome + 
                          (1|id), d, control=lmerControl(optimizer = "nloptwrap"))
    } else {
      ## selective
      # entropy 
      md_h <-  lmerTest::lmer(h ~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + 
                                v_entropy_wi  + 
                                outcome + 
                                (1|id), d)
      # entropy change 
      md_h_change <-  lme4::lmer(h ~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + 
                                   v_entropy_wi_change  + 
                                   outcome + 
                                   (1|id), d, control=lmerControl(optimizer = "nloptwrap"))
      # rt_vmax 
      md_h_rtvmax <-  lme4::lmer(h ~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + 
                                   rt_vmax_lag  + 
                                   outcome + 
                                   (1|id), d %>% filter(!is.na(d$rt_vmax_lag)), control=lmerControl(optimizer = "nloptwrap"))
      # pe 
      md_h_pe <-  lme4::lmer(h ~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + 
                               abs(pe_max)  + 
                               outcome + 
                               (1|id), d, control=lmerControl(optimizer = "nloptwrap"))
      
      ## full
      # entropy 
      mdf_h <-  lmerTest::lmer(h ~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + 
                                 v_entropy_wi_full + outcome + 
                                 (1|id), d) 
      # entropy change 
      mdf_h_change <-  lme4::lmer(h ~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + 
                                    v_entropy_wi_change_full  + 
                                    outcome + 
                                    (1|id), d, control=lmerControl(optimizer = "nloptwrap"))
      # rt_vmax 
      mdf_h_rtvmax <-  lme4::lmer(h ~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + 
                                    rt_vmax_lag_full  + 
                                    outcome + 
                                    (1|id), d %>% filter(!is.na(d$rt_vmax_lag)), control=lmerControl(optimizer = "nloptwrap"))
      # pe
      mdf_h_pe <-  lme4::lmer(h ~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + 
                                abs(pe_max_full)  + 
                                outcome + 
                                (1|id), d, control=lmerControl(optimizer = "nloptwrap"))
      
      effect <- c("entropy", "entropy.change", "RT.Vmax", "abs(PE)")
      AIC.diff <- c(AIC(mdf_h) - AIC(md_h), AIC(mdf_h_change) - AIC(md_h_change), 
                    AIC(mdf_h_rtvmax) - AIC(md_h_rtvmax), AIC(mdf_h_pe) - AIC(md_h_pe))
      dm <- as_tibble(cbind(effect, AIC.diff))
      dm$label <- label
      dm$t <- gsub(".*_", "\\1", label)
      # Akaike weights following Wagenmakers & Farrell, 2004
      dm$w.AIC.selective <- exp(AIC.diff/2)
      newlist[[label]]<-dm
    }
  }
  
  ddf <- do.call(rbind,newlist)
  ddf$t <- as.numeric(ddf$t)
  ddf$label <- as.factor(sub("_[^_]+$", "", ddf$label))
  ddf <- ddf %>% mutate(weight_order = as.factor(case_when(1/w.AIC.selective > 5 & 1/w.AIC.selective < 10  ~ 2,
                                                           1/w.AIC.selective > 10 & 1/w.AIC.selective < 100  ~ 3,
                                                           1/w.AIC.selective > 100   ~ 4,
                                                           w.AIC.selective < 5 & w.AIC.selective > .02  ~ 1, 
                                                           w.AIC.selective > 2 & w.AIC.selective < 10 ~ 2, 
                                                           w.AIC.selective > 10 & w.AIC.selective < 100 ~ 3, 
                                                           w.AIC.selective > 100  ~ 4)),
                        AIC.diff = as.numeric(AIC.diff))
  ddf$weight_order <- factor(ddf$weight_order, labels = c("indecisive", "5-10x", "10-100x", ">100x"))
  
  terms <- effect
  # no FDR correction ----
  ddf <- ddf  %>% mutate(side = substr(as.character(label), nchar(as.character(label)), nchar(as.character(label))),
                         region = substr(as.character(label), 1, nchar(as.character(label))-2))
  
  # plots ----
  
  if (online) {
    setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/dan/plots/clock_decode')
    epoch_label = "Time relative to clock onset, seconds"
  } else {setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/dan/plots/rt_decode/model_compare')
    epoch_label = "Time relative to outcome, seconds"}
  
  for (fe in terms) {
    edf <- ddf %>% filter(effect == paste(fe) & t < 8) 
    termstr <- str_replace_all(fe, "[^[:alnum:]]", "_")
    fname = paste(termstr, ".pdf", sep = "")
    pdf(fname, width = 11, height = 6)
    print(ggplot(edf, aes(t, region)) + geom_tile(aes(fill = AIC.diff, alpha = weight_order), size = 1) +  
            geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + facet_wrap(~side) +
            scale_fill_gradient2(midpoint = 0, low = "blue", mid = "white", high = "red", space = "Lab") + 
            scale_color_grey() + xlab(epoch_label) + ylab("Parcel") + 
            ggtitle(paste(termstr)) + labs(alpha = "Relative evidence \nfor winning \nmodel") +
            labs(fill = "AIC difference\nfavoring\nselective\nmaintenance"))
    dev.off()
    # save model stats ----
    save(file = "medusa_decode_model_comparison.output.Rdata", ddf)
  }
}

## RT prediction ----

# preliminary inspection of behavioral predictors:

# selective <-  lmer(scale(rt_next) ~  rt_csv_sc * outcome  +  scale(rt_vmax)  + scale(rt_vmax_full)  +
#                        rt_lag_sc + 
#                      (1|id), d %>% filter(!is.na(d$rt_vmax_lag) & !is.na(d$rt_vmax_lag_full)),
#                    control=lmerControl(optimizer = "nloptwrap"))
# full <-  lmer(scale(rt_next) ~  rt_csv_sc * outcome  + scale(rt_vmax_full)  +
#                          rt_lag_sc + 
#                         (1|id), d %>% filter(!is.na(d$rt_vmax_lag) & !is.na(d$rt_vmax_lag_full)),
#                control=lmerControl(optimizer = "nloptwrap"))


if (rt_predict) {
  newlist <- list()
  for (label in labels) {print(paste("Processing parcel", label,  sep = " "))
    d$h<-d[[label]]
    if (online) {
      md <-  lmer(scale(rt_next) ~ scale(h) * scale(rt_vmax)  + 
                    scale(h) * rt_csv_sc * last_outcome + scale(h) * rt_lag_sc + 
                    (1|id), d, control=lmerControl(optimizer = "nloptwrap"))
    } else {
      md_h_rtvmax <-  lmer(scale(rt_next) ~ scale(h) * rt_csv_sc * outcome  + scale(h) * scale(rt_vmax)  +
                             scale(h) * rt_lag_sc + scale(rt_vmax_full) +
                             (1|id), d %>% filter(!is.na(d$rt_vmax_lag) & !is.na(d$rt_vmax_lag_full)), control=lmerControl(optimizer = "nloptwrap"))
      mdf_h_rtvmax <-  lmer(scale(rt_next) ~ scale(h) * rt_csv_sc * outcome  + scale(h) * scale(rt_vmax_full)  +
                              scale(h) * rt_lag_sc + scale(rt_vmax) +
                              (1|id), d %>% filter(!is.na(d$rt_vmax_lag) & !is.na(d$rt_vmax_lag_full)), control=lmerControl(optimizer = "nloptwrap"))
    }
    
    # while (any(grepl("failed to converge", md@optinfo$conv$lme4$messages) )) {
    #   print(md@optinfo$conv$lme4$messages)
    #   ss <- getME(md,c("theta","fixef"))
    #   md <- update(md, start=ss)}
    effect <- c("RT.Vmax")
    AIC.diff <- c(AIC(mdf_h_rtvmax) - AIC(md_h_rtvmax))
    dm <- as_tibble(cbind(effect, AIC.diff))
    dm$label <- label
    dm$t <- gsub(".*_", "\\1", label)
    # Akaike weights following Wagenmakers & Farrell, 2004
    dm$w.AIC.selective <- exp(AIC.diff/2)
    newlist[[label]]<-dm
    # }
  }
  ddf <- do.call(rbind,newlist)
  ddf$t <- as.numeric(ddf$t)
  ddf$label <- as.factor(sub("_[^_]+$", "", ddf$label))
  ddf <- ddf %>% mutate(weight_order = as.factor(case_when(1/w.AIC.selective > 5 & 1/w.AIC.selective < 10  ~ 2,
                                                           1/w.AIC.selective > 10 & 1/w.AIC.selective < 100  ~ 3,
                                                           1/w.AIC.selective > 100   ~ 4,
                                                           w.AIC.selective < 5 & w.AIC.selective > .02  ~ 1, 
                                                           w.AIC.selective > 2 & w.AIC.selective < 10 ~ 2, 
                                                           w.AIC.selective > 10 & w.AIC.selective < 100 ~ 3, 
                                                           w.AIC.selective > 100  ~ 4)),
                        AIC.diff = as.numeric(AIC.diff))
  ddf$weight_order <- factor(ddf$weight_order, labels = c("indecisive", "5-10x", "10-100x", ">100x"))
  
  terms <- effect
  # no FDR correction ----
  ddf <- ddf  %>% mutate(side = substr(as.character(label), nchar(as.character(label)), nchar(as.character(label))),
                         region = substr(as.character(label), 1, nchar(as.character(label))-2))
  
  # plots ----
  if (online) {
    setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/dan/plots/clock_rt')
    epoch_label = "Time relative to clock onset, seconds"  
  } else {
    setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/dan/plots/rt_rt/model_compare')
    epoch_label = "Time relative to outcome, seconds"
  }
  for (fe in terms) {
    edf <- ddf %>% filter(effect == paste(fe) & t < 8) 
    termstr <- str_replace_all(fe, "[^[:alnum:]]", "_")
    fname = paste(termstr, ".pdf", sep = "")
    pdf(fname, width = 11, height = 6)
    print(ggplot(edf, aes(t, region)) + geom_tile(aes(fill = AIC.diff, alpha = weight_order), size = 1) +  
            geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + facet_wrap(~side) +
            scale_fill_gradient2(midpoint = 0, low = "blue", mid = "white", high = "red", space = "Lab") + scale_color_grey() + xlab(epoch_label) + ylab("Parcel") + 
            ggtitle(paste(termstr)) + labs(alpha = "Relative evidence \nfor winning \nmodel") +
            labs(fill = "AIC difference\nfavoring\nselective\nmaintenance"))
    dev.off()
    # save model stats ----
    save(file = "medusa_rt_predict_model_comparison.output.Rdata", ddf)
  }
}


