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


# select data
smooth_in_mask = T  # main analysis: data smoothed within mask
unsmoothed = F      # no smoothing whatsoever
newmask = F         # sensivitivy analysis: restrictive COBRA mask (default: Harvard-Oxford)

# what to run
plots = T
decode = T  # main analysis for Fig. 4 E-G
u = F       # exploratory analysis attempting to predict the uncertainty of the next choice
rt_predict = F # predicts next response based on signal and behavioral variables
exclude_first_run = T
online = F

# load data
setwd('~/Box/SCEPTIC_fMRI/dan_medusa/cache/')
load('feedback_dan_wide_ts.Rdata')
load('clock_dan_wide_ts.Rdata')

setwd('~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/')

# read in behavioral data
load('trial_df_and_vh_pe_clusters_u.Rdata')
cache_dir <- "~/Box/SCEPTIC_fMRI/dan_medusa/cache/"
repo_dir <- "~/code/clock_analysis"

# select relevant columns for compactness
df <- df %>% select(id, run, run_trial, rewFunc,emotion, last_outcome, rt_csv, score_csv, rt_next, pe_max, rt_vmax, rt_vmax_lag,
                    rt_vmax_change, v_max_wi, v_entropy_wi, v_entropy_b, v_entropy, v_max_b, u_chosen_quantile, u_chosen_quantile_lag, u_chosen_quantile_change, 
                    rt_vmax_lag_sc, rt_lag_sc, rt_csv_sc, trial_neg_inv_sc, Age, Female)

# add deconvolved timeseries
if (online) {
  d <- merge(df, clock_wide, by = c("id", "run", "run_trial"))
} else {
  d <- merge(df, fb_wide, by = c("id", "run", "run_trial")) }

d <- d %>% group_by(id, run) %>% arrange(id, run, run_trial) %>% mutate(rt_next = lead(rt_csv_sc),
                                                                        rt_change = rt_next - rt_csv_sc,
                                                                        v_entropy_wi_lead = lead(v_entropy_wi),
                                                                        v_entropy_wi_change = v_entropy_wi_lead-v_entropy_wi,
                                                                        u_chosen_quantile_next = lead(u_chosen_quantile),
                                                                        u_chosen_quantile_change_next = lead(u_chosen_quantile_change),
                                                                        outcome = case_when(
                                                                          score_csv>0 ~ 'Reward',
                                                                          score_csv==0 ~ "Omission"),
                                                                        abs_pe = abs(pe_max)
) %>% ungroup()
# remove first run
if (exclude_first_run) {
  d <- d %>% filter(run>1)
}



scale2 <- function(x, na.rm = FALSE) (x - mean(x, na.rm = na.rm)) / sd(x, na.rm)
# scale decon across subjects as a predictor
# choice uncertainty prediction analyses run on scaled 'ds' dataframe instead of 'd'
# ds <- d %>% mutate_at(vars(starts_with("dan")), scale2, na.rm = TRUE) %>% ungroup()


#######
# "decoding" analyses
# combined right and left hippocampus with side as a predictor
# if model does not converge, update with new starting values
labels <- names(d[grepl("_R_|_r_|_L_|_l_", names(d))])
if (decode) {
  newlist <- list()
  for (label in labels) {print(paste("Processing parcel", stringr::str_remove(pattern = gsub(".*_", "\\1", label), label),
                                     ", time", gsub(".*_", "\\1", label),  sep = " "))
    # for (side in c("l", "r")) {
    # for (t in -1:10) {
    d$h<-d[[label]]
    # d$h<-d[[paste(stringr::str_remove(pattern = gsub(".*_", "\\1", label), label), t, sep = "")]]
    md <-  lmer(h ~ trial_neg_inv_sc + scale(rt_csv) + scale(rt_vmax_lag)  + scale(rt_vmax_change)  + 
                  v_entropy_wi + v_entropy_wi_change  + v_max_wi  + pe_max + #u_chosen_quantile_change +
                  (1|id), d, control=lmerControl(optimizer = "nloptwrap"))
    while (any(grepl("failed to converge", md@optinfo$conv$lme4$messages) )) {
      print(md@optinfo$conv$lme4$conv)
      ss <- getME(md,c("theta","fixef"))
      md <- update(md, start=ss)}
    
    dm <- tidy(md)
    dm$label <- label
    # dm$side <- side
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
}
ddf <- do.call(rbind,newlist)
ddf$t <- as.numeric(ddf$t)
ddf$label <- as.factor(stringr::str_remove(pattern = paste0("_", gsub(".*_", "\\1", ddf$label), ""), ddf$label))
ddf$stat_order <- factor(ddf$stat_order, labels = c("NS", "|t| > 2", "|t| > 3"))
ddf$p_value <- factor(ddf$p_value, labels = c("NS", "p < .05", "p < .01", "p < .001"))

terms <- names(fixef(md))
# terms <- names(md$coefficients)
# FDR
ddf <- ddf  %>% group_by(term) %>% mutate(p_fdr = p.adjust(p.value, method = 'fdr'),
                                          p_level_fdr = as.factor(case_when(
                                            # p_fdr > .1 ~ '0',
                                            # p_fdr < .1 & p_fdr > .05 ~ '1',
                                            p_fdr > .05 ~ '1',
                                            p_fdr < .05 & p_fdr > .01 ~ '2',
                                            p_fdr < .01 & p_fdr > .001 ~ '3',
                                            p_fdr <.001 ~ '4'))#,
                                          # side_long = case_when(side=='l' ~ 'Left',
                                          #                       side=='r' ~ 'Right')
)

ddf$p_level_fdr <- factor(ddf$p_level_fdr, levels = c('1', '2', '3', '4'), labels = c("NS","p < .05", "p < .01", "p < .001"))

ddf$`p, FDR-corrected` = ddf$p_level_fdr

##############
# make plots
##############
if (online) {
  setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/dan/plots/online')
  epoch_label = "Time after clock onset, seconds"
} else {
  setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/dan/plots/')
  epoch_label = "Time after outcome, seconds"
}

for (fe in terms) {
  edf <- ddf %>% filter(term == paste(fe) & t < 8) 
  termstr <- str_replace_all(fe, "[^[:alnum:]]", "_")
  pdf(paste(termstr, ".pdf", sep = ""), width = 10, height = 12)
  print(ggplot(edf, aes(t, label)) + geom_tile(aes(fill = estimate, alpha = `p, FDR-corrected`), size = 1) +  
          scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) + ylab("Parcel") + 
          labs(alpha = expression(~italic(p)~', FDR-corrected'))) 
  dev.off()
  # save output for inspection
  save(file = "medusa_decode_output.Rdata", ddf)
}

#################
### RT prediction
#################

if (rt_predict) {
  labels <- names(d[grepl("_R_|_r_|_L_|_l_", names(d))])
  newlist <- list()
  for (label in labels) {print(paste("Processing parcel", stringr::str_remove(pattern = gsub(".*_", "\\1", label), label),
                                     ", time", gsub(".*_", "\\1", label),  sep = " "))
    # for (side in c("l", "r")) {
    # for (t in -1:10) {
    d$h<-d[[label]]
    # d$h<-d[[paste(stringr::str_remove(pattern = gsub(".*_", "\\1", label), label), t, sep = "")]]
    # md <-  lmer(scale(rt_next) ~ trial_neg_inv_sc + h*scale(rt_csv) + h*scale(rt_vmax)  + h*scale(rt_vmax_change)  + 
    #               (1|id), ds, control=lmerControl(optimizer = "nloptwrap"))
    # md <-  lmer(scale(rt_csv) ~ scale(h) * rt_lag_sc * last_outcome + scale(h) * scale(rt_vmax_lag) +
    md <-  lmer(scale(rt_csv) ~ scale(h) * rt_lag_sc * last_outcome + scale(h) * rewFunc +
                                
                  (1|id), d, control=lmerControl(optimizer = "nloptwrap"))
    # while (any(grepl("failed to converge", md@optinfo$conv$lme4$messages) )) {
    #   print(md@optinfo$conv$lme4$messages)
    #   ss <- getME(md,c("theta","fixef"))
    #   md <- update(md, start=ss)}
    # 
    dm <- tidy(md)
    dm$label <- label
    # dm$side <- side
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
}
ddf <- do.call(rbind,newlist)
ddf$t <- as.numeric(ddf$t)
ddf$label <- as.factor(stringr::str_remove(pattern = paste0("_", gsub(".*_", "\\1", ddf$label), ""), ddf$label))
ddf$stat_order <- factor(ddf$stat_order, labels = c("NS", "|t| > 2", "|t| > 3"))
ddf$p_value <- factor(ddf$p_value, labels = c("NS", "p < .05", "p < .01", "p < .001"))

terms <- names(fixef(md))
# terms <- names(md$coefficients)
# FDR
ddf <- ddf  %>% group_by(term) %>% mutate(p_fdr = p.adjust(p.value, method = 'fdr'),
                                          p_level_fdr = as.factor(case_when(
                                            # p_fdr > .1 ~ '0',
                                            # p_fdr < .1 & p_fdr > .05 ~ '1',
                                            p_fdr > .05 ~ '1',
                                            p_fdr < .05 & p_fdr > .01 ~ '2',
                                            p_fdr < .01 & p_fdr > .001 ~ '3',
                                            p_fdr <.001 ~ '4'))#,
                                          # side_long = case_when(side=='l' ~ 'Left',
                                          #                       side=='r' ~ 'Right')
)

ddf$p_level_fdr <- factor(ddf$p_level_fdr, levels = c('1', '2', '3', '4'), labels = c("NS","p < .05", "p < .01", "p < .001"))

ddf$`p, FDR-corrected` = ddf$p_level_fdr

##############
# make plots
##############
if (online) {
  setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/dan/plots/online')
  epoch_label = "Time after clock onset, seconds"
} else {
  setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/dan/plots/')
  epoch_label = "Time after outcome, seconds"
}
for (fe in terms) {
  edf <- ddf %>% filter(term == paste(fe) & t < 8) 
  termstr <- str_replace_all(fe, "[^[:alnum:]]", "_")
  pdf(paste(termstr, "cond_rt_predict.pdf", sep = ""), width = 10, height = 12)
  print(ggplot(edf, aes(t, label)) + geom_tile(aes(fill = estimate, alpha = `p, FDR-corrected`), size = 1) +  
          scale_fill_viridis(option = "plasma") + scale_color_grey() + 
          xlab(epoch_label) + 
          ylab("Parcel") + 
          labs(alpha = expression(~italic(p)~', FDR-corrected'))) 
  dev.off()
  # save output for inspection
  save(file = "medusa_rt_predict_output.Rdata", ddf)
}


########
# choice uncertainty prediction analyses
# diagnose trial/h collinearity
# ggplot(fb_comb, aes(run_trial, decon_interp, color = as.factor(bin_center))) + geom_smooth(method = 'gam', se = F)
# drop models without contingency and trial
# running with h scaled

# to think more about it, all behavioral variables should be represented in the brain
# what is represented is addressed in decoding analyses
# that said, it seems that uncertainty prediction analyses yield negative results over a range of models

# uncertainty prediction analyses also run on scaled data
if (u) {
  # for (trial_cont in c("TRUE", "FALSE")) {
  # for (trial_cont in c("FALSE")) {
  newlist <- list()
  for (slice in 1:9) {print(paste("Processing slice", slice, sep = " "))
    # for (side in c("l", "r")) {
    for (t in 0:10) {
      d$h<-d[[paste("vmPFC", slice, t, sep = "_")]]
      # if (trial_cont) {
      uf <- lmer(u_chosen_quantile_next ~  h * scale(run_trial) + u_chosen_quantile + (1|id/run), d)
      # else {
      # mf <-  lme4::lmer(rt_next ~ (scale(pe_max) + scale(rt_csv) + scale(rt_vmax_lag) + scale(rt_vmax_change) + scale(v_entropy_wi) + h)^2 + (1|id/run), ds)
      # uf <- lmer(u_chosen_next ~ scale(-1/run_trial)*scale(h) + scale(rt_csv)*scale(h) + scale(rt_vmax)*scale(h) + last_outcome*scale(h) + v_entropy_wi*scale(h) + scale(u_chosen) + (1|id/run), ds) 
      # }
      dm <- broom.mixed::tidy(uf,effects = "fixed") %>% mutate(term = str_remove(term, "TRUE")) # make betas compatible with ANOVA
      # run without anova
      # an <- broom.mixed::tidy(car::Anova(uf, '3')) %>% rename(anova_p = p.value, chisq = statistic)
      dm$slice <- slice
      # dm$side <- side
      dm$t <- t
      # dm <- inner_join(dm, an, by = "term") # this only works for continuous terms, unfortunately
      dm <- dm %>% mutate(stat_order = as.factor(case_when(abs(statistic) < 2 ~ '1', 
                                                           abs(statistic) > 2 & abs(statistic) < 3 ~ '2', 
                                                           abs(statistic) > 3 ~ '3')),
                          p_value = as.factor(case_when(p.value > .05 ~ '1',
                                                        p.value < .05 & p.value > .01 ~ '2',
                                                        p.value < .01 & p.value > .001 ~ '3',
                                                        p.value <.001 ~ '4')))
      # newlist[[paste("hipp", slice, side, t, sep = "_")]]<-dm
      newlist[[paste("hipp", slice, t, sep = "_")]]<-dm
    }
    # }
  }
  bdf <- do.call(rbind,newlist)
  bdf$slice <- as.factor(bdf$slice)
  bdf$stat_order <- factor(bdf$stat_order, labels = c("NS", "|t| > 2", "|t| > 3"))
  bdf <- bdf  %>% group_by(term) %>% mutate(p_fdr = p.adjust(p.value, method = 'fdr'),
                                            p_level_fdr = as.factor(case_when(
                                              # p_fdr > .1 ~ '0',
                                              # p_fdr < .1 & p_fdr > .05 ~ '1',
                                              p_fdr > .05 ~ '1',
                                              p_fdr < .05 & p_fdr > .01 ~ '2',
                                              p_fdr < .01 & p_fdr > .001 ~ '3',
                                              p_fdr <.001 ~ '4'))#,
                                            # side_long = case_when(side=='l' ~ 'Left',
                                            #                       side=='r' ~ 'Right')
  )
  
  # ddf$p_level_fdr <- factor(ddf$p_level_fdr, labels = c("NS","p < .1", "p < .05", "p < .01", "p < .001"))
  bdf$p_level_fdr <- factor(bdf$p_level_fdr, labels = c("NS","p < .05", "p < .01", "p < .001"))
  
  bdf$`p, FDR-corrected` = bdf$p_level_fdr
  
  # sanity check for FDR-corrected p value labels
  ggplot(bdf, aes(p_level_fdr, p_fdr)) + geom_point()
  
  terms <- unique(bdf$term)
  # if (trial_cont) {
  if (unsmoothed) {setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/hippo/figs/u_predict/unsmoothed')
    # } else {setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/hippo/figs/newmask/u_predict')}
  } else {setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/hippo/figs/u_predict')}
  
  # else {
  # if (unsmoothed) {setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/hippo/figs/u_predict/unsmoothed/no_trial_contingency')
  # } else {setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/hippo/figs/newmask/u_predict/no_trial_contingency/')}}
  setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/vmPFC/figs/u_predict')
  if (plots) {
    for (fe in terms) 
    {
      edf <- bdf %>% filter(term == paste(fe) & t < 8) 
      termstr <- str_replace_all(fe, "[^[:alnum:]]", "_")
      pdf(paste(termstr, "_vmPFC.pdf", sep = ""), width = 5, height = 3.5)
      print(ggplot(edf, aes(t, slice)) + geom_tile(aes(fill = estimate, alpha = `p, FDR-corrected`), size = 1) +  
              scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab("Time after outcome, seconds") + ylab("Posterior <-- Location --> Anterior\n (9 slices)") + 
              labs(alpha = expression(~italic(p)~', FDR-corrected'))) 
      dev.off()
    }
  }
}
# }
# change entropy to lead
# check v_max whole-brain 

