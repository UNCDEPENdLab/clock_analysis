# with 'decode = T' makes MEDUSA decoding plots for Fig. 4 E-G.
# loops over decoding, RT prediction and uncertainty prediction models for various hippocampal slices and post-feedback time points
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
unsmoothed = F
newmask = T
smooth_in_mask = F

# what to run
plots = T
rt = F
u = F
decode = T
centrality = F

# default is the more inclusive Harvard-Oxford mask
if (unsmoothed) {setwd("~/Box/SCEPTIC_fMRI/var/unsmoothed")
  # } else {setwd("/Users/localadmin/Box/SCEPTIC_fMRI/var/newmask")}
} else if (smooth_in_mask) {setwd("~/Box/SCEPTIC_fMRI/var/smooth_in_mask/")
  } else {setwd("~/Box/SCEPTIC_fMRI/var/")}

load('feedback_hipp_widest_by_timepoint_decon.Rdata')
load('feedback_hipp_wide_ts.Rdata')
setwd('~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/')
# load('trial_df_and_vhdkfpe_clusters.Rdata')
load('trial_df_and_vh_pe_clusters_u.Rdata')
if (unsmoothed) {cache_dir <- "~/Box/SCEPTIC_fMRI/var/unsmoothed"
# } else {cache_dir <- "~/Box/SCEPTIC_fMRI/var/newmask"}
} else {cache_dir <- "~/Box/SCEPTIC_fMRI/var"}
repo_dir <- "~/Data_Analysis/clock_analysis"

#super-wide variant used in lm analysis
# load(file.path(cache_dir, "feedback_hipp_tallest_by_timepoint_decon.Rdata"))
# 
# load(file.path(repo_dir, "/fmri/keuka_brain_behavior_analyses/trial_df_and_vhdkfpe_clusters.Rdata"))

attr(df, "labels") <- NULL #somehow this is holding a 560-row data.frame
df <- df %>% dplyr::ungroup()

# read in behavioral data
# select relevant columns for compactness
df <- df %>% select(id, run, run_trial, rewFunc,emotion, last_outcome, rt_csv, score_csv, rt_next, pe_max, rt_vmax, rt_vmax_lag,
                    rt_vmax_change, v_max_wi, v_entropy_wi, v_entropy_b, v_entropy, v_max_b, u_chosen_quantile, u_chosen_quantile_lag, u_chosen_quantile_change, 
                    rt_vmax_lag_sc, rt_lag_sc, rt_csv_sc, trial_neg_inv_sc, Age, Female)
# add deconvolved hippocampal timeseries
d <- merge(df, fb_wide_t, by = c("id", "run", "run_trial"))
d <- d %>% group_by(id, run) %>% arrange(id, run, run_trial) %>% mutate(rt_change = 100*rt_next - rt_csv, 
                                                                        v_entropy_wi_lead = lead(v_entropy_wi),
                                                                        v_entropy_wi_change = v_entropy_wi_lead-v_entropy_wi,
                                                                        u_chosen_quantile_next = lead(u_chosen_quantile),
                                                                        u_chosen_quantile_change_next = lead(u_chosen_quantile_change),
                                                                        outcome = lead(last_outcome)) %>% ungroup()
dbl <- merge(df, fb_wide_bl, by = c("id", "run", "run_trial"))
dbl <- dbl %>% group_by(id, run) %>% arrange(id, run, run_trial) %>% mutate(rt_change = 100*rt_next - rt_csv, 
                                                                        v_entropy_wi_lead = lead(v_entropy_wi),
                                                                        v_entropy_wi_change = v_entropy_wi_lead-v_entropy_wi,
                                                                        u_chosen_quantile_next = lead(u_chosen_quantile),
                                                                        u_chosen_quantile_change_next = lead(u_chosen_quantile_change),
                                                                        outcome = lead(last_outcome)) %>% ungroup()

scale2 <- function(x, na.rm = FALSE) (x - mean(x, na.rm = na.rm)) / sd(x, na.rm)
# scale decon across subjects as a predictor
# RT and U prediction analyses now run on 'ds' instead of 'd'
ds <- d %>% mutate_at(vars(starts_with("hipp")), scale2, na.rm = TRUE) %>% ungroup()


##########
# predict directional RT change
# try removing the entropy terms for simplicity
# strategically: pe*rt_t (or last_outcome*rt) - stay vs. shift, rt - explore in general, rt_vmax - exploit, prior, rt_vmax_change - exploit, posterior, that's it!
# for (trial_cont in c("TRUE", "FALSE")) {
if (rt) {
  for (trial_cont in c("TRUE")) {
    newlist <- list()
    for (slice in 1:12) {print(paste("Processing slice", slice, sep = " "))
      for (side in c("l", "r")) {
        for (t in 0:10) {
          d$h<-d[[paste("hipp", slice, side, t, sep = "_")]]
          # mf <-  lme4::lmer(rt_next ~ (scale(pe_max) + scale(rt_csv) + scale(rt_vmax_lag) + scale(rt_vmax_change) + scale(v_entropy_wi) + h)^2 + (1|id/run), ds)
          # mf <-  lmerTest::lmer(rt_next ~ scale(rt_csv) * last_outcome * h + scale(rt_vmax_lag) *  h + scale(rt_vmax_change) *  h + (1|id/run), ds)
          # a more detailed model
          mf <-  lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + rt_vmax_lag_sc + last_outcome + h)^2 + 
                        rt_lag_sc:last_outcome:h + 
                        rt_vmax_lag_sc:trial_neg_inv_sc:h + 
                        (1|id/run), ds)
          dm <- broom.mixed::tidy(mf,effects = "fixed")
          dm$slice <- slice
          dm$side <- side
          dm$t <- t
          dm <- dm %>% mutate(stat_order = as.factor(case_when(abs(statistic) < 2 ~ '1', 
                                                               abs(statistic) > 2 & abs(statistic) < 3 ~ '2', 
                                                               abs(statistic) > 3 ~ '3')),
                              p_value = as.factor(case_when(p.value > .05 ~ '1',
                                                            p.value < .05 & p.value > .01 ~ '2',
                                                            p.value < .01 & p.value > .001 ~ '3',
                                                            p.value <.001 ~ '4')))
          newlist[[paste("hipp", slice, side, t, sep = "_")]]<-dm
        }
      }
    }
    bdf <- do.call(rbind,newlist)
    bdf$slice <- as.factor(bdf$slice)
    bdf$stat_order <- factor(bdf$stat_order, labels = c("NS", "|t| > 2", "|t| > 3"))
    bdf <- bdf  %>% group_by(term,side) %>% mutate(p_fdr = p.adjust(p.value, method = 'fdr'),
                                                   p_level_fdr = as.factor(case_when(p_fdr > .05 ~ '1',
                                                                                     p_fdr < .05 & p_fdr > .01 ~ '2',
                                                                                     p_fdr < .01 & p_fdr > .001 ~ '3',
                                                                                     p_fdr <.001 ~ '4')))
    bdf$p_level_fdr <- factor(bdf$p_level_fdr, labels = c("NS", "p < .05", "p < .01", "p < .001"))
    
    terms <- names(fixef(mf))
    if (trial_cont) {
      if (unsmoothed) {setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/hippo/figs/rt_predict/unsmoothed')
        # } else {setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/hippo/figs/newmask/rt_predict')}}
      } else {setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/hippo/figs/rt_predict')}}
    
    else {
      if (unsmoothed) {setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/hippo/figs/rt_predict/unsmoothed/no_trial_contingency')
        # } else {setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/hippo/figs/newmask/rt_predict/no_trial_contingency/')}}
      } else {setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/hippo/figs/rt_predict/no_trial_contingency/')}}
    
    if (plots) {
      for (fe in terms)
      {edf <- bdf %>% filter(term == paste(fe) & t < 8)
      p1 <- ggplot(edf, aes(t, estimate, color = slice)) + geom_line() + 
        geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error), alpha = .5) + facet_wrap(~side) + 
        theme_dark() + scale_color_viridis_d() + geom_hline(yintercept = 0, lty = "dashed", color = "red") + labs(title = paste(fe))
      
      p2 <- ggplot(edf, aes(t, slice)) + geom_tile(aes(fill = estimate, alpha = p_level_fdr), size = 1) + facet_wrap(~side) + 
        scale_fill_viridis(option = "plasma") + scale_color_grey() + labs(title = paste(fe))
      
      termstr <- str_replace_all(fe, "[^[:alnum:]]", "_")
      pdf(paste(termstr, ".pdf", sep = ""), width = 12, height = 12)
      print(ggarrange(p1,p2,ncol = 1, nrow = 2))
      dev.off()
      }
    }
  }
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
  for (slice in 1:12) {print(paste("Processing slice", slice, sep = " "))
    for (side in c("l", "r")) {
      for (t in 0:10) {
        d$h<-d[[paste("hipp", slice, side, t, sep = "_")]]
        # if (trial_cont) {
        uf <- lmer(u_chosen_quantile_next ~  h * scale(run_trial) + u_chosen_quantile + (1|id/run), ds)
        # else {
        # mf <-  lme4::lmer(rt_next ~ (scale(pe_max) + scale(rt_csv) + scale(rt_vmax_lag) + scale(rt_vmax_change) + scale(v_entropy_wi) + h)^2 + (1|id/run), ds)
        # uf <- lmer(u_chosen_next ~ scale(-1/run_trial)*scale(h) + scale(rt_csv)*scale(h) + scale(rt_vmax)*scale(h) + last_outcome*scale(h) + v_entropy_wi*scale(h) + scale(u_chosen) + (1|id/run), ds) 
        # }
        dm <- broom.mixed::tidy(uf,effects = "fixed") %>% mutate(term = str_remove(term, "TRUE")) # make betas compatible with ANOVA
        # run without anova
        # an <- broom.mixed::tidy(car::Anova(uf, '3')) %>% rename(anova_p = p.value, chisq = statistic)
        dm$slice <- slice
        dm$side <- side
        dm$t <- t
        # dm <- inner_join(dm, an, by = "term") # this only works for continuous terms, unfortunately
        dm <- dm %>% mutate(stat_order = as.factor(case_when(abs(statistic) < 2 ~ '1', 
                                                             abs(statistic) > 2 & abs(statistic) < 3 ~ '2', 
                                                             abs(statistic) > 3 ~ '3')),
                            p_value = as.factor(case_when(p.value > .05 ~ '1',
                                                          p.value < .05 & p.value > .01 ~ '2',
                                                          p.value < .01 & p.value > .001 ~ '3',
                                                          p.value <.001 ~ '4')))
        newlist[[paste("hipp", slice, side, t, sep = "_")]]<-dm
      }
    }
  }
  bdf <- do.call(rbind,newlist)
  bdf$slice <- as.factor(bdf$slice)
  bdf$stat_order <- factor(bdf$stat_order, labels = c("NS", "|t| > 2", "|t| > 3"))
  bdf <- bdf  %>% group_by(term,side) %>% mutate(p_fdr = p.adjust(p.value, method = 'fdr'),
                                                 p_level_fdr = as.factor(case_when(p_fdr > .05 ~ '1',
                                                                                   p_fdr < .05 & p_fdr > .01 ~ '2',
                                                                                   p_fdr < .01 & p_fdr > .001 ~ '3',
                                                                                   p_fdr <.001 ~ '4'))
                                                 # p_anova_fdr = p.adjust(anova_p, method = 'fdr'),
                                                 # p_anova_level_fdr = as.factor(case_when(p_anova_fdr > .05 ~ '1',
                                                 # p_anova_fdr < .05 & p_anova_fdr > .01 ~ '2',
                                                 # p_anova_fdr < .01 & p_anova_fdr > .001 ~ '3',
                                                 # p_anova_fdr <.001 ~ '4')),
  )
  # bdf$p_level_fdr <- factor(bdf$p_level_fdr, labels = c("NS", "p < .05", "p < .01", "p < .001"))
  bdf$p_level_fdr <- factor(bdf$p_level_fdr, labels = c("NS", "p < .01", "p < .001"))
  # bdf$p_level_fdr <- factor(bdf$p_level_fdr, labels = c("NS", "p < .001"))
  
  # bdf$p_anova_level_fdr <- factor(bdf$p_anova_level_fdr, labels = c("NS", "p < .05", "p < .01", "p < .001"))
  
  terms <- unique(bdf$term)
  # if (trial_cont) {
  if (unsmoothed) {setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/hippo/figs/u_predict/unsmoothed')
    # } else {setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/hippo/figs/newmask/u_predict')}
  } else {setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/hippo/figs/u_predict')}
  
  # else {
  # if (unsmoothed) {setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/hippo/figs/u_predict/unsmoothed/no_trial_contingency')
  # } else {setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/hippo/figs/newmask/u_predict/no_trial_contingency/')}}
  if (plots) {
    for (fe in terms) 
    {
      edf <- bdf %>% filter(term == paste(fe) & t < 8)
      p1 <- ggplot(edf, aes(t, estimate, color = slice)) + geom_line() + 
        geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error), alpha = .5) + facet_wrap(~side) + 
        theme_dark() + scale_color_viridis_d() + geom_hline(yintercept = 0, lty = "dashed", color = "red") + labs(title = paste(fe))
      
      p2 <- ggplot(edf, aes(t, slice)) + geom_tile(aes(fill = estimate, alpha = p_level_fdr), size = 1) + facet_wrap(~side) + 
        scale_fill_viridis(option = "plasma") + scale_color_grey() + labs(title = paste(fe))
      # p3 <- ggplot(edf, aes(t, slice)) + geom_tile(aes(fill = chisq, alpha = p_anova_level_fdr), size = 1) + facet_wrap(~side) + 
      # scale_fill_viridis(option = "plasma") + scale_color_grey() + labs(title = paste(fe))
      termstr <- str_replace_all(fe, "[^[:alnum:]]", "_")
      pdf(paste(termstr, ".pdf", sep = ""), width = 12, height = 12)
      # pdf(paste(termstr, ".pdf", sep = ""), width = 12, height = 16)
      # print(ggarrange(p1,p2,p3,ncol = 1, nrow = 3))
      print(ggarrange(p1,p2,ncol = 1, nrow = 2))
      dev.off()
    }
  }
}
# }


#######


# "decoding" analyses
# this will always include trial and contingency

# getting rid of run nested within ID eliminates most singular fit warnings
# if (decode) {
#   newlist <- list()
#   for (slice in 1:12) {print(paste("Processing slice", slice, sep = " "))
#     for (side in c("l", "r")) {
#       for (t in -1:10) {
#         d$h<-d[[paste("hipp", slice, side, t, sep = "_")]]
#         md <-  lmer(h ~  scale(rt_csv)  + scale(rt_vmax_lag) + scale(rt_vmax_change) + v_entropy_wi + v_entropy_wi_change + #u_chosen_quantile_change +
#                       (1|id), d, control=lmerControl(optimizer = "nloptwrap"))
#         while (any(grepl("failed to converge", md@optinfo$conv$lme4$messages) )) {
#           print(md@optinfo$conv$lme4$conv)
#           ss <- getME(md,c("theta","fixef"))
#           md <- update(md, start=ss)}
#         
#         dm <- tidy(md)
#         dm$slice <- slice
#         dm$side <- side
#         dm$t <- t
#         dm <- dm %>% mutate(stat_order = as.factor(case_when(abs(statistic) < 2 ~ '1', 
#                                                              abs(statistic) > 2 & abs(statistic) < 3 ~ '2', 
#                                                              abs(statistic) > 3 ~ '3')),
#                             p_value = as.factor(case_when(p.value > .05 ~ '1',
#                                                           p.value < .05 & p.value > .01 ~ '2',
#                                                           p.value < .01 & p.value > .001 ~ '3',
#                                                           p.value <.001 ~ '4'))
#         )
#         newlist[[paste("hipp", slice, side, t, sep = "_")]]<-dm
#       }
#     }
#   }
#   ddf <- do.call(rbind,newlist)
#   ddf$slice <- as.factor(ddf$slice)
#   ddf$stat_order <- factor(ddf$stat_order, labels = c("NS", "|t| > 2", "|t| > 3"))
#   ddf$p_value <- factor(ddf$p_value, labels = c("NS", "p < .05", "p < .01", "p < .001"))
#   
#   terms <- names(fixef(md))
#   # terms <- names(md$coefficients)
#   # FDR
#   ddf <- ddf  %>% group_by(term,side) %>% mutate(p_fdr = p.adjust(p.value, method = 'fdr'),
#                                                  p_level_fdr = as.factor(case_when(
#                                                    p_fdr > .1 ~ '0',
#                                                    p_fdr < .1 & p_fdr > .05 ~ '1',
#                                                    p_fdr < .05 & p_fdr > .01 ~ '2',
#                                                    p_fdr < .01 & p_fdr > .001 ~ '3',
#                                                    p_fdr <.001 ~ '4')),
#                                                  side_long = case_when(side=='l' ~ 'Left',
#                                                                        side=='r' ~ 'Right')
#   )
#   
#   ddf$p_level_fdr <- factor(ddf$p_level_fdr, labels = c("NS","p < .1", "p < .05", "p < .01", "p < .001"))
#   ddf$`p, FDR-corrected` = ddf$p_level_fdr
#   if (unsmoothed) {setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/hippo/figs/decode/unsmoothed/lmer')
#   } else if (smooth_in_mask) {setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/hippo/figs/decode/smooth_in_mask/lmer')
#     } else {setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/hippo/figs/decode/lmer')} # manually indicate if this is the new COBRA lab mask
#   for (fe in terms) {
#     edf <- ddf %>% filter(term == paste(fe) & t < 8) 
#     termstr <- str_replace_all(fe, "[^[:alnum:]]", "_")
#     pdf(paste(termstr, ".pdf", sep = ""), width = 5, height = 3)
#     
#     # print(ggplot(edf, aes(t, slice)) + geom_tile(aes(fill = estimate, alpha = abs(statistic)>2), size = 1) + facet_wrap(~side) + 
#     # print(ggplot(edf, aes(t, slice)) + geom_tile(aes(fill = estimate, alpha = p_level_fdr), size = 1) + facet_wrap(~side) + 
#     print(ggplot(edf, aes(t, slice)) + geom_tile(aes(fill = estimate, alpha = `p, FDR-corrected`), size = 1, height = p_fdr>.05) + facet_wrap(~side_long) + 
#             scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab("Time after outcome, seconds") + ylab("Posterior <-- Location --> Anterior\n (12 slices)") + 
#             labs(alpha = expression(~italic(p)~', FDR-corrected'))) 
#     # print(ggarrange(p2,ncol = 1, labels = paste(fe), vjust = 4, font.label = list(color = "black", size = 16)))
#     dev.off()
#   }
# }

# combined right and left with side as a predictor
if (decode) {
  newlist <- list()
  for (slice in 1:12) {print(paste("Processing slice", slice, sep = " "))
    # for (side in c("l", "r")) {
      for (t in -1:10) {
        dbl$h<-dbl[[paste("hipp", slice, t, sep = "_")]]
        md <-  lmer(h ~  scale(rt_csv) + side  + scale(rt_vmax_lag)  + scale(rt_vmax_change)  + v_entropy_wi  + v_entropy_wi_change  + #u_chosen_quantile_change +
                      (1|id) + (1|side), dbl, control=lmerControl(optimizer = "nloptwrap"))
        while (any(grepl("failed to converge", md@optinfo$conv$lme4$messages) )) {
          print(md@optinfo$conv$lme4$conv)
          ss <- getME(md,c("theta","fixef"))
          md <- update(md, start=ss)}
        
        dm <- tidy(md)
        dm$slice <- slice
        # dm$side <- side
        dm$t <- t
        dm <- dm %>% mutate(stat_order = as.factor(case_when(abs(statistic) < 2 ~ '1', 
                                                             abs(statistic) > 2 & abs(statistic) < 3 ~ '2', 
                                                             abs(statistic) > 3 ~ '3')),
                            p_value = as.factor(case_when(p.value > .05 ~ '1',
                                                          p.value < .05 & p.value > .01 ~ '2',
                                                          p.value < .01 & p.value > .001 ~ '3',
                                                          p.value <.001 ~ '4'))
        )
        newlist[[paste("hipp", slice, t, sep = "_")]]<-dm
      # }
    }
  }
  ddf <- do.call(rbind,newlist)
  ddf$slice <- as.factor(ddf$slice)
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
  
  # ddf$p_level_fdr <- factor(ddf$p_level_fdr, labels = c("NS","p < .1", "p < .05", "p < .01", "p < .001"))
  ddf$p_level_fdr <- factor(ddf$p_level_fdr, labels = c("NS","p < .05", "p < .01", "p < .001"))
  
    ddf$`p, FDR-corrected` = ddf$p_level_fdr
  if (unsmoothed) {setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/hippo/figs/decode/unsmoothed/lmer')
  } else if (smooth_in_mask) {setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/hippo/figs/decode/smooth_in_mask/lmer')
  } else {setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/hippo/figs/decode/lmer')} # manually indicate if this is the new COBRA lab mask
  for (fe in terms) {
    edf <- ddf %>% filter(term == paste(fe) & t < 8) 
    termstr <- str_replace_all(fe, "[^[:alnum:]]", "_")
    pdf(paste(termstr, "_bl.pdf", sep = ""), width = 5, height = 3.5)
    print(ggplot(edf, aes(t, slice)) + geom_tile(aes(fill = estimate, alpha = `p, FDR-corrected`), size = 1) +  
            scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab("Time after outcome, seconds") + ylab("Posterior <-- Location --> Anterior\n (12 slices)") + 
            labs(alpha = expression(~italic(p)~', FDR-corrected'))) 
    dev.off()
  }
}

# centrality loop
if (centrality) {
  
  #########
  # add centrality for a similar RT prediction analysis
  load('~/OneDrive/collected_letters/papers/sceptic_fmri/hippo/brms_6slc_graphs.RData')
  l <- gobj_lh$all_nodal
  r <- gobj_rh$all_nodal
  # visualize -- same as Michael's plots
  # ggplot(l) + geom_boxplot(aes(node, strength_in), outlier.shape = NA) + geom_jitter(aes(node, strength_in, alpha = .1)) + 
  #   scale_x_discrete(limits = c("lh2", "lh4", "lh6", "lh8", "lh10", "lh12"))
  # ggplot(l) + geom_boxplot(aes(node, strength_out), outlier.shape = NA) + geom_jitter(aes(node, strength_out, alpha = .1)) + 
  #   scale_x_discrete(limits = c("lh2", "lh4", "lh6", "lh8", "lh10", "lh12"))
  l$id <- l$ID
  l <- l %>% select(id, node, strength_in, strength_out)
  lw <- as.tibble(dcast(setDT(l), id ~ node, value.var = c("strength_in", "strength_out")))
  r$id <- r$ID
  r <- r %>% select(id, node, strength_in, strength_out)
  rw <- as.tibble(dcast(setDT(r), id ~ node, value.var = c("strength_in", "strength_out")))
  dc <- merge(lw,rw, by = "id")
  dfc <- merge(dc,df, by = "id")
  dfc <- dfc %>% group_by(id, run) %>% arrange(id, run, run_trial) %>% mutate(last_outcome = score_csv > 0, 
                                                                              rt_change = 100*rt_next - rt_csv, 
                                                                              v_entropy_wi_lead = lead(v_entropy_wi),
                                                                              v_entropy_wi_change = v_entropy_wi_lead-v_entropy_wi) %>% ungroup()
  
  
  newlist <- list()
  for (slice in c(2,4,6,8,10,12)) {print(paste("Processing slice", slice, sep = " "))
    for (side in c("lh", "rh")) {
      for (dir in c("out", "in")) {
        dfc$c <- dfc[[paste(paste("strength", dir, side, sep = "_"),slice, sep = "")]]
        mf <-  lmerTest::lmer(rt_next ~ scale(-1/run_trial)*rewFunc + scale(rt_csv)  * last_outcome + scale(rt_csv)  * c + scale(rt_vmax_lag) *  c + scale(rt_vmax_change) *  c + (1|id), dfc)
        dm <- broom.mixed::tidy(mf,effects = "fixed")
        dm$slice <- slice
        dm$side <- side
        dm$direction <- dir
        dm <- dm %>% mutate(stat_order = as.factor(case_when(abs(statistic) < 2 ~ '1', 
                                                             abs(statistic) > 2 & abs(statistic) < 3 ~ '2', 
                                                             abs(statistic) > 3 ~ '3')),
                            p_value = as.factor(case_when(p.value > .05 ~ '1',
                                                          p.value < .05 & p.value > .01 ~ '2',
                                                          p.value < .01 & p.value > .001 ~ '3',
                                                          p.value <.001 ~ '4')),
                            p_fdr = p.adjust(p.value, method = 'fdr'),
                            p_level_fdr = as.factor(case_when(p_fdr > .05 ~ '1',
                                                              p_fdr < .05 & p_fdr > .01 ~ '2',
                                                              p_fdr < .01 & p_fdr > .001 ~ '3',
                                                              p_fdr <.001 ~ '4'))
        )
        newlist[[paste(paste("strength", dir, side, sep = "_"), slice, sep = "")]] <- dm
      }
    }
  }
  cdf <- do.call(rbind,newlist)
  cdf$slice <- as.factor(cdf$slice)
  cdf$stat_order <- factor(cdf$stat_order, labels = c("NS", "|t| > 2", "|t| > 3"))
  cdf$p_level_fdr <- factor(cdf$p_level_fdr, labels = c("NS", "p < .05", "p < .01", "p < .001"))
  
  terms <- names(fixef(mf))[c(15:17)]
  setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/hippo/figs/centrality')
  
  if (plots) {
    for (fe in terms)
    {edf <- cdf %>% filter(term == paste(fe))
    p2 <- ggplot(cdf[cdf$term==terms,], aes(direction, slice)) + geom_tile(aes(fill = estimate, alpha = p_level_fdr), size = 1) + facet_grid(term~side) + 
      scale_fill_distiller(palette = "Spectral") + scale_color_grey() + labs(title = paste(fe))
    
    termstr <- str_replace_all(fe, "[^[:alnum:]]", "_")
    pdf(paste(termstr, ".pdf", sep = ""), width = 6, height = 12)
    print(p2)
    dev.off()
    }
  }
  
  ggplot(cdf[cdf$term==terms,], aes(slice, estimate, color = direction)) + geom_line() +
    geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error), alpha = .5) + facet_wrap(term~side) +
    theme_dark() + scale_color_viridis_d() + geom_hline(yintercept = 0, lty = "dashed", color = "red")# + labs(title = paste(fe))
}
#for (trial_cont in c("TRUE", "FALSE")) {
#  newlist <- list()
#
#  # try to predict directional RT change
#
#  # loop over slices and timepoints
#  for (slice in 1:12) {
#    print(paste("Processing slice", slice, sep = " "))
#    for (side in c("l", "r")) {
#      for (t in -1:10) {
#        d$h<-d[[paste("hipp", slice, side, t, sep = "_")]]
#        if (trial_cont) {
#          mf <-  lmer(rt_next ~ scale(-1/run_trial)*rewFunc + (last_outcome + scale(rt_csv) + scale(rt_vmax_lag) + scale(rt_vmax) + h)^3 + (1|id/run), d)}
#        else {
#          mf <-  lmer(rt_next ~ (last_outcome + scale(rt_csv) + scale(rt_vmax_lag) + h)^3 + (1|id/run), d)}
#        dm <- tidy(mf)
#        dm$slice <- slice
#        dm$side <- side
#        dm$t <- t
#        newlist[[paste("hipp", slice, side, t, sep = "_")]]<-dm
#      }
#    }
#  }
#  
#  bdf <- do.call(rbind,newlist)
#  bdf$slice <- as.factor(bdf$slice)
#  terms <- names(fixef(mf))
#  
#  if (trial_cont) {
#    setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/hippo/figs')
#  }
#  else {
#    setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/hippo/figs/no_trial_contingency/')
#  }
#  
#  for (fe in terms) {
#    edf <- bdf %>% filter(term == paste(fe) & t < 8)
#    p1 <- ggplot(edf, aes(t, estimate, color = slice)) + geom_line() + 
#        geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error), alpha = .5) + facet_wrap(~side) + 
#        theme_dark() + scale_color_viridis_d() + geom_hline(yintercept = 0, lty = "dashed", color = "red") + labs(title = paste(fe))
#    
#    p2 <- ggplot(edf, aes(t, slice)) + geom_tile(aes(fill = estimate, alpha = abs(statistic)>2), size = 1) + facet_wrap(~side) + 
#        scale_fill_viridis(option = "plasma") + scale_color_grey() + labs(title = paste(fe))
#    
#    termstr <- str_replace_all(fe, "[^[:alnum:]]", "_")
#    pdf(paste(termstr, ".pdf", sep = ""), width = 12, height = 12)
#    print(ggarrange(p1,p2,ncol = 1, nrow = 2))
#    dev.off()
#  }
#}

# # "decoding" analyses
# # currently running lm
# for (trial_cont in c(TRUE, FALSE)) {
#   newlist <- list()
#   # loop over slices and timepoints
#   for (slice in 1:12) {
#     print(paste("Processing slice", slice, sep = " "))
#     for (side in c("l", "r")) {    
#       for (t in -1:10) {      
#         d$h<-d[[paste("hipp", slice, side, t, sep = "_")]]
#         #only run models with trial_cont=TRUE
#         if (trial_cont) {
#           mform <- as.formula(h ~ scale(-1/run_trial)*rewFunc + last_outcome + scale(rt_csv) + scale(rt_vmax_lag) + scale(rt_vmax_change) + v_entropy_wi + (1|id))
#           #TODO: Why does lmer return singular fit?
#           xtabs(~id + run_trial, d) #looks complete
#           d %>% group_by(id) %>% summarize(mean(h, na.rm=T), sd(h, na.rm=T))
#           
#           lattice::histogram(d$h)
#           md1 <-  lmer(h ~ scale(-1/run_trial)*rewFunc + last_outcome + scale(rt_csv) + scale(rt_vmax_lag) + scale(rt_vmax_change) + v_entropy_wi + (1|id/run), d)
#           md2 <-  lmer(h ~ scale(-1/run_trial)*rewFunc + last_outcome + scale(rt_csv) + scale(rt_vmax_lag) + scale(rt_vmax_change) + v_entropy_wi + (1|id), d)
#           md3 <-  lmer(h ~ scale(-1/run_trial)*rewFunc + last_outcome + scale(rt_csv) + scale(rt_vmax_lag) + scale(rt_vmax_change) + v_entropy_wi + (1|id:run), d)
#           
#           md2 <-  lmer(h ~ 1 + (1|id), d)
#           
#           md3 <-  lm(h ~ scale(-1/run_trial)*rewFunc + last_outcome + scale(rt_csv) + scale(rt_vmax_lag) + scale(rt_vmax_change) + v_entropy_wi, d)
#           
#           anova(md1, md2, md3)
#           anova(md1, md3)
#           library(RLRsim)
#           
#           md1_pure <- lme4::lmer(h ~ scale(-1/run_trial)*rewFunc + last_outcome + scale(rt_csv) + scale(rt_vmax_lag) + scale(rt_vmax_change) + v_entropy_wi + (1|id/run), d, REML=TRUE)
#           exactLRT(md1_pure, md3)
#           
#           malt1 <- update(md1_pure, . ~ . - (1|id:run))
#           malt2 <- update(md1_pure, . ~ . - (1|id))
#           exactRLRT(md1_pure, mA=malt1, m0=malt2)
#           exactRLRT(md1_pure, mA=malt2, m0=malt1)
#           
#           if (isSingular(md)) { browser() }
#         } else {
#           md <- lm(h ~ last_outcome + scale(rt_csv) + scale(rt_vmax_lag) + scale(rt_vmax_change) +  v_entropy_wi, d)
#         }
#         
#         
#         dm <- tidy(md) #need broom.mixed for this to work on lmer objects -- otherwise, base broom throws an error
#         dm$slice <- slice
#         dm$side <- side
#         dm$t <- t
#         dm <- dm %>% mutate(stat_order = as.factor(case_when(abs(statistic) < 2 ~ '1', 
#                                                              abs(statistic) > 2 & abs(statistic) < 3 ~ '2', 
#                                                              abs(statistic) > 3 ~ '3')),
#                             p_value = as.factor(case_when(p.value > .05 ~ '1',
#                                                           p.value < .05 & p.value > .01 ~ '2',
#                                                           p.value < .01 & p.value > .001 ~ '3',
#                                                           p.value <.001 ~ '4')))
#         newlist[[paste("hipp", slice, side, t, sep = "_")]]<-dm
#       }
#     }
#   }
#   ddf <- do.call(rbind,newlist)
#   ddf$slice <- as.factor(ddf$slice)
#   ddf$stat_order <- factor(ddf$stat_order, labels = c("NS", "|t| > 2", "|t| > 3"))
#   ddf$p_value <- factor(ddf$p_value, labels = c("NS", "p < .05", "p < .01", "p < .001"))
#   
#   # terms <- names(fixef(md))
#   terms <- names(md$coefficients)
#   if (trial_cont) {
#     setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/hippo/figs/decode')
#   } else {
#     setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/hippo/figs/decode/no_trial_contingency/')
#   }
#   
#   
#   for (fe in terms) {
#     edf <- ddf %>% filter(term == paste(fe) & t < 8)
#     # p1 <- ggplot(edf, aes(t, estimate, color = slice)) + geom_line() + 
#     #   geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error), alpha = .5) + facet_wrap(~side) + 
#     #   theme_dark() + scale_color_viridis_d() + geom_hline(yintercept = 0, lty = "dashed", color = "red")
#     
#     
#     termstr <- str_replace_all(fe, "[^[:alnum:]]", "_")
#     pdf(paste(termstr, ".pdf", sep = ""), width = 12, height = 7)
#     
#     # print(ggplot(edf, aes(t, slice)) + geom_tile(aes(fill = estimate, alpha = abs(statistic)>2), size = 1) + facet_wrap(~side) + 
#     print(ggplot(edf, aes(t, slice)) + geom_tile(aes(fill = estimate, alpha = p_value), size = 1) + facet_wrap(~side) + 
#             scale_fill_viridis(option = "plasma") + scale_color_grey() + labs(title = paste(fe)))
#     # print(ggarrange(p2,ncol = 1, labels = paste(fe), vjust = 4, font.label = list(color = "black", size = 16)))
#     dev.off()
#   }
# }
# 
# 