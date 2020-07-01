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
rt = T
u = F       # exploratory analysis attempting to predict the uncertainty of the next choice with hippocampal activity


# load data
if (unsmoothed) {setwd("~/Box/SCEPTIC_fMRI/var/unsmoothed")
  # } else {setwd("/Users/localadmin/Box/SCEPTIC_fMRI/var/newmask")}
} else if (smooth_in_mask) {setwd("~/Box/SCEPTIC_fMRI/var/smooth_in_mask/")
} else {setwd("~/Box/SCEPTIC_fMRI/var/")}

# load deconvolved signals
load('feedback_hipp_widest_by_timepoint_decon.Rdata')
load('feedback_hipp_wide_ts.Rdata')
setwd('~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/')
# load('trial_df_and_vhdkfpe_clusters.Rdata')

# load trial-level behavioral variables
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
                                                                        rt_swing_lead = abs(rt_change),
                                                                        v_entropy_wi_lead = lead(v_entropy_wi),
                                                                        v_entropy_wi_change = v_entropy_wi_lead-v_entropy_wi,
                                                                        u_chosen_quantile_next = lead(u_chosen_quantile),
                                                                        u_chosen_quantile_change_next = lead(u_chosen_quantile_change),
                                                                        outcome = lead(last_outcome)) %>% ungroup()
dbl <- merge(df, fb_wide_bl, by = c("id", "run", "run_trial"))
dbl <- dbl %>% group_by(id, run) %>% arrange(id, run, run_trial) %>% mutate(rt_change = 100*rt_next - rt_csv, 
                                                                            rt_swing_lead = abs(rt_change),
                                                                            v_entropy_wi_lead = lead(v_entropy_wi),
                                                                            v_entropy_wi_change = v_entropy_wi_lead-v_entropy_wi,
                                                                            u_chosen_quantile_next = lead(u_chosen_quantile),
                                                                            u_chosen_quantile_change_next = lead(u_chosen_quantile_change),
                                                                            outcome = lead(last_outcome),
                                                                            abs_pe_sc = scale(abs(pe_max))) %>% ungroup()

scale2 <- function(x, na.rm = FALSE) (x - mean(x, na.rm = na.rm)) / sd(x, na.rm)
# scale decon across subjects as a predictor
# choice uncertainty prediction analyses run on scaled 'ds' dataframe instead of 'd'
ds <- d %>% mutate_at(vars(starts_with("hipp")), scale2, na.rm = TRUE) %>% ungroup()


#######
# "decoding" analyses
# combined right and left hippocampus with side as a predictor
# if model does not converge, update with new starting values

# added RPE to the model: results are as expected, early posterior RPEs

if (decode) {
  newlist <- list()
  for (slice in 1:12) {print(paste("Processing slice", slice, sep = " "))
    # for (side in c("l", "r")) {
    for (t in -1:10) {
      dbl$h<-dbl[[paste("hipp", slice, t, sep = "_")]]
      md <-  lmer(h ~  scale(rt_csv) + side  + scale(rt_vmax_lag)  + scale(rt_vmax_change)  + v_entropy_wi  + v_entropy_wi_change  + pe_max +#u_chosen_quantile_change +
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
  if (unsmoothed) {setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/str_hippo/figs/decode/unsmoothed')
  } else if (smooth_in_mask) {setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/str_hippo/figs/decode')
  } else {setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/str_hippo/figs/decode/smoothed_whole_brain')} # manually indicate if this is the new COBRA lab mask
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

##########
# Predict RT swings

#########
# Diagnose multicollinearity
ivs <- dbl %>% select(rt_csv, rt_vmax, abs_pe_sc, outcome)

if (rt) {
  newlist <- list()
  for (slice in 1:12) {print(paste("Processing slice", slice, sep = " "))
    for (side in c("l", "r")) {
      for (t in 0:10) {
        dbl$h<-dbl[[paste("hipp", slice, t, sep = "_")]]
        # mf <-  lme4::lmer(rt_next ~ (scale(pe_max) + scale(rt_csv) + scale(rt_vmax_lag) + scale(rt_vmax_change) + scale(v_entropy_wi) + h)^2 + (1|id/run), ds)
        # mf <-  lmerTest::lmer(rt_next ~ scale(rt_csv) * last_outcome * h + scale(rt_vmax_lag) *  h + scale(rt_vmax_change) *  h + (1|id/run), ds)
        # a more detailed model
        # decompose pe into valence and magnitude
        mf <-  lmer(scale(rt_next) ~ #(trial_neg_inv_sc + rt_lag_sc + rt_vmax_lag_sc + pe_max + h)^2 + 
                      scale(rt_csv)*h + scale(rt_vmax)*h +
                      scale(rt_csv)*scale(abs_pe_sc)*h + 
                      scale(rt_csv)*outcome*h + 
                      # rt_vmax_lag_sc:trial_neg_inv_sc:h + 
                      (1|id/run), dbl)
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
  ddf <- do.call(rbind,newlist)
  ddf$slice <- as.factor(ddf$slice)
  ddf$stat_order <- factor(ddf$stat_order, labels = c("NS", "|t| > 2", "|t| > 3"))
  ddf$p_value <- factor(ddf$p_value, labels = c("NS", "p < .05", "p < .01", "p < .001"))
  
  terms <- names(fixef(mf))
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
  if (unsmoothed) {setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/str_hippo/figs/rt_predict/unsmoothed')
    # } else {setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/hippo/figs/newmask/rt_predict/no_trial_contingency/')}}
  } else {setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/str_hippo/figs/rt_predict/')}
  
  if (plots) {
    for (fe in terms) {
      edf <- ddf %>% filter(term == paste(fe) & t < 8) 
      termstr <- str_replace_all(fe, "[^[:alnum:]]", "_")
      pdf(paste(termstr, "_bl.pdf", sep = ""), width = 5, height = 3.5)
      print(ggplot(edf, aes(t, slice)) + geom_tile(aes(fill = estimate, alpha = `p, FDR-corrected`), size = 1) +  
              scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab("Time after outcome, seconds") + ylab("Posterior <-- Location --> Anterior\n (12 slices)") + 
              labs(alpha = expression(~italic(p)~', FDR-corrected'))) 
      dev.off()
      
      # for (fe in terms)
      # {edf <- bdf %>% filter(term == paste(fe) & t < 8)
      # p1 <- ggplot(edf, aes(t, estimate, color = slice)) + geom_line() + 
      #   geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error), alpha = .5) + facet_wrap(~side) + 
      #   theme_dark() + scale_color_viridis_d() + geom_hline(yintercept = 0, lty = "dashed", color = "red") + labs(title = paste(fe))
      # 
      # p2 <- ggplot(edf, aes(t, slice)) + geom_tile(aes(fill = estimate, alpha = p_level_fdr), size = 1) + facet_wrap(~side) + 
      #   scale_fill_viridis(option = "plasma") + scale_color_grey() + labs(title = paste(fe))
      # 
      # termstr <- str_replace_all(fe, "[^[:alnum:]]", "_")
      # pdf(paste(termstr, ".pdf", sep = ""), width = 12, height = 12)
      # print(ggarrange(p1,p2,ncol = 1, nrow = 2))
      # dev.off()
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


