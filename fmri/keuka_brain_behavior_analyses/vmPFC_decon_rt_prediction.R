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
decode = F  # main analysis for Fig. 4 E-G
u = T       # exploratory analysis attempting to predict the uncertainty of the next choice


# load data
# if (unsmoothed) {setwd("~/Box/SCEPTIC_fMRI/var/unsmoothed")
#   # } else {setwd("/Users/localadmin/Box/SCEPTIC_fMRI/var/newmask")}
# } else if (smooth_in_mask) {setwd("~/Box/SCEPTIC_fMRI/var/smooth_in_mask/")
#   } else {setwd("~/Box/SCEPTIC_fMRI/var/")}
setwd('~/Box/skinner/data/clock_task/vmPFC/')
# 
load('feedback_vmPFC_widest_by_timepoint_decon.Rdata')
# load('feedback_vmPFC_wide_ts.Rdata')
setwd('~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/')
# load('trial_df_and_vhdkfpe_clusters.Rdata')
load('trial_df_and_vh_pe_clusters_u.Rdata')
if (unsmoothed) {cache_dir <- "~/Box/SCEPTIC_fMRI/var/unsmoothed"
# } else {cache_dir <- "~/Box/SCEPTIC_fMRI/var/newmask"}
} else {cache_dir <- "~/Box/SCEPTIC_fMRI/var"}
repo_dir <- "~/code/clock_analysis"

#super-wide variant used in lm analysis
# load(file.path(cache_dir, "feedback_hipp_tallest_by_timepoint_decon.Rdata"))
# 
# load(file.path(repo_dir, "fmri/keuka_brain_behavior_analyses/trial_df_and_vhdkfpe_clusters.Rdata"))

# load(file = file.path(cache_dir,"sceptic_trial_df_for_medusa.RData"))
# df <- trial_df
# attr(df, "labels") <- NULL #somehow this is holding a 560-row data.frame
# df <- df %>% dplyr::ungroup()
# obtain within-subject v_max and entropy: correlated at -.37

# trial_df <- trial_df %>% group_by(id,run) %>% mutate(v_max_wi = scale(v_max),
#                                          v_max_wi_lag = lag(v_max_wi),
#                                          v_entropy_wi = scale(v_entropy),
#                                          v_max_b = mean(na.omit(v_max)),
#                                          v_entropy_b = mean(na.omit(v_entropy)),
#                                          rt_change = rt_csv - rt_lag,
#                                          pe_max_lag = lag(pe_max),
#                                          abs_pe_max_lag = abs(pe_max_lag),
#                                          rt_vmax_change = rt_vmax - rt_vmax_lag,
#                                          trial_neg_inv_sc = scale(-1/run_trial),
#                                          v_chosen_change = v_chosen - lag(v_chosen)) %>% ungroup() %>%
#   mutate(rt_lag_sc = scale(rt_lag),
#          rt_csv_sc = scale(rt_csv),
#          rt_vmax_lag_sc = scale(rt_vmax_lag))
# 
# u_df <- read_csv("~/Box/SCEPTIC_fMRI/sceptic_model_fits/mmclock_fmri_fixed_uv_ureset_fixedparams_fmri_ffx_trial_statistics.csv.gz")
# u_df <- u_df %>% select(id, run, trial, u_chosen, u_chosen_lag, u_chosen_change,
#                         u_chosen_quantile, u_chosen_quantile_lag, u_chosen_quantile_change,
#                         v_chosen_quantile, v_chosen_quantile_lag, v_chosen_quantile_change)
# 
# df <- inner_join(trial_df,u_df)


# read in behavioral data
# select relevant columns for compactness
df <- df %>% select(id, run, run_trial, rewFunc,emotion, last_outcome, rt_csv, score_csv, rt_next, pe_max, rt_vmax, rt_vmax_lag,
                    rt_vmax_change, v_max_wi, v_entropy_wi, v_entropy_b, v_entropy, v_max_b, u_chosen_quantile, u_chosen_quantile_lag, u_chosen_quantile_change, 
                    rt_vmax_lag_sc, rt_lag_sc, rt_csv_sc, trial_neg_inv_sc, Age, Female)
# add deconvolved timeseries
d <- merge(df, fb_wide_t, by = c("id", "run", "run_trial"))
d <- d %>% group_by(id, run) %>% arrange(id, run, run_trial) %>% mutate(rt_change = 100*rt_next - rt_csv, 
                                                                        v_entropy_wi_lead = lead(v_entropy_wi),
                                                                        v_entropy_wi_change = v_entropy_wi_lead-v_entropy_wi,
                                                                        u_chosen_quantile_next = lead(u_chosen_quantile),
                                                                        u_chosen_quantile_change_next = lead(u_chosen_quantile_change),
                                                                        outcome = lead(last_outcome)) %>% ungroup()
# dbl <- merge(df, fb_wide_bl, by = c("id", "run", "run_trial"))
# dbl <- dbl %>% group_by(id, run) %>% arrange(id, run, run_trial) %>% mutate(rt_change = 100*rt_next - rt_csv, 
#                                                                         v_entropy_wi_lead = lead(v_entropy_wi),
#                                                                         v_entropy_wi_change = v_entropy_wi_lead-v_entropy_wi,
#                                                                         u_chosen_quantile_next = lead(u_chosen_quantile),
#                                                                         u_chosen_quantile_change_next = lead(u_chosen_quantile_change),
#                                                                         outcome = lead(last_outcome)) %>% ungroup()

scale2 <- function(x, na.rm = FALSE) (x - mean(x, na.rm = na.rm)) / sd(x, na.rm)
# scale decon across subjects as a predictor
# choice uncertainty prediction analyses run on scaled 'ds' dataframe instead of 'd'
ds <- d %>% mutate_at(vars(starts_with("vmPFC")), scale2, na.rm = TRUE) %>% ungroup()


#######
# "decoding" analyses
# combined right and left hippocampus with side as a predictor
# if model does not converge, update with new starting values

if (decode) {
  newlist <- list()
  for (slice in 1:12) {print(paste("Processing slice", slice, sep = " "))
    # for (side in c("l", "r")) {
      for (t in -1:10) {
        d$h<-d[[paste("vmPFC", slice, t, sep = "_")]]
        md <-  lmer(h ~ trial_neg_inv_sc + scale(rt_csv) + scale(rt_vmax_lag)  + scale(rt_vmax_change)  + v_entropy_wi  + v_entropy_wi_change  + v_max_wi + #u_chosen_quantile_change +
                      (1|id), d, control=lmerControl(optimizer = "nloptwrap"))
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
        newlist[[paste("vmPFC", slice, t, sep = "_")]]<-dm
      }
    # }
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
  } else if (smooth_in_mask) {setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/vmPFC/figs/lmer')
  } else {setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/vmPFC/figs/')} # manually indicate if this is the new COBRA lab mask
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


