# loops over RT prediction models for various hippocampal slices and post-feedback time points
library(modelr)
library(tidyverse)
library(lme4)
library(afex)
library(broom)
library(broom.mixed) #plays will with afex p-values in lmer wrapper
library(ggpubr)
library(car)

setwd("/Users/localadmin/Box/SCEPTIC_fMRI/var")
load('feedback_hipp_tallest_by_timepoint_decon.Rdata')
setwd('~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/')
load('trial_df_and_vhdkfpe_clusters.Rdata')

cache_dir <- "~/Box/SCEPTIC_fMRI/var"
repo_dir <- "~/Data_Analysis/clock_analysis"

#super-wide variant used in lm analysis
# load(file.path(cache_dir, "feedback_hipp_tallest_by_timepoint_decon.Rdata"))
# 
# load(file.path(repo_dir, "/fmri/keuka_brain_behavior_analyses/trial_df_and_vhdkfpe_clusters.Rdata"))

attr(df, "labels") <- NULL #somehow this is holding a 560-row data.frame
df <- df %>% dplyr::ungroup()

# read in behavioral data
# select relevant columns for compactness
df <- df %>% select(id, run, run_trial, rewFunc,emotion, rt_csv, score_csv, rt_next, pe_max, rt_vmax, rt_vmax_lag,rt_vmax_change, v_max_wi, v_entropy_wi, v_entropy_b, v_entropy, v_max_b, Age, Female)
# add deconvolved hippocampal timeseries
d <- merge(df, fb_wide_t, by = c("id", "run", "run_trial"))
d <- d %>% group_by(id, run) %>% arrange(id, run, run_trial) %>% mutate(reward = score_csv > 0, 
                                                                        rt_change = 100*rt_next - rt_csv, 
                                                                        v_entropy_wi_lead = lead(v_entropy_wi),
                                                                        v_entropy_wi_change = v_entropy_wi_lead-v_entropy_wi) %>% ungroup()

# scale decons
# SKIP this step if running lmer on h

scale = T
if (scale) {
  scale2 <- function(x, na.rm = FALSE) (x - mean(x, na.rm = na.rm)) / sd(x, na.rm)
  # d <- d %>% group_by(id,run) %>%  mutate_at(vars(starts_with("hipp")), scale2, na.rm = TRUE) %>% ungroup()
  # try scaling across subjects
  d <- d %>% mutate_at(vars(starts_with("hipp")), scale2, na.rm = TRUE) %>% ungroup()
}

plots = T
rt = T
##########
# predict directional RT change
# try removing the entropy terms for simplicity
# strategically: pe*rt_t (or reward*rt) - stay vs. shift, rt - explore in general, rt_vmax - exploit, prior, rt_vmax_change - exploit, posterior, that's it!
# for (trial_cont in c("TRUE", "FALSE")) {
if (rt) {
  for (trial_cont in c("TRUE", "FALSE")) {
    newlist <- list()
    for (slice in 1:12) {print(paste("Processing slice", slice, sep = " "))
      for (side in c("l", "r")) {
        for (t in -1:10) {
          d$h<-d[[paste("hipp", slice, side, t, sep = "_")]]
          if (trial_cont) {
            mf <-  lmerTest::lmer(rt_next ~ scale(-1/run_trial)*rewFunc + scale(rt_csv) * scale(pe_max) * h + scale(rt_vmax_lag) *  h + scale(rt_vmax_change) * h + (1|id/run), d)
          }
          else {
            # mf <-  lme4::lmer(rt_next ~ (scale(pe_max) + scale(rt_csv) + scale(rt_vmax_lag) + scale(rt_vmax_change) + scale(v_entropy_wi) + h)^2 + (1|id/run), d)
            mf <-  lmerTest::lmer(rt_next ~ scale(rt_csv) * scale(pe_max) * h + scale(rt_vmax_lag) *  h + scale(rt_vmax_change) *  h + (1|id/run), d)
          }
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
      setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/hippo/figs/rt_predict')
    }
    else {
      setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/hippo/figs/rt_predict/no_trial_contingency/')}
    
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
# "decoding" analyses
for (trial_cont in c("F")) {
  newlist <- list()
  for (slice in 1:12) {print(paste("Processing slice", slice, sep = " "))
    for (side in c("l", "r")) {
      for (t in -1:10) {
        d$h<-d[[paste("hipp", slice, side, t, sep = "_")]]
        if (trial_cont) {
          md <-  lm(h ~ scale(-1/run_trial)*rewFunc + reward + scale(rt_csv) + scale(rt_vmax_lag) + scale(rt_vmax_change) + v_entropy_wi + v_entropy_wi_change, d)
        } else {
          md <-  lm(h ~ reward + scale(rt_csv) + scale(rt_vmax_lag) + scale(rt_vmax_change) +  v_entropy_wi + v_entropy_wi_change, d)}
        dm <- tidy(md)
        dm$slice <- slice
        dm$side <- side
        dm$t <- t
        dm <- dm %>% mutate(stat_order = as.factor(case_when(abs(statistic) < 2 ~ '1', 
                                                             abs(statistic) > 2 & abs(statistic) < 3 ~ '2', 
                                                             abs(statistic) > 3 ~ '3')),
                            p_value = as.factor(case_when(p.value > .05 ~ '1',
                                                          p.value < .05 & p.value > .01 ~ '2',
                                                          p.value < .01 & p.value > .001 ~ '3',
                                                          p.value <.001 ~ '4'))
        )
        newlist[[paste("hipp", slice, side, t, sep = "_")]]<-dm
      }
    }
  }
  ddf <- do.call(rbind,newlist)
  ddf$slice <- as.factor(ddf$slice)
  ddf$stat_order <- factor(ddf$stat_order, labels = c("NS", "|t| > 2", "|t| > 3"))
  ddf$p_value <- factor(ddf$p_value, labels = c("NS", "p < .05", "p < .01", "p < .001"))
  
  # terms <- names(fixef(md))
  terms <- names(md$coefficients)
  # FDR
  ddf <- ddf  %>% group_by(term,side) %>% mutate(p_fdr = p.adjust(p.value, method = 'fdr'),
                 p_level_fdr = as.factor(case_when(p_fdr > .05 ~ '1',
                                                   p_fdr < .05 & p_fdr > .01 ~ '2',
                                                   p_fdr < .01 & p_fdr > .001 ~ '3',
                                                   p_fdr <.001 ~ '4')))
  
  ddf$p_level_fdr <- factor(ddf$p_level_fdr, labels = c("NS", "p < .05", "p < .01", "p < .001"))
  if (trial_cont) {
    setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/hippo/figs/decode')
  } else {
    setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/hippo/figs/decode/no_trial_contingency/')
  }
  for (fe in terms) {
    edf <- ddf %>% filter(term == paste(fe) & t < 8) 
    termstr <- str_replace_all(fe, "[^[:alnum:]]", "_")
    pdf(paste(termstr, ".pdf", sep = ""), width = 12, height = 7)
    
    # print(ggplot(edf, aes(t, slice)) + geom_tile(aes(fill = estimate, alpha = abs(statistic)>2), size = 1) + facet_wrap(~side) + 
    print(ggplot(edf, aes(t, slice)) + geom_tile(aes(fill = estimate, alpha = p_level_fdr), size = 1) + facet_wrap(~side) + 
            scale_fill_viridis(option = "plasma") + scale_color_grey() + labs(title = paste(fe)))
    # print(ggarrange(p2,ncol = 1, labels = paste(fe), vjust = 4, font.label = list(color = "black", size = 16)))
    dev.off()
  }
}

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
dfc <- dfc %>% group_by(id, run) %>% arrange(id, run, run_trial) %>% mutate(reward = score_csv > 0, 
                                                                            rt_change = 100*rt_next - rt_csv, 
                                                                            v_entropy_wi_lead = lead(v_entropy_wi),
                                                                            v_entropy_wi_change = v_entropy_wi_lead-v_entropy_wi) %>% ungroup()

# God help me -- centrality loop
newlist <- list()
for (slice in c(2,4,6,8,10,12)) {print(paste("Processing slice", slice, sep = " "))
  for (side in c("lh", "rh")) {
    for (dir in c("out", "in")) {
      dfc$c <- dfc[[paste(paste("strength", dir, side, sep = "_"),slice, sep = "")]]
      mf <-  lmerTest::lmer(rt_next ~ scale(-1/run_trial)*rewFunc + scale(rt_csv)  * reward + scale(rt_csv)  * c + scale(rt_vmax_lag) *  c + scale(rt_vmax_change) *  c + (1|id/run), dfc)
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
#          mf <-  lmer(rt_next ~ scale(-1/run_trial)*rewFunc + (reward + scale(rt_csv) + scale(rt_vmax_lag) + scale(rt_vmax) + h)^3 + (1|id/run), d)}
#        else {
#          mf <-  lmer(rt_next ~ (reward + scale(rt_csv) + scale(rt_vmax_lag) + h)^3 + (1|id/run), d)}
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
#           mform <- as.formula(h ~ scale(-1/run_trial)*rewFunc + reward + scale(rt_csv) + scale(rt_vmax_lag) + scale(rt_vmax_change) + v_entropy_wi + (1|id))
#           #TODO: Why does lmer return singular fit?
#           xtabs(~id + run_trial, d) #looks complete
#           d %>% group_by(id) %>% summarize(mean(h, na.rm=T), sd(h, na.rm=T))
#           
#           lattice::histogram(d$h)
#           md1 <-  lmer(h ~ scale(-1/run_trial)*rewFunc + reward + scale(rt_csv) + scale(rt_vmax_lag) + scale(rt_vmax_change) + v_entropy_wi + (1|id/run), d)
#           md2 <-  lmer(h ~ scale(-1/run_trial)*rewFunc + reward + scale(rt_csv) + scale(rt_vmax_lag) + scale(rt_vmax_change) + v_entropy_wi + (1|id), d)
#           md3 <-  lmer(h ~ scale(-1/run_trial)*rewFunc + reward + scale(rt_csv) + scale(rt_vmax_lag) + scale(rt_vmax_change) + v_entropy_wi + (1|id:run), d)
#           
#           md2 <-  lmer(h ~ 1 + (1|id), d)
#           
#           md3 <-  lm(h ~ scale(-1/run_trial)*rewFunc + reward + scale(rt_csv) + scale(rt_vmax_lag) + scale(rt_vmax_change) + v_entropy_wi, d)
#           
#           anova(md1, md2, md3)
#           anova(md1, md3)
#           library(RLRsim)
#           
#           md1_pure <- lme4::lmer(h ~ scale(-1/run_trial)*rewFunc + reward + scale(rt_csv) + scale(rt_vmax_lag) + scale(rt_vmax_change) + v_entropy_wi + (1|id/run), d, REML=TRUE)
#           exactLRT(md1_pure, md3)
#           
#           malt1 <- update(md1_pure, . ~ . - (1|id:run))
#           malt2 <- update(md1_pure, . ~ . - (1|id))
#           exactRLRT(md1_pure, mA=malt1, m0=malt2)
#           exactRLRT(md1_pure, mA=malt2, m0=malt1)
#           
#           if (isSingular(md)) { browser() }
#         } else {
#           md <- lm(h ~ reward + scale(rt_csv) + scale(rt_vmax_lag) + scale(rt_vmax_change) +  v_entropy_wi, d)
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