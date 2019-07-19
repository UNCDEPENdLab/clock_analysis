# loops over RT prediction models for various hippocampal slices and post-feedback time points
library(modelr)
library(tidyverse)
library(lme4)
library(broom)
library(ggpubr)
setwd("/Users/localadmin/Box/SCEPTIC_fMRI/var")
load('feedback_hipp_tallest_by_timepoint_decon.Rdata')
setwd('~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/')
load('trial_df_and_vhdkfpe_clusters.Rdata')

# read in behavioral data
# select relevant columns for compactness
df <- df %>% select(id, run, run_trial, rewFunc,emotion, rt_csv, score_csv, rt_next, rt_vmax, rt_vmax_lag,rt_vmax_change, v_max_wi, v_entropy_wi, v_entropy_b, v_entropy, v_max_b, Age, Female)
# add deconvolved hippocampal timeseries
d <- merge(df, fb_wide_t, by = c("id", "run", "run_trial"))
d <- d %>% group_by(id, run) %>% arrange(id, run, run_trial) %>% mutate(reward = score_csv>0, 
                                                                        rt_change = 100*rt_next - rt_csv, 
                                                                        v_entropy_wi_lead = lead(v_entropy_wi),
                                                                        v_entropy_wi_change = v_entropy_wi_lead-v_entropy_wi) %>% ungroup()

# scale decons
# SKIP this step if running lmer on h
scale2 <- function(x, na.rm = FALSE) (x - mean(x, na.rm = na.rm)) / sd(x, na.rm)
d <- d %>% group_by(id,run) %>%  mutate_at(vars(starts_with("hipp")), scale2, na.rm = TRUE) %>% ungroup()
# test model
# assign("h",d$hipp_9_l_3)

for (trial_cont in c("TRUE", "FALSE")) {
newlist <- list()
# try to predict directional RT change
# loop over slices and timepoints
for (slice in 1:12) {print(paste("Processing slice", slice, sep = " "))
  for (side in c("l", "r")) {
    for (t in -1:10) {
      d$h<-d[[paste("hipp", slice, side, t, sep = "_")]]
      if (trial_cont) {
        mf <-  lmer(rt_next ~ scale(-1/run_trial)*rewFunc + (reward + scale(rt_csv) + scale(rt_vmax_lag) + scale(rt_vmax) + h)^3 + (1|id/run), d)}
      else {
        mf <-  lmer(rt_next ~ (reward + scale(rt_csv) + scale(rt_vmax_lag) + h)^3 + (1|id/run), d)}
      dm <- tidy(mf)
      dm$slice <- slice
      dm$side <- side
      dm$t <- t
      newlist[[paste("hipp", slice, side, t, sep = "_")]]<-dm
    }
  }
}
bdf <- do.call(rbind,newlist)
bdf$slice <- as.factor(bdf$slice)
terms <- names(fixef(mf))
if (trial_cont) {
  setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/hippo/figs')
}
else {
setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/hippo/figs/no_trial_contingency/')}

for (fe in terms)
{edf <- bdf %>% filter(term == paste(fe) & t < 8)
p1 <- ggplot(edf, aes(t, estimate, color = slice)) + geom_line() + 
  geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error), alpha = .5) + facet_wrap(~side) + 
  theme_dark() + scale_color_viridis_d() + geom_hline(yintercept = 0, lty = "dashed", color = "red") + labs(title = paste(fe))

p2 <- ggplot(edf, aes(t, slice)) + geom_tile(aes(fill = estimate, alpha = abs(statistic)>2), size = 1) + facet_wrap(~side) + 
  scale_fill_viridis(option = "plasma") + scale_color_grey() + labs(title = paste(fe))

termstr <- str_replace_all(fe, "[^[:alnum:]]", "_")
pdf(paste(termstr, ".pdf", sep = ""), width = 12, height = 12)
print(ggarrange(p1,p2,ncol = 1, nrow = 2))
dev.off()
}
}
# "decoding" analyses
# currently running lm
for (trial_cont in c("TRUE", "FALSE")) {
newlist <- list()
# loop over slices and timepoints
for (slice in 1:12) {print(paste("Processing slice", slice, sep = " "))
  for (side in c("l", "r")) {
    for (t in -1:10) {
      d$h<-d[[paste("hipp", slice, side, t, sep = "_")]]
      if (trial_cont) {
        md <-  lm(h ~ scale(-1/run_trial)*rewFunc + reward + scale(rt_csv) + scale(rt_vmax_lag) + scale(rt_vmax_change) + v_entropy_wi, d)
        } else {
        md <-  lm(h ~ reward + scale(rt_csv) + scale(rt_vmax_lag) + scale(rt_vmax_change) +  v_entropy_wi, d)}
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
                                                        p.value <.001 ~ '4')))
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
if (trial_cont) {
  setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/hippo/figs/decode')
} else {
  setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/hippo/figs/decode/no_trial_contingency/')}


for (fe in terms)
{edf <- ddf %>% filter(term == paste(fe) & t < 8)
# p1 <- ggplot(edf, aes(t, estimate, color = slice)) + geom_line() + 
#   geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error), alpha = .5) + facet_wrap(~side) + 
#   theme_dark() + scale_color_viridis_d() + geom_hline(yintercept = 0, lty = "dashed", color = "red")



termstr <- str_replace_all(fe, "[^[:alnum:]]", "_")
pdf(paste(termstr, ".pdf", sep = ""), width = 12, height = 7)
# print(ggplot(edf, aes(t, slice)) + geom_tile(aes(fill = estimate, alpha = abs(statistic)>2), size = 1) + facet_wrap(~side) + 
print(ggplot(edf, aes(t, slice)) + geom_tile(aes(fill = estimate, alpha = p_value), size = 1) + facet_wrap(~side) + 
  scale_fill_viridis(option = "plasma") + scale_color_grey() + labs(title = paste(fe)))
# print(ggarrange(p2,ncol = 1, labels = paste(fe), vjust = 4, font.label = list(color = "black", size = 16)))
dev.off()
}
}




# mf <-  lmer(rt_next ~ scale(-1/run_trial) * scale(rt_csv) +  scale(rt_csv)  + reward * scale(rt_csv)*hipp_1_r_5 + (1|id/run), d)
# summary(mf)
# car::Anova(mf)
# car::vif(mf)
# vif.lme(mf)
# mh <- lmer(hipp_9_r_4 ~ (scale(-1/run_trial) + scale(rt_csv) + reward + 
#                            v_entropy_wi)^2 + v_max_b + v_entropy_b + (1|id/run), d) 
# summary(mh)
# 
# 
# mb3 <-  lmer(rt_next ~ (scale(-1/run_trial) + scale(rt_csv) + scale(rt_vmax) + reward + 
#                           hipp_9_r_5)^3 +  (1|id/run), d)
# 
# summary(mb3)
# car::Anova(mb)
# car::vif(mb)
