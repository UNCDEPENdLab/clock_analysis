# loops over RT prediction models for various hippocampal slices and post-feedback time points
library(modelr)
library(tidyverse)
library(lme4)
library(broom)
setwd("/Users/localadmin/Box/SCEPTIC_fMRI/var")
load('feedback_hipp_tallest_by_timepoint_decon.Rdata')
setwd('~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/')
load('trial_df_and_vhdkfpe_clusters.Rdata')

# select relevant columns
df <- df %>% select(id, run, run_trial, rewFunc,emotion, rt_csv, score_csv, rt_next, rt_vmax, v_max_wi, v_entropy_wi, v_entropy_b, v_entropy, v_max_b, Age, Female)
d <- merge(df, fb_wide_t, by = c("id", "run", "run_trial"))
d <- d %>% group_by(id, run) %>% arrange(id, run, run_trial) %>% mutate(reward = score_csv>0) %>% ungroup()

# scale decons
scale2 <- function(x, na.rm = FALSE) (x - mean(x, na.rm = na.rm)) / sd(x, na.rm)
d <- d %>% group_by(id,run) %>%  mutate_at(vars(starts_with("hipp")), scale2, na.rm = TRUE) %>% ungroup()
# test model
# assign("h",d$hipp_9_l_3)

newlist <- list()
# loop over slices and timepoints
for (slice in 1:12) {print(paste("Processing slice", slice, sep = " "))
  for (side in c("l", "r")) {
    for (t in -1:10) {
      d$h<-d[[paste("hipp", slice, side, t, sep = "_")]]
      mf <-  lmer(rt_next ~ scale(-1/run_trial)*rewFunc + (reward + scale(rt_csv) + scale(rt_vmax) + h)^3 + (1|id/run), d)
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

# scale(-1/run_trial)                     -1.17674    7.12535  -0.165
# rewFuncCEVR                             -1.57670    3.09577  -0.509
# rewFuncDEV                              -1.12577    3.21316  -0.350
# rewFuncIEV                               0.08451    2.99803   0.028
# rewardTRUE                               0.05778    0.69452   0.083
# scale(rt_csv)                            0.74853    0.51851   1.444
# scale(rt_vmax)                           5.10485    0.53177   9.600
# h                                       -0.29986    0.57068  -0.525
# scale(-1/run_trial):rewFuncCEVR          3.63437    8.11871   0.448
# scale(-1/run_trial):rewFuncDEV           2.08506    8.06810   0.258
# scale(-1/run_trial):rewFuncIEV           1.61538    7.45955   0.217
# rewardTRUE:scale(rt_csv)                 0.89619    0.64330   1.393
# rewardTRUE:scale(rt_vmax)                1.41641    0.64746   2.188
# rewardTRUE:h                             0.56961    0.66931   0.851
# scale(rt_csv):scale(rt_vmax)             0.72352    0.35888   2.016
# scale(rt_csv):h                          0.46836    0.61410   0.763
# scale(rt_vmax):h                        -0.61496    0.60614  -1.015
# rewardTRUE:scale(rt_csv):scale(rt_vmax) -1.01659    0.46498  -2.186
# rewardTRUE:scale(rt_csv):h              -0.42836    0.78961  -0.542
# rewardTRUE:scale(rt_vmax):h              0.83799    0.77683   1.079
# scale(rt_csv):scale(rt_vmax):h          -0.07784    0.28449  -0.274

edf <- bdf %>% filter(term == "scale(rt_vmax):h" & t < 8)
ggplot(edf, aes(t, estimate, color = slice)) + geom_line() + 
  geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error), alpha = .5) + facet_wrap(~side) + 
  theme_dark() + scale_color_viridis_d() + geom_hline(yintercept = 0, lty = "dashed", color = "red")

ggplot(edf, aes(t, slice)) + geom_tile(aes(fill = estimate))


mf <-  lmer(rt_next ~ scale(-1/run_trial) * scale(rt_csv) +  scale(rt_csv)  + reward * scale(rt_csv)*hipp_1_r_5 + (1|id/run), d)
summary(mf)
car::Anova(mf)
car::vif(mf)
vif.lme(mf)
mh <- lmer(hipp_9_r_4 ~ (scale(-1/run_trial) + scale(rt_csv) + reward + 
                           v_entropy_wi)^2 + v_max_b + v_entropy_b + (1|id/run), d) 
summary(mh)


mb3 <-  lmer(rt_next ~ (scale(-1/run_trial) + scale(rt_csv) + scale(rt_vmax) + reward + 
                          hipp_9_r_5)^3 +  (1|id/run), d)

summary(mb3)
car::Anova(mb)
car::vif(mb)
