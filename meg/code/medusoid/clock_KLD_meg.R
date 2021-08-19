library(LaplacesDemon)
library(tidyverse)
library(lme4)
library(ggpubr)
library(car)


repo_directory <- "~/code/clock_analysis/meg"
setwd(repo_directory)
df <- readRDS("MEG_n63_behavioral_data_preprocessed_trial_df.rds")


get_kldsum <- function(v1, v2) {
  require(LaplacesDemon)
  stopifnot(length(v1) == length(v2))
  if (any(is.na(v1)) || any(is.na(v2))) { return(NA_real_) }
  kk <- KLD(v1, v2)
  return(kk$sum.KLD.px.py)
}
df <- df %>% group_by(id, run) %>% arrange(id, run, run_trial) %>% mutate(
  rt_lag2 = lag(rt_lag),
  rt_lag3 = lag(rt_lag2),
  rt_lag4 = lag(rt_lag3),
  rt_lag5 = lag(rt_lag4),
  rt_swing = abs(rt_csv_sc - rt_lag_sc),
  rt_swing_lag2 = lag(rt_swing_lag)
  ) %>% ungroup() %>%
  rowwise() %>% mutate(
    kld4 = get_kldsum(c(rt_lag4, rt_lag3, rt_lag2, rt_lag), c(rt_lag5, rt_lag4, rt_lag3, rt_lag2)),
    kld3 = get_kldsum(c(rt_lag3, rt_lag2, rt_lag), c(rt_lag4, rt_lag3, rt_lag2))) %>%
  ungroup() %>% mutate(rt_swing_lag_sc = scale(rt_swing_lag),
                       rt_swing_lag2_sc = scale(rt_swing_lag2))

# # inspect: r(KLD3, KLD4) = .88, similar timecourses; r(kld_rew3, kld_rew4) = 0.86
# # r(kld3, rt_swing_lag) = .53
# p1 <- ggplot(df, aes(run_trial, v_entropy, color = rewFunc)) + geom_smooth(method = 'loess') + facet_wrap(~total_earnings>median(total_earnings))
# p2 <- ggplot(df, aes(run_trial, kld_rew3, color = rewFunc)) + geom_smooth(method = 'loess') + facet_wrap(~total_earnings>median(total_earnings))
# p3 <- ggplot(df, aes(run_trial, kld3, color = rewFunc)) + geom_smooth(method = 'loess') + facet_wrap(~total_earnings>median(total_earnings))
# p4 <- ggplot(df, aes(run_trial, rt_swing, color = rewFunc)) + geom_smooth(method = 'loess') + facet_wrap(~total_earnings>median(total_earnings))
# pdf("kld_diagnostics.pdf", height = 6, width = 8)
# ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2)
# dev.off()
p1 <- ggplot(df, aes(run_trial, v_entropy, color = rewFunc)) + geom_smooth(method = 'loess') + facet_wrap(~total_earnings>median(total_earnings))
p2 <- ggplot(df, aes(run_trial, kld3, color = rewFunc)) + geom_smooth(method = 'loess') + facet_wrap(~total_earnings>median(total_earnings))
p3 <- ggplot(df, aes(run_trial, kld3_cum2, color = rewFunc)) + geom_smooth(method = 'loess') + facet_wrap(~total_earnings>median(total_earnings))
p4 <- ggplot(df, aes(run_trial, rt_swing, color = rewFunc)) + geom_smooth(method = 'loess') + facet_wrap(~total_earnings>median(total_earnings))
pdf("kld_cum_diagnostics.pdf", height = 6, width = 8)
ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2)
dev.off()

# does entropy predict RT above and beyond KLD?
# their correlation is only 0.15
# just selection history KLD
m1 <- lmer(rt_csv_sc ~ rt_lag_sc*kld3 + rt_lag_sc*omission_lag + (1|ID), df )
summary(m1)
# reinforcement history KLD: not much going on
m2 <- lmer(rt_csv_sc ~ rt_lag_sc*kld3*kld_rew3 + rt_lag_sc*omission_lag*kld3 + rt_lag_sc*omission_lag*kld_rew3 + (1|ID), df )
summary(m2)
# add entropy
m3 <- lmer(rt_csv_sc ~ rt_lag_sc*omission_lag*kld3 +  rt_lag_sc*omission_lag*kld_rew3  + rt_lag_sc*v_entropy_wi +  (1|ID), df )
summary(m3)
# add RT_vmax, remove reward KLD
m4 <- lmer(rt_csv_sc ~ rt_lag_sc*omission_lag*kld3 + rt_vmax_lag_sc*kld3 + 
             rt_lag_sc*v_entropy_wi +  rt_vmax_lag_sc*v_entropy_wi + (1|ID), df )
summary(m4)
# add lagged RT swing
m5 <- lmer(rt_csv_sc ~ rt_swing_lag_sc*rt_lag_sc + rt_lag_sc*omission_lag*kld3 + rt_vmax_lag_sc*kld3 + 
             rt_lag_sc*v_entropy_wi +  rt_vmax_lag_sc*v_entropy_wi + rt_vmax_lag_sc*rt_swing_lag_sc + (1|ID), df )
summary(m5)
Anova(m5)
# remove 3-way interaction rt_lag*reward*KLD
m6 <- lmer(rt_csv_sc ~ rt_lag_sc*omission_lag + rt_swing_lag_sc*rt_lag_sc + rt_lag_sc*kld3 + rt_vmax_lag_sc*kld3 + 
             rt_lag_sc*v_entropy_wi +  rt_vmax_lag_sc*v_entropy_wi + rt_vmax_lag_sc*rt_swing_lag_sc + (1|ID), df )
summary(m6)
Anova(m6, '3')

# cumulative KLD of the last two updates, first vs. rt_swing_lag
# check correlations: 
mc1 <- lmer(rt_csv_sc ~ rt_lag_sc*omission_lag + rt_swing_lag_sc*rt_lag_sc + rt_lag_sc*kld3_cum2 + rt_vmax_lag_sc*kld3_cum2 + 
              rt_lag_sc*v_entropy_wi +  rt_vmax_lag_sc*v_entropy_wi + rt_vmax_lag_sc*rt_swing_lag_sc + (1|ID), df )
summary(mc1)
Anova(mc1, '3')
# add rt_swing_lag2 -- a lot to ask!
# r(kld_cum2, rt_swing_lag) = .39
# r(kld_cum2, rt_swing_lag2) = .61
mk <- lmer(kld3_cum2 ~ rt_swing_lag_sc + rt_swing_lag2_sc + (1|ID), df)
summary(mk)
mc2 <- lmer(rt_csv_sc ~ rt_lag_sc*omission_lag + rt_swing_lag_sc*rt_lag_sc + rt_swing_lag2_sc*rt_lag_sc + rt_lag_sc*kld3_cum2 + rt_vmax_lag_sc*kld3_cum2 + 
              rt_lag_sc*v_entropy_wi +  rt_vmax_lag_sc*v_entropy_wi + rt_vmax_lag_sc*rt_swing_lag_sc + (1|ID), df )
summary(mc2)
Anova(mc2, '3')

