#signal correlation examination
library(dplyr)
library(tidyverse)
global_file <- read.csv("/Users/mnh5174/Data_Analysis/clock_analysis/fmri/data/mmclock_fmri_decay_factorize_selective_psequate_mfx_sceptic_global_statistics.csv")

#group fixed parameters (looks good)
#global_file <- read.csv("/Users/mnh5174/Data_Analysis/clock_analysis/fmri/data/mmclock_fmri_decay_factorize_selective_psequate_fixedparams_ffx_sceptic_global_statistics.csv")

hist(global_file$alpha)
hist(global_file$gamma)
hist(global_file$beta)
mean(global_file$alpha)
mean(global_file$gamma)
mean(global_file$beta)

trial_df <- read.csv("/Users/mnh5174/Data_Analysis/clock_analysis/fmri/data/mmclock_fmri_decay_factorize_selective_psequate_mfx_trial_statistics.csv.gz") %>%
#trial_df <- read.csv("/Users/mnh5174/Data_Analysis/clock_analysis/fmri/data/mmclock_fmri_decay_factorize_mfx_trial_statistics.csv.gz") %>%
  filter(!id %in% c(11335, 11332, 11282, 11246, 10662)) %>% mutate(d_auc = -1*d_auc)

trial_df <- trial_df %>%
  group_by(id, run) %>%  dplyr::mutate(rt_swing = abs( c(NA, diff(rt_csv)) ), run_trial=1:50) %>% ungroup() #compute rt_swing within run and subject


sum_df <- trial_df %>% group_by(id) %>% summarize(total_earnings = sum(score_csv)) %>% arrange(total_earnings)
ids_to_select <- sum_df %>% top_n(5) %>% bind_rows(sum_df %>% top_n(-5)) %>% pull(id)
head(sum_df, n=5)
tail(sum_df, n=5)

trial_df_subset <- trial_df %>% filter(id %in% ids_to_select)

library(ggplot2)
trial_df_long <- trial_df_subset %>% mutate_at(vars(v_auc, d_auc, v_entropy), scale) %>% gather(key="signal", value="sig_val", v_auc, d_auc, v_entropy) %>%
  mutate(id_type=case_when(
    id %in% c(11243, 11162, 11322, 11262, 11298) ~ paste("high", id, sep=":"),
    TRUE ~ paste("low", id, sep=":")
  ))

trial_df_long$run_trial <- 1:50

trial_df_long$run_condition <- with(trial_df_long, paste(rewFunc, emotion, sep=":"))
ggplot(trial_df_long, aes(x=run_trial, y=sig_val, color=signal)) + geom_point() + stat_smooth() + facet_wrap(~id_type, nrow=2) #geom_line() +


ggplot(trial_df_long %>% filter(id %in% c(11243, 11162, 11322, 11262, 11298)), aes(x=run_trial, y=sig_val, color=run_condition)) +
  geom_point() + geom_line() + stat_smooth() + facet_grid(signal~id_type)

xx <- split(trial_df_subset, trial_df_subset$id)
lapply(xx, function(subdf) {
  cor(subset(subdf, select=c(v_auc, v_entropy, d_auc)), use="pairwise.complete.obs")
})

#d_auc driven by rt_swing?

library(lme4)
#trial_df_subset <- trial_df_subset %>%

m1 <- lmer(rt_swing ~ 1 + run_trial + d_auc + v_entropy + (1 | id/run), trial_df_subset )
summary(m1)

m2 <- lmer(rt_swing ~ 1 + run_trial + d_auc + (1 | id/run), trial_df_subset )
summary(m2)

m3 <- lmer(d_auc ~ 1 + run_trial + v_entropy + rt_swing + (1 | id/run), trial_df_subset )
summary(m3)

summary(lmer(rt_swing ~ 1 + run_trial + d_auc + (1 | id/run), filter(trial_df, run_trial > 2 & run_trial <= 10)))
summary(lmer(rt_swing ~ 1 + run_trial + d_auc + (1 | id/run), filter(trial_df, run_trial > 40)))

summary(m3)

summary(m4 <- lmer(rt_swing ~ 1 + run_trial + d_auc + v_entropy + (1 | id/run), trial_df))
vif.lme(m4)

summary(m5 <- lmer(rt_swing ~ 1 + run_trial + v_auc + v_chosen + v_entropy + (1 | id/run), trial_df))
vif.lme(m5)

summary(m5 <- lmer(rt_swing ~ 1 + run_trial + v_auc + v_max + v_entropy + (1 | id/run), trial_df))
vif.lme(m5)

summary(m5 <- lmer(rt_swing ~ 1 + run_trial + v_auc + v_max + v_entropy + (1 | id/run), trial_df))
vif.lme(m5)

summary(m6 <- lmer(v_entropy ~ 1 + run_trial + v_auc +  + (1 | id/run), trial_df))
vif.lme(m6)

summary(m6 <- lmer(v_entropy ~ 1 + run_trial * v_auc +  + (1 | id/run), trial_df))
vif.lme(m6)



vif.lme <- function (fit) {
  ## adapted from rms::vif
  v <- vcov(fit)
  nam <- names(fixef(fit))
  ## exclude intercepts
  ns <- sum(1 * (nam == "Intercept" | nam == "(Intercept)"))
  if (ns > 0) {
    v <- v[-(1:ns), -(1:ns), drop = FALSE]
    nam <- nam[-(1:ns)] }
  d <- diag(v)^0.5
  v <- diag(solve(v/(d %o% d)))
  names(v) <- nam
  v }

one_run <- trial_df_subset %>% filter(id==11243 & run==2)
ccf(na.omit(subset(one_run, select=c(rt_swing, d_auc))))
