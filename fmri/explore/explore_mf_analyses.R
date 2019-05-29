# analyzes behavioral preliminary data from clock reversal task in explore
# first run run_fsl_pipeline_explore.R
library(tidyverse)
library(lme4)
library(car)
library(emmeans)
df <- vba_output
subject_df$id <- as.character(subject_df$redcapid)
df <- inner_join(df, subject_df)
df <- df %>% group_by(id) %>% mutate(rev_trial = case_when(
  trial < 41 ~ trial,
  trial > 40 & trial < 81 ~ trial - 40,
  trial > 80 & trial < 121 ~ trial - 80,
  trial > 120 & trial < 161 ~ trial - 120,
  trial > 160 & trial < 201 ~ trial - 160,
  trial > 200 ~ trial - 200
  ),
  rt_lag  = lag(rt_csv),
  rt_lag_scale  = scale(lag(rt_csv)),
  rev_trial_neginv_scale = scale( - 100/rev_trial),
  trial_neginv_sc = scale( - 100/trial),
  rew_lag = lag(score_csv>0),
  pe_max_lag = lag(pe_max)) %>% ungroup() %>% mutate(rt_sc = scale(rt_csv))
ggplot(df, aes(rev_trial, rt_csv, color = Group, lty = rewFunc)) + geom_smooth()
ggplot(df, aes(rev_trial, rt_csv, color = female, lty = rewFunc)) + geom_smooth()
ggplot(df, aes(trial, rt_csv, color = Group, lty = rewFunc)) + geom_smooth(method = 'gam')

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
summary(m1 <- lmer(rt_sc ~ rt_lag_scale + rewFunc *  rev_trial_neginv_scale + (1|id), df))
Anova(m1, '3')
summary(m2 <- lmer(rt_sc ~ rt_lag_scale + rewFunc *  rev_trial_neginv_scale * Group + (1|id), df))
Anova(m2, '3')
summary(m3 <- lmer(rt_sc ~ rt_lag_scale + rewFunc *  rev_trial_neginv_scale * Group + rewFunc * rev_trial_neginv_scale * scale(age) + (1|id), df))
Anova(m3, '3')
summary(m4 <- lmer(rt_sc ~ rt_lag_scale + rewFunc *  rev_trial_neginv_scale * Group + rewFunc * rev_trial_neginv_scale * scale(age) + rewFunc * rev_trial_neginv_scale * female + (1|id), df))
Anova(m4, '3')
anova(m1,m2,m3,m4)

summary(m1a <- lmer(rt_sc ~ rt_lag_scale + rewFunc *  rev_trial_neginv_scale * trial_neginv_sc + (1|id), df))
Anova(m1a,'3')
anova(m1, m1a)
summary(m4a <- lmer(rt_sc ~ rt_lag_scale + rewFunc *  rev_trial_neginv_scale * trial_neginv_sc * Group + rewFunc * rev_trial_neginv_scale * scale(age) + rewFunc * rev_trial_neginv_scale * female + (1|id), df))
summary(m4a)
Anova(m4a, '3')
anova(m1,m2,m3,m4,m1a, m4a)

em4 <- emmeans(m4a, Group ~ )

summary(mt4a <- lmer(rt_sc ~ rt_lag_scale + rewFunc *  rev_trial_neginv_scale * trial_neginv_sc * Group + rewFunc * rev_trial_neginv_scale * scale(age) + rewFunc * rev_trial_neginv_scale * female + 
                       rew_lag * Group + rew_lag * scale(age) + rew_lag * female + (1|id), df)))
Anova(mt4a)

summary(mt5a <- lmer(rt_sc ~ rt_lag_scale + (rewFunc +  rev_trial_neginv_scale + trial_neginv_sc + Group) ^3 + 
                       (rewFunc +  rev_trial_neginv_scale + trial_neginv_sc + scale(age)) ^3 + 
                     (rewFunc +  rev_trial_neginv_scale + trial_neginv_sc + female) ^3 +
                       rew_lag * rewFunc * Group + rew_lag * rewFunc * scale(age) + rew_lag * rewFunc * female + (1|id), df))
Anova(mt5a, '3')

summary(ms1 <- lmer(rt_sc ~ (rt_lag_scale + rev_trial_neginv_scale + trial_neginv_sc + rt_lag_scale + rewFunc)^2 + (1|id), df))
summary(ms2 <- lmer(rt_sc ~ (rt_lag_scale + rev_trial_neginv_scale + trial_neginv_sc + rt_lag_scale + rewFunc)^3 + (1|id), df))
anova(ms1,ms2)

ggplot(df, aes(trial, abs(rt_sc-rt_lag_scale), lty = v_entropy>0, color = Group)) + geom_smooth(method = 'gam')
ggplot(df, aes(v_entropy, abs(rt_sc-rt_lag_scale), color = Group)) + geom_smooth(method = 'gam') + facet_wrap(~rewFunc)


summary(mb1 <- lmer(rt_sc ~ (rt_lag_scale + v_entropy + rev_trial_neginv_scale)^2 + (1|id), df))
summary(mb2 <- lmer(rt_sc ~ (rt_lag_scale + v_entropy + rev_trial_neginv_scale + rewFunc)^2 + (1|id), df))
summary(mb3 <- lmer(rt_sc ~ (rt_lag_scale + v_entropy + rev_trial_neginv_scale + rewFunc + Group)^2 + (1|id), df))
vif.lme(mb3)
vif(mb3)
summary(mb4 <- lmer(rt_sc ~ (rt_lag_scale + v_entropy + rev_trial_neginv_scale + rewFunc + Group)^3 + (1|id), df))
Anova(mb4, '3')

