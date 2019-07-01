# analyzes behavioral preliminary data from clock reversal task in explore
# first run run_fsl_pipeline_explore.R
library(tidyverse)
library(lme4)
library(car)
library(emmeans)
library(ggpubr)
setwd('~/OneDrive/grants/explore_renewal_drafts/A1/data')

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
  block = case_when(
    trial < 41 ~ 1,
    trial > 40 & trial < 81 ~ 2,
    trial > 80 & trial < 121 ~ 3,
    trial > 120 & trial < 161 ~ 4,
    trial > 160 & trial < 201 ~ 5,
    trial > 200 ~ 6
  ),
  group_full = case_when(
    Group == 'HC' ~ 'Controls',
    Group == 'DEP' ~ 'Depressed',
    Group == 'IDE' ~ 'Ideators',
    Group == 'ATT' ~ 'Attempters'),
  rt_lag  = lag(rt_csv),
  rt_lag_scale  = scale(lag(rt_csv)),
  rev_trial_neginv_scale = scale( - 100/rev_trial),
  trial_neginv_sc = scale( - 100/trial),
  rew_lag = lag(score_csv>0),
  rt_swing = abs(rt_csv - rt_lag),
  pe_max_lag = lag(pe_max)) %>% ungroup() %>% mutate(rt_sc = scale(rt_csv)) %>% filter(rt_csv<5000)
save(list = c('df', 'subject_df'), file = 'explore_clock_05_2019.RData')
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
pdf('rt_swings_by_entropy_and_group.pdf', width = 4, height = 3)
ggplot(df, aes(v_entropy, abs(rt_sc-rt_lag_scale), color = group_full)) + geom_smooth(method = 'gam') + 
  xlim(-3,3) + xlab('Entropy, scaled') + ylab('Response time swing') + theme_bw() + theme(legend.title=element_blank()) +
  annotate("text", x = 0, y = .25, label = "Entropy * Group\n Chisq = 12.00, p = .012") + 
  annotate("text", x = -1, y = 1.1, label = "Group\n Chisq = 9.65, p = .022") 

dev.off()
ggplot(df, aes(rev_trial, abs(rt_sc-rt_lag_scale), lty = v_entropy>0, color = Group)) + geom_smooth(method = 'gam') + facet_wrap(~rewFunc)



summary(mb1 <- lmer(rt_sc ~ (rt_lag_scale + v_entropy + rev_trial_neginv_scale)^2 + (1|id), df))
summary(mb2 <- lmer(rt_sc ~ (rt_lag_scale + v_entropy + rev_trial_neginv_scale + rewFunc)^2 + (1|id), df))
summary(mb3 <- lmer(rt_sc ~ (rt_lag_scale + v_entropy + rev_trial_neginv_scale + rewFunc + Group)^2 + (1|id), df))
vif.lme(mb3)
vif(mb3)
summary(mb4 <- lmer(rt_sc ~ (rt_lag_scale + v_entropy + rev_trial_neginv_scale + rewFunc + Group)^3 + (1|id), df))
Anova(mb4, '3')

summary(mb5 <- lmer(ev ~ v_entropy * Group * rewFunc + (1|id),df))

summary(mbs1 <- lmer(rt_swing ~ v_entropy * Group + (1|id),df))
Anova(mbs1)

# understand entropy timecourses as a function of rewFunc
pdf('h_by_rewFunc.pdf', width = 6, height = 8)
ggplot(df, aes(rev_trial, v_entropy, color = rewFunc)) + geom_smooth()
dev.off()

# just for Dr. Chen: by id
pdf('h_by_rewFunc_by_id.pdf', width = 30, height = 30)
ggplot(df, aes(rev_trial, v_entropy, color = rewFunc)) + geom_smooth() + facet_wrap(~id)
dev.off()

# h is confounded with RT
pdf('h_by_rt.pdf', width = 6, height = 8)
ggplot(df, aes(rt_csv, v_entropy, color = rewFunc)) + geom_smooth()
dev.off()
pdf('rt_swing_by_rewFunc.pdf', width = 6, height = 6)
ggplot(df, aes(rev_trial, rt_swing, color = rewFunc)) + geom_smooth()
dev.off()
pdf('rt_by_rewFunc.pdf', width = 6, height = 6)
ggplot(df, aes(rev_trial, rt_csv, color = rewFunc)) + geom_smooth()
dev.off()
pdf('rt_by_rewFunc_by_id.pdf', width = 30, height = 30)
ggplot(df, aes(rev_trial, rt_csv, color = rewFunc)) + geom_smooth() + facet_wrap(~id)
dev.off()

g1 <- ggplot(df, aes(rt_csv,probability, color = rewFunc)) + geom_smooth()
g2 <- ggplot(df, aes(rt_csv,magnitude, color = rewFunc)) + geom_smooth()
g3 <- ggplot(df, aes(rt_csv,ev, color = rewFunc)) + geom_smooth()
pdf('contingencies_clock_reversal.pdf', width = 6, height = 6)
ggarrange(g1,g2,g3)
dev.off()

g1 <- ggplot(df, aes(rt_csv,probability, color = rewFunc)) + geom_point() + facet_wrap(~id)
g2 <- ggplot(df, aes(rt_csv,magnitude, color = rewFunc)) + geom_point() + facet_wrap(~id)
g3 <- ggplot(df, aes(rt_csv,ev, color = rewFunc)) + geom_point() + facet_wrap(~id)
pdf('contingencies_clock_reversal_by_id.pdf', width = 30, height = 30)
ggarrange(g1,g2,g3)
dev.off()

# exclude very short responses
ldf <- df %>% filter(rt_csv>500) %>% filter(!is.na(v_entropy))
# rt swings by condition and group
pdf('rt_swing_by_rewFunc_group.pdf', width = 6, height = 6)
ggplot(ldf, aes(rev_trial, rt_swing, lty = rewFunc, color = Group)) + geom_smooth(method = 'gam')
dev.off()

pdf('rt_swing_by_rewFunc_h.pdf', width = 6, height = 6)
ggplot(ldf, aes(rev_trial, rt_swing, lty = rewFunc, color = v_entropy>0)) + geom_smooth(method = 'gam')
dev.off()

