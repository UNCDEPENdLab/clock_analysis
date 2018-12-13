library(dplyr)
library(tidyverse)
library(psych)
library(corrplot)
library(lme4)

load('trial_df_and_vhd_clusters.Rdata')
######
# end of preprocessing

#####
## "model-free analyses"
# RT by condition, trial and cluster activation (may need to recode rewFunc)
summary(r1 <- lmer(rt_csv ~ (scale(-1/run_trial) + rewFunc)^3 + (1|id/run), df))
summary(r2 <- lmer(rt_csv ~ (scale(-1/run_trial) + rewFunc +h_f1_fp + h_f2_neg_paralimb)^3 + (1|id/run), df))
car::Anova(r2,'3')
summary(r3 <- lmer(rt_csv ~ (scale(-1/run_trial) + rewFunc + v_f1_neg_cog + v_f2_paralimb)^3 + (1|id/run), df))
car::Anova(r3,'3')
anova(r1,r2,r3)
# these simple RT swing prediction models only reveal that both networks catalyze convergence regardless of condition
summary(m01 <- lmer(log(rt_swing) ~ scale(-1/run_trial) * rewFunc + (1|id/run), df[df$rt_swing>0,]))
summary(m02 <- lmer(log(rt_swing) ~ (scale(-1/run_trial) + rewFunc + h_f1_fp + h_f2_neg_paralimb)^2 + (1|id/run), df[df$rt_swing>0,]))
# add V to check for dissociation
summary(m03 <- lmer(log(rt_swing) ~ (scale(-1/run_trial) + rewFunc + h_f1_fp + h_f2_neg_paralimb + low_v_fp_acc_vlpfc + v_paralimbic)^2 + (1|id/run), df[df$rt_swing>0,]))
car::Anova(m03,'3')
anova(m01,m02,m03)
# over-engineered, does not add much
# summary(m03 <- lmer(rt_csv ~ (scale(rt_lag) + scale(-1/run_trial) + rewFunc + h_f1_fp + h_f2_neg_paralimb)^4 + (1|id/run), df[df$rt_swing>0,]))

#####
## "model-based analyses"
# RT swings analyses: "exploration"
summary(m1 <- lmer(log(rt_swing) ~ scale(-1/run_trial) * scale(v_max) + rewFunc + (1|id/run), df[df$rt_swing>0,]))
# keeps flipping between log and untransformed...
summary(m2 <- lmer(log(rt_swing) ~ (scale(-1/run_trial) + scale(v_max) + h_f1_fp + h_f2_neg_paralimb) ^3 + rewFunc + (1|id/run), df[df$rt_swing>0,]))
car::Anova(m2, '3')
summary(m3 <- lmer(log(rt_swing) ~ (scale(-1/run_trial) + scale(v_max) + v_f1_neg_cog + v_f2_paralimb) ^3 + rewFunc + (1|id/run), df[df$rt_swing>0,]))
car::Anova(m3, '3')
summary(m4 <- lmer(log(rt_swing) ~ (scale(-1/run_trial) + scale(v_max) + h_f1_fp + h_f2_neg_paralimb + v_f1_neg_cog + v_f2_paralimb) ^3 + rewFunc + (1|id/run), df[df$rt_swing>0,]))
anova(m4)
car::Anova(m4, '3')
# adding d introduces multi-collinearity
summary(m5 <- lmer(log(rt_swing) ~ (scale(-1/run_trial) + scale(v_max) + h_f1_fp + h_f2_neg_paralimb + v_f1_neg_cog + v_f2_paralimb + d_f1 + d_f2 + d_f3) ^3 + rewFunc + (1|id/run), df[df$rt_swing>0,]))

anova(m1,m2,m3,m4,m5)
# value clusters explain more than entropy (diffAIC = 31), but each set explains unique variance; d_auc adds little

###########
# EV analyses: "exploitation" -- run removed from RE to avoid singular fit
# Differences emerge in IEV (explains why total earnings are not informative)
summary(ev1 <- lmer(ev ~ scale(-1/run_trial) * rewFunc + (1|id), df))
summary(ev2 <- lmer(ev ~ (I(-scale(1/run_trial)) + rewFunc + h_f1_fp) ^3 + (1|id), df))
summary(ev3 <- lmer(ev ~ (I(-scale(1/run_trial)) + rewFunc + I(-h_f2_neg_paralimb)) ^3 + (1|id), df))
summary(ev4 <- lmer(ev ~ (I(-scale(1/run_trial)) + rewFunc + h_f1_fp + I(-h_f2_neg_paralimb)) ^4 + (1|id), df))
anova(ev1,ev2,ev3,ev4)

pdf('ev_by_condition_and_h_betas.pdf', height = 8, width = 8)
ggplot(df, aes(run_trial, ev, color = low_h_paralimbic, lty = h_fp)) + geom_smooth(method = "loess") + facet_wrap (~rewFunc)
dev.off()
pdf('ev_by_condition_and_v_betas.pdf', height = 8, width = 8)
ggplot(df, aes(run_trial, ev, color = v_paralimbic, lty = low_v_fp_acc_vlpfc)) + geom_smooth(method = "loess") + facet_wrap (~rewFunc)
dev.off()
#Hmmm...  this is very puzzling.  The v_max positive factor (vmPFC et al.) looks like a maladaptive response, even
# after removing cerebellum and fusiform!
# sanity check
# ggplot(df, aes(v_paralimbic,v_f2_paralimb)) + geom_boxplot()
# ggplot(df, aes(low_v_fp_acc_vlpfc,v_f1_neg_cog)) + geom_boxplot()
# check H timecourses
pdf('h_timecourse_by_condition_and_v_betas.pdf', height = 8, width = 8)
ggplot(df, aes(run_trial, v_entropy, color = v_paralimbic, lty = low_v_fp_acc_vlpfc)) + geom_smooth(method = "loess") + facet_wrap (~rewFunc)
dev.off()
pdf('h_timecourse_by_condition_and_h_betas.pdf', height = 8, width = 8)
ggplot(df, aes(run_trial, v_entropy, color = low_h_paralimbic, lty = h_fp)) + geom_smooth(method = "loess") + facet_wrap (~rewFunc)
dev.off()


pdf('rt_swings_by_condition_and_h_betas.pdf', height = 8, width = 8)
ggplot(df, aes(run_trial, log(rt_swing), color = low_h_paralimbic, lty = h_fp)) + geom_smooth(method = "loess") + facet_wrap(~rewFunc)
dev.off()

pdf('rt_swings_by_condition_and_v_betas.pdf', height = 8, width = 8)
ggplot(df, aes(run_trial, log(rt_swing), color = v_paralimbic, lty = low_v_fp_acc_vlpfc)) + geom_smooth(method = "loess") + facet_wrap(~rewFunc)
dev.off()

pdf('rt_by_condition_and_h_betas.pdf', height = 8, width = 8)
ggplot(df, aes(run_trial, rt_csv, color = low_h_paralimbic, lty = h_fp)) + geom_smooth(method = "loess") + facet_wrap(~rewFunc)
dev.off()

pdf('rt_by_condition_and_v_betas.pdf', height = 8, width = 8)
ggplot(df, aes(run_trial, rt_csv, color = v_paralimbic, lty = low_v_fp_acc_vlpfc)) + geom_smooth(method = "loess") + facet_wrap(~rewFunc)
dev.off()


# coupling between entropy and RT swings -- betas don't moderate effects of entropy.  Rather, they influence entropy as it unfolds.
# summary(m7 <- lmer(log(rt_swing) ~ (scale(-1/run_trial) + scale(v_max) + scale(v_entropy))^2 + rewFunc + (1|id/run), df[df$rt_swing>0,]))
# summary(m8 <- lmer(log(rt_swing) ~ (scale(-1/run_trial) + scale(v_max) + scale(v_entropy) + h_f1_fp) ^2 + rewFunc + (1|id/run), df[df$rt_swing>0,]))
# car::Anova(m8, '3')
# summary(m9 <- lmer(log(rt_swing) ~ (scale(-1/run_trial) + scale(v_max) + scale(v_entropy) + h_f2_neg_paralimb) ^2 + rewFunc + (1|id/run), df[df$rt_swing>0,]))
# car::Anova(m9, '3')
# summary(m10 <- lmer(log(rt_swing) ~ (scale(-1/run_trial) + scale(v_max) + scale(v_entropy) +h_f1_fp + h_f2_neg_paralimb) ^2 + rewFunc + (1|id/run), df[df$rt_swing>0,]))
# car::Anova(m10, '3')
# # a bunch of relatively weak interactions, set a side for now:
# summary(m11 <- lmer(log(rt_swing) ~ (scale(-1/run_trial) + scale(v_max) + scale(v_entropy) +h_f1_fp + h_f2_neg_paralimb + v_f1_neg_cog + v_f2_paralimb) ^3 + rewFunc + (1|id/run), df[df$rt_swing>0,]))
# car::Anova(m11, '3')
# anova(m7,m8,m9, m10)

# plot out the relationships
ggplot(df, aes(v_entropy,rt_swing, color = h_f1_fp>0)) + geom_smooth(method = "gam")
ggplot(df, aes(v_max,rt_swing,  color = low_h_paralimbic, lty = h_fp)) + geom_smooth(method = "gam")

# just entropy -- VIFs ok
summary(rt0 <- lmer(scale(rt_csv) ~ (scale(-1/run_trial) + scale(rt_vmax) + scale(rt_lag))^3 + rewFunc + (1|id/run), df))

summary(rth <- lmer(scale(rt_csv) ~ (scale(-1/run_trial) + scale(rt_vmax) + scale(rt_lag) + h_f1_fp + h_f2_neg_paralimb)^3 + rewFunc + (1|id/run), df))
car::Anova(rt1)
ggplot(df,aes(rt_vmax,rt_csv, color = low_h_paralimbic, lty = h_fp)) + facet_wrap(~run_trial>20) + geom_smooth(method = "glm")
ggplot(df,aes(rt_lag,rt_csv, color = low_h_paralimbic, lty = h_fp)) + facet_wrap(~run_trial>20) + geom_smooth(method = "glm")

# just value -- VIFs ok
summary(rtv <- lmer(scale(rt_csv) ~ (scale(-1/run_trial) + scale(rt_vmax) + scale(rt_lag) + v_f1_neg_cog + v_f2_paralimb)^3 + rewFunc + (1|id/run), df))
# + d_auc
summary(rtd <- lmer(scale(rt_csv) ~ (scale(-1/run_trial) + scale(rt_vmax) + scale(rt_lag) + d_f1 + d_f2 + d_f3)^3 + rewFunc + (1|id/run), df))
# only H and V
summary(rtvh <- lmer(scale(rt_csv) ~ (scale(-1/run_trial) + scale(rt_vmax) + scale(rt_lag) + h_f1_fp + h_f2_neg_paralimb + v_f1_neg_cog +v_f2_paralimb)^3 + rewFunc + (1|id/run), df))
# all
summary(rtvhd <- lmer(scale(rt_csv) ~ (scale(-1/run_trial) + scale(rt_vmax) + scale(rt_lag) + h_f1_fp + h_f2_neg_paralimb  + d_f1 + d_f2 + d_f3)^3 + rewFunc + (1|id/run), df))
anova(rt0,rtv,rtd,rth,rtvh,rtvhd) # H predicts best, each addition (V, D_AUC) improves the prediction further (by >300 AIC points)
# but careful with cluster*cluster interactions -- high VIFs!


# the model to rule them all?
summary(rtwh <- lmer(scale(rt_csv) ~ (scale(run_trial) + scale(rt_vmax) + scale(rt_lag) + v_entropy_wi + v_max_wi + h_f1_fp + h_f2_neg_paralimb)^3 + rewFunc + (1|id/run), df))

# plot vs. raw data
ggplot(df, aes(rt_vmax, rt_csv, color = h_f1_fp>0)) + geom_smooth(method = "glm") + facet_wrap(~run_trial>20)
ggplot(df, aes(rt_vmax, rt_csv, color = h_f2_neg_paralimb>0)) + geom_smooth(method = "glm") + facet_wrap(~run_trial>20)

ggplot(df, aes(rt_lag, rt_csv, color = h_f1_fp>0)) + geom_smooth(method = "glm") + facet_wrap(~run_trial>20)

ggplot(df, aes(rt_lag, rt_csv, color = low_h_paralimbic, lty = h_fp)) + geom_smooth(method = "glm")+ facet_wrap(~run_trial>20)
ggplot(df, aes(rt_vmax, rt_csv, color = low_h_paralimbic, lty = h_fp)) + geom_smooth(method = "glm")+ facet_wrap(~run_trial>20)

pdf("h_timecourse_brain_fixed.pdf", width = 8, height = 8)
ggplot(df, aes(run_trial, v_entropy,color = low_h_paralimbic, lty = h_fp)) + geom_smooth(method = "loess") + facet_wrap(~gamma>0)
dev.off()


ggplot(df, aes(run_trial, v_entropy, color = rewFunc)) + geom_smooth()
summary(m0 <- lmer(log(rt_swing) ~ scale(-1/run_trial) * scale(v_max) * scale(v_entropy) * rewFunc + (1|id/run), df[df$rt_swing>0,]))

# plot out

# toward the goal of dissociating exploitative behavioral adjustment from H-driven exploration
# dm = dissociation model
# unfortunately does not work all that well due to collinearity
summary(dm1 <- lmer(log(rt_swing) ~ (scale(-1/run_trial) + omission_lag + scale(v_max))^3 + scale(rt_lag) + (1|id/run), df[df$rt_swing>0,]))
ggplot(na.omit(df), aes(v_entropy, rt_swing, color = v_max>25, lty = omission_lag)) + geom_smooth(method = "gam")

summary(dm2 <- lmer(rt_csv ~ (scale(rt_lag) + omission_lag + scale(rt_vmax))^3 + rewFunc + (1|id/run), df[df$rt_swing>0,]))

summary(dm3 <- lmer(rt_csv ~ (scale(rt_lag) + omission_lag + scale(rt_vmax) + h_f1_fp + h_f2_neg_paralimb)^3 + rewFunc + (1|id/run), df[df$rt_swing>0,]))

# more focused test of dislodgment
summary(dm3a <- lmer(rt_csv ~ (scale(rt_lag) + omission_lag + scale(rt_vmax) + h_f1_fp)^3 + rewFunc + (1|id/run), df[df$rt_swing>0,]))
summary(dm3b <- lmer(rt_csv ~ (scale(rt_lag) + omission_lag + scale(rt_vmax) + h_f2_neg_paralimb)^3 + rewFunc + (1|id/run), df[df$rt_swing>0,]))

# plot 3-way interactions
ggplot(df, aes(rt_lag, rt_csv, color = omission_lag, lty = h_f1_fp>0)) + geom_smooth(method = "glm")
ggplot(df, aes(rt_lag, rt_csv, color = omission_lag, lty = h_f2_neg_paralimb<0)) + geom_smooth(method = "glm")
# the interesting one:
ggplot(na.omit(df), aes(rt_vmax, rt_csv, color = omission_lag, lty = h_f1_fp>0)) + geom_smooth(method = "glm")
ggplot(na.omit(df), aes(rt_vmax, rt_csv, color = omission_lag, lty = h_f2_neg_paralimb<0)) + geom_smooth(method = "glm")

anova(dm1,dm2,dm3)

summary(dm2 <- lmer(log(rt_swing) ~ (scale(-1/run_trial) +  scale(v_entropy) + omission_lag + scale(v_max) + h_f1_fp + h_f2_neg_paralimb)^3 + (1|id/run), df[df$rt_swing>0,]))
anova(dm1,dm2)

summary(w1 <- lmer(log(rt_swing) ~ (scale(run_trial) + v_max_wi + v_entropy_wi)^2 + rewFunc + (1|id/run), df[df$rt_swing>0,]))
summary(w2 <- lmer(log(rt_swing) ~ (scale(run_trial) + v_max_wi + v_entropy_wi + h_f1 + h_f2)^3 + rewFunc + (1|id/run), df[df$rt_swing>0,]))
car::Anova(w2)

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
