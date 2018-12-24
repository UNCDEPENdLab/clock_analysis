library(dplyr)
library(tidyverse)
library(psych)
library(corrplot)
library(lme4)
library(lmerTest)
library(ggpubr)
library(grid)
# source('~/code/Rhelpers/')
setwd('~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/')
load('trial_df_and_vhd_bs.Rdata')
######
# end of preprocessing

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

# check VIFs of significant effects
screen.lmerTest <- function (mod,p=NULL) {
  if (is.null(p)) {p <- .05}
  c1 <- as.data.frame(coef(summary(mod))[,4:5])
  dd <- cbind(c1[2:nrow(c1),],as.data.frame(vif.lme(mod)))
  names(dd)[3] <- 'VIF'
  dd$`Pr(>|t|)` <- as.numeric(dd$`Pr(>|t|)`)
  print(dd[dd$`Pr(>|t|)`<p,c(1,3)], digits = 3)}
#####

# for further analyses, we want the following regions:
# hb_f1_DAN (alt. vb_f1_lo_DAN_vlPFC)
# hb_f2_paralimbic (alt vb_f2_hi_paralimbic)
# db_f4_ACC_ins -- compression
# db_f1_rIFG_SMA -- ?overload

dfc <- na.omit(df[df$rt_swing>0,])

# parallel to the beta models

## "model-free analyses"
mf1 <- lmer(log(rt_swing) ~ scale(-1/run_trial) +  rewFuncIEVsum + (1|id/run), dfc)
screen.lmerTest(mf1)
## favorite simple model ##
mf2 <- lmer(log(rt_swing) ~ (scale(-1/run_trial) + rewFuncIEVsum + hb_f1_DAN_vlPFC + hb_f2_neg_paralimb + gamma)^3 + 
              rt_lag + (1|id/run), dfc)
screen.lmerTest(mf2)
##
mf3 <- lmer(log(rt_swing) ~ scale(-1/run_trial) + rewFuncIEVsum + hb_f1_DAN_vlPFC + hb_f2_neg_paralimb + db_f1_rIFG_rSMA + db_f4_ACC_ins + (1|id/run), dfc)
screen.lmerTest(mf3)
mf4 <- lmer(log(rt_swing) ~ scale(-1/run_trial) + rewFuncIEVsum + hb_f1_DAN_vlPFC + hb_f2_neg_paralimb + 
              vb_f3_hi_blITG + vb_f4_lo_cerebell_crus + vb_f5_lo_ACC + 
              db_f1_rIFG_rSMA + db_f2_VS + db_f3_occ_parietal + db_f4_ACC_ins + (1|id/run), dfc)
screen.lmerTest(mf4)
# model-free RT
mfr1 <- lmer(rt_csv ~ (scale(-1/run_trial) +  rewFuncIEVsum)^2 + (1|id/run), dfc)
screen.lmerTest(mfr1)
## favorite simple model ##
mfr2 <- lmer(rt_csv ~ (scale(-1/run_trial) + rewFuncIEVsum + hb_f1_DAN_vlPFC + hb_f2_neg_paralimb + gamma)^3 + (1|id/run), dfc)
screen.lmerTest(mfr2)
##
mfr3 <- lmer(rt_csv ~ (scale(-1/run_trial) + rewFuncIEVsum + hb_f1_DAN_vlPFC + hb_f2_neg_paralimb + db_f1_rIFG_rSMA + db_f4_ACC_ins)^3 + (1|id/run), dfc)
screen.lmerTest(mfr3)
# just paralimbic and ACC
mfr3a <- lmer(rt_csv ~ (scale(-1/run_trial) + rewFuncIEVsum + hb_f2_neg_paralimb + db_f4_ACC_ins)^4 + (1|id/run), dfc)
screen.lmerTest(mfr3a)
# 
ggplot(dfc, aes(hb_f2_neg_paralimb,rt_csv, color = db_f4_ACC_ins>0)) + geom_smooth(method = "loess")
ggplot(dfc, aes(hb_f1_DAN_vlPFC,log(rt_swing), color = gamma>0)) + geom_smooth(method = "loess")
ggplot(dfc, aes(hb_f2_neg_paralimb,log(rt_swing),color = gamma>0)) + geom_smooth(method = "loess")

mfr4 <- lmer(rt_csv ~ scale(-1/run_trial) + rewFuncIEVsum + hb_f1_DAN_vlPFC + hb_f2_neg_paralimb + 
              vb_f3_hi_blITG + vb_f4_lo_cerebell_crus + vb_f5_lo_ACC + 
              db_f1_rIFG_rSMA + db_f2_VS + db_f3_occ_parietal + db_f4_ACC_ins + (1|id/run), dfc)
screen.lmerTest(mf4)


# model-based

# first, does DAN improve prediction of RTs and EXPLAIN effects of v_entropy_wi?
# it improves prediction, but seems to contain information independent of entropy_wi...
w0 <- lmer(log(rt_swing) ~ (scale(-1/run_trial) + omission_lag + v_max_wi + v_entropy_wi + gamma)^2 + 
             v_max_b + v_entropy_b + scale(rt_swing_lag) + scale(rt_lag) + (1|id/run), dfc)
screen.lmerTest(w0)

w0dan <- lmer(log(rt_swing) ~ (scale(-1/run_trial) + omission_lag + v_max_wi + v_entropy_wi + hb_f1_DAN_vlPFC + gamma)^2 + 
             v_max_b + v_entropy_b + scale(rt_swing_lag) + scale(rt_lag) + (1|id/run), dfc)
screen.lmerTest(w0dan)
anova(w0, w0dan)


w1 <- lmer(log(rt_swing) ~ (scale(-1/run_trial) + omission_lag + v_max_wi + v_entropy_wi + gamma)^2 + v_max_b + v_entropy_b +scale(rt_swing_lag) + scale(rt_lag) + (1|id/run), dfc)
screen.lmerTest(w1)
w2h <- lmer(log(rt_swing) ~ (scale(-1/run_trial) + omission_lag + v_max_wi + v_entropy_wi + hb_f1_DAN_vlPFC + gamma)^3 + 
              (scale(-1/run_trial) + omission_lag + v_max_wi + v_entropy_wi + hb_f2_neg_paralimb + gamma)^2  + v_max_b + v_entropy_b + scale(rt_swing_lag) + scale(rt_lag) + (1|id/run), dfc)
screen.lmerTest(w2h)
# reduced model without entropy -- now h1 interacts with v_max, but just barely
w1h <- lmer(log(rt_swing) ~ (scale(-1/run_trial) + omission_lag + v_max_wi +  hb_f1_DAN_vlPFC + gamma)^2 + 
              (scale(-1/run_trial) + omission_lag + v_max_wi +  hb_f2_neg_paralimb + gamma)^2  + v_max_b + v_entropy_b + scale(rt_swing_lag) + scale(rt_lag) + (1|id/run), dfc)
screen.lmerTest(w1h)

# the same as full
w1ha <- lmer(log(rt_swing) ~ (scale(-1/run_trial) + omission_lag + v_entropy_wi +  hb_f1_DAN_vlPFC + gamma)^2 + 
              (scale(-1/run_trial) + omission_lag + v_entropy_wi +  hb_f2_neg_paralimb + gamma)^2  + v_max_b + v_entropy_b + scale(rt_swing_lag) + scale(rt_lag) + (1|id/run), dfc)
screen.lmerTest(w1ha)

# screen factors in separate models b/c collinearity
w2d1 <- lmer(log(rt_swing) ~ (scale(-1/run_trial) + omission_lag + v_max_wi + v_entropy_wi + db_f1_rIFG_rSMA + gamma)^2 +scale(rt_swing_lag) +scale(rt_swing_lag) +  v_max_b + v_entropy_b + scale(rt_lag) + (1|id/run), dfc)
screen.lmerTest(w2d1)
# VS seems to accentuate the effect of v_max
w2d2 <- lmer(log(rt_swing) ~ (scale(-1/run_trial) + omission_lag + v_max_wi + v_entropy_wi + db_f2_VS + gamma)^2 + v_max_b + v_entropy_b + scale(rt_swing_lag) + scale(rt_lag) + (1|id/run), dfc)
screen.lmerTest(w2d2)
# effect of d3
w2d3 <- lmer(log(rt_swing) ~ (scale(-1/run_trial) + omission_lag + v_max_wi + v_entropy_wi + db_f3_occ_parietal + gamma)^2 + v_max_b + v_entropy_b + scale(rt_swing_lag) + scale(rt_lag) + (1|id/run), dfc)
screen.lmerTest(w2d3)
# effect of d4 relatively subtle
w2d4 <- lmer(log(rt_swing) ~ (scale(-1/run_trial) + omission_lag + v_max_wi + v_entropy_wi + db_f4_ACC_ins + gamma)^2 + v_max_b + v_entropy_b + scale(rt_swing_lag) + scale(rt_lag) + (1|id/run), dfc)
screen.lmerTest(w2d4,.1)
# no effect of ITG/v3
w2v3 <- lmer(log(rt_swing) ~ (scale(-1/run_trial) + omission_lag + v_max_wi + v_entropy_wi + vb_f3_hi_blITG + gamma)^2 + v_max_b + v_entropy_b + scale(rt_swing_lag) + scale(rt_lag) + (1|id/run), dfc)
screen.lmerTest(w2v3)
# no effect of cerebellar v4
w2v4 <- lmer(log(rt_swing) ~ (scale(-1/run_trial) + omission_lag + v_max_wi + v_entropy_wi + vb_f4_lo_cerebell_crus + gamma)^2 + v_max_b + v_entropy_b + scale(rt_swing_lag) + scale(rt_lag) + (1|id/run), dfc)
screen.lmerTest(w2v4)
# the more properly anterior v5 ACC predicts RT and interacts with gamma
w2v5 <- lmer(log(rt_swing) ~ (scale(-1/run_trial) + omission_lag + v_max_wi + v_entropy_wi + vb_f5_lo_ACC + gamma)^2 + v_max_b + v_entropy_b + scale(rt_swing_lag) + scale(rt_lag) + (1|id/run), dfc)
screen.lmerTest(w2v5)

# conclusions: focus on d1, d4 (despite moderate effect), v5


w3hd <- lmer(log(rt_swing) ~ (scale(-1/run_trial) + omission_lag + v_max_wi + v_entropy_wi + hb_f1_DAN_vlPFC + gamma)^2 + 
                (scale(-1/run_trial) + omission_lag + v_max_wi + v_entropy_wi + db_f1_rIFG_rSMA + gamma)^2  + 
                (scale(-1/run_trial) + omission_lag + v_max_wi + v_entropy_wi + db_f4_ACC_ins + gamma)^2  + 
                v_max_b + v_entropy_b + scale(rt_swing_lag) + scale(rt_lag) + (1|id/run), dfc)
screen.lmerTest(w3hd)
w4vhd <- lmer(log(rt_swing) ~ (scale(-1/run_trial) + omission_lag + v_max_wi + v_entropy_wi + hb_f1_DAN_vlPFC + gamma)^2 + 
                (scale(-1/run_trial) + omission_lag + v_max_wi + v_entropy_wi + hb_f2_neg_paralimb + gamma)^2  + 
                (scale(-1/run_trial) + omission_lag + v_max_wi + v_entropy_wi + vb_f3_hi_blITG + gamma)^2  + 
                (scale(-1/run_trial) + omission_lag + v_max_wi + v_entropy_wi + vb_f4_lo_cerebell_crus + gamma)^2  + 
                (scale(-1/run_trial) + omission_lag + v_max_wi + v_entropy_wi + vb_f5_lo_ACC + gamma)^2  + 
                (scale(-1/run_trial) + omission_lag + v_max_wi + v_entropy_wi + db_f1_rIFG_rSMA + gamma)^2  + 
                (scale(-1/run_trial) + omission_lag + v_max_wi + v_entropy_wi + db_f2_VS + gamma)^2  + 
                (scale(-1/run_trial) + omission_lag + v_max_wi + v_entropy_wi + db_f3_occ_parietal + gamma)^2  + 
                (scale(-1/run_trial) + omission_lag + v_max_wi + v_entropy_wi + db_f4_ACC_ins + gamma)^2  + 
                 v_max_b + v_entropy_b + scale(rt_swing_lag) + scale(rt_lag) + (1|id/run), dfc)
screen.lmerTest(w4vhd)
anova(w1,w2h,w3hd,w4vhd)
# these other regions explain a lot of variance, but the multi-collinearity is tough

# RTs
wr1 <- lmer(rt_csv ~ (scale(-1/run_trial) + scale(rt_lag) + scale(rt_vmax_lag) + omission_lag + v_max_wi_lag + v_entropy_wi + gamma)^2 + v_max_b + v_entropy_b + scale(rt_lag) + (1|id/run), dfc)
screen.lmerTest(wr1)
wr2h <- lmer(rt_csv ~ (scale(-1/run_trial) + scale(rt_lag) + scale(rt_vmax_lag) + omission_lag + v_max_wi_lag + v_entropy_wi + hb_f1_DAN_vlPFC + gamma)^3 
             +(scale(-1/run_trial) + scale(rt_lag) + scale(rt_vmax_lag) + omission_lag + v_max_wi_lag + v_entropy_wi + hb_f2_neg_paralimb + gamma)^3 
             + v_max_b + v_entropy_b +  (1|id/run), dfc)
screen.lmerTest(wr2h)
# reduced
wr2h <- lmer(rt_csv ~ (scale(-1/run_trial) + scale(rt_lag) + scale(rt_vmax_lag) + omission_lag + v_max_wi_lag + hb_f1_DAN_vlPFC + gamma)^3 
             +(scale(-1/run_trial) + scale(rt_lag) + scale(rt_vmax_lag) + omission_lag + v_max_wi_lag +  hb_f2_neg_paralimb + gamma)^3 
             + v_max_b + v_entropy_b +  (1|id/run), dfc)


wr3hd <- lmer(rt_csv ~ (scale(-1/run_trial) + scale(rt_lag) + scale(rt_vmax_lag) + omission_lag + v_max_wi_lag + v_entropy_wi + hb_f1_DAN_vlPFC + gamma)^3 
             +(scale(-1/run_trial) + scale(rt_lag) + scale(rt_vmax_lag) + omission_lag + v_max_wi_lag + v_entropy_wi + hb_f2_neg_paralimb + gamma)^3 
             +(scale(-1/run_trial) + scale(rt_lag) + scale(rt_vmax_lag) + omission_lag + v_max_wi_lag + v_entropy_wi + db_f1_rIFG_rSMA + gamma)^3 
             +(scale(-1/run_trial) + scale(rt_lag) + scale(rt_vmax_lag) + omission_lag + v_max_wi_lag + v_entropy_wi + db_f4_ACC_ins + gamma)^3 
             + v_max_b + v_entropy_b +  (1|id/run), dfc)
wr4hvd <- lmer(rt_csv ~ (scale(-1/run_trial) + scale(rt_lag) + scale(rt_vmax_lag) + omission_lag + v_max_wi_lag + v_entropy_wi + hb_f1_DAN_vlPFC + gamma)^3 
              +(scale(-1/run_trial) + scale(rt_lag) + scale(rt_vmax_lag) + omission_lag + v_max_wi_lag + v_entropy_wi + hb_f2_neg_paralimb + gamma)^3 
              +(scale(-1/run_trial) + scale(rt_lag) + scale(rt_vmax_lag) + omission_lag + v_max_wi_lag + v_entropy_wi + vb_f3_hi_blITG + gamma)^3 
              +(scale(-1/run_trial) + scale(rt_lag) + scale(rt_vmax_lag) + omission_lag + v_max_wi_lag + v_entropy_wi + vb_f4_lo_cerebell_crus + gamma)^3 
              +(scale(-1/run_trial) + scale(rt_lag) + scale(rt_vmax_lag) + omission_lag + v_max_wi_lag + v_entropy_wi + vb_f5_lo_ACC + gamma)^3 
              +(scale(-1/run_trial) + scale(rt_lag) + scale(rt_vmax_lag) + omission_lag + v_max_wi_lag + v_entropy_wi + db_f1_rIFG_rSMA + gamma)^3 
              +(scale(-1/run_trial) + scale(rt_lag) + scale(rt_vmax_lag) + omission_lag + v_max_wi_lag + v_entropy_wi + db_f2_VS + gamma)^3 
              +(scale(-1/run_trial) + scale(rt_lag) + scale(rt_vmax_lag) + omission_lag + v_max_wi_lag + v_entropy_wi + db_f3_occ_parietal + gamma)^3 
              +(scale(-1/run_trial) + scale(rt_lag) + scale(rt_vmax_lag) + omission_lag + v_max_wi_lag + v_entropy_wi + db_f4_ACC_ins + gamma)^3 
              + v_max_b + v_entropy_b +  (1|id/run), dfc)

screen.lmerTest(wr4hvd)
anova(wr1,wr2h,wr3hd,wr4hvd)

#####
## prediction of BS
# lagged relationship between ACC and DAN
b1 <- lmer(hb_f1_DAN_vlPFC ~ (scale(-1/run_trial) + omission_lag + v_max_wi + v_entropy_wi + h1_lag + gamma)^2 
              + v_max_b + v_entropy_b +  (1|id) + (1|run), df)
screen.lmerTest(b1)
# shockingly, it is present in high-gamma subjects
b2 <- lmer(hb_f1_DAN_vlPFC ~ (scale(-1/run_trial) + omission_lag + v_max_wi + v_entropy_wi + gamma + h1_lag + d4_lag)^2 
           + v_max_b + v_entropy_b + (1|run), df)
screen.lmerTest(b2)
# does this happen with any region?
b3 <- lmer(hb_f1_DAN_vlPFC ~ (scale(-1/run_trial) + omission_lag + v_max_wi + v_entropy_wi + gamma + h1_lag + h2_lag)^2 
           + v_max_b + v_entropy_b + (1|run), df)
screen.lmerTest(b3)
# it does for h2_paralimbic, but not for VS
b4 <- lmer(hb_f1_DAN_vlPFC ~ (scale(-1/run_trial) + omission_lag + v_max_wi + v_entropy_wi + gamma + h1_lag + db_f1_rIFG_rSMA + d1_lag)^2 
           + v_max_b + v_entropy_b + (1|run), df)
screen.lmerTest(b4)
b5 <- lmer(hb_f1_DAN_vlPFC ~ (scale(-1/run_trial) + omission_lag + v_max_wi + v_entropy_wi + gamma + h1_lag + d2_lag)^2 
           + v_max_b + v_entropy_b + (1|run), df)
screen.lmerTest(b5)
# is it really a lagged relationship?
# yes, although VIF = 2.73
b5 <- lmer(hb_f1_DAN_vlPFC ~ (scale(-1/run_trial) + omission_lag + v_max_wi + v_entropy_wi + gamma + h1_lag + db_f4_ACC_ins + d4_lag)^2 
           + v_max_b + v_entropy_b + (1|run), df)
screen.lmerTest(b5)
# seems to hold with h2, although VIFs are atrocious
b6 <- lmer(hb_f1_DAN_vlPFC ~ (scale(-1/run_trial) + omission_lag + v_max_wi + v_entropy_wi + gamma + 
                                h1_lag + db_f4_ACC_ins + d4_lag + hb_f2_neg_paralimb + h2_lag)^2 
           + v_max_b + v_entropy_b + (1|run), df)
screen.lmerTest(b6)

#####
## value or entropy in DAN?
## ENTROPY even in V-negative region!
v1 <- lmer(vb_f1_lo_DAN ~ (scale(-1/run_trial) + omission_lag + v_max_wi + v_entropy_wi + v1_lag + gamma)^2 
           + v_max_b + v_entropy_b +  (1|id) + (1|run), dfc)
screen.lmerTest(v1)
v2 <- lmer(vb_f2_hi_vmPFC_cOFC ~ (scale(-1/run_trial) + omission_lag + v_max_wi + v_entropy_wi + v2_lag + gamma)^2 
           + v_max_b + v_entropy_b +  (1|id) + (1|run), dfc)
screen.lmerTest(v2)


