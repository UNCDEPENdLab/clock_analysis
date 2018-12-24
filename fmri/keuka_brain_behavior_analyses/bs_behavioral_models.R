library(dplyr)
library(tidyverse)
library(psych)
library(corrplot)
library(lme4)
library(lmerTest)
library(ggpubr)
library(grid)
source('~/code/Rhelpers/')
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
mf2 <- lmer(log(rt_swing) ~ scale(-1/run_trial) + rewFuncIEVsum + hb_f1_DAN_vlPFC + hb_f2_neg_paralimb + (1|id/run), dfc)
screen.lmerTest(mf2)
##
mf3 <- lmer(log(rt_swing) ~ scale(-1/run_trial) + rewFuncIEVsum + hb_f1_DAN_vlPFC + hb_f2_neg_paralimb + db_f1_rIFG_rSMA + db_f4_ACC_ins + (1|id/run), dfc)
screen.lmerTest(mf3)
mf4 <- lmer(log(rt_swing) ~ scale(-1/run_trial) + rewFuncIEVsum + hb_f1_DAN_vlPFC + hb_f2_neg_paralimb + 
              vb_f3_hi_blITG + vb_f4_lo_cerebell_crus + vb_f5_lo_ACC + 
              db_f1_rIFG_rSMA + db_f2_VS + db_f3_occ_parietal + db_f4_ACC_ins + (1|id/run), dfc)
screen.lmerTest(mf4,1)

# model-based

w1 <- lmer(log(rt_swing) ~ (scale(-1/run_trial) + omission_lag + v_max_wi + v_entropy_wi + gamma)^2 + v_max_b + v_entropy_b + scale(rt_lag) + (1|id/run), dfc)
screen.lmerTest(w1)
w2h <- lmer(log(rt_swing) ~ (scale(-1/run_trial) + omission_lag + v_max_wi + v_entropy_wi + hb_f1_DAN_vlPFC + gamma)^2 + 
              (scale(-1/run_trial) + omission_lag + v_max_wi + v_entropy_wi + hb_f2_neg_paralimb + gamma)^2  + v_max_b + v_entropy_b + scale(rt_lag) + (1|id/run), dfc)
screen.lmerTest(w2h)
# reduced model without entropy -- now h1 interacts with v_max
w1h <- lmer(log(rt_swing) ~ (scale(-1/run_trial) + omission_lag + v_max_wi +  hb_f1_DAN_vlPFC + gamma)^2 + 
              (scale(-1/run_trial) + omission_lag + v_max_wi +  hb_f2_neg_paralimb + gamma)^2  + v_max_b + v_entropy_b + scale(rt_lag) + (1|id/run), dfc)
screen.lmerTest(w1h)

# the same as full
w1ha <- lmer(log(rt_swing) ~ (scale(-1/run_trial) + omission_lag + v_entropy_wi +  hb_f1_DAN_vlPFC + gamma)^2 + 
              (scale(-1/run_trial) + omission_lag + v_entropy_wi +  hb_f2_neg_paralimb + gamma)^2  + v_max_b + v_entropy_b + scale(rt_lag) + (1|id/run), dfc)
screen.lmerTest(w1ha)

# screen factors in separate models b/c collinearity
w2d1 <- lmer(log(rt_swing) ~ (scale(-1/run_trial) + omission_lag + v_max_wi + v_entropy_wi + db_f1_rIFG_rSMA + gamma)^2 + v_max_b + v_entropy_b + scale(rt_lag) + (1|id/run), dfc)
screen.lmerTest(w2d1)
w2d2 <- lmer(log(rt_swing) ~ (scale(-1/run_trial) + omission_lag + v_max_wi + v_entropy_wi + db_f2_VS + gamma)^2 + v_max_b + v_entropy_b + scale(rt_lag) + (1|id/run), dfc)
screen.lmerTest(w2d2)
w2d3 <- lmer(log(rt_swing) ~ (scale(-1/run_trial) + omission_lag + v_max_wi + v_entropy_wi + db_f3_occ_parietal + gamma)^2 + v_max_b + v_entropy_b + scale(rt_lag) + (1|id/run), dfc)
screen.lmerTest(w2d3)
w2d4 <- lmer(log(rt_swing) ~ (scale(-1/run_trial) + omission_lag + v_max_wi + v_entropy_wi + db_f4_ACC_ins + gamma)^2 + v_max_b + v_entropy_b + scale(rt_lag) + (1|id/run), dfc)
screen.lmerTest(w2d4,.1)
w2v3 <- lmer(log(rt_swing) ~ (scale(-1/run_trial) + omission_lag + v_max_wi + v_entropy_wi + vb_f3_hi_blITG + gamma)^2 + v_max_b + v_entropy_b + scale(rt_lag) + (1|id/run), dfc)
screen.lmerTest(w2v3)
w2v4 <- lmer(log(rt_swing) ~ (scale(-1/run_trial) + omission_lag + v_max_wi + v_entropy_wi + vb_f4_lo_cerebell_crus + gamma)^2 + v_max_b + v_entropy_b + scale(rt_lag) + (1|id/run), dfc)
screen.lmerTest(w2v4)
# predicts and interacts with gamma
w2v5 <- lmer(log(rt_swing) ~ (scale(-1/run_trial) + omission_lag + v_max_wi + v_entropy_wi + vb_f5_lo_ACC + gamma)^2 + v_max_b + v_entropy_b + scale(rt_lag) + (1|id/run), dfc)
screen.lmerTest(w2v5)




w3hd <- lmer(log(rt_swing) ~ (scale(-1/run_trial) + omission_lag + v_max_wi + v_entropy_wi + hb_f1_DAN_vlPFC + gamma)^2 + 
                (scale(-1/run_trial) + omission_lag + v_max_wi + v_entropy_wi + db_f1_rIFG_rSMA + gamma)^2  + 
                (scale(-1/run_trial) + omission_lag + v_max_wi + v_entropy_wi + db_f4_ACC_ins + gamma)^2  + 
                v_max_b + v_entropy_b + scale(rt_lag) + (1|id/run), dfc)
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
                 v_max_b + v_entropy_b + scale(rt_lag) + (1|id/run), dfc)
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
