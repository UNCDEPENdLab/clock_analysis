library(dplyr)
library(tidyverse)
library(psych)
library(corrplot)
library(lme4)
library(lmerTest)
# source('~/code/Rhelpers/')
setwd('~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/')
load('trial_df_and_vhdkfpe_clusters.Rdata')
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
## "model-free analyses"


# deal with multicollinearity in RT models
vv <- df[,c("rt_lag", "rt_vmax", "rt_vmax_lag", "run_trial", "omission_lag", "v_max_wi", "v_max_wi_lag", "v_entropy_wi")]
corr.test(vv)
# NB: omissions have two effects: RT shortening and non-directional RT swing
# collinearity remediated by taking lags of rt_vmax and v_max, conceptually the same variables, just discard the effect of last outcome
wr1 <- lmer(rt_csv ~ (scale(-1/run_trial) + scale(rt_lag) + scale(rt_vmax_lag) + omission_lag + v_max_wi_lag + v_entropy_wi)^2 + v_max_b + v_entropy_b + scale(rt_lag) + (1|id/run), df)
screen.lmerTest(wr1)
wr2h <- lmer(rt_csv ~ (scale(-1/run_trial) + scale(rt_lag) + scale(rt_vmax_lag) + omission_lag + v_max_wi_lag + v_entropy_wi + h_f1_fp + I(-h_f2_neg_paralimb))^2 + v_max_b + v_entropy_b +  (1|id/run), df)
screen.lmerTest(wr2h)
wr3hpe <- lmer(rt_csv ~ (scale(-1/run_trial) + scale(rt_lag) + scale(rt_vmax_lag) + omission_lag + v_max_wi_lag + v_entropy_wi + h_f1_fp + I(-h_f2_neg_paralimb) + pe_f1_cort_str + pe_f2_hipp)^2 + v_max_b + v_entropy_b +   (1|id/run), df)
screen.lmerTest(wr3hpe)
# control for performance
wr3hpe_perf <- lmer(rt_csv ~ (scale(-1/run_trial) + scale(rt_lag) + scale(rt_vmax_lag) + omission_lag + v_max_wi_lag + v_entropy_wi + h_f1_fp + I(-h_f2_neg_paralimb) + pe_f1_cort_str + pe_f2_hipp + scale(total_earnings))^2 + 
                      scale(rt_lag)*pe_f2_hipp*scale(total_earnings) + v_max_b + v_entropy_b +   (1|id/run), df)
screen.lmerTest(wr3hpe_perf)

# purely hippocampal version:
# for plotting -- calculate negative AH and other scaled variables
wr3hpe_perf_hipp <- lmer(rt_csv ~ (scale(-1/run_trial) + scale(rt_lag) + scale(rt_vmax_lag) + omission_lag + v_max_wi_lag + v_entropy_wi + 
                                     I(-h_HippAntL) + pe_f2_hipp + scale(total_earnings))^2 + 
                      scale(rt_lag)*pe_f2_hipp*scale(total_earnings) + 
                        scale(rt_lag)*I(-h_HippAntL)*scale(total_earnings) +  
                        scale(rt_vmax_lag)*I(-h_HippAntL)*scale(total_earnings) +  
                        v_max_b + v_entropy_b +   (1|id/run), df)
screen.lmerTest(wr3hpe_perf_hipp)
# unpack interaction
####### AD stopped here 7/22/19
r1 <- emtrends(wr3hpe_perf_hipp, var = 'rt_lag', specs = c('I(-h_HippAntL)','scale(total_earnings)'))
r1 <- as.data.frame(r1)
ggplot(r1, aes(evt_time_f, bin_num_f, color = rt_vmax_change.trend)) + geom_tile()



##### BEST MODEL
wr3hpe3 <-  update(wr3hpe, . ~ . + 
                     scale(rt_lag):omission_lag:h_f1_fp + scale(rt_lag):omission_lag:I(-h_f2_neg_paralimb) + 
                     scale(rt_lag):omission_lag:pe_f1_cort_str +scale(rt_lag):omission_lag:pe_f2_hipp + 
                     scale(rt_vmax_lag):scale(-1/run_trial):h_f1_fp + scale(rt_vmax_lag):scale(-1/run_trial):I(-h_f2_neg_paralimb) + 
                     scale(rt_vmax_lag):scale(-1/run_trial):pe_f1_cort_str + scale(rt_vmax_lag):scale(-1/run_trial):pe_f2_hipp  +
                     (1|id/run), df)
screen.lmerTest(wr3hpe3, .01)




# wr4hvd3a <-  update(wr4hvd3, . ~ . + 
#                       scale(rt_vmax_lag):v_max_wi_lag:h_f1_fp + scale(rt_vmax_lag):v_max_wi_lag:I(-h_f2_neg_paralimb) + 
#                       scale(rt_vmax_lag):v_max_wi_lag:I(-v_f1_neg_cog) +scale(rt_vmax_lag):v_max_wi_lag:v_f2_paralimb + 
#                       scale(rt_vmax_lag):v_max_wi_lag:d_f1_FP_SMA +scale(rt_vmax_lag):v_max_wi_lag:d_f2_VS + scale(rt_vmax_lag):v_max_wi_lag:d_f3_ACC_ins + (1|id/run), df)
# screen.lmerTest(wr4hvd3a, .01)

anova(wr1,wr2h,wr2v,wr3hv,wr2d,wr2pe, wr2dh, wr2dhp,wr4hvd, wr4hvd3,wr3hpe3)
# at the end of the day, the h clusters explain the most
# pe is better than KLD or entropy change
###########
## PLOTS ##
#####
## "model-based analyses": plots

# behavioral exegisis on Dauc betas
ggplot(df,aes(run_trial,rt_csv, lty = pe_f2_hipp_resp)) + geom_smooth() + facet_wrap(~rewFunc)

# not a huge impact of KLD betas on RT swings...

