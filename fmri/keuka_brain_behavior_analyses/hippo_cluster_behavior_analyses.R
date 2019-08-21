library(dplyr)
library(tidyverse)
library(psych)
library(corrplot)
library(lme4)
library(lmerTest)
library(emmeans)
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

# for plotting -- calculate negative AH and other scaled variables
df <- df %>% group_by(id, run) %>% arrange(id, run, run_trial)  %>% 
  mutate(trial_neg_inv_sc = scale(-1/run_trial),
         rt_lag_sc = scale(rt_lag),
         rt_vmax_lag_sc = scale(rt_vmax_lag),
         v_max_wi_lag_sc = scale(v_max_wi_lag),
         v_entropy_wi_sc = scale(v_entropy_wi)) %>% ungroup() %>% mutate(AH_neg_sc = scale(-h_HippAntL),total_earnings_sc = scale(total_earnings))

# deal with multicollinearity in RT models
vv <- df[,c("rt_lag", "rt_vmax", "rt_vmax_lag", "run_trial", "omission_lag", "v_max_wi", "v_max_wi_lag", "v_entropy_wi")]
corr.test(vv)
# NB: omissions have two effects: RT shortening and non-directional RT swing
# collinearity remediated by taking lags of rt_vmax and v_max, conceptually the same variables, just discard the effect of last outcome
wr1 <- lmer(rt_csv ~ (trial_neg_inv_sc + rt_lag_sc + rt_vmax_lag_sc + omission_lag + v_entropy_wi)^2 + v_max_b + v_entropy_b + rt_lag_sc + (1|id/run), df)
screen.lmerTest(wr1)
wr2h <- lmer(rt_csv ~ (trial_neg_inv_sc + rt_lag_sc + rt_vmax_lag_sc + omission_lag + v_entropy_wi + h_f1_fp + I(-h_f2_neg_paralimb))^2 + v_max_b + v_entropy_b +  (1|id/run), df)
screen.lmerTest(wr2h)
wr3hpe <- lmer(rt_csv ~ (trial_neg_inv_sc + rt_lag_sc + rt_vmax_lag_sc + omission_lag + v_entropy_wi + AH_neg_sc + pe_f2_hipp)^2 + v_max_b + v_entropy_b +   (1|id/run), df)
screen.lmerTest(wr3hpe)
# control for performance
wr3hpe_perf <- lmer(rt_csv ~ (trial_neg_inv_sc + rt_lag_sc + rt_vmax_lag_sc + omission_lag + v_max_wi_lag + v_entropy_wi + h_f1_fp + I(-h_f2_neg_paralimb) + pe_f1_cort_str + pe_f2_hipp + scale(total_earnings))^2 + 
                      rt_lag_sc*pe_f2_hipp*scale(total_earnings) + v_max_b + v_entropy_b +   (1|id/run), df)
screen.lmerTest(wr3hpe_perf)

# hippocampal model for publication:
wr3hpe_hipp <- lmer(rt_csv ~ (trial_neg_inv_sc + rt_lag_sc + rt_vmax_lag_sc + omission_lag + v_entropy_wi_sc + 
                                     AH_neg_sc + pe_f2_hipp)^2 + 
                           v_max_b + v_entropy_b + trial_neg_inv_sc*rewFunc + (1|id/run), df)
screen.lmerTest(wr3hpe_hipp, .01)
summary(wr3hpe_hipp)

# interactions with performance
wr3hpe_perf_hipp <- lmer(rt_csv ~ (trial_neg_inv_sc + rt_lag_sc + rt_vmax_lag_sc + omission_lag + v_entropy_wi_sc + 
                                     AH_neg_sc + pe_f2_hipp + total_earnings_sc)^2 + 
                      rt_lag_sc*pe_f2_hipp*total_earnings_sc + 
                        rt_lag_sc*AH_neg_sc*total_earnings_sc +  
                        v_max_b + v_entropy_b + trial_neg_inv_sc*rewFunc + (1|id/run), df)
screen.lmerTest(wr3hpe_perf_hipp)
# unpack AH*performance(earnings) interaction: the increased stickiness with higher AH is seen mostly in poorly performing subjects
# interpretation: they 
r1 <- emtrends(wr3hpe_perf_hipp, var = "rt_lag_sc", specs = c("AH_neg_sc","total_earnings_sc"), at = (list(rt_lag_sc = c(-2,0,2), AH_neg_sc  = c(-2,0,2), total_earnings_sc =c(-2,0,2))))
r1 <- as.data.frame(r1)
ggplot(r1, aes(total_earnings_sc, rt_lag_sc.trend, color = AH_neg_sc)) + geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL)) + geom_point(size = 5)

ggplot(r1, aes(rt_lag_sc, emmean, color = AH_neg_sc)) + geom_point(aes(rt_lag_sc, emmean, color = AH_neg_sc)) + geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL)) + facet_wrap(~total_earnings_sc)



##### more complete model, not significantly better than wr3hpe (2 AIC points)
wr3hpe3 <-  update(wr3hpe, . ~ . + 
                     rt_lag_sc:v_entropy_wi_sc:AH_neg_sc + rt_lag_sc:v_entropy_wi_sc:pe_f2_hipp + 
                     rt_vmax_lag_sc:trial_neg_inv_sc:AH_neg_sc + rt_vmax_lag_sc:trial_neg_inv_sc:pe_f2_hipp  +
                     (1|id/run), df)
screen.lmerTest(wr3hpe3)
car::Anova(wr3hpe3)
r1 <- emtrends(wr3hpe3, var = "rt_vmax_lag_sc", specs = c("AH_neg_sc"), at = (list(rt_vmax_lag_sc = c(-2,0,2), AH_neg_sc  = c(-2,0,2))))
r1 <- as.data.frame(r1)
ggplot(r1, aes(AH_neg_sc, rt_vmax_lag_sc.trend)) + geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL)) + geom_point(size = 5)
r1 <- emmeans(wr3hpe3, specs = c("rt_vmax_lag_sc","AH_neg_sc"), at = (list(rt_vmax_lag_sc = c(-2,0,2), AH_neg_sc  = c(-2,0,2))))
r1 <- as.data.frame(r1)
r1$`Predicted RT` <- r1$emmean
ggplot(r1, aes(rt_vmax_lag_sc, `Predicted RT`, group = AH_neg_sc, color = AH_neg_sc)) + geom_line() + geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL)) + geom_point(size = 5)



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

