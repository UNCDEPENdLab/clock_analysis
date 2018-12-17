library(dplyr)
library(tidyverse)
library(psych)
library(corrplot)
library(lme4)
library(lmerTest)
source('~/code/Rhelpers/')
setwd('~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/')
load('trial_df_and_vhd_clusters.Rdata')
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

# RT by condition, trial and cluster activation (may need to recode rewFunc)
# these models suffer from multicollinearity, set aside
r1 <- lmer(rt_csv ~ (scale(-1/run_trial) + rewFuncIEVsum)^3 + (1|id/run), df)
screen.lmerTest(r1)
r2 <- lmer(rt_csv ~ (scale(-1/run_trial) + rewFuncIEVsum +h_f1_fp + I(-h_f2_neg_paralimb))^3 + (1|id/run), df)
screen.lmerTest(r2)
car::Anova(r2,'3')
r3 <- lmer(rt_csv ~ (scale(-1/run_trial) + rewFunc + I(-v_f1_neg_cog) + v_f2_paralimb)^3 + (1|id/run), df)
screen.lmerTest(r3)
car::Anova(r3,'3')
anova(r1,r2,r3)

# these simple RT swing prediction models only reveal that both networks catalyze convergence regardless of condition
mf1 <- lmer(log(rt_swing) ~ scale(-1/run_trial) * rewFuncIEVsum + (1|id/run), df[df$rt_swing>0,])
screen.lmerTest(mf1)
## favorite simple model ##
mf2 <- lmer(log(rt_swing) ~ (scale(-1/run_trial) + rewFuncIEVsum + h_f1_fp + I(-h_f2_neg_paralimb))^2 + (1|id/run), df[df$rt_swing>0,])
screen.lmerTest(mf2)
##
# add V to check for dissociation
mf3 <- lmer(log(rt_swing) ~ (scale(-1/run_trial) + rewFuncIEVsum + h_f1_fp + I(-h_f2_neg_paralimb) + I(-v_f1_neg_cog) + v_f2_paralimb)^2 + (1|id/run), df[df$rt_swing>0,])
screen.lmerTest(mf3)
# H effects stand, V does not add to fit
mf4 <- lmer(log(rt_swing) ~ (scale(-1/run_trial) + rewFuncIEVsum + h_f1_fp + I(-h_f2_neg_paralimb) + I(-v_f1_neg_cog) + v_f2_paralimb + d_f1_FP_SMA + d_f2_VS + d_f3_ACC_ins)^2 + (1|id/run), df[df$rt_swing>0,])
screen.lmerTest(mf4)
# addition of decay does not help beyond H alone
anova(mf1,mf2,mf3, mf4)


# continuing with mb* models, add within-run entropy
# within-subject models -- VIF<2.96
w1 <- lmer(log(rt_swing) ~ (scale(-1/run_trial) + omission_lag + v_max_wi + v_entropy_wi)^2 + v_max_b + v_entropy_b + scale(rt_lag) + (1|id/run), df[df$rt_swing>0,])
screen.lmerTest(w1)
w2h <- lmer(log(rt_swing) ~ (scale(-1/run_trial) + omission_lag + v_max_wi + v_entropy_wi + h_f1_fp + I(-h_f2_neg_paralimb))^2 + v_max_b + v_entropy_b + scale(rt_lag) + (1|id/run), df[df$rt_swing>0,])
screen.lmerTest(w2h)
w2v <- lmer(log(rt_swing) ~ (scale(-1/run_trial) + omission_lag + v_max_wi + v_entropy_wi + I(-v_f1_neg_cog) + v_f2_paralimb)^2 + v_max_b + v_entropy_b + scale(rt_lag) + (1|id/run), df[df$rt_swing>0,])
screen.lmerTest(w2v)
w3hv <- lmer(log(rt_swing) ~ (scale(-1/run_trial) + omission_lag + v_max_wi + v_entropy_wi + h_f1_fp + I(-h_f2_neg_paralimb) + I(-v_f1_neg_cog) + v_f2_paralimb)^2 + v_max_b + v_entropy_b + scale(rt_lag) + (1|id/run), df[df$rt_swing>0,])
screen.lmerTest(w3hv)
w2d <- lmer(log(rt_swing) ~ (scale(-1/run_trial) + omission_lag + v_max_wi + v_entropy_wi + d_f1_FP_SMA + d_f2_VS + d_f3_ACC_ins)^2 + v_max_b + v_entropy_b + scale(rt_lag) + (1|id/run), df[df$rt_swing>0,])
screen.lmerTest(w2d)
w4hvd <- lmer(log(rt_swing) ~ (scale(-1/run_trial) + omission_lag + v_max_wi + v_entropy_wi + h_f1_fp + I(-h_f2_neg_paralimb) + 
                                 I(-v_f1_neg_cog) + v_f2_paralimb + 
                                 d_f1_FP_SMA + d_f2_VS + d_f3_ACC_ins)^2 + v_max_b + v_entropy_b + scale(rt_lag) + (1|id/run), df[df$rt_swing>0,])
screen.lmerTest(w4hvd)
anova(w1,w2h,w2v,w3hv,w2d,w4hvd)
# plot out scale(run_trial):v_entropy_wi:h_f2_neg_paralimb
ggplot(df, aes(v_entropy_wi, log(rt_swing),color = low_h_paralimbic, lty = h_fp)) +
  geom_smooth(method = 'gam') + facet_wrap(~run_trial>20)


# same with RT instead of RT swing
# deal with multicollinearity in RT models
vv <- df[,c("rt_lag", "rt_vmax", "rt_vmax_lag", "run_trial", "omission_lag", "v_max_wi", "v_max_wi_lag", "v_entropy_wi")]
corr.test(vv)
# NB: omissions have two effects: RT shortening and non-directional RT swing
# collinearity remediated by taking lags of rt_vmax and v_max, conceptually the same variables, just discard the effect of last outcome
wr1 <- lmer(rt_csv ~ (scale(-1/run_trial) + scale(rt_lag) + scale(rt_vmax_lag) + omission_lag + v_max_wi_lag + v_entropy_wi)^2 + v_max_b + v_entropy_b + scale(rt_lag) + (1|id/run), df)
screen.lmerTest(wr1)
wr2h <- lmer(rt_csv ~ (scale(-1/run_trial) + scale(rt_lag) + scale(rt_vmax_lag) + omission_lag + v_max_wi_lag + v_entropy_wi + h_f1_fp + I(-h_f2_neg_paralimb))^2 + v_max_b + v_entropy_b +  (1|id/run), df)
screen.lmerTest(wr2h)
wr2v <- lmer(rt_csv ~ (scale(-1/run_trial) + scale(rt_lag) + scale(rt_vmax_lag) + omission_lag + v_max_wi_lag + v_entropy_wi + I(-v_f1_neg_cog) + v_f2_paralimb)^2 + v_max_b + v_entropy_b +   (1|id/run), df)
screen.lmerTest(wr2v)
wr3hv <- lmer(rt_csv ~ (scale(-1/run_trial) + scale(rt_lag) + scale(rt_vmax_lag) + omission_lag + v_max_wi_lag  + h_f1_fp + I(-h_f2_neg_paralimb) + I(-v_f1_neg_cog) + v_f2_paralimb)^2 + v_max_b + v_entropy_b +  (1|id/run), df)
screen.lmerTest(wr3hv)
wr2d <- lmer(rt_csv ~ (scale(-1/run_trial) + scale(rt_lag) + scale(rt_vmax_lag) + omission_lag + v_max_wi_lag  + d_f1_FP_SMA + d_f2_VS + d_f3_ACC_ins)^2 + v_max_b + v_entropy_b +  (1|id/run), df)
screen.lmerTest(wr2d)
wr4hvd <- lmer(rt_csv ~ (scale(-1/run_trial) + scale(rt_lag) + scale(rt_vmax_lag) + omission_lag + v_max_wi_lag  + h_f1_fp + I(-h_f2_neg_paralimb) + 
                                 I(-v_f1_neg_cog) + v_f2_paralimb + 
                                 d_f1_FP_SMA + d_f2_VS + d_f3_ACC_ins)^2 + v_max_b + v_entropy_b +  (1|id/run), df)
screen.lmerTest(wr4hvd, .01)
wr4hvd3 <-  update(wr4hvd, . ~ . + 
                     scale(rt_lag):omission_lag:h_f1_fp + scale(rt_lag):omission_lag:I(-h_f2_neg_paralimb) + 
                  scale(rt_lag):omission_lag:I(-v_f1_neg_cog) +scale(rt_lag):omission_lag:v_f2_paralimb + 
                    scale(rt_lag):omission_lag:d_f1_FP_SMA +scale(rt_lag):omission_lag:d_f2_VS + scale(rt_lag):omission_lag:d_f3_ACC_ins + (1|id/run), df)
screen.lmerTest(wr4hvd3, .01)

wr4hvd3a <-  update(wr4hvd3, . ~ . + 
                      scale(rt_vmax_lag):v_max_wi_lag:h_f1_fp + scale(rt_vmax_lag):v_max_wi_lag:I(-h_f2_neg_paralimb) + 
                      scale(rt_vmax_lag):v_max_wi_lag:I(-v_f1_neg_cog) +scale(rt_vmax_lag):v_max_wi_lag:v_f2_paralimb + 
                      scale(rt_vmax_lag):v_max_wi_lag:d_f1_FP_SMA +scale(rt_vmax_lag):v_max_wi_lag:d_f2_VS + scale(rt_vmax_lag):v_max_wi_lag:d_f3_ACC_ins + (1|id/run), df)
screen.lmerTest(wr4hvd3a, .01)

anova(wr1,wr2h,wr2v,wr3hv,wr2d,wr4hvd)
anova(wr1,wr2h,wr4hvd)
# v and d add a bit to H, but H is the best among single-signal models

###########
## PLOTS ##
#####
## "model-based analyses": plots

# behavioral exegisis on Dauc betas
# crescendo contributions of decay factors to RT swings
ggplot(df,aes(run_trial,log(rt_swing), lty = d_f1_FP_SMAresp)) + geom_smooth() + facet_wrap(~rewFunc)
ggplot(df,aes(run_trial,log(rt_swing), size = d_f2_VSresp)) + geom_smooth() + facet_wrap(~rewFunc)
ggplot(df,aes(run_trial,log(rt_swing), color = d_f3_ACC_ins_resp)) + geom_smooth() + facet_wrap(~rewFunc)
ggplot(df,aes(run_trial,log(rt_swing), color = d_f3_ACC_ins_resp, lty = d_f1_FP_SMAresp, size = d_f2_VSresp)) + geom_smooth() + facet_wrap(~rewFunc)

ggplot(df,aes(v_max_wi,log(rt_swing), lty = d_f1_FP_SMAresp)) + geom_smooth() + facet_wrap(~rewFunc)
ggplot(df,aes(v_max_wi,log(rt_swing), size = d_f2_VSresp)) + geom_smooth() + facet_wrap(~rewFunc)
ggplot(df,aes(v_max_wi,log(rt_swing), color = d_f3_ACC_ins_resp)) + geom_smooth() + facet_wrap(~rewFunc)
# this may mean that Dauc indexes heterogeneous cognitive processes
# does decay control H?
p1 <- ggplot(df[df$run>1,],aes(run_trial,v_entropy, color = d_f1_FP_SMAresp)) + geom_smooth() + facet_wrap(~rewFunc) + scale_color_brewer(palette = 'Set3', direction = -1)+ scale_x_continuous(breaks = c(1,50)) + guides(color=FALSE) + theme(axis.title.x=element_blank())
p2 <- ggplot(df[df$run>1,],aes(run_trial,v_entropy, color = d_f2_VSresp)) + geom_smooth() + facet_wrap(~rewFunc) + scale_color_brewer(palette = 'Dark2', direction = 1)+ scale_x_continuous(breaks = c(1,50)) + guides(color=FALSE)+ theme(axis.title.x=element_blank())
p3 <- ggplot(df[df$run>1,],aes(run_trial,v_entropy, color = d_f3_ACC_ins_resp)) + geom_smooth() + facet_wrap(~rewFunc) + scale_color_brewer(palette = 'Set1', direction = 1)+ scale_x_continuous(breaks = c(1,50)) + guides(color=FALSE)
p4 <- ggplot(df[df$run>1,],aes(run_trial,v_entropy, color = low_v_fp_acc_vlpfc)) + geom_smooth() + facet_wrap(~rewFunc) + scale_color_brewer(palette = 'Paired', direction = -1)+ scale_x_continuous(breaks = c(1,50)) + guides(color=FALSE)+ theme(axis.title.x=element_blank())
p5 <- ggplot(df[df$run>1,],aes(run_trial,v_entropy, color = v_paralimbic)) + geom_smooth() + facet_wrap(~rewFunc) + scale_color_brewer(palette = 'Pastel1', direction = 1)+ scale_x_continuous(breaks = c(1,50)) + guides(color=FALSE)+ theme(axis.title.x=element_blank())
p6 <- ggplot(df[df$run>1,],aes(run_trial,v_entropy, color = h_fp)) + geom_smooth() + facet_wrap(~rewFunc) + scale_color_brewer(palette = 'Pastel2', direction = -1)+ scale_x_continuous(breaks = c(1,50)) + guides(color=FALSE)+ theme(axis.title.x=element_blank())
p7 <- ggplot(df[df$run>1,],aes(run_trial,v_entropy, color = low_h_paralimbic)) + geom_smooth() + facet_wrap(~rewFunc) + scale_color_brewer(palette = 'Accent', direction = -1)+ scale_x_continuous(breaks = c(1,50)) + guides(color=FALSE)+ theme(axis.title.x=element_blank())
pdf("all_betas_H_timecourse.pdf", width = 12, height = 18)
ggarrange(p4,p5,p6,p7,p1,p2,p3, ncol = 2,nrow = 7, common.legend = F, align = c('v'))
dev.off()

t1 <- ggplot(df[df$run>1,],aes(run_trial,log(rt_swing), color = d_f1_FP_SMAresp)) + geom_smooth() + facet_wrap(~rewFunc) + scale_color_brewer(palette = 'Set3', direction = -1) + scale_x_continuous(breaks = c(1,50)) + guides(color=FALSE)+ theme(axis.title.x=element_blank())
t2 <- ggplot(df[df$run>1,],aes(run_trial,log(rt_swing), color = d_f2_VSresp)) + geom_smooth() + facet_wrap(~rewFunc) + scale_color_brewer(palette = 'Dark2', direction = 1)+ scale_x_continuous(breaks = c(1,50)) + guides(color=FALSE)+ theme(axis.title.x=element_blank())
t3 <- ggplot(df[df$run>1,],aes(run_trial,log(rt_swing), color = d_f3_ACC_ins_resp)) + geom_smooth() + facet_wrap(~rewFunc) + scale_color_brewer(palette = 'Set1', direction = 1)+ scale_x_continuous(breaks = c(1,50)) + guides(color=FALSE)
t4 <- ggplot(df[df$run>1,],aes(run_trial,log(rt_swing), color = low_v_fp_acc_vlpfc)) + geom_smooth() + facet_wrap(~rewFunc) + scale_color_brewer(palette = 'Paired', direction = -1)+ scale_x_continuous(breaks = c(1,50)) + guides(color=FALSE)+ theme(axis.title.x=element_blank())
t5 <- ggplot(df[df$run>1,],aes(run_trial,log(rt_swing), color = v_paralimbic)) + geom_smooth() + facet_wrap(~rewFunc) + scale_color_brewer(palette = 'Pastel1', direction = 1)+ scale_x_continuous(breaks = c(1,50)) + guides(color=FALSE)+ theme(axis.title.x=element_blank())
t6 <- ggplot(df[df$run>1,],aes(run_trial,log(rt_swing), color = h_fp)) + geom_smooth() + facet_wrap(~rewFunc) + scale_color_brewer(palette = 'Pastel2', direction = -1)+ scale_x_continuous(breaks = c(1,50)) + guides(color=FALSE)+ theme(axis.title.x=element_blank())
t7 <- ggplot(df[df$run>1,],aes(run_trial,log(rt_swing), color = low_h_paralimbic)) + geom_smooth() + facet_wrap(~rewFunc) + scale_color_brewer(palette = 'Accent', direction = -1)+ scale_x_continuous(breaks = c(1,50)) + guides(color=FALSE)+ theme(axis.title.x=element_blank())
pdf("all_betas_RT_swing_timecourse.pdf", width = 12, height = 18)
ggarrange(t4,t5,t6,t7,t1,t2,t3, ncol = 2,nrow = 7, common.legend = F, align = c('v'))
dev.off()

r1 <- ggplot(df[df$run>1,],aes(run_trial,rt_csv, color = d_f1_FP_SMAresp)) + geom_smooth() + facet_wrap(~rewFunc) + scale_color_brewer(palette = 'Set3', direction = -1)+ scale_x_continuous(breaks = c(1,50))+ theme(axis.title.x=element_blank())
r2 <- ggplot(df[df$run>1,],aes(run_trial,rt_csv, color = d_f2_VSresp)) + geom_smooth() + facet_wrap(~rewFunc) + scale_color_brewer(palette = 'Dark2', direction = 1)+ scale_x_continuous(breaks = c(1,50))+ theme(axis.title.x=element_blank())
r3 <- ggplot(df[df$run>1,],aes(run_trial,rt_csv, color = d_f3_ACC_ins_resp)) + geom_smooth() + facet_wrap(~rewFunc) + scale_color_brewer(palette = 'Set1', direction = 1)+ scale_x_continuous(breaks = c(1,50))
r4 <- ggplot(df[df$run>1,],aes(run_trial,rt_csv, color = low_v_fp_acc_vlpfc)) + geom_smooth() + facet_wrap(~rewFunc) + scale_color_brewer(palette = 'Paired', direction = -1)+ scale_x_continuous(breaks = c(1,50))+ theme(axis.title.x=element_blank())
r5 <- ggplot(df[df$run>1,],aes(run_trial,rt_csv, color = v_paralimbic)) + geom_smooth() + facet_wrap(~rewFunc) + scale_color_brewer(palette = 'Pastel1', direction = 1)+ scale_x_continuous(breaks = c(1,50))+ theme(axis.title.x=element_blank())
r6 <- ggplot(df[df$run>1,],aes(run_trial,rt_csv, color = h_fp)) + geom_smooth() + facet_wrap(~rewFunc) + scale_color_brewer(palette = 'Pastel2', direction = -1)+ scale_x_continuous(breaks = c(1,50))+ theme(axis.title.x=element_blank())
r7 <- ggplot(df[df$run>1,],aes(run_trial,rt_csv, color = low_h_paralimbic)) + geom_smooth() + facet_wrap(~rewFunc) + scale_color_brewer(palette = 'Accent', direction = -1)+ scale_x_continuous(breaks = c(1,50))+ theme(axis.title.x=element_blank())
pdf("all_betas_RT_timecourse.pdf", width = 12, height = 18)
ggarrange(r4,r5,r6,r7,r1,r2,r3, ncol = 2,nrow = 7, common.legend = F, align = c('v'))
dev.off()

allplots <- ggarrange(p4, t4, r4,
          p5, t5, r5,
          p6, t6, r6,
          p7, t7, r7,
          p1, t1, r1,
          p2, t2, r2,
          p3, t3, r3, ncol = 3,nrow = 7, common.legend = F)


pdf("all_betas_H_rt_swing_RT_timecourse.pdf", width = 12, height = 16)
ggarrange(p4, t4, r4,
          p5, t5, r5,
          p6, t6, r6,
          p7, t7, r7,
          p1, t1, r1,
          p2, t2, r2,
          p3, t3, r3, ncol = 3,nrow = 7, common.legend = F)
dev.off()


ggplot(df,aes(run_trial,v_entropy_wi, lty = d_f1_FP_SMAresp)) + geom_smooth() + facet_wrap(~rewFunc)
ggplot(df,aes(run_trial,v_entropy_wi, size = d_f2_VSresp)) + geom_smooth() + facet_wrap(~rewFunc)
ggplot(df,aes(run_trial,v_entropy_wi, color = d_f3_ACC_ins_resp)) + geom_smooth() + facet_wrap(~rewFunc) + scale_color_brewer(palette = 'BrBG', direction = 1)


# reward omission
ggplot(na.omit(df),aes(d_f1_FP_SMAresp,log(rt_swing), lty = last_outcome)) + geom_boxplot()
ggplot(na.omit(df),aes(d_f2_VSresp,log(rt_swing), lty = last_outcome)) + geom_boxplot()
ggplot(na.omit(df),aes(d_f3_ACC_ins_resp,log(rt_swing), lty = last_outcome)) + geom_boxplot()

ggplot(na.omit(df),aes(last_outcome,log(rt_swing), lty = d_f1_FP_SMAresp, size = d_f2_VSresp, color = d_f3_ACC_ins_resp)) + geom_boxplot()
ggplot(na.omit(df),aes(d_f2_VSresp,log(rt_swing), lty = last_outcome)) + geom_boxplot()
ggplot(na.omit(df),aes(d_f3_ACC_ins_resp,log(rt_swing), lty = last_outcome)) + geom_boxplot()


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

ggplot(df, aes(v_max,rt_swing,  color = low_h_paralimbic, lty = h_fp)) + geom_smooth(method = "gam")

# plot vs. raw data
ggplot(df, aes(rt_lag, rt_csv, color = low_h_paralimbic, lty = h_fp)) + geom_smooth(method = "glm")+ facet_wrap(~run_trial>20)
ggplot(df, aes(rt_vmax, rt_csv, color = low_h_paralimbic, lty = h_fp)) + geom_smooth(method = "glm")+ facet_wrap(~run_trial>20)
pdf("h_timecourse_brain_fixed.pdf", width = 8, height = 8)
ggplot(df, aes(run_trial, v_entropy,color = low_h_paralimbic, lty = h_fp)) + geom_smooth(method = "loess") + facet_wrap(~gamma>0)
dev.off()


ggplot(df, aes(run_trial, v_entropy, color = rewFunc)) + geom_smooth()
summary(m0 <- lmer(log(rt_swing) ~ scale(-1/run_trial) * scale(v_max) * scale(v_entropy) * rewFunc + (1|id/run), df[df$rt_swing>0,]))

# plot out


save(file = 'vhd_models.Rdata', list = ls(all.names = TRUE))
 # load('vhd_models.Rdata')
# 
# # also ran:
# # RT swings analyses: "exploration"
# # NB: covarying for contingency does not change anything, keep as sensitivity analysis for the future
# mb1 <- lmer(log(rt_swing) ~ scale(-1/run_trial) * scale(v_max_wi) + v_max_b + v_entropy_b + (1|id/run), df[df$rt_swing>0,])
# screen.lmerTest(mb1)
# # effect of v_max:beta was suspiciously changing sign depending on whether we log-transformed RT swings.
# # # plot to understand the relationship
# # ggplot(df, aes(v_max_wi, log(rt_swing), color = low_h_paralimbic, lty = h_fp)) + geom_smooth()
# # remediated by decomposing within-run vs. between-r v_max
# # NB cannot allow between-run Vmax and H to interact with betas, multicolllinearity explodes (VIFs>300)
# mb2h <- lmer(log(rt_swing) ~ (scale(-1/run_trial) + scale(v_max_wi) + h_f1_fp + I(-h_f2_neg_paralimb)) ^2 + v_max_b + v_entropy_b + (1|id/run), df[df$rt_swing>0,])
# screen.lmerTest(mb2h)
# mb2v <- lmer(log(rt_swing) ~ (scale(-1/run_trial) + scale(v_max_wi) + I(-v_f1_neg_cog) + v_f2_paralimb) ^3  +v_max_b + v_entropy_b + (1|id/run), df[df$rt_swing>0,])
# screen.lmerTest(mb2v)
# mb2d <- lmer(log(rt_swing) ~ (scale(-1/run_trial) + scale(v_max_wi) + d_f1_FP_SMA + d_f2_VS + d_f3_ACC_ins) ^2 + v_max_b + v_entropy_b  + (1|id/run), df[df$rt_swing>0,])
# screen.lmerTest(mb2d)
# mb3hv <- lmer(log(rt_swing) ~ (scale(-1/run_trial) + scale(v_max_wi) + h_f1_fp + I(-h_f2_neg_paralimb) + I(-v_f1_neg_cog) + v_f2_paralimb) ^2 + v_max_b + v_entropy_b  + (1|id/run), df[df$rt_swing>0,])
# screen.lmerTest(mb3hv)
# # adding d introduces multi-collinearity
# mb4vhd <- lmer(log(rt_swing) ~ (scale(-1/run_trial) + scale(v_max_wi) + h_f1_fp + I(-h_f2_neg_paralimb) + I(-v_f1_neg_cog) + v_f2_paralimb + d_f1_FP_SMA + d_f2_VS + d_f3_ACC_ins) ^2 + v_max_b + v_entropy_b  + (1|id/run), df[df$rt_swing>0,])
# screen.lmerTest(mb4vhd)
# anova(mb1,mb2h,mb2v,mb2d,mb3hv,mb4vhd)
# # H, V, Dauc all explain unique variance.  The effect of ACC/mid-insula dAUC factor is striking -- should we resurrect the decay story?

# ###########
# # EV analyses: "exploitation" -- run removed from RE to avoid singular fit
# # Differences emerge in IEV (explains why total earnings are not informative)
# # don't love these models because of multicollinearity, keep in reserve
# # summary(ev1 <- lmer(ev ~ scale(-1/run_trial) * rewFunc + (1|id), df))
# # summary(ev2 <- lmer(ev ~ (I(-scale(1/run_trial)) + rewFunc + h_f1_fp) ^3 + (1|id), df))
# # summary(ev3 <- lmer(ev ~ (I(-scale(1/run_trial)) + rewFunc + I(-h_f2_neg_paralimb)) ^3 + (1|id), df))
# ev4 <- lmerTest::lmer(ev ~ (scale(-1/run_trial) + rewFuncIEVsum + h_f1_fp + I(-h_f2_neg_paralimb)) ^3 + (1|id), df)
# screen.lmerTest(ev4)
# # anova(ev1,ev2,ev3,ev4)
# 
# pdf('ev_by_condition_and_h_betas.pdf', height = 8, width = 8)
# ggplot(df, aes(run_trial, ev, color = low_h_paralimbic, lty = h_fp)) + geom_smooth(method = "loess") + facet_wrap (~rewFunc)
# dev.off()
# pdf('ev_by_condition_and_v_betas.pdf', height = 8, width = 8)
# ggplot(df, aes(run_trial, ev, color = v_paralimbic, lty = low_v_fp_acc_vlpfc)) + geom_smooth(method = "loess") + facet_wrap (~rewFunc)
# dev.off()
# pdf('ev_by_condition_v_maxWi_and_h_betas.pdf', height = 8, width = 8)
# ggplot(df, aes(v_max_wi, ev, color = low_h_paralimbic, lty = h_fp)) + geom_smooth(method = "gam") + facet_wrap (~rewFunc)
# dev.off()
# 
# # just entropy -- VIFs ok
# rt0 <- lmer(scale(rt_csv) ~ (scale(-1/run_trial) + scale(rt_vmax) + scale(rt_lag))^3 + rewFunc + (1|id/run), df)
# screen.lmerTest(rt0)
# rth <- lmerTest::lmer(scale(rt_csv) ~ (scale(-1/run_trial) + scale(rt_vmax) + scale(rt_lag) + h_f1_fp + I(-h_f2_neg_paralimb))^3 + rewFunc + (1|id/run), df)
# screen.lmerTest(rth)
# car::Anova(rt1)
# ggplot(df,aes(rt_vmax,rt_csv, color = low_h_paralimbic, lty = h_fp)) + facet_wrap(~run_trial>20) + geom_smooth(method = "glm")
# ggplot(df,aes(rt_lag,rt_csv, color = low_h_paralimbic, lty = h_fp)) + facet_wrap(~run_trial>20) + geom_smooth(method = "glm")
# 
# # just value -- VIFs ok
# rtv <- lmer(scale(rt_csv) ~ (scale(-1/run_trial) + scale(rt_vmax) + scale(rt_lag) + I(-v_f1_neg_cog) + v_f2_paralimb)^3 + rewFunc + (1|id/run), df)
# screen.lmerTest(rtv)
# # + d_auc
# rtd <- lmer(scale(rt_csv) ~ (scale(-1/run_trial) + scale(rt_vmax) + scale(rt_lag) + d_f1_FP_SMA + d_f2_VS + d_f3_ACC_ins)^3 + rewFunc + (1|id/run), df)
# screen.lmerTest(rtd)
# 
# # only H and V
# rtvh <- lmer(scale(rt_csv) ~ (scale(-1/run_trial) + scale(rt_vmax) + scale(rt_lag) + h_f1_fp + I(-h_f2_neg_paralimb) + I(-v_f1_neg_cog) +v_f2_paralimb)^3 + rewFunc + (1|id/run), df)
# # all
# rtvhd <- lmer(scale(rt_csv) ~ (scale(-1/run_trial) + scale(rt_vmax) + scale(rt_lag) + h_f1_fp + I(-h_f2_neg_paralimb)  + d_f1_FP_SMA + d_f2_VS + d_f3_ACC_ins)^3 + rewFunc + (1|id/run), df)
# screen.lmerTest(rtvhd,.05)
# anova(rt0,rtv,rtd,rth,rtvh,rtvhd) # H predicts best, each addition (V, D_AUC) improves the prediction further (by >300 AIC points)
# # but careful with cluster*cluster interactions -- high VIFs!
# 
# 
# # coupling between entropy and RT swings -- betas don't moderate effects of entropy.  Rather, they influence entropy as it unfolds.
# # summary(m7 <- lmer(log(rt_swing) ~ (scale(-1/run_trial) + scale(v_max) + scale(v_entropy))^2 + rewFunc + (1|id/run), df[df$rt_swing>0,]))
# # summary(m8 <- lmer(log(rt_swing) ~ (scale(-1/run_trial) + scale(v_max) + scale(v_entropy) + h_f1_fp) ^2 + rewFunc + (1|id/run), df[df$rt_swing>0,]))
# # car::Anova(m8, '3')
# # summary(m9 <- lmer(log(rt_swing) ~ (scale(-1/run_trial) + scale(v_max) + scale(v_entropy) + h_f2_neg_paralimb) ^2 + rewFunc + (1|id/run), df[df$rt_swing>0,]))
# # car::Anova(m9, '3')
# # summary(m10 <- lmer(log(rt_swing) ~ (scale(-1/run_trial) + scale(v_max) + scale(v_entropy) +h_f1_fp + h_f2_neg_paralimb) ^2 + rewFunc + (1|id/run), df[df$rt_swing>0,]))
# # car::Anova(m10, '3')
# # # a bunch of relatively weak interactions, set a side for now:
# # summary(m11 <- lmer(log(rt_swing) ~ (scale(-1/run_trial) + scale(v_max) + scale(v_entropy) +h_f1_fp + h_f2_neg_paralimb + v_f1_neg_cog + v_f2_paralimb) ^3 + rewFunc + (1|id/run), df[df$rt_swing>0,]))
# # car::Anova(m11, '3')
# # anova(m7,m8,m9, m10)
# 
# # plot out the relationships
# ggplot(df, aes(v_entropy,rt_swing, color = h_f1_fp>0)) + geom_smooth(method = "gam")
# 
# # toward the goal of dissociating exploitative behavioral adjustment from H-driven exploration
# # dm = dissociation model
# # unfortunately does not work all that well due to collinearity
# summary(dm1 <- lmer(log(rt_swing) ~ (scale(-1/run_trial) + omission_lag + scale(v_max))^3 + scale(rt_lag) + (1|id/run), df[df$rt_swing>0,]))
# ggplot(na.omit(df), aes(v_entropy, rt_swing, color = v_max>25, lty = omission_lag)) + geom_smooth(method = "gam")
# 
# summary(dm2 <- lmer(rt_csv ~ (scale(rt_lag) + omission_lag + scale(rt_vmax))^3 + rewFunc + (1|id/run), df[df$rt_swing>0,]))
# 
# summary(dm3 <- lmer(rt_csv ~ (scale(rt_lag) + omission_lag + scale(rt_vmax) + h_f1_fp + I(-h_f2_neg_paralimb))^3 + rewFunc + (1|id/run), df[df$rt_swing>0,]))
# 
# # more focused test of dislodgment
# summary(dm3a <- lmer(rt_csv ~ (scale(rt_lag) + omission_lag + scale(rt_vmax) + h_f1_fp)^3 + rewFunc + (1|id/run), df[df$rt_swing>0,]))
# summary(dm3b <- lmer(rt_csv ~ (scale(rt_lag) + omission_lag + scale(rt_vmax) + I(-h_f2_neg_paralimb))^3 + rewFunc + (1|id/run), df[df$rt_swing>0,]))
# 
# # plot 3-way interactions
# ggplot(df, aes(rt_lag, rt_csv, color = omission_lag, lty = h_f1_fp>0)) + geom_smooth(method = "glm")
# ggplot(df, aes(rt_lag, rt_csv, color = omission_lag, lty = h_f2_neg_paralimb<0)) + geom_smooth(method = "glm")
# # the interesting one:
# ggplot(na.omit(df), aes(rt_vmax, rt_csv, color = omission_lag, lty = h_f1_fp>0)) + geom_smooth(method = "glm")
# ggplot(na.omit(df), aes(rt_vmax, rt_csv, color = omission_lag, lty = h_f2_neg_paralimb<0)) + geom_smooth(method = "glm")
# 
# anova(dm1,dm2,dm3)
# 
# summary(dm2 <- lmer(log(rt_swing) ~ (scale(-1/run_trial) +  scale(v_entropy) + omission_lag + scale(v_max) + h_f1_fp + I(-h_f2_neg_paralimb))^3 + (1|id/run), df[df$rt_swing>0,]))
# anova(dm1,dm2)
# 
# 
