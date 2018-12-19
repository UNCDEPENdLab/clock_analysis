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

###########
## PLOTS ##
#####

# BS by trial and condition
df$decay <- NA
df$decay[df$gamma>0] <- 'high'
df$decay[df$gamma<0] <- 'low'

h <-  ggplot(df,aes(run_trial, v_entropy, color = decay)) + geom_smooth() + facet_wrap(~rewFunc)
v <-  ggplot(df,aes(run_trial, v_max, color = decay)) + geom_smooth() + facet_wrap(~rewFunc)
d <-  ggplot(df,aes(run_trial, d_auc, color = decay)) + geom_smooth() + facet_wrap(~rewFunc)

h1 <- ggplot(df,aes(run_trial, hb_f1_DAN_vlPFC, color = decay)) + geom_smooth() + facet_wrap(~rewFunc)
h2 <- ggplot(df,aes(run_trial, hb_f2_neg_paralimb, color = decay)) + geom_smooth() + facet_wrap(~rewFunc)
v1 <- ggplot(df,aes(run_trial, vb_f1_lo_DAN, color = decay)) + geom_smooth() + facet_wrap(~rewFunc)
v2 <- ggplot(df,aes(run_trial, vb_f2_hi_vmPFC_cOFC, color = decay)) + geom_smooth() + facet_wrap(~rewFunc)
v3 <- ggplot(df,aes(run_trial, vb_f3_lo_ACC, color = decay)) + geom_smooth() + facet_wrap(~rewFunc)
v4 <- ggplot(df,aes(run_trial, vb_f4_lo_cerebell_crus, color = decay)) + geom_smooth() + facet_wrap(~rewFunc)
v5 <- ggplot(df,aes(run_trial, vb_f5_hi_blITG, color = decay)) + geom_smooth() + facet_wrap(~rewFunc)
d1 <- ggplot(df,aes(run_trial, db_f1_rIFG_rSMA, color = decay)) + geom_smooth() + facet_wrap(~rewFunc)
d2 <- ggplot(df,aes(run_trial, db_f2_VS, color = decay)) + geom_smooth() + facet_wrap(~rewFunc)
d3 <- ggplot(df,aes(run_trial, db_f3_occ_parietal, color = decay)) + geom_smooth() + facet_wrap(~rewFunc)
d4 <- ggplot(df,aes(run_trial, db_f4_ACC_ins, color = decay)) + geom_smooth() + facet_wrap(~rewFunc)

pdf("h_bs_timecourse_by_condition_gamma.pdf", width = 12, height = 10)
ggarrange(h,h1,h2, ncol = 3)
dev.off()
pdf("v_bs_timecourse_by_condition_gamma.pdf", width = 20, height = 10)
ggarrange(v,v1,v2,v3,v4,v5,ncol = 6)
dev.off()
pdf("d_bs_timecourse_by_condition_gamma.pdf", width = 16, height = 10)
ggarrange(d, d1,d2,d3,d4, ncol = 5)
dev.off()
h <-  ggplot(df,aes(run_trial, v_entropy, color = performance)) + geom_smooth() + facet_wrap(~rewFunc)
v <-  ggplot(df,aes(run_trial, v_max, color = performance)) + geom_smooth() + facet_wrap(~rewFunc)
d <-  ggplot(df,aes(run_trial, d_auc, color = performance)) + geom_smooth() + facet_wrap(~rewFunc)

h1 <- ggplot(df,aes(run_trial, hb_f1_DAN_vlPFC, color = performance)) + geom_smooth() + facet_wrap(~rewFunc)
h2 <- ggplot(df,aes(run_trial, hb_f2_neg_paralimb, color = performance)) + geom_smooth() + facet_wrap(~rewFunc)
v1 <- ggplot(df,aes(run_trial, vb_f1_lo_DAN, color = performance)) + geom_smooth() + facet_wrap(~rewFunc)
v2 <- ggplot(df,aes(run_trial, vb_f2_hi_vmPFC_cOFC, color = performance)) + geom_smooth() + facet_wrap(~rewFunc)
v3 <- ggplot(df,aes(run_trial, vb_f3_lo_ACC, color = performance)) + geom_smooth() + facet_wrap(~rewFunc)
v4 <- ggplot(df,aes(run_trial, vb_f4_lo_cerebell_crus, color = performance)) + geom_smooth() + facet_wrap(~rewFunc)
v5 <- ggplot(df,aes(run_trial, vb_f5_hi_blITG, color = performance)) + geom_smooth() + facet_wrap(~rewFunc)
d1 <- ggplot(df,aes(run_trial, db_f1_rIFG_rSMA, color = performance)) + geom_smooth() + facet_wrap(~rewFunc)
d2 <- ggplot(df,aes(run_trial, db_f2_VS, color = performance)) + geom_smooth() + facet_wrap(~rewFunc)
d3 <- ggplot(df,aes(run_trial, db_f3_occ_parietal, color = performance)) + geom_smooth() + facet_wrap(~rewFunc)
d4 <- ggplot(df,aes(run_trial, db_f4_ACC_ins, color = performance)) + geom_smooth() + facet_wrap(~rewFunc)
pdf("h_bs_timecourse_by_condition_perf.pdf", width = 12, height = 10)
ggarrange(h,h1,h2, ncol = 3)
dev.off()
pdf("v_bs_timecourse_by_condition_perf.pdf", width = 20, height = 10)
ggarrange(v,v1,v2,v3,v4,v5,ncol = 6)
dev.off()
pdf("d_bs_timecourse_by_condition_perf.pdf", width = 16, height = 10)
ggarrange(d, d1,d2,d3,d4, ncol = 5)
dev.off()




pdf("all_bs_timecourse_by_condition.pdf", width = 20, height = 20)
ggarrange(v1,v2,v3,v4,v5,d1,d2,d3,d4,h1,h2, ncol = 5, nrow = 3)
dev.off()
