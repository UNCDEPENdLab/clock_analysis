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
## "model-free analyses"

###########
## PLOTS ##
#####

# for further analyses, we want the following regions:
# hb_f1_DAN (alt. vb_f1_lo_DAN)
# hb_f2_paralimbic (alt vb_f2_hi_paralimbic)
# db_f4_ACC_ins -- compression
# db_f1_rIFG_SMA -- ?overload


h <-  ggplot(df,aes(run_trial, v_entropy, color = decay)) + geom_smooth() + facet_wrap(~rewFunc) + scale_x_continuous(breaks = c(1,50))
hneg <-  ggplot(df,aes(run_trial, -v_entropy, color = decay)) + geom_smooth() + facet_wrap(~rewFunc)+ scale_x_continuous(breaks = c(1,50)) + guides(color = F)
v <-  ggplot(df,aes(run_trial, v_max, color = decay)) + geom_smooth() + facet_wrap(~rewFunc)+ scale_x_continuous(breaks = c(1,50))+ guides(color = F)
vneg <-  ggplot(df,aes(run_trial, -v_max, color = decay)) + geom_smooth() + facet_wrap(~rewFunc)+ scale_x_continuous(breaks = c(1,50))+ guides(color = F)
d <-  ggplot(df,aes(run_trial, d_auc, color = decay)) + geom_smooth() + facet_wrap(~rewFunc)+ scale_x_continuous(breaks = c(1,50))+ guides(color = F)

h1 <- ggplot(df,aes(run_trial, hb_f1_DAN_vlPFC, color = decay)) + geom_smooth() + facet_wrap(~rewFunc)+ scale_x_continuous(breaks = c(1,50))+ guides(color = F)
h2 <- ggplot(df,aes(run_trial, hb_f2_neg_paralimb, color = decay)) + geom_smooth() + facet_wrap(~rewFunc)+ scale_x_continuous(breaks = c(1,50))+ guides(color = F)
v1 <- ggplot(df,aes(run_trial, vb_f1_lo_DAN, color = decay)) + geom_smooth() + facet_wrap(~rewFunc)+ scale_x_continuous(breaks = c(1,50))+ guides(color = F)
v2 <- ggplot(df,aes(run_trial, vb_f2_hi_vmPFC_cOFC, color = decay)) + geom_smooth() + facet_wrap(~rewFunc)+ scale_x_continuous(breaks = c(1,50))+ guides(color = F)
d1 <- ggplot(df,aes(run_trial, db_f1_rIFG_rSMA, color = decay)) + geom_smooth() + facet_wrap(~rewFunc)+ scale_x_continuous(breaks = c(1,50))+ guides(color = F)
d4 <- ggplot(df,aes(run_trial, db_f4_ACC_ins, color = decay)) + geom_smooth() + facet_wrap(~rewFunc)+ scale_x_continuous(breaks = c(1,50))+ guides(color = F)
pdf("bs_timecourse_by_condition_gamma.pdf", width = 12, height = 10)
ggpubr::ggarrange(h,vneg,h1,v,hneg,h2,d,h, d1,d4, ncol = 3, nrow = 3)
dev.off()

# not entirely convinced that we need to rescale within-subject, but let's move on with analyses

h <-  ggplot(df,aes(run_trial, v_entropy, color = performance)) + geom_smooth() + facet_wrap(~rewFunc) + scale_x_continuous(breaks = c(1,50))
hneg <-  ggplot(df,aes(run_trial, -v_entropy, color = performance)) + geom_smooth() + facet_wrap(~rewFunc)+ scale_x_continuous(breaks = c(1,50))
v <-  ggplot(df,aes(run_trial, v_max, color = performance)) + geom_smooth() + facet_wrap(~rewFunc)+ scale_x_continuous(breaks = c(1,50))
vneg <-  ggplot(df,aes(run_trial, -v_max, color = performance)) + geom_smooth() + facet_wrap(~rewFunc)+ scale_x_continuous(breaks = c(1,50))
d <-  ggplot(df,aes(run_trial, d_auc, color = performance)) + geom_smooth() + facet_wrap(~rewFunc)+ scale_x_continuous(breaks = c(1,50))

h1 <- ggplot(df,aes(run_trial, hb_f1_DAN_vlPFC, color = performance)) + geom_smooth() + facet_wrap(~rewFunc)+ scale_x_continuous(breaks = c(1,50))
h2 <- ggplot(df,aes(run_trial, hb_f2_neg_paralimb, color = performance)) + geom_smooth() + facet_wrap(~rewFunc)+ scale_x_continuous(breaks = c(1,50))
v1 <- ggplot(df,aes(run_trial, vb_f1_lo_DAN, color = performance)) + geom_smooth() + facet_wrap(~rewFunc)+ scale_x_continuous(breaks = c(1,50))
v2 <- ggplot(df,aes(run_trial, vb_f2_hi_vmPFC_cOFC, color = performance)) + geom_smooth() + facet_wrap(~rewFunc)+ scale_x_continuous(breaks = c(1,50))
d1 <- ggplot(df,aes(run_trial, db_f1_rIFG_rSMA, color = performance)) + geom_smooth() + facet_wrap(~rewFunc)+ scale_x_continuous(breaks = c(1,50))
d4 <- ggplot(df,aes(run_trial, db_f4_ACC_ins, color = performance)) + geom_smooth() + facet_wrap(~rewFunc)+ scale_x_continuous(breaks = c(1,50))
pdf("bs_timecourse_by_condition_performance.pdf", width = 12, height = 10)
ggarrange(h,vneg,h1,v1,v,hneg,h2,v2,d,h, d1,d4, ncol = 4, nrow = 3)
dev.off()
'

# pdf("all_bs_timecourse_by_condition.pdf", width = 20, height = 20)
# ggarrange(v1,v2,v3,v4,v5,d1,d2,d3,d4,h1,h2, ncol = 5, nrow = 3)
# dev.off()
