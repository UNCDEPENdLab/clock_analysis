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

ggplot(dfc[dfc$run>1,],aes(run_trial,hb_f1_DAN_vlPFC, color = rewFunc, lty = performance)) + geom_smooth(method = "loess") + facet_wrap(~rewFunc)
ggplot(dfc[dfc$run>1,],aes(run_trial,hb_f2_neg_paralimb, color = rewFunc, lty = performance)) + geom_smooth(method = "loess") + facet_wrap(~rewFunc)

ggplot(dfc[dfc$run>1,],aes(trial,hb_f1_DAN_vlPFC, lty = performance)) + geom_smooth(method = "loess") 


h <-  ggplot(dfc[dfc$run>1,],aes(run_trial, v_entropy, color = decay)) + geom_smooth() + facet_wrap(~rewFunc) + scale_x_continuous(breaks = c(1,50))
hneg <-  ggplot(dfc[dfc$run>1,],aes(run_trial, -v_entropy, color = decay)) + geom_smooth() + facet_wrap(~rewFunc)+ scale_x_continuous(breaks = c(1,50)) + guides(color = F)
v <-  ggplot(dfc[dfc$run>1,],aes(run_trial, v_max, color = decay)) + geom_smooth() + facet_wrap(~rewFunc)+ scale_x_continuous(breaks = c(1,50))+ guides(color = F)
vneg <-  ggplot(dfc[dfc$run>1,],aes(run_trial, -v_max, color = decay)) + geom_smooth() + facet_wrap(~rewFunc)+ scale_x_continuous(breaks = c(1,50))+ guides(color = F)
d <-  ggplot(dfc[dfc$run>1,],aes(run_trial, d_auc, color = decay)) + geom_smooth() + facet_wrap(~rewFunc)+ scale_x_continuous(breaks = c(1,50))+ guides(color = F)

h1 <- ggplot(dfc[dfc$run>1,],aes(run_trial, hb_f1_DAN_vlPFC, color = decay)) + geom_smooth() + facet_wrap(~rewFunc)+ scale_x_continuous(breaks = c(1,50))+ guides(color = F)
h2 <- ggplot(dfc[dfc$run>1,],aes(run_trial, hb_f2_neg_paralimb, color = decay)) + geom_smooth() + facet_wrap(~rewFunc)+ scale_x_continuous(breaks = c(1,50))+ guides(color = F)
v1 <- ggplot(dfc[dfc$run>1,],aes(run_trial, vb_f1_lo_DAN, color = decay)) + geom_smooth() + facet_wrap(~rewFunc)+ scale_x_continuous(breaks = c(1,50))+ guides(color = F)
v2 <- ggplot(dfc[dfc$run>1,],aes(run_trial, vb_f2_hi_vmPFC_cOFC, color = decay)) + geom_smooth() + facet_wrap(~rewFunc)+ scale_x_continuous(breaks = c(1,50))+ guides(color = F)
d1 <- ggplot(dfc[dfc$run>1,],aes(run_trial, db_f1_rIFG_rSMA, color = decay)) + geom_smooth() + facet_wrap(~rewFunc)+ scale_x_continuous(breaks = c(1,50))+ guides(color = F)
d4 <- ggplot(dfc[dfc$run>1,],aes(run_trial, db_f4_ACC_ins, color = decay)) + geom_smooth() + facet_wrap(~rewFunc)+ scale_x_continuous(breaks = c(1,50))+ guides(color = F)
pdf("bs_timecourse_by_condition_gamma.pdf", width = 12, height = 10)
ggpubr::ggarrange(h,vneg,h1,v,hneg,h2,d,h, d1,d4, ncol = 3, nrow = 3)
dev.off()

# not entirely convinced that we need to rescale within-subject, but let's move on with analyses

h <-  ggplot(dfc[dfc$run>1,],aes(run_trial, v_entropy, color = performance)) + geom_smooth(method = "loess") + facet_wrap(~rewFunc) + scale_x_continuous(breaks = c(1,50))
hneg <-  ggplot(dfc[dfc$run>1,],aes(run_trial, -v_entropy, color = performance)) + geom_smooth(method = "loess") + facet_wrap(~rewFunc)+ scale_x_continuous(breaks = c(1,50)) + guides(color = F)
v <-  ggplot(dfc[dfc$run>1,],aes(run_trial, v_max, color = performance)) + geom_smooth(method = "loess") + facet_wrap(~rewFunc)+ scale_x_continuous(breaks = c(1,50))+ guides(color = F)
vneg <-  ggplot(dfc[dfc$run>1,],aes(run_trial, -v_max, color = performance)) + geom_smooth(method = "loess") + facet_wrap(~rewFunc)+ scale_x_continuous(breaks = c(1,50))+ guides(color = F)
d <-  ggplot(dfc[dfc$run>1,],aes(run_trial, d_auc, color = performance)) + geom_smooth(method = "loess") + facet_wrap(~rewFunc)+ scale_x_continuous(breaks = c(1,50))+ guides(color = F)

h1 <- ggplot(dfc[dfc$run>1,],aes(run_trial, hb_f1_DAN_vlPFC, color = performance)) + geom_smooth(method = "loess") + facet_wrap(~rewFunc)+ scale_x_continuous(breaks = c(1,50))+ guides(color = F)
h2 <- ggplot(dfc[dfc$run>1,],aes(run_trial, hb_f2_neg_paralimb, color = performance)) + geom_smooth(method = "loess") + facet_wrap(~rewFunc)+ scale_x_continuous(breaks = c(1,50))+ guides(color = F)
v1 <- ggplot(dfc[dfc$run>1,],aes(run_trial, vb_f1_lo_DAN, color = performance)) + geom_smooth(method = "loess") + facet_wrap(~rewFunc)+ scale_x_continuous(breaks = c(1,50))+ guides(color = F)
v2 <- ggplot(dfc[dfc$run>1,],aes(run_trial, vb_f2_hi_vmPFC_cOFC, color = performance)) + geom_smooth(method = "loess") + facet_wrap(~rewFunc)+ scale_x_continuous(breaks = c(1,50))+ guides(color = F)
d1 <- ggplot(dfc[dfc$run>1,],aes(run_trial, db_f1_rIFG_rSMA, color = performance)) + geom_smooth(method = "loess") + facet_wrap(~rewFunc)+ scale_x_continuous(breaks = c(1,50))+ guides(color = F)
d4 <- ggplot(dfc[dfc$run>1,],aes(run_trial, db_f4_ACC_ins, color = performance)) + geom_smooth(method = "loess") + facet_wrap(~rewFunc)+ scale_x_continuous(breaks = c(1,50))+ guides(color = F)
pdf("bs_timecourse_by_condition_performance.pdf", width = 12, height = 10)
ggarrange(h,vneg,h1,v,hneg,h2,d,d1,d4, ncol = 3, nrow = 3)
dev.off()
# d4 looks more like paralimbic in good subjects and more like DAN in bad
# review bs by ID
dflh <- gather(dfc, network, signal, hb_f1_DAN_vlPFC:hb_f2_neg_paralimb, factor_key = T)
pdf("h1_bs_timecourse_by_condition_subject.pdf", width = 20, height = 20)
ggplot(dflh,aes(run_trial, signal, color = network)) + geom_line() + facet_wrap(ID~rewFunc)+ scale_x_continuous(breaks = c(1,50))
dev.off()


# pdf("all_bs_timecourse_by_condition.pdf", width = 20, height = 20)
# ggarrange(v1,v2,v3,v4,v5,d1,d2,d3,d4,h1,h2, ncol = 5, nrow = 3)
# dev.off()
