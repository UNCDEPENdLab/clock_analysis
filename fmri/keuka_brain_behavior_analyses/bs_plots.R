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
setwd('~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/plots')

######
# end of preprocessing
#####
## "model-free analyses"

###########
## PLOTS ##
#####

dfc <- na.omit(df[df$rt_swing>0,])


# for further analyses, we want the following regions:
# hb_f1_DAN (alt. vb_f1_lo_DAN)
# hb_f2_paralimbic (alt vb_f2_hi_paralimbic)
# db_f4_ACC_ins -- compression
# db_f1_rIFG_SMA -- ?overload

# sanity check spaghetti
pdf("sanity_check_pe_spaghetti.pdf", height = 20, width = 20)
ggplot(dfc,aes(run_trial,peb_f2_p_hipp, color = rewFunc)) + geom_smooth(method = 'glm') + facet_grid(~id)
dev.off()

pdf("sanity_check_pe.pdf", height = 10, width = 10)
ggplot(dfc,aes(run_trial,peb_f2_p_hipp, lty = performance, color = rewFunc)) + geom_smooth(method = 'glm') 
dev.off()


# compare ant. vs. post hippocampus
rts <- ggplot(dfc[dfc$run>1,],aes(run_trial, rt_csv, color = performance)) + geom_smooth(method = 'gam',  formula = y ~ s(x, bs = "ad")) + facet_wrap(~rewFunc)+ scale_x_continuous(breaks = c(1,50))#+ guides(color = F)
swings <- ggplot(dfc[dfc$run>1,],aes(run_trial, rt_swing, color = performance)) + geom_smooth(method = 'gam',  formula = y ~ s(x, bs = "ad")) + facet_wrap(~rewFunc)+ scale_x_continuous(breaks = c(1,50))#+ guides(color = F)
ah <- ggplot(dfc[dfc$run>1,],aes(run_trial, h_ant_hipp_b_f, color = performance)) + geom_smooth(method = 'gam',  formula = y ~ s(x, bs = "ad")) + facet_wrap(~rewFunc)+ scale_x_continuous(breaks = c(1,50))#+ guides(color = F)
ph <- ggplot(dfc[dfc$run>1,],aes(run_trial, peb_f2_p_hipp, color = performance)) + geom_smooth(method = 'gam',  formula = y ~ s(x, bs = "ad")) + facet_wrap(~rewFunc)+ scale_x_continuous(breaks = c(1,50))#+ guides(color = F)
h <-  ggplot(dfc[dfc$run>1,],aes(run_trial, v_entropy, color = performance)) + geom_smooth(method = 'gam',  formula = y ~ s(x, bs = "ad")) + facet_wrap(~rewFunc) + scale_x_continuous(breaks = c(1,50))
peabs <-  ggplot(dfc[dfc$run>1,],aes(run_trial, abs(pe_max), color = performance)) + geom_smooth(method = 'gam',  formula = y ~ s(x, bs = "ad")) + facet_wrap(~rewFunc) + scale_x_continuous(breaks = c(1,50))
pe <-  ggplot(dfc[dfc$run>1,],aes(run_trial, pe_max, color = performance)) + geom_smooth(method = 'gam',  formula = y ~ s(x, bs = "ad")) + facet_wrap(~rewFunc) + scale_x_continuous(breaks = c(1,50))

pdf("hipp_bs_beh_timecourse_by_condition_gamma.pdf", width = 12, height = 10)
ggpubr::ggarrange(rts, swings, h, ah,ph, pe, ncol = 3, nrow = 2)
dev.off()

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
hneg <-  ggplot(dfc[dfc$run>1,],aes(run_trial, -v_entropy, color = performance)) + geom_smooth() + facet_wrap(~rewFunc)+ scale_x_continuous(breaks = c(1,50))
v <-  ggplot(dfc[dfc$run>1,],aes(run_trial, v_max, color = performance)) + geom_smooth() + facet_wrap(~rewFunc)+ scale_x_continuous(breaks = c(1,50))
vneg <-  ggplot(dfc[dfc$run>1,],aes(run_trial, -v_max, color = performance)) + geom_smooth() + facet_wrap(~rewFunc)+ scale_x_continuous(breaks = c(1,50))
d <-  ggplot(dfc[dfc$run>1,],aes(run_trial, d_auc, color = performance)) + geom_smooth() + facet_wrap(~rewFunc)+ scale_x_continuous(breaks = c(1,50))

h1 <- ggplot(dfc[dfc$run>1,],aes(run_trial, hb_f1_DAN_vlPFC, color = performance)) + geom_smooth(method = "loess") + facet_wrap(~rewFunc)+ scale_x_continuous(breaks = c(1,50))
h2 <- ggplot(dfc[dfc$run>1,],aes(run_trial, hb_f2_neg_paralimb, color = performance)) + geom_smooth(method = "loess") + facet_wrap(~rewFunc)+ scale_x_continuous(breaks = c(1,50))
v1 <- ggplot(dfc[dfc$run>1,],aes(run_trial, vb_f1_lo_DAN, color = performance)) + geom_smooth() + facet_wrap(~rewFunc)+ scale_x_continuous(breaks = c(1,50))
v2 <- ggplot(dfc[dfc$run>1,],aes(run_trial, vb_f2_hi_vmPFC_cOFC, color = performance)) + geom_smooth() + facet_wrap(~rewFunc)+ scale_x_continuous(breaks = c(1,50))
d1 <- ggplot(dfc[dfc$run>1,],aes(run_trial, db_f1_rIFG_rSMA, color = performance)) + geom_smooth() + facet_wrap(~rewFunc)+ scale_x_continuous(breaks = c(1,50))
d4 <- ggplot(dfc[dfc$run>1,],aes(run_trial, db_f4_ACC_ins, color = performance)) + geom_smooth() + facet_wrap(~rewFunc)+ scale_x_continuous(breaks = c(1,50))
pdf("bs_timecourse_by_condition_performance.pdf", width = 12, height = 10)
ggarrange(h,vneg,h1,v1,v,hneg,h2,v2,d,h, d1,d4, ncol = 4, nrow = 3)
dev.off()

pdf("h_bs_timecourse_by_condition_performance.pdf", width = 12, height = 10)
ggarrange(h,h1,h2, ncol = 2, nrow = 2)
dev.off()



# pdf("all_bs_timecourse_by_condition.pdf", width = 20, height = 20)
# ggarrange(v1,v2,v3,v4,v5,d1,d2,d3,d4,h1,h2, ncol = 5, nrow = 3)
# dev.off()
