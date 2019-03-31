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
# load('trial_df_and_vhd_bs.Rdata')
load('clusters_and_beta_series.Rdata')
setwd('~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/plots')

######
# end of preprocessing
#####
## "model-free analyses"

###########
## PLOTS ##
#####

dfc <- (df)


# for further analyses, we want the following regions:
# hb_f1_DAN (alt. vb_f1_lo_DAN)
# hb_f2_paralimbic (alt vb_f2_hi_paralimbic)
# db_f4_ACC_ins -- compression
# db_f1_rIFG_SMA -- ?overload

# sanity check spaghetti
pdf("sanity_check_ph_spaghetti.pdf", height = 20, width = 20)
ggplot(dfc,aes(run_trial,peb_f2_p_hipp, color = rewFunc)) + geom_smooth(method = 'glm') + facet_grid(~id)
dev.off()

pdf("sanity_check_ph.pdf", height = 10, width = 10)
ggplot(dfc,aes(run_trial,peb_f2_p_hipp, lty = performance, color = rewFunc)) + geom_smooth(method = 'gam', formula = y ~ s(x, bs = "ad")) + facet_grid(~rewFunc)
dev.off()

ph <- ggplot(dfc,aes(run_trial,peb_f2_p_hipp, lty = performance, color = rewFunc)) + geom_smooth(method = 'gam', formula = y ~ s(x, bs = "ad"), formula = y ~ s(x, bs = "ad")) + facet_grid(~rewFunc)
ah <- ggplot(dfc,aes(run_trial,h_ant_hipp_b_f, lty = performance, color = rewFunc)) + geom_smooth(method = 'gam', formula = y ~ s(x, bs = "ad")) + facet_grid(~rewFunc)
rts <- ggplot(dfc[dfc$run>1,],aes(run_trial, rt_csv, lty = performance, color = rewFunc)) + geom_smooth(method = 'gam', formula = y ~ s(x, bs = "ad")) + facet_grid(~rewFunc)
swings <- ggplot(dfc[dfc$run>1,],aes(run_trial, rt_swing, lty = performance, color = rewFunc)) + geom_smooth(method = 'gam', formula = y ~ s(x, bs = "ad")) + facet_grid(~rewFunc)
h <-  ggplot(dfc[dfc$run>1,],aes(run_trial, v_entropy, lty = performance, color = rewFunc)) + geom_smooth(method = 'gam', formula = y ~ s(x, bs = "ad")) + facet_grid(~rewFunc)
peabs <-  ggplot(dfc[dfc$run>1,],aes(run_trial, abs(pe_max), lty = performance, color = rewFunc)) + geom_smooth(method = 'gam', formula = y ~ s(x, bs = "ad")) + facet_grid(~rewFunc)
pe <-  ggplot(dfc[dfc$run>1,],aes(run_trial, abs(pe_max), lty = performance, color = rewFunc)) + geom_smooth(method = 'gam', formula = y ~ s(x, bs = "ad")) + facet_grid(~rewFunc)

pdf("sanity_check_aph.pdf", height = 10, width = 20)
ggarrange(ah,ph,nrow = 2)
dev.off()

# individual time courses for hippo bs, censored
pdf("ph_ind_timecourse_by_condition.pdf", height = 20, width = 20)
ggplot(df[abs(df$peb_f2_p_hipp)<1,],aes(run_trial,peb_f2_p_hipp, color = rewFunc)) + geom_smooth() + facet_wrap(~id)
dev.off()

pt <- ggplot(df[abs(df$peb_f2_p_hipp)<2,],aes(run_trial,peb_f2_p_hipp/rt_csv, color = rewFunc)) + geom_smooth()
at <- ggplot(df[abs(df$h_ant_hipp_b_f)<2,],aes(run_trial,h_ant_hipp_b_f/rt_csv, color = rewFunc)) + geom_smooth()
pdf("h_timecourse_by_condition.pdf", height = 6, width = 12)
ggarrange(at,pt, ncol = 2, nrow = 1)
dev.off()


pdf("ah_ind_timecourse_by_condition.pdf", height = 20, width = 20)
ggplot(df,aes(run_trial,h_ant_hipp_b_f, color = rewFunc)) + geom_smooth() + facet_wrap(~id)
dev.off()




pdf("rts_signals_aph.pdf", height = 20, width = 20)
ggarrange(rts, swings, h, pe, ah,ph,nrow = 6)
dev.off()


pdf("ph_bs_by_pe_resp.pdf", height = 10, width = 10)
ggplot(dfc,aes(run_trial,peb_f2_p_hipp, lty = pe_f2_hipp_resp, color = rewFunc)) + geom_smooth(method = 'gam', formula = y ~ s(x, bs = "ad")) + facet_grid(~rewFunc)
dev.off()


pdf("sanity_check_ah.pdf", height = 10, width = 10)
ggplot(dfc,aes(run_trial,h_ant_hipp_b_f, lty = performance, color = rewFunc)) + geom_smooth(method = 'gam', formula = y ~ s(x, bs = "ad")) + facet_grid(~rewFunc)
dev.off()

pdf("ah_bs_by_h_resp.pdf", height = 10, width = 10)
ggplot(dfc,aes(run_trial,h_ant_hipp_b_f, lty = h_HippAntL_resp, color = rewFunc)) + geom_smooth(method = 'gam', formula = y ~ s(x, bs = "ad")) + facet_grid(~rewFunc)
dev.off()

# how does Hipp scale with PEs?
p2 <- ggplot(dfc,aes(pe_max_lag,peb_f2_p_hipp, color = rewFunc)) + geom_smooth(method = 'loess') + facet_grid(~rewFunc)
p1 <- ggplot(dfc,aes(pe_max_lag,peb_f1_cort_str, color = rewFunc)) + geom_smooth(method = 'loess') + facet_grid(~rewFunc)

pdf("peb_vs_pe.pdf", height = 10, width = 20)
ggarrange(p1,p2,ncol = 2,nrow = 1)
dev.off()


# plot the between-regions interaction on RTs
pdf("ant_by_post_on_rt.pdf", width = 10, height = 10)
ggplot(df[!is.na(df$h_ant_hipp_b_f) & !is.na(df$peb_f2_p_hipp),], aes(rt_lag, rt_csv, lty = h_ant_hipp_b_f>0, color = peb_f2_p_hipp>0)) + 
  # geom_smooth(method = 'gam', formula = y ~ s(x, bs = "tp")) + facet_grid(~rewFunc)
  geom_smooth(method = 'loess') + facet_grid(~rewFunc)

dev.off()

pdf("ant_by_post_rt_timecourse.pdf", width = 10, height = 10)
ggplot(df[!is.na(df$h_ant_hipp_b_f) & !is.na(df$peb_f2_p_hipp),], aes(run_trial, rt_csv, lty = h_ant_hipp_b_f>0, color = peb_f2_p_hipp>0)) + 
  # geom_smooth(method = 'gam', formula = y ~ s(x, bs = "tp")) + facet_grid(~rewFunc)
  geom_smooth(method = 'loess') + facet_grid(~rewFunc)
dev.off()

# looks like anterior decreases and posterior increases RT swings
as <- ggplot(df[!is.na(df$h_ant_hipp_b_f) & !is.na(df$peb_f2_p_hipp),], aes(h_ant_hipp_b_f, rt_swing)) + geom_smooth(method = 'gam') 
ps <- ggplot(df[!is.na(df$h_ant_hipp_b_f) & !is.na(df$peb_f2_p_hipp),], aes(peb_f2_p_hipp, rt_swing)) + geom_smooth(method = 'gam') 
pdf("ant_post_rt_swing.pdf", width = 16, height = 10)
ggpubr::ggarrange(as,ps, ncol = 2)
dev.off()

# check within subjects -- not day and night...
asi <- ggplot(df[!is.na(df$h_ant_hipp_b_f) & !is.na(df$peb_f2_p_hipp),], aes(h_ant_hipp_b_f, rt_swing)) + geom_point() + facet_wrap(~id)
psi <- ggplot(df[!is.na(df$h_ant_hipp_b_f) & !is.na(df$peb_f2_p_hipp),], aes(peb_f2_p_hipp, rt_swing)) + geom_point() + facet_wrap(~id)
pdf("ant_post_rtswings_subjects.pdf", width = 16, height = 10)
ggpubr::ggarrange(asi,psi, ncol = 2)
dev.off()

# subtle, but anterior causes convergence and posterior, exploration
am <- ggplot(df[!is.na(df$h_ant_hipp_b_f) & !is.na(df$peb_f2_p_hipp),], aes(rt_vmax,rt_csv, color = h_ant_hipp_b_f>0)) + geom_smooth(method = 'gam') 
pm <- ggplot(df[!is.na(df$h_ant_hipp_b_f) & !is.na(df$peb_f2_p_hipp),], aes(rt_vmax,rt_csv, color = peb_f2_p_hipp>0)) + geom_smooth(method = 'gam') 
pdf("ant_post_rt_vmax.pdf", width = 16, height = 10)
ggpubr::ggarrange(am,pm, ncol = 2)
dev.off()

# RT convolution artifact on BS
pdf("rt_hipp.pdf", width = 16, height = 10)
ggplot(df[!is.na(df$h_ant_hipp_b_f) & !is.na(df$peb_f2_p_hipp),], aes(peb_f2_p_hipp, rt_csv, color = h_ant_hipp_b_f>0)) + geom_smooth(method = 'loess') + facet_wrap(~rewFunc)
dev.off()

# what about next RT?
pdf("rtnext_hipp_ind.pdf", width = 20, height = 20)
ggplot(df[!is.na(df$h_ant_hipp_b_f) & !is.na(df$peb_f2_p_hipp),], aes(peb_f2_p_hipp, rt_next, color = h_ant_hipp_b_f>0)) + geom_smooth(method = 'loess') + facet_wrap(~id)
dev.off()


# further investigation of RT convolution effect on beta series
rt_bs <- ggplot(df[!is.na(df$h_ant_hipp_b_f) & !is.na(df$peb_f2_p_hipp),], aes(rt_csv, peb_f2_p_hipp)) + geom_smooth() + facet_wrap(~rewFunc)
pdf("post_rt.pdf", width = 16, height = 10)
rt_bs
dev.off()

ph <- ggplot(df[!is.na(df$h_ant_hipp_b_f) & !is.na(df$peb_f2_p_hipp),], aes(peb_f2_p_hipp,rt_csv, color = rewFunc)) + geom_smooth()
ah <- ggplot(df[!is.na(df$h_ant_hipp_b_f) & !is.na(df$peb_f2_p_hipp),], aes(h_ant_hipp_b_f,rt_csv, color = rewFunc)) + geom_smooth()
dan <- ggplot(df[!is.na(df$h_ant_hipp_b_f) & !is.na(df$peb_f2_p_hipp),], aes(hb_f1_DAN_vlPFC,rt_csv, color = rewFunc)) + geom_smooth()
pe1 <- ggplot(df[!is.na(df$h_ant_hipp_b_f) & !is.na(df$peb_f2_p_hipp),], aes(peb_f1_cort_str,rt_csv, color = rewFunc)) + geom_smooth()

pdf("rt_bs.pdf", width = 20, height = 20)
ggpubr::ggarrange(ah, dan, ph, pe1, ncol = 2, nrow = 2)
dev.off()


# compare ant. vs. post hippocampus
rts <- ggplot(dfc[dfc$run>1,],aes(run_trial, rt_csv, color = performance)) + geom_smooth(method = 'gam', formula = y ~ s(x, bs = "ad")) + facet_wrap(~rewFunc)+ scale_x_continuous(breaks = c(1,50))#+ guides(color = F)
swings <- ggplot(dfc[dfc$run>1,],aes(run_trial, rt_swing, color = performance)) + geom_smooth(method = 'gam', formula = y ~ s(x, bs = "ad")) + facet_wrap(~rewFunc)+ scale_x_continuous(breaks = c(1,50))#+ guides(color = F)
ah <- ggplot(dfc[dfc$run>1,],aes(run_trial, h_ant_hipp_b_f, color = performance)) + geom_smooth(method = 'gam', formula = y ~ s(x, bs = "ad")) + facet_wrap(~rewFunc)+ scale_x_continuous(breaks = c(1,50))#+ guides(color = F)
ph <- ggplot(dfc[dfc$run>1,],aes(run_trial, peb_f2_p_hipp, color = performance)) + geom_smooth(method = 'gam', formula = y ~ s(x, bs = "ad")) + facet_wrap(~rewFunc)+ scale_x_continuous(breaks = c(1,50))#+ guides(color = F)
h <-  ggplot(dfc[dfc$run>1,],aes(run_trial, v_entropy, color = performance)) + geom_smooth(method = 'gam', formula = y ~ s(x, bs = "ad")) + facet_wrap(~rewFunc)+ scale_x_continuous(breaks = c(1,50))#+ guides(color = F)
peabs <-  ggplot(dfc[dfc$run>1,],aes(run_trial, abs(pe_max), color = performance)) + geom_smooth(method = 'gam', formula = y ~ s(x, bs = "ad")) + facet_wrap(~rewFunc)+ scale_x_continuous(breaks = c(1,50))#+ guides(color = F)
pe <-  ggplot(dfc[dfc$run>1,],aes(run_trial, pe_max, color = performance)) + geom_smooth(method = 'gam', formula = y ~ s(x, bs = "ad")) + facet_wrap(~rewFunc)+ scale_x_continuous(breaks = c(1,50))#+ guides(color = F)

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
