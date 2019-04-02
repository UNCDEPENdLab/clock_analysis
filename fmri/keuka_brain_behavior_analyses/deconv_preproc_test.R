library(dplyr)
library(tidyverse)
library(psych)
library(corrplot)
library(lme4)

setwd('~/Box Sync/skinner/projects_analyses/SCEPTIC/fMRI_paper/signals_review/MMClock_aroma_preconvolve_fse_groupfixed/')
d <- read_csv("10637_run4_long_axis_l_2.3mm.nii.gz_deconvolved.csv")
pdf('decon_sample.pdf', width = 20, height = 10)
ggplot(d,aes(time, decon, color = atlas_value>.5)) + geom_smooth()
dev.off()
# get trial_level data
trial_df <- read_csv("~/code/clock_analysis/fmri/data/mmclock_fmri_decay_factorize_selective_psequate_mfx_trial_statistics.csv.gz")
trial_df <- trial_df %>%
  group_by(id, run) %>%  dplyr::mutate(rt_swing = abs(c(NA, diff(rt_csv))),
                                       rt_swing_lr = abs(log(rt_csv/lag(rt_csv))),
                                       rt_lag = lag(rt_csv) ,
                                       omission_lag = lag(score_csv==0),
                                       rt_vmax_lag = lag(rt_vmax),
                                       v_entropy_wi = scale(v_entropy),
                                       run_trial=1:50) %>% ungroup() #compute rt_swing within run and subject
# remove trial variable to avoid confusion
trial_df <- trial_df[,c(1:4,6:41)]
# performance
sum_df <- trial_df %>% group_by(id) %>% dplyr::summarize(total_earnings = sum(score_csv)) %>% arrange(total_earnings)
beta_fscores$id <- beta_fscores$ID
beta_fscores <- inner_join(beta_fscores,sum_df)

# model parameters
params <- read_csv("~/code/clock_analysis/fmri/data/mmclock_fmri_decay_factorize_selective_psequate_mfx_sceptic_global_statistics.csv")
bdf <- inner_join(beta_fscores,params)

# params_beta <- bdf[,c(names(dvh), "total_earnings", "LL", "alpha", "gamma", "beta")]
# 
# param_cor <- corr.test(params_beta,method = 'pearson', adjust = 'none')
# 
# setwd('~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/')
# pdf("dvh_bs_param_corr_fixed.pdf", width=12, height=12)
# corrplot(param_cor$r, cl.lim=c(-1,1),
#          method = "circle", tl.cex = 1.5, type = "upper", tl.col = 'black',
#          order = "hclust", diag = FALSE,
#          addCoef.col="black", addCoefasPercent = FALSE,
#          p.mat = param_cor$p, sig.level=0.05, insig = "blank")
# dev.off()
# # they don't correlate

# merge into trial-level data
bdf <- bdf[,c(names(beta_fscores), "LL", "alpha", "gamma", "beta")]

df <- inner_join(trial_df,bdf)

# # check correlation with RTs -- at least not uniform
# dvh_rt <- df[,c(names(dvh), "rt_csv")]
# rt_cor <- corr.test(dvh_rt,method = 'pearson', adjust = 'none')
# 
# setwd('~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/')
# pdf("dvh_bs_rt_corr.pdf", width=12, height=12)
# corrplot(rt_cor$r, cl.lim=c(-1,1),
#          method = "circle", tl.cex = 1.5, type = "upper", tl.col = 'black',
#          order = "hclust", diag = FALSE,
#          addCoef.col="black", addCoefasPercent = FALSE,
#          p.mat = rt_cor$p, sig.level=0.05, insig = "blank")
# dev.off()


df$rewFunc <- relevel(as.factor(df$rewFunc),ref = "CEV")
df$rewFuncIEVsum <- df$rewFunc
contrasts(df$rewFuncIEVsum) <- contr.sum
colnames(contrasts(df$rewFuncIEVsum)) <- c('CEV','CEVR', 'DEV')
# # dichotomize betas for plotting
# df$h_fp <- 'low'
# df$h_fp[df$h_f1_fp>0] <- 'high'
# df$low_h_paralimbic <- 'low'
# df$low_h_paralimbic[df$h_f2_neg_paralimb<0] <- 'high'
# 
# df$v_paralimbic <- 'low'
# df$v_paralimbic[df$v_f2_paralimb>0] <- 'high'
# #NB: v_f2_paralimb is mostly vlPFC driven now...
# df$low_v_fp_acc_vlpfc <- 'low'
# df$low_v_fp_acc_vlpfc[df$v_f1_neg_cog<0] <- 'high'

df$learning_epoch <- 'trials 1-10'
df$learning_epoch[df$run_trial>10] <- 'trials 11-50'




# obtain within-subject v_max and entropy: correlated at -.37
df <- df %>% group_by(id,run) %>% mutate(v_max_wi = scale(v_max),
                                         v_max_wi_lag = lag(v_max_wi),
                                         v_entropy_wi = scale(v_entropy),
                                         v_max_b = mean(na.omit(v_max)),
                                         v_entropy_b = mean(na.omit(v_entropy)),
                                         rt_change = rt_csv - rt_lag,
                                         pe_max_lag = lag(pe_max), 
                                         abs_pe_max_lag = abs(pe_max_lag), 
                                         rt_vmax_change = rt_vmax - rt_vmax_lag,
                                         # lag  bs
                                         peb_f1_cort_str_lag = lag(peb_f1_cort_str),
                                         peb_f2_p_hipp_lag = lag(peb_f2_p_hipp),
                                         h_ant_hipp_b_f_lag = lag(h_ant_hipp_b_f)
)


df$performance <- cut_number(df$total_earnings,2)
levels(df$performance) <- c( "below median", "above median")
#sanity check
# ggplot(df,aes(performance, total_earnings)) + geom_boxplot()
df$last_outcome <- NA
df$last_outcome[df$omission_lag] <- 'Omission'
df$last_outcome[!df$omission_lag] <- 'Reward'
# BS by trial and condition
df$decay <- NA
df$decay[df$gamma>0] <- 'high'
df$decay[df$gamma<0] <- 'low'

# dichotomize BS for plotting

save(file = 'trial_df_and_vhd_bs.Rdata', df)

