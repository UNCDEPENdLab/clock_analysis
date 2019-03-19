library(dplyr)
library(tidyverse)
library(psych)
library(corrplot)
library(lme4)

remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}

scale_within <- T

# get cluster meta-data, first for H
setwd('~/Box Sync/skinner/projects_analyses/SCEPTIC/fMRI_paper/signals_review/beta_series_dec2018/')
metaH <- read_csv("v_entropy_cluster_metadata.csv")
metaH$label <- substr(metaH$label,22,100)
meta_overall <- metaH[metaH$l2_contrast == 'overall' & metaH$l3_contrast == 'Intercept' & metaH$model == 'Intercept-Age',]
hbs <- read_csv("v_entropy_roi_beta_series.csv")
# filter out clusters <150 voxels
hb <- hbs[hbs$l2_contrast == 'overall' & hbs$l3_contrast == 'Intercept' & hbs$model == 'Intercept-Age',1:6] %>% filter(cluster_number<10)

# head(merge(h,meta))
rois <- distinct(meta_overall[,c(5,12)])

# inspect distributions
hb <- inner_join(hb,meta_overall)
hb$labeled_cluster <- paste(hb$cluster_number,hb$label)
# ggplot(hb,aes(scale(bs_value))) + geom_histogram() + facet_wrap(~label)
# outliers -- 2.5% looks best by eye
hb <- hb %>% mutate_at(.vars = 'bs_value', remove_outliers)
# ggplot(hb,aes(scale(bs_value))) + geom_histogram() + facet_wrap(~labeled_cluster)

ggplot(hb,aes(bs_value)) + geom_histogram() + facet_wrap(~labeled_cluster)

# scale each region within subject
if (scale_within)
{hb <- hb %>% dplyr::group_by(feat_input_id,labeled_cluster) %>% mutate(bs_value = scale(bs_value)) %>% ungroup}
ggplot(hb,aes(bs_value)) + geom_histogram() + facet_wrap(~labeled_cluster)

#just the bs and 

hb_wide <- spread(hb[,c("feat_input_id", "run", "trial", "bs_value", "labeled_cluster")],labeled_cluster,bs_value)
# head(h_wide)
# with group-fixed parameters we don't seem to need the winsorization step!!!
# h_wide <- h_wide %>% mutate_if(is.double, winsor,trim = .075)

just_bs <- hb_wide[,3:ncol(hb_wide)]
# winsorize to deal with beta ouliers

# non-parametric correlations to deal with outliers
bs_cor <- corr.test(just_bs,method = 'pearson')
# parametric correlations on winsorised betas
# clust_cor <- cor(just_rois_w,method = 'pearson')

setwd('~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/')
pdf("h_bs_wi_corr.pdf", width=12, height=12)
corrplot(bs_cor$r, cl.lim=c(-1,1),
         method = "circle", tl.cex = 1.5, type = "upper", tl.col = 'black',
         order = "hclust", diag = FALSE,
         addCoef.col="black", addCoefasPercent = FALSE,
         p.mat = bs_cor$p, sig.level=0.01, insig = "blank")
dev.off()
# remove trial
just_bs <- just_bs[,2:ncol(just_bs)]

# 2-factor solution has a left post-central loading on the paralimbic hb_f2
# may refine later, let's prototype
hbs.fa = psych::fa(just_bs, nfactors=2)
hfscores <- factor.scores(x = just_bs, hbs.fa)$scores
hb_wide$hb_f1_DAN_vlPFC <- hfscores[,1]
hb_wide$hb_f2_neg_paralimb <- hfscores[,2]


######## Vmax
# use v_max as theoretically interesting and unbiased by choice
setwd('~/Box Sync/skinner/projects_analyses/SCEPTIC/fMRI_paper/signals_review/beta_series_dec2018/')
metaV <- read_csv('v_max_cluster_metadata.csv')
metaV$label <- substr(metaV$label,22,100)
metaV_overall <- metaV[metaV$l2_contrast == 'overall' & metaV$l3_contrast == 'Intercept' & metaV$model == 'Intercept-Age',]
vbs <- read_csv("v_max_roi_beta_series.csv")
# filter out clusters <150 voxels
vb <- vbs[vbs$l2_contrast == 'overall' & vbs$l3_contrast == 'Intercept' & vbs$model == 'Intercept-Age',1:6] %>% filter(cluster_number<15)

# head(merge(h,meta))
vroi_list <- distinct(metaV_overall[,c(5,12)])

# inspect distributions
vb <- inner_join(vb,metaV_overall)
vb$labeled_cluster <- paste(vb$cluster_number,vb$label)
ggplot(vb,aes(bs_value)) + geom_histogram() + facet_wrap(~label)
# outliers -- start with 2.5% as in H; looks pretty good to me, tails just a bit long
# vb <- vb %>% mutate_at(.vars = 'bs_value', winsor, trim = .025)
# ggplot(hb,aes(scale(bs_value))) + geom_histogram() + facet_wrap(~labeled_cluster)

# alternatively, remove outliers
vb <- vb %>% mutate_at(.vars = 'bs_value', remove_outliers)
ggplot(vb,aes(scale(bs_value))) + geom_histogram() + facet_wrap(~labeled_cluster)
#just the bs and 
if (scale_within)
  
{vb <- vb %>% dplyr::group_by(feat_input_id,labeled_cluster) %>% mutate(bs_value = scale(bs_value)) %>% ungroup}


vb_wide <- spread(vb[,c("feat_input_id", "run", "trial", "bs_value", "labeled_cluster")],labeled_cluster,bs_value)
# head(h_wide)
# with group-fixed parameters we don't seem to need the winsorization step!!!
# h_wide <- h_wide %>% mutate_if(is.double, winsor,trim = .075)

just_bsv <- vb_wide[,3:ncol(vb_wide)]
# winsorize to deal with beta ouliers

# non-parametric correlations to deal with outliers
bs_cor <- corr.test(just_bsv,method = 'pearson')
# parametric correlations on winsorised betas
# clust_cor <- cor(just_rois_w,method = 'pearson')

setwd('~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/')
pdf("v_bs_censored_corr.pdf", width=12, height=12)
corrplot(bs_cor$r, cl.lim=c(-1,1),
         method = "circle", tl.cex = 1.5, type = "upper", tl.col = 'black',
         order = "hclust", diag = FALSE,
         addCoef.col="black", addCoefasPercent = FALSE,
         p.mat = bs_cor$p, sig.level=0.01, insig = "blank")
dev.off()

# try within-subject version
# vb_wide <- vb_wide %>% dplyr::group_by(id,run) %>% mutate(v_max_wi = scale(v_max),
#                                          v_max_wi_lag = lag(v_max_wi),
#                                          v_entropy_wi = scale(v_entropy),
#                                          v_max_b = mean(na.omit(v_max)),
#                                          v_entropy_b = mean(na.omit(v_entropy))
# )
# 

# remove trial
just_bsv <- just_bsv[,2:ncol(just_bsv)]

# more factors than the betas, but the structure is similar
# may refine later
vbs.fa = psych::fa(just_bsv, nfactors=5)
vfscores <- factor.scores(just_bsv, vbs.fa)$scores
vb_wide$vb_f1_lo_DAN <- vfscores[,1]
vb_wide$vb_f2_hi_vmPFC_cOFC <- vfscores[,2]
vb_wide$vb_f5_lo_ACC <- vfscores[,5]
vb_wide$vb_f4_lo_cerebell_crus <- vfscores[,4]
vb_wide$vb_f3_hi_blITG <- vfscores[,3]


######### 
# DECAY
setwd('~/Box Sync/skinner/projects_analyses/SCEPTIC/fMRI_paper/signals_review/beta_series_dec2018/')
metaD <- read_csv('d_auc_cluster_metadata.csv')
metaD$label <- substr(metaD$label,22,100)
metaD_overall <- metaD[metaD$l2_contrast == 'overall' & metaD$l3_contrast == 'Intercept' & metaD$model == 'Intercept-Age',]
dbs <- read_csv("d_auc_roi_beta_series.csv")
# filter out clusters <150 voxels
db <- dbs[dbs$l2_contrast == 'overall' & dbs$l3_contrast == 'Intercept' & dbs$model == 'Intercept-Age',1:6] %>% filter(cluster_number<15)

# head(merge(h,meta))
droi_list <- distinct(metaD_overall[,c(5,12)])

# inspect distributions
db <- inner_join(db,metaD_overall)
db$labeled_cluster <- paste(db$cluster_number,db$label)
ggplot(db,aes(scale(bs_value))) + geom_histogram() + facet_wrap(~label)
# outliers -- start with 2.5% as in H; looks pretty good to me, tails just a bit long like with Vmax
# db <- db %>% mutate_at(.vars = 'bs_value', winsor, trim = .025)
db <- db %>% mutate_at(.vars = 'bs_value', remove_outliers)
ggplot(db,aes(bs_value)) + geom_histogram() + facet_wrap(~labeled_cluster)
if (scale_within)
{db <- db %>% dplyr::group_by(feat_input_id,labeled_cluster) %>% mutate(bs_value = scale(bs_value)) %>% ungroup}


#just the bs and 

db_wide <- spread(db[,c("feat_input_id", "run", "trial", "bs_value", "labeled_cluster")],labeled_cluster,bs_value)
# head(h_wide)
# with group-fixed parameters we don't seem to need the winsorization step!!!
# h_wide <- h_wide %>% mutate_if(is.double, winsor,trim = .075)

just_bsd <- db_wide[,3:ncol(db_wide)]
# winsorize to deal with beta ouliers

# non-parametric correlations to deal with outliers
bs_cor <- corr.test(just_bsd,method = 'pearson')
# parametric correlations on winsorised betas
# clust_cor <- cor(just_rois_w,method = 'pearson')

setwd('~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/')
pdf("d_bs_corr.pdf", width=12, height=12)
corrplot(bs_cor$r, cl.lim=c(-1,1),
         method = "circle", tl.cex = 1.5, type = "upper", tl.col = 'black',
         order = "hclust", diag = FALSE,
         addCoef.col="black", addCoefasPercent = FALSE,
         p.mat = bs_cor$p, sig.level=0.01, insig = "blank")
dev.off()
# remove trial
just_bsd <- just_bsd[,2:ncol(just_bsd)]

# clear 4-factor solution
# very nice that the same factors are recovered, esp f3 and f4!
dbs.fa = psych::fa(just_bsd, nfactors=4)
dfscores <- factor.scores(just_bsd, dbs.fa)$scores
db_wide$db_f1_rIFG_rSMA <- dfscores[,1]
db_wide$db_f2_VS <- dfscores[,2]
db_wide$db_f3_occ_parietal <- dfscores[,3]
db_wide$db_f4_ACC_ins <- dfscores[,4]

vh_b_wide <- inner_join(vb_wide,hb_wide, by = c("feat_input_id", "run", "trial"))
dvh_b_wide <- inner_join(vh_b_wide,db_wide, by = c("feat_input_id", "run", "trial"))

# check beta correlations across signals at factor levels

dvh <-  dvh_b_wide[,grepl("b_f", names(dvh_b_wide))]
dvh_cor <- corr.test(dvh,method = 'pearson', adjust = 'none')

setwd('~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/')
pdf("dvh_b_corr_fixed_wi.pdf", width=20, height=20)
corrplot(dvh_cor$r, cl.lim=c(-.1,1),
         method = "circle", tl.cex = 1.5, type = "upper", tl.col = 'black',
         order = "hclust", diag = FALSE,
         addCoef.col="black", addCoefasPercent = FALSE,
         p.mat = dvh_cor$p, sig.level=0.05, insig = "blank")
dev.off()
# within-subject scaling makes no difference
dvh_b_wide$run_trial <- dvh_b_wide$trial

# factor analysis of factors!
dvh.fa = psych::fa(dvh, nfactors=1)
# without rescaling, only a single-factor solution works

# this was instructive
# for further analyses, we want the following regions:
# hb_f1_DAN (alt. vb_f1_lo_DAN)
# hb_f2_paralimbic (alt vb_f2_hi_paralimbic)
# db_f4_ACC_ins -- compression
# db_f1_rIFG_SMA -- ?overload

# # try once again with al clusters
# bs_clusters <-  dvh_b_wide[,!grepl("b_f", names(dvh_b_wide))]
# 
# bs_clusters <- bs_clusters %>% select_if(is.double)
# clust.fa <- psych::fa(bs_clusters,nfactors = 6)
# # similar solution, won't bother

dvh_bs_factors_wide <- dvh_b_wide[,c(names(dvh), "feat_input_id", "run", "run_trial")]

# add ids

############# Add new beta series, starting with PE
setwd('~/Box Sync/skinner/projects_analyses/SCEPTIC/fMRI_paper/signals_review/beta_series_dec2018/')

metaPE <- read_csv("PE_max_cluster_metadata.csv")
metaPE$label <- substr(metaPE$label,22,100)
meta_PE_overall <- metaPE[metaPE$l2_contrast == 'overall' & metaPE$l3_contrast == 'Intercept' & metaPE$model == 'Intercept-Age',]
pebs <- read_csv("pe_max_roi_beta_series.csv.gz")
# filter out clusters <150 voxels
peb <- pebs[pebs$l2_contrast == 'overall' & pebs$l3_contrast == 'Intercept' & pebs$model == 'Intercept-Age',1:6] %>% filter(cluster_number<8 | cluster_number == 10 | cluster_number == 11)

# head(merge(h,meta))
rois <- distinct(meta_PE_overall[,c(5,12)])

# inspect distributions
peb <- inner_join(peb,meta_PE_overall)
peb$labeled_cluster <- paste(peb$cluster_number,peb$label)
# ggplot(peb,aes(scale(bs_value))) + geom_histogram() + facet_wrap(~label)
# outliers -- 2.5% looks best by eye
peb <- peb %>% mutate_at(.vars = 'bs_value', remove_outliers)
# ggplot(peb,aes(scale(bs_value))) + geom_histogram() + facet_wrap(~labeled_cluster)

ggplot(peb,aes(bs_value)) + geom_histogram() + facet_wrap(~labeled_cluster)

# scale each region within subject
if (scale_within)
{peb <- peb %>% dplyr::group_by(feat_input_id,labeled_cluster) %>% mutate(bs_value = scale(bs_value)) %>% ungroup}
ggplot(peb,aes(labeled_cluster,bs_value)) + geom_boxplot()

#just the bs and 

peb_wide <- spread(peb[,c("feat_input_id", "run", "trial", "bs_value", "labeled_cluster")],labeled_cluster,bs_value)
# head(h_wide)
# with group-fixed parameters we don't seem to need the winsorization step!!!
# h_wide <- h_wide %>% mutate_if(is.double, winsor,trim = .075)

just_bs <- peb_wide[,3:ncol(peb_wide)]
# winsorize to deal with beta ouliers

# non-parametric correlations to deal with outliers
bs_cor <- corr.test(just_bs,method = 'pearson')
# parametric correlations on winsorised betas
# clust_cor <- cor(just_rois_w,method = 'pearson')

setwd('~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/')
pdf("pe_bs_wi_corr.pdf", width=12, height=12)
corrplot(bs_cor$r, cl.lim=c(-1,1),
         method = "circle", tl.cex = 1.5, type = "upper", tl.col = 'black',
         order = "hclust", diag = FALSE,
         addCoef.col="black", addCoefasPercent = FALSE,
         p.mat = bs_cor$p, sig.level=0.01, insig = "blank")
dev.off()
# remove trial
just_bs <- just_bs[,2:ncol(just_bs)]

# 2-factor solution has a left post-central loading on the paralimbic peb_f2
# may refine later, let's prototype
mpebs <- nfactors(bs_cor$r, n=5, rotate = "oblimin", diagonal = FALSE,fm = "pa", n.obs = 70, SMC = FALSE)
# seems to suggest a 3-factor solution, but the cortico-striatal factors then split right-left; went with 2 factors
pebs.fa = psych::fa(just_bs, nfactors=2)
pefscores <- factor.scores(x = just_bs, pebs.fa)$scores
peb_wide$peb_f1_cort_str <- pefscores[,1]
peb_wide$peb_f2_p_hipp <- pefscores[,2]

# just the factors
peb_wide <-  peb_wide[,c(1:3,13:14)]
peb_wide$run_trial <- peb_wide$trial
# merge with other BS
dvhpe_b_wide <- inner_join(dvh_bs_factors_wide,peb_wide, by = c("feat_input_id", "run", "run_trial"))
# add in a. hipp
dvhpe_b_wide <- inner_join(dvhpe_b_wide,hb_wide[,c(1:3,12)])
dvhpe_b_wide$h_ant_hipp_b_f <- dvhpe_b_wide$`9 Left Hippocampus`

dvhpe <-  dvhpe_b_wide[,grepl("b_f", names(dvhpe_b_wide))]
dvhpe_cor <- corr.test(dvhpe,method = 'pearson', adjust = 'none')

setwd('~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/')
pdf("dvhpe_b_corr_fixed_wi.pdf", width=20, height=20)
corrplot(dvhpe_cor$r, cl.lim=c(0,1),
         method = "circle", tl.cex = 1.5, type = "upper", tl.col = 'black',
         order = "hclust", diag = FALSE,
         addCoef.col="black", addCoefasPercent = FALSE,
         p.mat = dvhpe_cor$p, sig.level=0.05, insig = "blank")
dev.off()

map_df  <- as.tibble(read.csv("~/Box Sync/skinner/projects_analyses/SCEPTIC/fMRI_paper/signals_review/MMClock_aroma_preconvolve_fse_groupfixed/sceptic-clock-feedback-v_entropy-preconvolve_fse_groupfixed/v_entropy/v_entropy-Intercept_design.txt", sep=""))

beta_fscores <- inner_join(dvhpe_b_wide,map_df[,c(1:2,4:15)])
beta_fscores$id <- beta_fscores$ID

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
beta_fscores <- inner_join(beta_fscores,sum_df)

# model parameters
params <- read_csv("~/code/clock_analysis/fmri/data/mmclock_fmri_decay_factorize_selective_psequate_mfx_sceptic_global_statistics.csv")
bdf <- inner_join(beta_fscores,params)
params_beta <- bdf[,c(names(dvh), "total_earnings", "LL", "alpha", "gamma", "beta")]
param_cor <- corr.test(params_beta,method = 'pearson', adjust = 'none')

setwd('~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/')
pdf("dvh_bs_param_corr_fixed.pdf", width=12, height=12)
corrplot(param_cor$r, cl.lim=c(-1,1),
         method = "circle", tl.cex = 1.5, type = "upper", tl.col = 'black',
         order = "hclust", diag = FALSE,
         addCoef.col="black", addCoefasPercent = FALSE,
         p.mat = param_cor$p, sig.level=0.05, insig = "blank")
dev.off()
# they don't correlate

# merge into trial-level data
bdf <- bdf[,c(names(beta_fscores), "LL", "alpha", "gamma", "beta")]

df <- inner_join(trial_df,bdf)

# check correlation with RTs -- at least not uniform
dvh_rt <- df[,c(names(dvh), "rt_csv")]
rt_cor <- corr.test(dvh_rt,method = 'pearson', adjust = 'none')

setwd('~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/')
pdf("dvh_bs_rt_corr.pdf", width=12, height=12)
corrplot(rt_cor$r, cl.lim=c(-1,1),
         method = "circle", tl.cex = 1.5, type = "upper", tl.col = 'black',
         order = "hclust", diag = FALSE,
         addCoef.col="black", addCoefasPercent = FALSE,
         p.mat = rt_cor$p, sig.level=0.05, insig = "blank")
dev.off()


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
                                         rt_vmax_change = rt_vmax - rt_vmax_lag
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


save(file = 'trial_df_and_vhd_bs.Rdata', df)
