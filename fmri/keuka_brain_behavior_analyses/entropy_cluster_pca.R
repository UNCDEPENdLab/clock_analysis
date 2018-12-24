library(dplyr)
library(tidyverse)
library(psych)
library(corrplot)
library(lme4)

# get H betas
setwd('~/Box Sync/skinner/projects_analyses/SCEPTIC/fMRI_paper/signals_review/MMClock_aroma_preconvolve_fse_groupfixed/sceptic-v_entropy-preconvolve_fse_groupfixed/v_entropy/')
meta <- read_csv("~/Box Sync/skinner/projects_analyses/SCEPTIC/fMRI_paper/signals_review/MMClock_aroma_preconvolve_fse_groupfixed/sceptic-v_entropy-preconvolve_fse_groupfixed/v_entropy/v_entropy_cluster_metadata.csv")
meta$label <- substr(meta$label,14,100)
meta_overall <- meta[meta$l2_contrast == 'overall' & meta$l3_contrast == 'Intercept' & meta$model == 'Intercept-Age',]
Hbetas <- read_csv("v_entropy_roi_betas.csv")
h <- as.tibble(Hbetas[Hbetas$l2_contrast == 'overall' & Hbetas$l3_contrast == 'Intercept' & Hbetas$model == 'Intercept-Age',1:3]) %>% filter(cluster_number<11)
# head(merge(h,meta))
rois <- distinct(meta_overall[,c(5,12)])

# inspect distributions
hrois <- inner_join(h,meta_overall)
hrois$labeled_cluster <- paste(hrois$cluster_number,hrois$label)
ggplot(hrois,aes(scale(cope_value))) + geom_histogram() + facet_wrap(~label)
h_wide <- spread(h,cluster_number,cope_value)
# head(h_wide)
# with group-fixed parameters we don't seem to need the winsorization step!!!
# h_wide <- h_wide %>% mutate_if(is.double, winsor,trim = .075)

just_rois <- h_wide[,2:11]
# winsorize to deal with beta ouliers

# non-parametric correlations to deal with outliers
clust_cor <- cor(just_rois,method = 'pearson')
# parametric correlations on winsorised betas
# clust_cor <- cor(just_rois_w,method = 'pearson')

setwd('~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/')
pdf("h_cluster_corr_fixed.pdf", width=12, height=12)
corrplot(clust_cor, cl.lim=c(-1,1),
         method = "circle", tl.cex = 1.5, type = "upper", tl.col = 'black',
         order = "hclust", diag = FALSE,
         addCoef.col="black", addCoefasPercent = FALSE,
         p.mat = 1-clust_cor, sig.level=0.75, insig = "blank")
dev.off()

h.pca = prcomp((just_rois),scale = TRUE, center = TRUE)
h_wide$h_pc1 <- h.pca$x[,1]
h_wide$h_pc2 <- h.pca$x[,2]

h.fa = psych::fa(just_rois, nfactors=2)
fscores <- factor.scores(just_rois, h.fa)$scores
h_wide$h_f1_fp <- fscores[,1]
h_wide$h_f2_neg_paralimb <- fscores[,2]
h_wide$ParaHippL <- h_wide$`9`
h_wide$vmPFC <- h_wide$`6`

# add value betas -- use v_max as theoretically interesting and unbiased by choice
setwd('~/Box Sync/skinner/projects_analyses/SCEPTIC/fMRI_paper/signals_review/MMClock_aroma_preconvolve_fse_groupfixed/sceptic-v_max-preconvolve_fse_groupfixed/v_max/')
vmeta <- read_csv("v_max_cluster_metadata.csv")
vmeta$label <- substr(vmeta$label,14,100)
vmeta_overall <- vmeta[vmeta$l2_contrast == 'overall' & vmeta$l3_contrast == 'Intercept' & vmeta$model == 'Intercept-Age',]
vbetas <- read_csv("v_max_roi_betas.csv")
v <- as.tibble(vbetas[vbetas$l2_contrast == 'overall' & vbetas$l3_contrast == 'Intercept' & vbetas$model == 'Intercept-Age',1:3]) %>% filter(cluster_number<15)
vrois_list <- distinct(vmeta_overall[,c(5,12)])


# inspect distributions
vrois <- inner_join(v,vmeta_overall)
vrois$labeled_cluster <- paste(vrois$cluster_number,vrois$label)
ggplot(vrois,aes(scale(cope_value))) + geom_histogram() + facet_wrap(~labeled_cluster)
v_wide <- spread(v,cluster_number,cope_value)
# try v winsorization given the left outliers
# v_wide <- v_wide %>% mutate_if(is.double, winsor,trim = .075)

# few outliers, let's hold off on winsorizing for now
# h_wide <- h_wide %>% mutate_if(is.double, winsor,trim = .075)

# I am worried about these scores: a lot of the V+ f2 variance comes from fusiform gyri (4,11)
# let's redo using only fronto-striatal clusters (vmPFC - 10, BG - 2, vlPFC - 9,14) removing fusiform (4,11)
# and cerebellar (8) clusters
v_wide <- subset(v_wide, select = -c(`4`, `11`,`8`))
v_wide <- v_wide %>% mutate_if(is.double, winsor,trim = .075)


vjust_rois <- v_wide[,2:ncol(v_wide)]
# winsorize to deal with beta ouliers

# non-parametric correlations to deal with outliers
# parametric correlations on winsorised betas
# clust_cor <- cor(just_rois_w,method = 'pearson')
v.pca = prcomp((vjust_rois),scale = TRUE, center = TRUE)
v_wide$v_pc1 <- v.pca$x[,1]
v_wide$v_pc2 <- v.pca$x[,2]

v.fa = psych::fa(vjust_rois, nfactors=3)
vfscores <- factor.scores(vjust_rois, v.fa)$scores
v_wide$v_f1_neg_cog <- vfscores[,1]
v_wide$v_f2_paralimb <- vfscores[,2]

vclust_cor <- corr.test(v_wide %>% select_if(is.double),method = 'pearson', adjust = 'none')

setwd('~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/')
pdf("v_cluster_corr_fixed.pdf", width=12, height=12)
corrplot(vclust_cor$r, cl.lim=c(-1,1),
         method = "circle", tl.cex = 1.5, type = "upper", tl.col = 'black',
         order = "hclust", diag = FALSE,
         addCoef.col="black", addCoefasPercent = FALSE,
         p.mat = vclust_cor$p, sig.level=0.05, insig = "blank")
dev.off()

# a bit crazy, but let's put v and h into a single factor analysis
vh_wide <- inner_join(subset(v_wide, select = c("feat_input_id","v_f1_neg_cog","v_f2_paralimb")),
                      subset(h_wide, select = c("feat_input_id","h_f1_fp","h_f2_neg_paralimb")), by = "feat_input_id")
vhjust_rois <- vh_wide %>% select_if(is.double)
vhclust_cor <- corr.test(vhjust_rois,method = 'pearson', adjust = 'none')
pdf("vh_cluster_corr_fixed.pdf", width=12, height=12)
corrplot(vhclust_cor$r, cl.lim=c(-1,1),
         method = "circle", tl.cex = 1.5, type = "upper", tl.col = 'black',
         order = "hclust", diag = FALSE,
         addCoef.col="black", addCoefasPercent = FALSE,
         p.mat = vhclust_cor$p, sig.level=0.05, insig = "blank")
dev.off()
# the first value (V+) factor is correlated with the entropy hippocampal factor (.55)

# OK, that went well -- let's add dAUC
setwd('~/Box Sync/skinner/projects_analyses/SCEPTIC/fMRI_paper/signals_review/MMClock_aroma_preconvolve_fse_groupfixed/sceptic-d_auc-preconvolve_fse_groupfixed/d_auc/')
dmeta <- read_csv("d_auc_cluster_metadata.csv")
dmeta$label <- substr(dmeta$label,14,100)
dmeta_overall <- dmeta[dmeta$l2_contrast == 'overall' & dmeta$l3_contrast == 'Intercept' & dmeta$model == 'Intercept-Age',]
dbetas <- read_csv("d_auc_roi_betas.csv")
d <- as.tibble(dbetas[dbetas$l2_contrast == 'overall' & dbetas$l3_contrast == 'Intercept' & dbetas$model == 'Intercept-Age',1:3]) %>% filter(cluster_number<15)
# head(merge(h,meta))
drois_list <- distinct(dmeta_overall[,c(5,12)])

# inspect distributions
# for some reason, the clusters do not seem to fit the map -- check with MNH (big clusters broken down?) !!
drois <- inner_join(d,dmeta_overall)
drois$labeled_cluster <- paste(drois$cluster_number,drois$label)
ggplot(drois,aes(scale(cope_value))) + geom_histogram() + facet_wrap(~labeled_cluster)
d_wide <- spread(d,cluster_number,cope_value)

# few outliers, let's hold off on winsorizing for now
# h_wide <- h_wide %>% mutate_if(is.double, winsor,trim = .075)

djust_rois <- d_wide[,2:ncol(d_wide)]
# winsorize to deal with beta ouliers

# non-parametric correlations to deal with outliers
dclust_cor <- cor(djust_rois,method = 'pearson')
# parametric correlations on winsorised betas
# clust_cor <- cor(just_rois_w,method = 'pearson')

setwd('~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/')
pdf("d_cluster_corr_fixed.pdf", width=12, height=12)
corrplot(dclust_cor, cl.lim=c(-1,1),
         method = "circle", tl.cex = 1.5, type = "upper", tl.col = 'black',
         order = "hclust", diag = FALSE,
         addCoef.col="black", addCoefasPercent = FALSE,
         p.mat = 1-dclust_cor, sig.level=0.75, insig = "blank")
dev.off()
# looks like a single dimension

d.pca = prcomp((djust_rois),scale = TRUE, center = TRUE)
d_wide$d_pc1 <- d.pca$x[,1]
d_wide$d_pc2 <- d.pca$x[,2]

d.fa = psych::fa(djust_rois, nfactors=3)
dfscores <- factor.scores(djust_rois, d.fa)$scores
d_wide$d_f1_FP_SMA <- dfscores[,1]
d_wide$d_f2_VS <- dfscores[,2]
d_wide$d_f3_ACC_ins <- dfscores[,3]

# a bit crazy, but let's put v and h into a single factor analysis
dvh_wide <- inner_join(vh_wide,d_wide[,c(1,18:20)], by = "feat_input_id")
dvhjust_rois <- dvh_wide[,2:ncol(dvh_wide)]
dvhclust_cor <- corr.test(dvhjust_rois,method = 'pearson', adjust = 'none')
pdf("dvh_cluster_corr_fixed.pdf", width=12, height=12)
corrplot(dvhclust_cor$r, cl.lim=c(-1,1),
         method = "circle", tl.cex = 1.5, type = "upper", tl.col = 'black',
         order = "hclust", diag = FALSE,
         addCoef.col="black", addCoefasPercent = FALSE,
         p.mat = dvhclust_cor$p, sig.level=0.05, insig = "blank")
dev.off()

# add ids
map_df  <- as.tibble(read.csv("~/Box Sync/skinner/projects_analyses/SCEPTIC/fMRI_paper/signals_review/MMClock_aroma_preconvolve_fse_groupfixed/sceptic-v_entropy-preconvolve_fse_groupfixed/v_entropy/v_entropy-Intercept_design.txt", sep=""))

pc_scores <- inner_join(dvh_wide,map_df[,c(1:2,4:15)])
pc_scores$id <- pc_scores$ID


# get trial_level data
trial_df <- read_csv("~/code/clock_analysis/fmri/data/mmclock_fmri_decay_factorize_selective_psequate_mfx_trial_statistics.csv.gz")
trial_df <- trial_df %>%
  group_by(id, run) %>%  dplyr::mutate(rt_swing = abs(c(NA, diff(rt_csv))),
                                       rt_swing_lr = abs(log(rt_csv/lag(rt_csv))),
                                       rt_lag = lag(rt_csv) ,
                                       rt_swing_lag = lag(rt_swing),
                                       omission_lag = lag(score_csv==0),
                                       rt_vmax_lag = lag(rt_vmax),
                                       run_trial=1:50) %>% ungroup() #compute rt_swing within run and subject

# performance
sum_df <- trial_df %>% group_by(id) %>% dplyr::summarize(total_earnings = sum(score_csv)) %>% arrange(total_earnings)
beta_sum <- inner_join(pc_scores,sum_df)
plot(beta_sum$h_f1_fp,beta_sum$total_earnings)
plot(beta_sum$h_f2_neg_paralimb,beta_sum$total_earnings)
plot(beta_sum$v_f1_neg_cog,beta_sum$total_earnings)
plot(beta_sum$v_f2_paralimb,beta_sum$total_earnings)

# model parameters
params <- read_csv("~/code/clock_analysis/fmri/data/mmclock_fmri_decay_factorize_selective_psequate_mfx_sceptic_global_statistics.csv")
sub_df <- inner_join(beta_sum,params)
params_beta <- sub_df[,c("v_f1_neg_cog","v_f2_paralimb","h_f1_fp", "h_f2_neg_paralimb","d_f1_FP_SMA","d_f2_VS","d_f3_ACC_ins", "total_earnings", "LL", "alpha", "gamma", "beta")]
param_cor <- corr.test(params_beta,method = 'pearson', adjust = 'none')

setwd('~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/')
pdf("dvhbeta_param_corr_fixed.pdf", width=12, height=12)
corrplot(param_cor$r, cl.lim=c(-1,1),
         method = "circle", tl.cex = 1.5, type = "upper", tl.col = 'black',
         order = "hclust", diag = FALSE,
         addCoef.col="black", addCoefasPercent = FALSE,
         p.mat = param_cor$p, sig.level=0.05, insig = "blank")
dev.off()

# merge into trial-level data
df <- inner_join(trial_df,sub_df)
df$rewFunc <- relevel(as.factor(df$rewFunc),ref = "CEV")
df$rewFuncIEVsum <- df$rewFunc
contrasts(df$rewFuncIEVsum) <- contr.sum
colnames(contrasts(df$rewFuncIEVsum)) <- c('CEV','CEVR', 'DEV')
# dichotomize betas for plotting
df$h_fp <- 'low'
df$h_fp[df$h_f1_fp>0] <- 'high'
df$low_h_paralimbic <- 'low'
df$low_h_paralimbic[df$h_f2_neg_paralimb<0] <- 'high'

df$v_paralimbic <- 'low'
df$v_paralimbic[df$v_f2_paralimb>0] <- 'high'
#NB: v_f2_paralimb is mostly vlPFC driven now...
df$low_v_fp_acc_vlpfc <- 'low'
df$low_v_fp_acc_vlpfc[df$v_f1_neg_cog<0] <- 'high'

df$learning_epoch <- 'trials 1-10'
df$learning_epoch[df$run_trial>10] <- 'trials 11-50'




# obtain within-subject v_max and entropy: correlated at -.37

df <- df %>% group_by(id,run) %>% mutate(v_max_wi = scale(v_max),
                                         v_max_wi_lag = lag(v_max_wi),
                                         v_entropy_wi = scale(v_entropy),
                                         v_max_b = mean(na.omit(v_max)),
                                         v_entropy_b = mean(na.omit(v_entropy))
)

# correlate between-subject V and H with clusters
b_df <- df %>% group_by(id) %>% dplyr::summarise(v_maxB = mean(v_max, na.rm = T),
                                          v_entropyB = mean(v_entropy, na.rm = T))

sub_df <- inner_join(sub_df, b_df, by = 'id')
bdf <- sub_df[,c("v_f1_neg_cog","v_f2_paralimb","h_f1_fp", "h_f2_neg_paralimb","d_f1_FP_SMA","d_f2_VS","d_f3_ACC_ins", 
                 "total_earnings", "LL", "alpha", "gamma", "beta", "v_maxB", "v_entropyB")]
b_cor <- corr.test(bdf,method = 'pearson', adjust = 'none')

setwd('~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/')
pdf("between_subject_v_h_beh_corr_fixed.pdf", width=12, height=12)
corrplot(b_cor$r, cl.lim=c(-1,1),
         method = "circle", tl.cex = 1.5, type = "upper", tl.col = 'black',
         order = "hclust", diag = FALSE,
         addCoef.col="black", addCoefasPercent = FALSE,
         p.mat = b_cor$p, sig.level=0.05, insig = "blank")
dev.off()
df$d_f1_FP_SMAresp <- 'low'
df$d_f1_FP_SMAresp[df$d_f1_FP_SMA>0] <- 'high'
df$d_f2_VSresp <- 'low'
df$d_f2_VSresp[df$d_f2_VS>0] <- 'high'
df$d_f3_ACC_ins_resp <- 'low'
df$d_f3_ACC_ins_resp[df$d_f3_ACC_ins>0] <- 'high'
df$d_f2_VSresp <- as.factor(df$d_f2_VSresp)
df$d_f2_VSresp <- relevel(df$d_f2_VSresp, ref = 'low') 

df$last_outcome <- NA
df$last_outcome[df$omission_lag] <- 'Omission'
df$last_outcome[!df$omission_lag] <- 'Reward'

# circular, but just check to what extent each area conforms to SCEPTIC-SM: looks like there is an interaction, both need to be involved again
ggplot(df, aes(run_trial, v_entropy_wi, color = low_h_paralimbic, lty = h_fp)) + geom_smooth(method = "loess") #+ facet_wrap(~gamma>0)
# I think this means again that H estimates are more precise for high-paralimbic people, esp. late in learning
ggplot(df, aes( v_entropy_wi,log(rt_swing), color = low_h_paralimbic, lty = h_fp)) + geom_smooth(method = "gam") + facet_wrap(~learning_epoch)
ggplot(df, aes( v_max_wi,log(rt_swing), color = low_h_paralimbic, lty = h_fp)) + geom_smooth(method = "gam") + facet_wrap(~run_trial > 10)

save(file = 'trial_df_and_vhd_clusters.Rdata', df)
