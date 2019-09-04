library(dplyr)
library(tidyverse)
library(psych)
library(corrplot)
library(lme4)

# get H betas
setwd('~/Box/skinner/projects_analyses/SCEPTIC/fMRI_paper/signals_review/MMClock_aroma_preconvolve_fse_groupfixed/sceptic-clock-feedback-v_entropy-preconvolve_fse_groupfixed/v_entropy/')
meta <- read_csv("~/Box/skinner/projects_analyses/SCEPTIC/fMRI_paper/signals_review/MMClock_aroma_preconvolve_fse_groupfixed/sceptic-clock-feedback-v_entropy-preconvolve_fse_groupfixed/v_entropy/v_entropy_cluster_metadata.csv")
meta$label <- substr(meta$label,22,100)
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
mh <- nfactors(clust_cor, n=5, rotate = "oblimin", diagonal = FALSE,fm = "pa", n.obs = 70, SMC = FALSE)
h.fa = psych::fa(just_rois, nfactors=2, rotate = "oblimin", fm = "pa")
fscores <- factor.scores(just_rois, h.fa)$scores
h_wide$h_f1_fp <- fscores[,1]
h_wide$h_f2_neg_paralimb <- fscores[,2]
h_wide$h_HippAntL <- h_wide$`9`
h_wide$h_vmPFC <- h_wide$`6`

h_wide <- subset(h_wide, select = c("feat_input_id","h_f1_fp","h_f2_neg_paralimb", "h_HippAntL"))

#####
# add PE
setwd('~/Box/skinner/projects_analyses/SCEPTIC/fMRI_paper/signals_review/MMClock_aroma_preconvolve_fse_groupfixed/sceptic-clock-feedback-pe_max-preconvolve_fse_groupfixed/pe_max/')
pemeta <- read_csv("pe_max_cluster_metadata.csv")
pemeta$label <- substr(pemeta$label,22,100)
pemeta_overall <- pemeta[pemeta$l2_contrast == 'overall' & pemeta$l3_contrast == 'Intercept' & pemeta$model == 'Intercept-Age',]
pebetas <- read_csv("pe_max_roi_betas.csv")
pe <- as.tibble(pebetas[pebetas$l2_contrast == 'overall' & pebetas$l3_contrast == 'Intercept' & pebetas$model == 'Intercept-Age',1:3]) %>% filter(cluster_number<8 | cluster_number == 10 | cluster_number == 11)
# head(merge(h,meta))
perois_list <- distinct(pemeta_overall[c(1:7, 10:11),c(5,12)])

# inspect distributions
# for some reason, the clusters do not seem to fit the map -- check with MNH (big clusters broken down?) !!
perois <- inner_join(pe,pemeta_overall)
perois$labeled_cluster <- paste(perois$cluster_number,perois$label)
pemeta_overall$labeled_cluster <- paste(pemeta_overall$cluster_number,pemeta_overall$label)
ggplot(perois,aes(scale(cope_value))) + geom_histogram() + facet_wrap(~labeled_cluster)

pe_labeled <- inner_join(pe,perois_list)
pe_labeled$labeled_cluster <- paste(pe_labeled$cluster_number,pe_labeled$label)
pe_num <- select(pe_labeled,c(1,2,3))
pe_labeled <- select(pe_labeled,c(1,3,5))

pe_wide <- spread(pe_labeled,labeled_cluster,cope_value)
pe_wide_num <- spread(pe_num,cluster_number,cope_value)

# some outliers, let's winsorize for now
## THIS WAS A MINOR ERROR -- we are no longer winsorizing betas
# pe_wide <- pe_wide %>% mutate_if(is.double, winsor,trim = .075)
# pe_wide_num <- pe_wide_num %>% mutate_if(is.double, winsor,trim = .075)

pejust_rois <- pe_wide[,2:ncol(pe_wide)]
pejust_rois_num <- pe_wide_num[,2:ncol(pe_wide_num)]

# winsorize to deal with beta ouliers

# non-parametric correlations to deal with outliers
peclust_cor <- corr.test(pejust_rois,method = 'pearson', adjust = 'none')
# parametric correlations on winsorised betas
# clust_cor <- cor(just_rois_w,method = 'pearson')

# setwd('~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/')
# pdf("pe_cluster_corr_fixed.pdf", width=12, height=12)
# corrplot(peclust_cor$r, cl.lim=c(-1,1),
#          method = "circle", tl.cex = 1.5, type = "upper", tl.col = 'black',
#          order = "hclust", diag = FALSE,
#          addCoef.col="black", addCoefasPercent = FALSE,
#          p.mat = peclust_cor$p, sig.level=0.05, insig = "blank")
# dev.off()

mpe <- nfactors(peclust_cor$r, n=5, rotate = "oblimin", diagonal = FALSE,fm = "pa", n.obs = 70, SMC = FALSE)
pe.fa = psych::fa(pejust_rois_num, nfactors=2, rotate = "varimax", fm = "pa")

# library("lavaan")
# msyn <- '
# cort_str =~ 1*`1` + `11` + `2` + `3` +
#             `4` + `5` +
#             `6` #putting the 1* here for clarity, but it is the lavaan default
# hipp =~ 1*`10` + `7`
# cort_str ~~ 0*hipp #force orthogonal (uncorrelated) factors
# cort_str ~~ cort_str #explicitly indicate that lavaan should estimate a free parameters for factor variances
# hipp ~~ hipp
# '
# msyn_cor <- '
# cort_str =~ 1*`1` + `11` + `2` + `3` +
#             `4` + `5` +
# `6` #putting the 1* here for clarity, but it is the lavaan default
# hipp =~ 1*`10` + `7`
# cort_str ~~ hipp #force orthogonal (uncorrelated) factors
# cort_str ~~ cort_str #explicitly indicate that lavaan should estimate a free parameters for factor variances
# hipp ~~ hipp
# '
# 
# 
# mcfa <- cfa(msyn, pejust_rois_num)
# summary(mcfa, standardized=TRUE)
# mcfa_cor <- cfa(msyn_cor, pejust_rois_num)
# summary(mcfa_cor, standardized=TRUE)
# anova(mcfa,mcfa_cor)

pe.fa = psych::fa(pejust_rois, nfactors=2)
pefscores <- factor.scores(pejust_rois, pe.fa)$scores
pe_wide$pe_f1_cort_str <- pefscores[,1]
pe_wide$pe_f2_hipp <- pefscores[,2]

# dvhkpe_wide <- inner_join(dvhkf_wide,pe_wide[,c("feat_input_id","pe_f1_cort_str", "pe_f2_hipp")])
hpe_wide <- inner_join(h_wide,pe_wide[,c("feat_input_id","pe_f1_cort_str", "pe_f2_hipp")])

#####
# add ids
map_df  <- as.tibble(read.csv("~/Box/skinner/projects_analyses/SCEPTIC/fMRI_paper/signals_review/MMClock_aroma_preconvolve_fse_groupfixed/sceptic-clock-feedback-v_entropy-preconvolve_fse_groupfixed/v_entropy/v_entropy-Intercept_design.txt", sep=""))

pc_scores <- inner_join(hpe_wide,map_df[,c(1:2,4:15)])
pc_scores$id <- pc_scores$ID


# get trial_level data
trial_df <- read_csv("~/code/clock_analysis/fmri/data/mmclock_fmri_decay_factorize_selective_psequate_mfx_trial_statistics.csv.gz")
# trial_df <- read_csv("~/Box/SCEPTIC_fMRI/mmclock_fmri_fixed_uv_mfx_trial_statistics.csv.gz")
u_df <- read_csv("~/Box/SCEPTIC_fMRI/mmclock_fmri_fixed_uv_ureset_mfx_trial_statistics.csv.gz")
trial_df <- trial_df %>%
  group_by(id, run) %>%  dplyr::mutate(rt_swing = abs(c(NA, diff(rt_csv))),
                                       rt_swing_lr = abs(log(rt_csv/lag(rt_csv))),
                                       rt_lag = lag(rt_csv) ,
                                       rt_swing_lag = lag(rt_swing),
                                       omission_lag = lag(score_csv==0),
                                       rt_vmax_lag = lag(rt_vmax),
                                       run_trial=1:50) %>% ungroup() #compute rt_swing within run and subject
u_df <- u_df %>% select(id, run, trial, u_chosen, u_chosen_lag, u_chosen_change)

trial_df <- inner_join(trial_df,u_df)

# performance
sum_df <- trial_df %>% group_by(id) %>% dplyr::summarize(total_earnings = sum(score_csv)) %>% arrange(total_earnings)
beta_sum <- inner_join(pc_scores,sum_df)

# model parameters
params <- read_csv("~/code/clock_analysis/fmri/data/mmclock_fmri_decay_factorize_selective_psequate_mfx_sceptic_global_statistics.csv")
sub_df <- inner_join(beta_sum,params)
sub_df$id <- sub_df$ID
params_beta <- sub_df[,c("h_f1_fp", "h_f2_neg_paralimb","h_HippAntL",
                         "pe_f1_cort_str", "pe_f2_hipp", 
                         "total_earnings", "LL", "alpha", "gamma", "beta")]
param_cor <- corr.test(params_beta,method = 'pearson', adjust = 'none')

setwd('~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/')

# merge into trial-level data
df <- inner_join(trial_df,sub_df, by = "id")
df$rewFunc <- relevel(as.factor(df$rewFunc),ref = "CEV")
df$rewFuncIEVsum <- df$rewFunc
contrasts(df$rewFuncIEVsum) <- contr.sum
colnames(contrasts(df$rewFuncIEVsum)) <- c('CEV','CEVR', 'DEV')
# dichotomize betas for plotting
df$h_fp <- 'low'
df$h_fp[df$h_f1_fp>0] <- 'high'
df$low_h_paralimbic <- 'low'
df$low_h_paralimbic[df$h_f2_neg_paralimb<0] <- 'high'

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

# correlate between-subject V and H with clusters
b_df <- df %>% group_by(id) %>% dplyr::summarise(v_maxB = mean(v_max, na.rm = T),
                                          v_entropyB = mean(v_entropy, na.rm = T))

sub_df <- inner_join(sub_df, b_df, by = 'id')
bdf <- sub_df[,c("h_f1_fp", "h_f2_neg_paralimb",
                 "pe_f1_cort_str", "pe_f2_hipp",
                 "total_earnings", "LL", "alpha", "gamma", "beta", "v_maxB", "v_entropyB")]
b_cor <- corr.test(bdf,method = 'pearson', adjust = 'none')

setwd('~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/')
pdf("between_subject_v_h_kf_pe_beh_corr_fixed.pdf", width=12, height=12)
corrplot(b_cor$r, cl.lim=c(-1,1),
         method = "circle", tl.cex = 1.5, type = "upper", tl.col = 'black',
         order = "hclust", diag = FALSE,
         addCoef.col="black", addCoefasPercent = FALSE,
         p.mat = b_cor$p, sig.level=0.05, insig = "blank")
dev.off()

df$h_f1_fp_resp <- 'low'
df$h_f1_fp_resp[df$h_f1_fp>0] <- 'high'
df$pe_f1_cort_str_resp <- 'low'
df$pe_f1_cort_str_resp[df$pe_f1_cort_str>0] <- 'high'
df$pe_f2_hipp_resp <- 'low'
df$pe_f2_hipp_resp[df$pe_f2_hipp>0] <- 'high'

df$h_HippAntL_resp <- 'low'
df$h_HippAntL_resp[df$h_HippAntL<mean(df$h_HippAntL)] <- 'high'
df$last_outcome <- NA
df$last_outcome[df$omission_lag] <- 'Omission'
df$last_outcome[!df$omission_lag] <- 'Reward'

# circular, but just check to what extent each area conforms to SCEPTIC-SM: looks like there is an interaction, both need to be involved again
# ggplot(df, aes(run_trial, v_entropy_wi, color = low_h_paralimbic, lty = h_fp)) + geom_smooth(method = "loess") #+ facet_wrap(~gamma>0)
# # I think this means again that H estimates are more precise for high-paralimbic people, esp. late in learning
# ggplot(df, aes( v_entropy_wi,log(rt_swing), color = low_h_paralimbic, lty = h_fp)) + geom_smooth(method = "gam") + facet_wrap(~learning_epoch)
# ggplot(df, aes( v_max_wi,log(rt_swing), color = low_h_paralimbic, lty = h_fp)) + geom_smooth(method = "gam") + facet_wrap(~run_trial > 10)
# 
# Okay, some behavioral relevance of KLD
# ggplot(df, aes(run_trial, v_entropy_wi, color = k_f1_all_pos_resp)) + geom_smooth(method = "loess")


save(file = 'trial_df_and_vh_pe_clusters_u.Rdata', df)
