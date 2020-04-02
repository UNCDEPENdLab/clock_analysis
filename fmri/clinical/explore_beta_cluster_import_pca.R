# preprocess voxel-wise fMRI data for brain-to-behavior analyses
# extracts hippocampal activation coefficients for prediction error and entropy analyses
# merges with behavioral data from the fMRI and replication (MEG sessions)

library(dplyr)
library(tidyverse)
library(psych)
library(corrplot)
library(lme4)
library(stargazer)

clock_folder <- "~/code/clock_analysis" #alex
# clock_folder <- "~/Data_Analysis/clock_analysis" #michael


unsmoothed = F # unsmoothed data used only for sensitivity analyses

# get H (entropy) betas
sort_names <- function(data) {
  name  <- names(data)
  chars <- keep(name, grepl, pattern = "[^0-9]") %>% sort()
  nums  <- discard(name, grepl, pattern = "[^0-9]") %>% 
    as.numeric() %>% 
    sort() %>% 
    sprintf("`%s`", .)
  
  select_(data, .dots = c(chars, nums))
}

if (unsmoothed) {
  setwd('~/Box/SCEPTIC_fMRI/')
  meta <- read_csv("~/Box/SCEPTIC_fMRI/v_entropy_unsmoothed_cluster_metadata.csv")
  Hbetas <- read_csv("v_entropy_unsmoothed_subj_betas.csv")
} else {
  setwd('~/Box/skinner/data/MRI/clock_explore/')
  meta <- read_csv("v_entropy_z_index.csv")
  load("explore_extracted_roivalue.rdata")
}

betas <- as_tibble(extracted_roi_df)
# AH entropy cluster: #9
# PH PE clusters: #7 and #10
ggplot(betas,aes(scale(pe_max_cluster_7_3mm))) + geom_histogram()
ggplot(betas,aes(scale(pe_max_cluster_10_3mm))) + geom_histogram()
ggplot(betas,aes(scale(v_entropy_cluster_9_3mm))) + geom_histogram()


just_rois <- betas %>% select(-c(ID, pe_max_left_hipp_3mm,pe_max_right_hipp_3mm,v_entropy_left_hipp_3mm,v_entropy_right_hipp_3mm))
# winsorize to deal with beta ouliers
# check missingness
library(VIM)
df_aggr = aggr(just_rois, col=mdc(1:2), numbers=TRUE, sortVars=TRUE, labels=names(just_rois), cex.axis=.7, gap=3, ylab=c("Proportion of missingness","Missingness Pattern"))

# non-parametric correlations to deal with outliers
clust_cor <- cor(just_rois,method = 'pearson', use = "complete.obs")
# parametric correlations on winsorised betas
# clust_cor <- cor(just_rois_w,method = 'pearson')

setwd(file.path(clock_folder, 'fmri/clinical/'))
if (unsmoothed) {
  pdf("h_cluster_corr_fixed_unsmoothed.pdf", width=12, height=12)  
} else {
  pdf("cluster_corr_fixed.pdf", width=24, height=24)  
}
corrplot(clust_cor, cl.lim=c(-1,1),
         method = "circle", tl.cex = 1.5, type = "upper", tl.col = 'black',
         order = "hclust", diag = FALSE,
         addCoef.col="black", addCoefasPercent = FALSE,
         p.mat = 1-clust_cor, sig.level=0.75, insig = "blank")
dev.off()

###########
# go for the jugular and replicate the brain-behavior effects from Nat Comm paper
betas <- betas %>% rename(id = ID)
allbetas <- betas
betas <- betas %>% mutate(
  PH_pe = (pe_max_cluster_7_3mm + pe_max_cluster_10_3mm)/2,
  AH_h_neg = - v_entropy_cluster_9_3mm) %>%  
  select(PH_pe, AH_h_neg, id)

save(betas, file = "~/Box/skinner/data/MRI/clock_explore/explore_betas.rdata")
save(allbetas, file = "~/Box/skinner/data/MRI/clock_explore/explore_allbetas.rdata")

#########


hwide <- betas %>% select(contains('v_entropy_clu'))
# not really recovering the factor structure from Nat Comm
h.fa = fa.sort(psych::fa(hwide, nfactors=3, rotate = "varimax", fm = "pa"))
# fscores <- factor.scores(just_rois, h.fa)$scores
# h_wide$h_f1_fp <- fscores[,1]
# h_wide$h_f2_neg_paralimb <- fscores[,2]
# h_wide$h_HippAntL <- h_wide$`9 Left Hippocampus`
# if (unsmoothed) {
#   h_wide$h_vmPFC <- h_wide$`5 Left Mid Orbital Gyrus`
# } else {h_wide$h_vmPFC <- h_wide$`6`}
# hf <- as_tibble(cbind(rownames(unclass(round(h.fa$Structure, digits = 3))),unclass(round(h.fa$Structure, digits = 3)))) %>% arrange(desc(PA1) ) %>%  
#   rename(region = V1, `Factor 1` = PA1, `Factor 2` = PA2)
# hf$region <- gsub('[0-9;]+', '', hf$region)
# stargazer(hf, type = 'html', out = '~/OneDrive/collected_letters/papers/sceptic_fmri/hippo/supp/h_fa_structure.html', summary = F)



#####
# add PE
if (unsmoothed) {
  setwd('~/Box/SCEPTIC_fMRI/')
  pemeta <- read_csv("pe_max_unsmoothed_cluster_metadata.csv")
  pebetas <- read_csv("pe_max_unsmoothed_subj_betas.csv")
} else {
  setwd('~/Box/skinner/projects_analyses/SCEPTIC/fMRI_paper/signals_review/MMClock_aroma_preconvolve_fse_groupfixed/sceptic-clock-feedback-pe_max-preconvolve_fse_groupfixed/pe_max/')
  pemeta <- read_csv("pe_max_cluster_metadata.csv")
  pebetas <- read_csv("pe_max_roi_betas.csv")
}
pemeta$label <- substr(pemeta$label,22,100)
pemeta_overall <- pemeta[pemeta$l2_contrast == 'overall' & pemeta$l3_contrast == 'Intercept' & pemeta$model == 'Intercept-Age',]
pe <- as.tibble(pebetas[pebetas$l2_contrast == 'overall' & pebetas$l3_contrast == 'Intercept' & pebetas$model == 'Intercept-Age',1:3]) %>% filter(cluster_number<8 | cluster_number == 10 | cluster_number == 11)
# head(merge(h,meta))
perois_list <- distinct(pemeta_overall[c(1:7, 10:11),c(5,12)])

# inspect distributions
# for some reason, the clusters do not seem to fit the map -- check with MNH (big clusters broken down?) !!
perois <- inner_join(pe,pemeta_overall)
perois$labeled_cluster <- paste(perois$cluster_number,perois$label)
pemeta_overall$labeled_cluster <- paste(pemeta_overall$cluster_number,pemeta_overall$label)
# ggplot(perois,aes(scale(cope_value))) + geom_histogram() + facet_wrap(~labeled_cluster)

pe_labeled <- inner_join(pe,perois_list)
pe_labeled$labeled_cluster <- paste(pe_labeled$cluster_number,pe_labeled$label)
pe_num <- select(pe_labeled,c(1,2,3))
pe_labeled <- select(pe_labeled,c(1,3,5))

pe_wide <- spread(pe_labeled,labeled_cluster,cope_value)
pe_wide_num <- spread(pe_num,cluster_number,cope_value)

## we are no longer winsorizing betas
# pe_wide <- pe_wide %>% mutate_if(is.double, winsor,trim = .075)
# pe_wide_num <- pe_wide_num %>% mutate_if(is.double, winsor,trim = .075)

pejust_rois <- pe_wide[,2:ncol(pe_wide)]
pejust_rois_num <- pe_wide_num[,2:ncol(pe_wide_num)]


# non-parametric correlations to deal with outliers
peclust_cor <- corr.test(pejust_rois,method = 'pearson', adjust = 'none')
# parametric correlations on winsorised betas
# clust_cor <- cor(just_rois_w,method = 'pearson')

setwd(file.path(clock_folder, 'fmri/keuka_brain_behavior_analyses/'))
if (unsmoothed) {
  pdf("pe_cluster_corr_fixed_unsmoothed.pdf", width=12, height=12)
} else {pdf("pe_cluster_corr_fixed.pdf", width=12, height=12)}
corrplot(peclust_cor$r, cl.lim=c(-1,1),
         method = "circle", tl.cex = 1.5, type = "upper", tl.col = 'black',
         order = "hclust", diag = FALSE,
         addCoef.col="black", addCoefasPercent = FALSE,
         p.mat = peclust_cor$p, sig.level=0.05, insig = "blank")
dev.off()

mpe <- nfactors(peclust_cor$r, n=5, rotate = "oblimin", diagonal = FALSE,fm = "pa", n.obs = 70, SMC = FALSE)
pe.fa = psych::fa(pejust_rois, nfactors=2, rotate = "varimax", fm = "pa")


pe.fa = fa.sort(psych::fa(pejust_rois, nfactors=2))
pefscores <- factor.scores(pejust_rois, pe.fa)$scores
pe_wide$pe_f1_cort_str <- pefscores[,1]
pe_wide$pe_f2_hipp <- pefscores[,2]
pe_wide$pe_PH <- scale(rowMeans(cbind(pe_wide$`10 Left Hippocampus`, pe_wide$`7 Right Hippocampus`)))
pe_wide$pe_PH_l <- scale(pe_wide$`10 Left Hippocampus`)
pe_wide$pe_PH_r <- scale(pe_wide$`7 Right Hippocampus`)

# save loadings
setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/hippo/supp')
pf <- as_tibble(cbind(rownames(unclass(round(pe.fa$Structure, digits = 3))),unclass(round(pe.fa$Structure, digits = 3)))) %>% arrange(desc(MR1) ) %>%  
  rename(region = V1, `Factor 1` = MR1, `Factor 2` = MR2)
pf$region <- gsub('[0-9;]+', '', pf$region)
stargazer(pf, type = 'html', out = 'pe_fa_structure.html', summary = F)

 if (unsmoothed) {
  pe_wide$pe_PH <- (pe_wide$`7 Right Hippocampus` + pe_wide$`10 Left Hippocampus`)/2
  hpe_wide <- inner_join(h_wide,pe_wide[,c("feat_input_id","pe_f1_cort_str", "pe_f2_hipp", "pe_PH")])
} else {  hpe_wide <- inner_join(h_wide,pe_wide[,c("feat_input_id","pe_f1_cort_str", "pe_f2_hipp", "pe_PH", "pe_PH_l", "pe_PH_r")])  }


#####
# add ids
map_df  <- as.tibble(read.csv("~/Box/skinner/projects_analyses/SCEPTIC/fMRI_paper/signals_review/MMClock_aroma_preconvolve_fse_groupfixed/sceptic-clock-feedback-v_entropy-preconvolve_fse_groupfixed/v_entropy/v_entropy-Intercept_design.txt", sep=""))
pc_scores <- inner_join(hpe_wide,map_df[,c(1:2,4:15)])
pc_scores$id <- pc_scores$ID


# get trial_level data
trial_df <- read_csv(file.path(clock_folder, "fmri/data/mmclock_fmri_decay_factorize_selective_psequate_mfx_trial_statistics.csv.gz"))
# trial_df <- read_csv("~/Box/SCEPTIC_fMRI/mmclock_fmri_fixed_uv_mfx_trial_statistics.csv.gz")
# u_df <- read_csv("~/Box/SCEPTIC_fMRI/mmclock_fmri_fixed_uv_ureset_mfx_trial_statistics.csv.gz")
# use fixed-parameter U
u_df <- read_csv("~/Box/SCEPTIC_fMRI/sceptic_model_fits/mmclock_fmri_fixed_uv_ureset_fixedparams_fmri_ffx_trial_statistics.csv.gz")

trial_df <- trial_df %>%
  group_by(id, run) %>%  dplyr::mutate(rt_swing = abs(c(NA, diff(rt_csv))),
                                       rt_swing_lr = abs(log(rt_csv/lag(rt_csv))),
                                       rt_lag = lag(rt_csv) ,
                                       rt_swing_lag = lag(rt_swing),
                                       omission_lag = lag(score_csv==0),
                                       rt_vmax_lag = lag(rt_vmax),
                                       v_chosen_lag = lag(v_chosen),
                                       run_trial=1:50) %>% ungroup() #compute rt_swing within run and subject
u_df <- u_df %>% select(id, run, trial, u_chosen, u_chosen_lag, u_chosen_change, 
                        u_chosen_quantile, u_chosen_quantile_lag, u_chosen_quantile_change,
                        v_chosen_quantile, v_chosen_quantile_lag, v_chosen_quantile_change)

trial_df <- inner_join(trial_df,u_df)

# performance
sum_df <- trial_df %>% group_by(id) %>% dplyr::summarize(total_earnings = sum(score_csv)) %>% arrange(total_earnings)
beta_sum <- inner_join(pc_scores,sum_df)

# model parameters
params <- read_csv(file.path(clock_folder, "fmri/data/mmclock_fmri_decay_factorize_selective_psequate_mfx_sceptic_global_statistics.csv"))
sub_df <- inner_join(beta_sum,params)
sub_df$id <- sub_df$ID
if (unsmoothed) {
  params_beta <- sub_df[,c("h_f1_fp", "h_f2_neg_paralimb","h_HippAntL",
                           "pe_f1_cort_str", "pe_f2_hipp", "pe_PH", "pe_PH_l", "pe_PH_r",
                           "total_earnings", "LL", "alpha", "gamma", "beta")]
} else {
params_beta <- sub_df[,c("h_f1_fp", "h_f2_neg_paralimb","h_HippAntL",
                         "pe_f1_cort_str", "pe_f2_hipp", "pe_PH",
                         "total_earnings", "LL", "alpha", "gamma", "beta")]}
param_cor <- corr.test(params_beta,method = 'pearson', adjust = 'none')

setwd(file.path(clock_folder, 'fmri/keuka_brain_behavior_analyses/'))

# merge into trial-level data
df <- inner_join(trial_df,sub_df, by = "id")
df$rewFunc <- relevel(as.factor(df$rewFunc),ref = "CEV")
df$rewFuncIEVsum <- df$rewFunc
contrasts(df$rewFuncIEVsum) <- contr.sum
colnames(contrasts(df$rewFuncIEVsum)) <- c('CEV','CEVR', 'DEV')

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
                                         trial_neg_inv_sc = scale(-1/run_trial),
                                         v_chosen_change = v_chosen - lag(v_chosen)) %>% ungroup() %>% 
  mutate(rt_lag_sc = scale(rt_lag),
         rt_csv_sc = scale(rt_csv),
         rt_vmax_lag_sc = scale(rt_vmax_lag))

# correlate between-subject V and H with clusters
b_df <- df %>% group_by(id) %>% dplyr::summarise(v_maxB = mean(v_max, na.rm = T),
                                                 v_entropyB = mean(v_entropy, na.rm = T))

sub_df <- inner_join(sub_df, b_df, by = 'id')
if (unsmoothed) {
  bdf <- sub_df[,c("h_f1_fp", "h_f2_neg_paralimb",
                   "pe_f1_cort_str", "pe_f2_hipp", "pe_PH", "pe_PH_l", "pe_PH_r",
                   "total_earnings", "LL", "alpha", "gamma", "beta", "v_maxB", "v_entropyB")]  
} else {
bdf <- sub_df[,c("h_f1_fp", "h_f2_neg_paralimb",
                 "pe_f1_cort_str", "pe_f2_hipp", "pe_PH", "pe_PH_l", "pe_PH_r",
                 "total_earnings", "LL", "alpha", "gamma", "beta", "v_maxB", "v_entropyB")]}
b_cor <- corr.test(bdf,method = 'pearson', adjust = 'none')

setwd(file.path(clock_folder, 'fmri/keuka_brain_behavior_analyses/'))
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
df$h_HippAntL_neg <- -df$h_HippAntL
df$h_HippAntL_resp <- 'low'
df$h_HippAntL_resp[df$h_HippAntL<mean(df$h_HippAntL)] <- 'high'
df$last_outcome <- NA
df$last_outcome[df$omission_lag] <- 'Omission'
df$last_outcome[!df$omission_lag] <- 'Reward'
df$last_outcome <- relevel(as.factor(df$last_outcome), ref = "Reward")

# circular, but just check to what extent each area conforms to SCEPTIC-SM: looks like there is an interaction, both need to be involved again
# ggplot(df, aes(run_trial, v_entropy_wi, color = low_h_paralimbic, lty = h_fp)) + geom_smooth(method = "loess") #+ facet_wrap(~gamma>0)
# # I think this means again that H estimates are more precise for high-paralimbic people, esp. late in learning
# ggplot(df, aes( v_entropy_wi,log(rt_swing), color = low_h_paralimbic, lty = h_fp)) + geom_smooth(method = "gam") + facet_wrap(~learning_epoch)
# ggplot(df, aes( v_max_wi,log(rt_swing), color = low_h_paralimbic, lty = h_fp)) + geom_smooth(method = "gam") + facet_wrap(~run_trial > 10)
# 
# Okay, some behavioral relevance of KLD
# ggplot(df, aes(run_trial, v_entropy_wi, color = k_f1_all_pos_resp)) + geom_smooth(method = "loess")

# add MEG behavioral data
# mtdf = MEG trial df
mtdf <- read_csv("~/Box/SCEPTIC_fMRI/sceptic_model_fits/mmclock_meg_decay_factorize_selective_psequate_fixedparams_meg_ffx_trial_statistics.csv.gz")
mtdf <- mtdf %>%
  group_by(id, run) %>%  dplyr::mutate(rt_swing = abs(c(NA, diff(rt_csv))),
                                       rt_swing_lr = abs(log(rt_csv/lag(rt_csv))),
                                       rt_lag = lag(rt_csv) ,
                                       rt_swing_lag = lag(rt_swing),
                                       omission_lag = lag(score_csv==0),
                                       rt_vmax_lag = lag(rt_vmax),
                                       run_trial=1:63, 
                                       v_max_wi = scale(v_max),
                                       v_max_wi_lag = lag(v_max_wi),
                                       v_entropy_wi = scale(v_entropy),
                                       v_max_b = mean(na.omit(v_max)),
                                       v_entropy_b = mean(na.omit(v_entropy)),
                                       rt_change = rt_csv - rt_lag,
                                       pe_max_lag = lag(pe_max), 
                                       abs_pe_max_lag = abs(pe_max_lag), 
                                       rt_vmax_change = rt_vmax - rt_vmax_lag,
                                       v_chosen_change = v_chosen - lag(v_chosen),
                                       trial_neg_inv_sc = scale(-1/run_trial)) %>% ungroup() %>% 
  mutate(rt_lag_sc = scale(rt_lag),
         rt_csv_sc = scale(rt_csv),
         rt_vmax_lag_sc = scale(rt_vmax_lag),
         id = as.integer(substr(id, 1, 5)))#compute rt_swing within run and subject
# add fMRI betas
mdf <- inner_join(mtdf,sub_df, by = "id")
mdf$rewFunc <- relevel(as.factor(mdf$rewFunc),ref = "CEV")
mdf$rewFuncIEVsum <- mdf$rewFunc
contrasts(mdf$rewFuncIEVsum) <- contr.sum
colnames(contrasts(mdf$rewFuncIEVsum)) <- c('CEV','CEVR', 'DEV')
mdf$h_f1_fp_resp <- 'low'
mdf$h_f1_fp_resp[mdf$h_f1_fp>0] <- 'high'
mdf$pe_f1_cort_str_resp <- 'low'
mdf$pe_f1_cort_str_resp[mdf$pe_f1_cort_str>0] <- 'high'
mdf$pe_f2_hipp_resp <- 'low'
mdf$pe_f2_hipp_resp[mdf$pe_f2_hipp>0] <- 'high'
mdf$h_HippAntL_resp <- 'low'
mdf$h_HippAntL_resp[mdf$h_HippAntL<mean(mdf$h_HippAntL)] <- 'high'
mdf$last_outcome <- NA
mdf$last_outcome[mdf$omission_lag] <- 'Omission'
mdf$last_outcome[!mdf$omission_lag] <- 'Reward'
mdf$last_outcome <- relevel(as.factor(mdf$last_outcome), ref = "Reward")

mdf$learning_epoch <- 'trials 1-10'
mdf$learning_epoch[df$run_trial>10] <- 'trials 11-50'
mdf$h_HippAntL_neg <- -mdf$h_HippAntL

mu_df <- read_csv("~/Box/SCEPTIC_fMRI/sceptic_model_fits/mmclock_meg_fixed_uv_ureset_fixedparams_meg_ffx_trial_statistics.csv.gz")
mu_df <- mu_df %>% select(id, run, trial, u_chosen, u_chosen_lag, u_chosen_change,
                        u_chosen_quantile, u_chosen_quantile_lag, u_chosen_quantile_change) %>% mutate(id = as.integer(substr(id, 1, 5)))
# mu_df <- mu_df %>% select(id, run, trial, u_chosen, u_chosen_lag, u_chosen_change) %>% mutate(
#     id = as.integer(substr(id, 1, 5))) 

mdf <- inner_join(mdf,mu_df, by = c("id", "run", "trial"))


if (unsmoothed) {
  save(file = 'trial_df_and_vh_pe_clusters_u_unsmoothed.Rdata', df, mdf)
} else {save(file = 'trial_df_and_vh_pe_clusters_u.Rdata', df, mdf)}


