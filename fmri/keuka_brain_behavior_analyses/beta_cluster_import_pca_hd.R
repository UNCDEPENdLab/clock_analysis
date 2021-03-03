# preprocess voxel-wise fMRI data for brain-to-behavior analyses
# extracts hippocampal activation coefficients for prediction error and entropy analyses
# merges with behavioral data from the fMRI and replication (MEG sessions)

library(dplyr)
library(tidyverse)
library(psych)
library(corrplot)
library(lme4)
library(stargazer)
# detach(package:MASS)
# detach(package:plyr)
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
  setwd('~/Box/SCEPTIC_fMRI/MMClock_aroma_preconvolve_fse_groupfixed/sceptic-clock-feedback-v_entropy-preconvolve_fse_groupfixed/v_entropy/')
  meta <- read_csv("~/Box/SCEPTIC_fMRI/MMClock_aroma_preconvolve_fse_groupfixed/sceptic-clock-feedback-v_entropy-preconvolve_fse_groupfixed/v_entropy/v_entropy_cluster_metadata.csv")
  Hbetas <- read_csv("v_entropy_roi_betas.csv")
}
meta$label <- substr(meta$label,22,100)
meta_overall <- meta[meta$l2_contrast == 'overall' & meta$l3_contrast == 'Intercept' & meta$model == 'Intercept-Age',]
h <- as_tibble(Hbetas[Hbetas$l2_contrast == 'overall' & Hbetas$l3_contrast == 'Intercept' & Hbetas$model == 'Intercept-Age',1:3]) %>% filter(cluster_number<12)

# head(merge(h,meta))
rois <- distinct(meta_overall[,c(5,6,12)])
# DAN: ROIs 1 (BL parietal), 2 (R SFG), 6 (L SFG)

# inspect distributions
hrois <- inner_join(h,meta_overall)
hrois$labeled_cluster <- paste(hrois$cluster_number,hrois$label)
hrois <- hrois %>% select(-c(model, l1_contrast, l2_contrast, l3_contrast, cluster_size, cluster_number, cluster_threshold,
                             z_threshold, x, y, z, label))
# ggplot(hrois,aes(scale(cope_value))) + geom_histogram() + facet_wrap(~labeled_cluster)
h_wide <- spread(hrois,labeled_cluster,cope_value) 
# head(h_wide)
# with group-fixed parameters we don't seem to need the winsorization step!!!
# h_wide <- h_wide %>% mutate_if(is.double, winsor,trim = .075)

just_rois <- h_wide[,2:12]
# winsorize to deal with beta ouliers

# non-parametric correlations to deal with outliers
clust_cor <- cor(just_rois,method = 'pearson')
# parametric correlations on winsorised betas
# clust_cor <- cor(just_rois_w,method = 'pearson')

setwd(file.path(clock_folder, 'fmri/keuka_brain_behavior_analyses/dan'))
if (unsmoothed) {
  pdf("h_cluster_corr_fixed_unsmoothed.pdf", width=12, height=12)  
} else {
  pdf("h_cluster_corr_fixed.pdf", width=12, height=12)  
}
corrplot(clust_cor, cl.lim=c(-1,1),
         method = "circle", tl.cex = 1.5, type = "upper", tl.col = 'black',
         order = "hclust", diag = FALSE,
         addCoef.col="black", addCoefasPercent = FALSE,
         p.mat = 1-clust_cor, sig.level=0.75, insig = "blank")
dev.off()
# mh <- nfactors(clust_cor, n=5, rotate = "oblimin", diagonal = FALSE,fm = "pa", n.obs = 70, SMC = FALSE)
h.fa = fa.sort(psych::fa(just_rois, nfactors=2, rotate = "varimax", fm = "pa"))
fscores <- factor.scores(just_rois, h.fa)$scores
h_wide$h_f1_fp <- fscores[,1]
h_wide$h_f2_neg_paralimb <- fscores[,2]
h_wide$h_HippAntL <- h_wide$`9 Left Hippocampus`

# high-entropy factor = DAN + cerebellum + rlPFC/dlPFC

h_wide <- h_wide %>% mutate(DAN = 
                              rowMeans(cbind(`1 Right Precuneus`, `2 Right Superior Frontal Gyrus`, `6 Left Precentral Gyrus`)),
                            dan_parietal = `1 Right Precuneus`,
                            dan_r_sfg = `2 Right Superior Frontal Gyrus`,
                            dan_l_sfg = `6 Left Precentral Gyrus`)
if (unsmoothed) {
  h_wide$h_vmPFC <- h_wide$`5 Left Mid Orbital Gyrus`
} else {h_wide$h_vmPFC <- h_wide$`5 Left Mid Orbital Gyrus`}
hf <- as_tibble(cbind(rownames(unclass(round(h.fa$Structure, digits = 3))),unclass(round(h.fa$Structure, digits = 3)))) %>% arrange(desc(PA1) ) %>%  
  rename(region = V1, `Factor 1` = PA1, `Factor 2` = PA2)
hf$region <- gsub('[0-9;]+', '', hf$region)
stargazer(hf, type = 'html', out = '~/OneDrive/collected_letters/papers/sceptic_fmri/hippo/supp/h_fa_structure.html', summary = F)


h_wide <- subset(h_wide, select = c("feat_input_id","h_f1_fp","h_f2_neg_paralimb", "h_HippAntL", "DAN", "dan_parietal", "dan_l_sfg", "dan_r_sfg", "7 Right Middle Frontal Gyrus"))


#####
# add Schaefer parcel-based DAN entropy betas
schaefer <- as_tibble(read.csv('~/Box/SCEPTIC_fMRI/MMClock_aroma_preconvolve_fse_groupfixed/sceptic-clock-feedback-v_entropy-preconvolve_fse_groupfixed/v_entropy/entropy_beta_bifactor_fscores.csv'))
schaefer$feat_input_id <- schaefer$numid
h_wide <- inner_join(h_wide, schaefer %>% select(-numid)) %>% mutate(entropy_vlPFC = `7 Right Middle Frontal Gyrus`) %>% dplyr::select(-`7 Right Middle Frontal Gyrus`)
#####
# add PE
if (unsmoothed) {
  setwd('~/Box/SCEPTIC_fMRI/')
  pemeta <- read_csv("pe_max_unsmoothed_cluster_metadata.csv")
  pebetas <- read_csv("pe_max_unsmoothed_subj_betas.csv")
} else {
  setwd('~/Box/SCEPTIC_fMRI/MMClock_aroma_preconvolve_fse_groupfixed/sceptic-clock-feedback-pe_max-preconvolve_fse_groupfixed/pe_max/')
  # pemeta <- read_csv("pe_max_cluster_metadata.csv")
  pemeta <- read_delim("pe_clusters_z4p41_labels.txt", 
                       "\t", escape_double = FALSE, trim_ws = TRUE)
  pebetas <- read_csv("pe_max_betas.csv.gz")
}
# pemeta$label <- substr(pemeta$label,22,100)
# pemeta_overall <- pemeta[pemeta$l2_contrast == 'overall' & pemeta$l3_contrast == 'Intercept' & pemeta$model == 'Intercept-Age',]
pemeta$atlas_value <- pemeta$roi
pemeta_overall <- pemeta
pe <- as_tibble(pebetas[pebetas$l2_contrast == 'overall' ,1:5]) 
# pe <- as_tibble(pebetas[pebetas$l2_contrast == 'overall' & pebetas$l3_contrast == 'Intercept' & pebetas$model == 'Intercept-Age',1:3]) %>% filter(cluster_number<8 | cluster_number == 10 | cluster_number == 11)

# head(merge(h,meta))
# perois_list <- distinct(pemeta_overall[c(1:7, 10:11),c(5,12)])
perois_list <- distinct(pemeta_overall)
pemeta_overall$atlas_value <- pemeta_overall$roi

# inspect distributions
perois <- inner_join(pe,pemeta_overall)
# perois$labeled_cluster <- paste(perois$cluster_number,perois$label)
perois$labeled_cluster <- paste(perois$roi,perois$label)
# pemeta_overall$labeled_cluster <- paste(pemeta_overall$cluster_number,pemeta_overall$label)
pemeta_overall$labeled_cluster <- paste(pemeta_overall$roi,pemeta_overall$label)
# ggplot(perois,aes(scale(cope_value))) + geom_histogram() + facet_wrap(~labeled_cluster)

pe_labeled <- inner_join(pe,perois_list)
pe_labeled$labeled_cluster <- paste(pe_labeled$roi,pe_labeled$label)
# pe_labeled$labeled_cluster <- paste(pe_labeled$cluster_number,pe_labeled$label)
pe_num <- select(pe_labeled,c("id", "beta", "numid", "roi"))
pe_labeled <- select(pe_labeled,c("numid", "id", "beta", "labeled_cluster"))

# pe_wide <- spread(pe_labeled,labeled_cluster,cope_value)
pe_wide <- spread(pe_labeled,labeled_cluster,beta)

# pe_wide_num <- spread(pe_num,cluster_number,cope_value)
pe_wide_num <- spread(pe_num,roi,beta)

## we are no longer winsorizing betas
# pe_wide <- pe_wide %>% mutate_if(is.double, winsor,trim = .075)
# pe_wide_num <- pe_wide_num %>% mutate_if(is.double, winsor,trim = .075)

pejust_rois <- pe_wide[,3:ncol(pe_wide)]
# pejust_rois <- pe_wide[,2:ncol(pe_wide)]
pejust_rois_num <- pe_wide_num[,3:ncol(pe_wide_num)]


# non-parametric correlations to deal with outliers
peclust_cor <- corr.test(pejust_rois,method = 'pearson', adjust = 'none')
# parametric correlations on winsorised betas
# clust_cor <- cor(just_rois_w,method = 'pearson')

setwd(file.path(clock_folder, 'fmri/keuka_brain_behavior_analyses/'))
if (unsmoothed) {
  pdf("pe_cluster_corr_fixed_unsmoothed_broken_up.pdf", width=12, height=12)
} else {pdf("pe_cluster_corr_fixed_broken_up.pdf", width=12, height=12)}
corrplot(peclust_cor$r, cl.lim=c(-1,1),
         method = "circle", tl.cex = 1.5, type = "upper", tl.col = 'black',
         order = "hclust", diag = FALSE,
         addCoef.col="black", addCoefasPercent = FALSE,
         p.mat = peclust_cor$p, sig.level=0.05, insig = "blank")
dev.off()

mpe <- nfactors(peclust_cor$r, n=5, rotate = "oblimin", diagonal = FALSE,fm = "pa", n.obs = 70, SMC = FALSE)
pe.fa = psych::fa(pejust_rois, nfactors=3, rotate = "oblimin", fm = "pa")
pe.faba = psych::bassAckward(pejust_rois, nfactors=6, rotate = "oblimin", fm = "pa")


pe.fa = fa.sort(psych::fa(pejust_rois, nfactors=3))
pefscores <- factor.scores(pejust_rois, pe.fa)$scores
pe_wide$pe_f1_cort_hipp <- pefscores[,1]
pe_wide$pe_f2_cerebell <- pefscores[,2]
pe_wide$pe_f3_str <- pefscores[,3]
# pe_wide$pe_PH <- scale(rowMeans(cbind(pe_wide$`10 Left Hippocampus`, pe_wide$`7 Right Hippocampus`)))
# pe_wide$pe_PH_l <- scale(pe_wide$`10 Left Hippocampus`)
pe_wide$pe_PH_r <- scale(pe_wide$`18 Right Hippocampus`)
pe_wide$pe_ips <- scale(rowMeans(cbind(pe_wide$`1 Right Angular Gyrus`, pe_wide$`3 Left Inferior Parietal Lobule` )))
# pe_wide$pe_str <- scale(pe_wide$`2 ; Left Caudate Nucleus`)
# add 
pe_wide$feat_input_id <- pe_wide$numid
# save loadings
setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/hippo/supp')
pf <- as_tibble(cbind(rownames(unclass(round(pe.fa$Structure, digits = 3))),unclass(round(pe.fa$Structure, digits = 3)))) %>% arrange(desc(MR1) ) %>%  
  rename(region = V1, `Factor 1` = MR1, `Factor 2` = MR2)
pf$region <- gsub('[0-9;]+', '', pf$region)
stargazer(pf, type = 'html', out = 'pe_fa_structure.html', summary = F)

if (unsmoothed) {
  pe_wide$pe_PH <- (pe_wide$`7 Right Hippocampus` + pe_wide$`10 Left Hippocampus`)/2
  hpe_wide <- inner_join(h_wide,pe_wide[,c("feat_input_id","pe_f1_cort_hipp", "pe_f2_cerebell", "pe_f3_str", "pe_PH_r", "pe_ips")])
} else {  hpe_wide <- inner_join(h_wide,pe_wide[,c("feat_input_id","pe_f1_cort_hipp", "pe_f2_cerebell", "pe_f3_str", "pe_PH_r", "pe_ips")])  }

#####
# add entropy change betas, dh, for "delta H"
setwd('~/Box/skinner/projects_analyses/SCEPTIC/fMRI_paper/mmc/signals_review/compiled_outputs/change_betas/sceptic-clock-feedback-v_entropy_change-preconvolve_fse_groupfixed/v_entropy_change')
hdmeta <- read_csv("v_entropy_change_cluster_metadata.csv")
hdbetas <- read_csv("v_entropy_change_roi_betas.csv")

hdmeta$label <- substr(hdmeta$label,22,100)
hdmeta_overall <- hdmeta[hdmeta$l2_contrast == 'overall' & hdmeta$l3_contrast == 'Intercept' & hdmeta$model == 'Intercept-Age',]
hd <- as_tibble(hdbetas[hdbetas$l2_contrast == 'overall' & hdbetas$l3_contrast == 'Intercept' & hdbetas$model == 'Intercept-Age',1:3]) %>% 
  filter(cluster_number<12 & cluster_number != 7 & cluster_number != 10)

# head(merge(h,hdmeta))
rois <- distinct(hdmeta_overall[,c(5,6,12)])
# DAN: ROIs 1 (BL parietal), 2 (R SFG), 6 (L SFG)

# inspect distributions
hdrois <- inner_join(hd,hdmeta_overall)
hdrois$labeled_cluster <- paste(hdrois$cluster_number,hdrois$label)
hdrois <- hdrois %>% select(-c(model, l1_contrast, l2_contrast, l3_contrast, cluster_size, cluster_number, cluster_threshold,
                               z_threshold, x, y, z, label))
# ggplot(hdrois,aes(scale(cope_value))) + geom_histogram() + facet_wrap(~labeled_cluster)
hd_wide <- spread(hdrois,labeled_cluster,cope_value) 
# head(h_wide)
# with group-fixed parameters we don't seem to need the winsorization step!!!
# h_wide <- h_wide %>% mutate_if(is.double, winsor,trim = .075)

just_rois <- hd_wide[,2:10]
# winsorize to deal with beta ouliers

# non-parametric correlations to deal with outliers
clust_cor <- cor(just_rois,method = 'pearson')
# parametric correlations on winsorised betas
# clust_cor <- cor(just_rois_w,method = 'pearson')

setwd(file.path(clock_folder, 'fmri/keuka_brain_behavior_analyses/'))
pdf("hd_cluster_corr_fixed.pdf", width=12, height=12)  
corrplot(clust_cor, cl.lim=c(-1,1),
         method = "circle", tl.cex = 1.5, type = "upper", tl.col = 'black',
         order = "hclust", diag = FALSE,
         addCoef.col="black", addCoefasPercent = FALSE,
         p.mat = 1-clust_cor, sig.level=0.75, insig = "blank")
dev.off()
mhd <- nfactors(clust_cor, n=5, rotate = "oblimin", diagonal = FALSE,fm = "pa", n.obs = 70, SMC = FALSE)
hd.fa = fa.sort(psych::fa(just_rois, nfactors=3, rotate = "varimax", fm = "pa"))
fscores <- factor.scores(just_rois, hd.fa)$scores

hd_wide$hd_f1_salience <- pefscores[,1]
hd_wide$hd_f2_DAN <- pefscores[,2]
hd_wide$hd_f3_rlPFC <- pefscores[,3]
hd_wide$hd_DAN_mean <- scale(rowMeans(cbind(hd_wide$`1 Right Precuneus`, hd_wide$`3 Right Middle Frontal Gyrus`, hd_wide$`5 Left Superior Frontal Gyrus`)))
hpehd_wide <- inner_join(hpe_wide,hd_wide[,c("feat_input_id","hd_f1_salience", "hd_f2_DAN", "hd_f3_rlPFC", "hd_DAN_mean")])  
#####
# add ids
map_df  <- as_tibble(read.csv("~/Box/SCEPTIC_fMRI/MMClock_aroma_preconvolve_fse_groupfixed/sceptic-clock-feedback-v_entropy-preconvolve_fse_groupfixed/v_entropy/v_entropy-Intercept_design.txt", sep=""))
pc_scores <- inner_join(hpehd_wide,map_df[,c(1:2,4:15)])
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
                                       rt_lag = lag(rt_csv),
                                       rt_lag2 = lag(rt_csv,2),
                                       rt_lag3 = lag(rt_csv,3),
                                       rt_swing_lag = lag(rt_swing),
                                       omission_lag = lag(score_csv==0),
                                       omission_lag2 = lag(score_csv==0, 2),
                                       omission_lag3 = lag(score_csv==0, 3),
                                       rt_vmax_lag = lag(rt_vmax),
                                       rt_vmax_lag2 = lag(rt_vmax_lag),
                                       v_chosen_lag = lag(v_chosen),
                                       run_trial=1:50) %>% ungroup() %>%
  dplyr::mutate(rt_lag2_sc = scale(rt_lag2),
                rt_lag3_sc = scale(rt_lag3),
                rt_vmax_lag2_sc = scale(rt_vmax_lag2))
#compute rt_swing within run and subject
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
                           "DAN", "dan_r_sfg", "dan_l_sfg", "dan_parietal",
                           "pe_f1_cort_hipp", "pe_f2_cerebell", "pe_f3_str", "pe_PH_r", "pe_ips",
                           "total_earnings", "LL", "alpha", "gamma", "beta")]
} else {
  params_beta <- sub_df[,c("h_f1_fp", "h_f2_neg_paralimb","h_HippAntL",
                           "DAN", "dan_r_sfg", "dan_l_sfg", "dan_parietal",
                           "pe_f1_cort_hipp", "pe_f2_cerebell", "pe_f3_str", "pe_PH_r", "pe_ips",
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

get_kldsum <- function(v1, v2) {
  require(LaplacesDemon)
  stopifnot(length(v1) == length(v2))
  if (any(is.na(v1)) || any(is.na(v2))) { return(NA_real_) }
  kk <- KLD(v1, v2)
  return(kk$sum.KLD.px.py)
}
df <- df %>% group_by(ID, run) %>% arrange(ID, run, run_trial) %>% mutate(
  rt_lag2 = lag(rt_lag),
  rt_lag3 = lag(rt_lag2),
  rt_lag4 = lag(rt_lag3),
  rt_lag5 = lag(rt_lag4),
  rt_swing_lag2 = lag(rt_swing_lag),
  omission_lag4 = lag(omission_lag3),
  omission_lag5 = lag(omission_lag4)) %>% ungroup() %>%
  rowwise() %>% mutate(
    kld4 = get_kldsum(c(rt_lag4, rt_lag3, rt_lag2, rt_lag), c(rt_lag5, rt_lag4, rt_lag3, rt_lag2)),
    kld3 = get_kldsum(c(rt_lag3, rt_lag2, rt_lag), c(rt_lag4, rt_lag3, rt_lag2)),
    kld_rew3 = get_kldsum(c(omission_lag3, omission_lag2, omission_lag), c(omission_lag4, omission_lag3, omission_lag2)),
    kld_rew4 = get_kldsum(c(omission_lag4, omission_lag3, omission_lag2, omission_lag), c(omission_lag5, omission_lag4, omission_lag3, omission_lag2))) %>%
  ungroup() %>% group_by(ID, run) %>% mutate(kld3_lag = lag(kld3),
                                             kld4_lag = lag(kld4),
                                             kld3_cum2 = kld3 + kld3_lag,
                                             kld4_cum2 = kld4 + kld4_lag
  ) %>%
  ungroup() %>% mutate(rt_swing_lag_sc = scale(rt_swing_lag),
                       rt_swing_lag2_sc = scale(rt_swing_lag2))

df$last_outcome <- NA
df$last_outcome[df$omission_lag] <- 'Omission'
df$last_outcome[!df$omission_lag] <- 'Reward'
df$last_outcome <- relevel(as.factor(df$last_outcome), ref = "Reward")
df <- df %>% group_by(id, run) %>%  dplyr::mutate(last_outcome_lag2 = lag(last_outcome),
                                                  last_outcome_lag3 = lag(last_outcome,2)) %>% ungroup()

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
                                       rt_lag2 = lag(rt_csv,2),
                                       rt_lag3 = lag(rt_csv,3),
                                       rt_swing_lag = lag(rt_swing),
                                       omission_lag = lag(score_csv==0),
                                       omission_lag2 = lag(score_csv==0, 2),
                                       omission_lag3 = lag(score_csv==0, 3),
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
         id = as.integer(substr(id, 1, 5)),
         rt_lag2_sc = scale(rt_lag2),
         rt_lag3_sc = scale(rt_lag3))#compute rt_swing within run and subject
# add fMRI betas
mdf <- inner_join(mtdf,sub_df, by = "id")
mdf$rewFunc <- relevel(as.factor(mdf$rewFunc),ref = "CEV")
mdf$rewFuncIEVsum <- mdf$rewFunc
contrasts(mdf$rewFuncIEVsum) <- contr.sum
colnames(contrasts(mdf$rewFuncIEVsum)) <- c('CEV','CEVR', 'DEV')
mdf$last_outcome <- NA
mdf$last_outcome[mdf$omission_lag] <- 'Omission'
mdf$last_outcome[!mdf$omission_lag] <- 'Reward'
mdf$last_outcome <- relevel(as.factor(mdf$last_outcome), ref = "Reward")
mdf <- mdf %>% group_by(id, run) %>%  dplyr::mutate(last_outcome_lag2 = lag(last_outcome),
                                                    last_outcome_lag3 = lag(last_outcome,2)) %>% ungroup()

mdf$learning_epoch <- 'trials 1-10'
mdf$learning_epoch[df$run_trial>10] <- 'trials 11-50'
mdf$h_HippAntL_neg <- -mdf$h_HippAntL

mu_df <- read_csv("~/Box/SCEPTIC_fMRI/sceptic_model_fits/mmclock_meg_fixed_uv_ureset_fixedparams_meg_ffx_trial_statistics.csv.gz")
mu_df <- mu_df %>% select(id, run, trial, u_chosen, u_chosen_lag, u_chosen_change,
                          u_chosen_quantile, u_chosen_quantile_lag, u_chosen_quantile_change) %>% mutate(id = as.integer(substr(id, 1, 5)))
# mu_df <- mu_df %>% select(id, run, trial, u_chosen, u_chosen_lag, u_chosen_change) %>% mutate(
#     id = as.integer(substr(id, 1, 5))) 

mdf <- inner_join(mdf,mu_df, by = c("id", "run", "trial"))
## dichotomize betas
df <- df %>% mutate(pe_ips_resp = case_when(
  pe_ips < median(pe_ips) ~ 'low_pe_ips',
  pe_ips > median(pe_ips) ~ 'high_pe_ips'),
  dan_h_resp = case_when(
    DAN < median(DAN) ~ 'low_h_dan',
    DAN > median(DAN) ~ 'high_h_dan'),
  dan_general_entropy_resp = case_when(
    general_entropy < median(general_entropy) ~ 'low_general_entropy_dan',
    general_entropy > median(general_entropy) ~ 'high_general_entropy_dan'),
  pe_f1_cort_hipp_resp = case_when(
    pe_f1_cort_hipp < median(pe_f1_cort_hipp) ~ 'low_pe_f1_cort_hipp',
    pe_f1_cort_hipp > median(pe_f1_cort_hipp) ~ 'high_pe_f1_cort_hipp'),
  pe_f3_str_resp  = case_when(
    pe_f3_str < median(pe_f3_str) ~ 'low_pe_f3_str',
    pe_f3_str > median(pe_f3_str) ~ 'high_pe_f3_str'),
)
mdf <- mdf %>% mutate(pe_ips_resp = case_when(
  pe_ips < median(pe_ips) ~ 'low_pe_ips',
  pe_ips > median(pe_ips) ~ 'high_pe_ips'),
  dan_h_resp = case_when(
    DAN < median(DAN) ~ 'low_h_dan',
    DAN > median(DAN) ~ 'high_h_dan'),
  dan_general_entropy_resp = case_when(
    general_entropy < median(general_entropy) ~ 'low_general_entropy_dan',
    general_entropy > median(general_entropy) ~ 'high_general_entropy_dan'),
  pe_f1_cort_hipp_resp = case_when(
    pe_f1_cort_hipp < median(pe_f1_cort_hipp) ~ 'low_pe_f1_cort_hipp',
    pe_f1_cort_hipp > median(pe_f1_cort_hipp) ~ 'high_pe_f1_cort_hipp'),
  pe_f3_str_resp  = case_when(
    pe_f3_str < median(pe_f3_str) ~ 'low_pe_f3_str',
    pe_f3_str > median(pe_f3_str) ~ 'high_pe_f3_str'),
)
# checked box plots -- looks good


if (unsmoothed) {
  save(file = 'trial_df_and_vh_pe_hd_clusters_u_unsmoothed.Rdata', df, mdf)
} else {save(file = 'trial_df_and_vh_pe_hd_clusters_u.Rdata', df, mdf)}


