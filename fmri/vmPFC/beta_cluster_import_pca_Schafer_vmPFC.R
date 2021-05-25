# preprocess voxel-wise fMRI data for brain-to-behavior analyses
# extracts vmPFC activation coefficients for prediction error and entropy analyses
# merges with behavioral data from the fMRI and replication (MEG sessions)

library(dplyr)
library(tidyverse)
library(psych)
library(corrplot)
library(lme4)
library(stargazer)

# clock_folder <- "~/code/clock_analysis" #alex
# clock_folder <- "~/Data_Analysis/clock_analysis" #michael
clock_folder <- "/Users/andypapale/vmPFC_clock_analysis/clock_analysis" # Andrew

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

# 
# if (unsmoothed) {
#   setwd('~/Box/SCEPTIC_fMRI/')
#   meta <- read_csv("~/Box/SCEPTIC_fMRI/v_entropy_unsmoothed_cluster_metadata.csv")
#   Hbetas <- read_csv("v_entropy_unsmoothed_subj_betas.csv")
# } else {
#setwd('/Users/andypapale/vmPFC_clock_analysis/atlas_betas_sep2020/sceptic-clock-feedback-v_chosen-preconvolve_fse_groupfixed/clock')
#setwd("~/vmPFC/DropOut-Weighted-Betas-2021-01-29/sceptic-clock-feedback-v_chosen-preconvolve_fse_groupfixed/feedback/")
#setwd("~/vmPFC/DropOut-Weighted-Betas-2021-01-29/sceptic-clock-feedback-v_entropy-preconvolve_fse_groupfixed/feedback/")
#setwd("~/vmPFC/DropOut-Weighted-Betas-2021-01-29/sceptic-clock-feedback-pe_max-preconvolve_fse_groupfixed/feedback/")
#   meta <- read_csv("~/Box/SCEPTIC_fMRI/MMClock_aroma_preconvolve_fse_groupfixed/sceptic-clock-feedback-v_chosen-preconvolve_fse_groupfixed/v_entropy/")

#setwd('~/vmPFC/DropOut-Weighted-Betas-2021-03-04/sceptic-clock-feedback-pe_max-preconvolve_fse_groupfixed/pe_max/')
#setwd('~/vmPFC/DropOut-Weighted-Betas-2021-03-04/sceptic-clock-feedback-rt_vmax_change-preconvolve_fse_groupfixed/rt_vmax_change/')
#setwd('~/vmPFC/DropOut-Weighted-Betas-2021-03-04/sceptic-clock-feedback-v_chosen-preconvolve_fse_groupfixed/v_chosen/')
setwd('~/vmPFC/DropOut-Weighted-Betas-2021-03-04/sceptic-clock-feedback-v_entropy-preconvolve_fse_groupfixed/v_entropy/')
Hbetas <- read_csv("v_entropy_Schaeffer_2018_vmPFC_mask_betas.csv.gz")
# }
# meta$label <- substr(meta$label,22,100)
# meta_overall <- meta[meta$l2_contrast == 'overall' & meta$l3_contrast == 'Intercept' & meta$model == 'Intercept-Age',]
#h <- as_tibble(Hbetas %>% filter(atlas_name=="vmpfc_ba_max_2.3mm_filled.nii.gz" & l2_contrast=="overall"))
#h <- as_tibble(Hbetas %>% filter(atlas_name=="Schaefer_136_2.3mm.nii.gz" & l2_contrast=="overall") %>% filter(atlas_value %in% c(55,56,65,66,84,86,88,89,159,160,161,170,171,191,192,194)))
h <- as_tibble(Hbetas %>% filter(l2_contrast=="overall") %>% filter(atlas_value %in% c(55,56,65,66,67,84,86,88,89,159,160,161,170,171,191,192,194)))

   
# head(merge(h,meta))
# rois <- distinct(meta_overall[,c(5,6,12)])
# DAN: ROIs 1 (BL parietal), 2 (R SFG), 6 (L SFG)

# inspect distributions
#hrois <- inner_join(h,meta_overall)
#hrois$labeled_cluster <- paste(hrois$cluster_number,hrois$label)
#hrois <- hrois %>% select(-c(model, l1_contrast, l2_contrast, l3_contrast, cluster_size, cluster_number, cluster_threshold,
#                             z_threshold, x, y, z, label))
# ggplot(hrois,aes(scale(cope_value))) + geom_histogram() + facet_wrap(~labeled_cluster)
#h_wide <- spread(hrois,labeled_cluster,cope_value) 
region0 <- c("11mL","11mR","14cL","14cR","14mL","14mR","14rL","14rR","14rrL","14rrR","24L","24R","25L","25R","32L","32R")

# 1 11mL
# 2 11mR
# 3 14cL
# 4 14cR
# 5 14mL
# 6 14mR
# 7 14rL
# 8 14rR
# 9 14rrL
# 10 14rrR
# 11 24L
# 12 24R
# 13 25L
# 14 25R
# 15 32L
# 16 32R


h <- h %>% mutate(Schaefer_label = case_when(
  atlas_value  == 159 ~ "pfc14r14c11mR",
  atlas_value  == 191 ~ "pfc14m3225R",
  atlas_value  == 160 ~ "pfc11aR",
  atlas_value  == 170 ~ "pfc11b47R",
  atlas_value  == 161 ~ "pfcfp10R",
  atlas_value  == 192 ~ "pfc2432R",
  atlas_value  == 194 ~ "pfcd10R",
  atlas_value  == 171 ~ "pfcl10R",
  atlas_value  == 56 ~ "pfc14r14c11mL",  # limbic
  atlas_value  == 84 ~ "pfc14m3225L",
  atlas_value  == 86 ~ "pfcfp10L",
  atlas_value  == 88 ~ "pfc2432L",
  atlas_value  == 89 ~ "pfcrdb10L",
  atlas_value  == 55 ~ "pfc11m13L",   # limbic
  atlas_value  == 65 ~ "pfc11L",
  atlas_value  == 66 ~ "pfc47L",
  atlas_value  == 67 ~ "pfcpfc14rc2L"
  #atlas_value  == 15 ~ region0[15],
  #atlas_value  == 16 ~ region0[16],
  ), side = substr(Schaefer_label, nchar(Schaefer_label), nchar(Schaefer_label)),
  region_S = substr(Schaefer_label, 1, nchar(Schaefer_label)-1)
  )

# winsorize all betas across regions to remove susceptibility outliers
#h <- h %>% mutate_at(c("beta"), winsor,trim = .075)

#region0 <- c("11mL","11mR","14cL","14cR","14mL","14mR","14rL","14rR","14rrL","14rrR","24L","24R","25L","25R","32L","32R")
#region <- matrix(rep(t(region0),70))
#region <- as.character((region))

#h <- cbind(h,region)

h_wide <- pivot_wider(h, names_from=Schaefer_label, values_from=beta, id_cols = id) # 2020-09-30 AndyP

# h_wide <- h_wide %>% select(-c(atlas_name,fsldir,l1_contrast,l2_contrast))
# head(h_wide)
# with group-fixed parameters we don't seem to need the winsorization step!!!
# h_wide <- h_wide %>% mutate_if(is.double, winsor,trim = .075)

just_rois <- h_wide[,3:ncol(h_wide)]
# winsorize to deal with beta outliers


ggplot(h, aes(Schaefer_label,beta))+geom_boxplot()
pdf("boxplot1.pdf", width = 12, height = 30)
dev.off()
ggplot(h, aes(beta))+geom_histogram()+facet_wrap(~Schaefer_label)

# non-parametric correlations to deal with outliers
clust_cor <- cor(just_rois,method = 'pearson',use="complete.obs")
# parametric correlations on winsorised betas
# clust_cor <- cor(just_rois_w,method = 'pearson')

#setwd(file.path(clock_folder, 'fmri/keuka_brain_behavior_analyses/'))
if (unsmoothed) {
  pdf("h_cluster_corr_fixed_unsmoothed.pdf", width=12, height=12)  
} else {
  pdf("h_cluster_corr_fixed_vmPFC.pdf", width=12, height=12)  
}
# corrplot(clust_cor, cl.lim=c(-1,1),
#          method = "circle", tl.cex = 1.5, type = "upper", tl.col = 'black',
#          order = "hclust", diag = FALSE,
#          addCoef.col="black", addCoefasPercent = FALSE,
#          p.mat = 1-clust_cor, sig.level=0.75, insig = "blank")
 
#colnames(clust_cor) <- region0
#rownames(clust_cor) <- region0
corrplot(clust_cor, cl.lim=c(-1,1),
          method = "circle", tl.cex = 1.5, type = "upper", tl.col = 'black', diag = FALSE,
          order = 'hclust',addCoef.col="black", addCoefasPercent = FALSE,
          p.mat = 1-clust_cor, sig.level=0.75, insig = "blank")



dev.off()

# test against 0
# first visualize
ggplot(h, aes(beta, Schaefer_label)) + geom_boxplot() + geom_vline(xintercept = 0)

summary(m1 <- lmer(beta ~ region_S * side + (1|id), h))
# plot model-predicted means by region
em1 <- as_tibble(emmeans::emmeans(m1, c("region_S", "side")))
em1$beta <- em1$emmean
ggplot(em1, aes(region_S, beta, color = side)) + geom_point() + geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL))
car::Anova(m1, '3')
 mh <- nfactors(clust_cor, n=5, rotate = "oblimin", diagonal = FALSE,fm = "pa", n.obs = 70, SMC = FALSE)
h.fa = fa.sort(psych::fa(just_rois, nfactors=2, rotate = "oblimin", fm = "gls"))  # two factors are vmPFC and central OFC
#$fscores <- factor.scores(just_rois, h.fa)$scores
#h_wide$h_f1_fp <- fscores[,1]
#h_wide$h_f2_neg_paralimb <- fscores[,2]
#h_wide$h_HippAntL <- h_wide$`9 Left Hippocampus`

# high-entropy factor = DAN + cerebellum + rlPFC/dlPFC

# rename atlas_vals to vmPFC MP areas here
# h_wide <- h_wide %>% mutate(DAN = 
#                               rowMeans(cbind(`1 Right Precuneus`, `2 Right Superior Frontal Gyrus`, `6 Left Precentral Gyrus`)),
#                             dan_parietal = `1 Right Precuneus`,
#                             dan_r_sfg = `2 Right Superior Frontal Gyrus`,
#                             dan_l_sfg = `6 Left Precentral Gyrus`)
# if (unsmoothed) {
#   h_wide$h_vmPFC <- h_wide$`5 Left Mid Orbital Gyrus`
# } else {h_wide$h_vmPFC <- h_wide$`6`}
# hf <- as_tibble(cbind(rownames(unclass(round(h.fa$Structure, digits = 3))),unclass(round(h.fa$Structure, digits = 3)))) %>% arrange(desc(PA1) ) %>%  
#   rename(region = V1, `Factor 1` = PA1, `Factor 2` = PA2)
# hf$region <- gsub('[0-9;]+', '', hf$region)
# stargazer(hf, type = 'html', out = '~/OneDrive/collected_letters/papers/sceptic_fmri/hippo/supp/h_fa_structure.html', summary = F)
# 
# 
# h_wide <- subset(h_wide, select = c("feat_input_id","h_f1_fp","h_f2_neg_paralimb", "h_HippAntL", "DAN", "dan_parietal", "dan_l_sfg", "dan_r_sfg"))

#####
# add PE
# if (unsmoothed) {
#   setwd('~/Box/SCEPTIC_fMRI/')
#   pemeta <- read_csv("pe_max_unsmoothed_cluster_metadata.csv")
#   pebetas <- read_csv("pe_max_unsmoothed_subj_betas.csv")
# } else {
#   setwd('/Users/andypapale/Box/SCEPTIC_fMRI/MMClock_aroma_preconvolve_fse_groupfixed/sceptic-clock-feedback-pe_max-preconvolve_fse_groupfixed/pe_max/')
#   pemeta <- read_csv("pe_max_cluster_metadata.csv")
#   pebetas <- read_csv("pe_max_roi_betas.csv")
# }
# pemeta$label <- substr(pemeta$label,22,100)
# pemeta_overall <- pemeta[pemeta$l2_contrast == 'overall' & pemeta$l3_contrast == 'Intercept' & pemeta$model == 'Intercept-Age',]
# pe <- as_tibble(pebetas[pebetas$l2_contrast == 'overall' & pebetas$l3_contrast == 'Intercept' & pebetas$model == 'Intercept-Age',1:3]) %>% filter(cluster_number<8 | cluster_number == 10 | cluster_number == 11)
# # head(merge(h,meta))
# perois_list <- distinct(pemeta_overall[c(1:7, 10:11),c(5,12)])
# 
# # inspect distributions
# # for some reason, the clusters do not seem to fit the map -- check with MNH (big clusters broken down?) !!
# perois <- inner_join(pe,pemeta_overall)
# perois$labeled_cluster <- paste(perois$cluster_number,perois$label)
# pemeta_overall$labeled_cluster <- paste(pemeta_overall$cluster_number,pemeta_overall$label)
# # ggplot(perois,aes(scale(cope_value))) + geom_histogram() + facet_wrap(~labeled_cluster)
# 
# pe_labeled <- inner_join(pe,perois_list)
# pe_labeled$labeled_cluster <- paste(pe_labeled$cluster_number,pe_labeled$label)
# pe_num <- select(pe_labeled,c(1,2,3))
# pe_labeled <- select(pe_labeled,c(1,3,5))
# 
# pe_wide <- spread(pe_labeled,labeled_cluster,cope_value)
# pe_wide_num <- spread(pe_num,cluster_number,cope_value)
# 
# ## we are no longer winsorizing betas
# # pe_wide <- pe_wide %>% mutate_if(is.double, winsor,trim = .075)
# # pe_wide_num <- pe_wide_num %>% mutate_if(is.double, winsor,trim = .075)
# 
# pejust_rois <- pe_wide[,2:ncol(pe_wide)]
# pejust_rois_num <- pe_wide_num[,2:ncol(pe_wide_num)]
# 
# 
# # non-parametric correlations to deal with outliers
# peclust_cor <- corr.test(pejust_rois,method = 'pearson', adjust = 'none')
# # parametric correlations on winsorised betas
# # clust_cor <- cor(just_rois_w,method = 'pearson')

# setwd(file.path(clock_folder, 'fmri/keuka_brain_behavior_analyses/'))
# if (unsmoothed) {
#   pdf("pe_cluster_corr_fixed_unsmoothed.pdf", width=12, height=12)
# } else {pdf("pe_cluster_corr_fixed.pdf", width=12, height=12)}
# corrplot(peclust_cor$r, cl.lim=c(-1,1),
#          method = "circle", tl.cex = 1.5, type = "upper", tl.col = 'black',
#          order = "hclust", diag = FALSE,
#          addCoef.col="black", addCoefasPercent = FALSE,
#          p.mat = peclust_cor$p, sig.level=0.05, insig = "blank")
# dev.off()

# mpe <- nfactors(peclust_cor$r, n=5, rotate = "oblimin", diagonal = FALSE,fm = "pa", n.obs = 70, SMC = FALSE)
# pe.fa = psych::fa(pejust_rois, nfactors=2, rotate = "varimax", fm = "pa")
# 
# 
# pe.fa = fa.sort(psych::fa(pejust_rois, nfactors=2))
# pefscores <- factor.scores(pejust_rois, pe.fa)$scores
# pe_wide$pe_f1_cort_str <- pefscores[,1]
# pe_wide$pe_f2_hipp <- pefscores[,2]
# pe_wide$pe_PH <- scale(rowMeans(cbind(pe_wide$`10 Left Hippocampus`, pe_wide$`7 Right Hippocampus`)))
# pe_wide$pe_PH_l <- scale(pe_wide$`10 Left Hippocampus`)
# pe_wide$pe_PH_r <- scale(pe_wide$`7 Right Hippocampus`)
# 
# # save loadings
# setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/hippo/supp')
# pf <- as_tibble(cbind(rownames(unclass(round(pe.fa$Structure, digits = 3))),unclass(round(pe.fa$Structure, digits = 3)))) %>% arrange(desc(MR1) ) %>%  
#   rename(region = V1, `Factor 1` = MR1, `Factor 2` = MR2)
# pf$region <- gsub('[0-9;]+', '', pf$region)
# stargazer(pf, type = 'html', out = 'pe_fa_structure.html', summary = F)
# 
# if (unsmoothed) {
#   pe_wide$pe_PH <- (pe_wide$`7 Right Hippocampus` + pe_wide$`10 Left Hippocampus`)/2
#   hpe_wide <- inner_join(h_wide,pe_wide[,c("feat_input_id","pe_f1_cort_str", "pe_f2_hipp", "pe_PH")])
# } else {  hpe_wide <- inner_join(h_wide,pe_wide[,c("feat_input_id","pe_f1_cort_str", "pe_f2_hipp", "pe_PH", "pe_PH_l", "pe_PH_r")])  }
# 

#####
# add ids
map_df  <- as.tibble(read.csv("~/Box/SCEPTIC_fMRI/MMClock_aroma_preconvolve_fse_groupfixed/sceptic-clock-feedback-v_chosen-preconvolve_fse_groupfixed/v_chosen/v_chosen-Intercept_design.txt", sep=""))
h_wide$ID <- h_wide$id
pc_scores <- inner_join(h_wide,map_df[,c(1:2,4:15)])
#pc_scores$id <- pc_scores$ID


# get trial_level data
trial_df <- read_csv(file.path(clock_folder, "fmri/data/mmclock_fmri_decay_factorize_selective_psequate_mfx_trial_statistics.csv.gz"))
# trial_df <- read_csv("~/Box/SCEPTIC_fMRI/mmclock_fmri_fixed_uv_mfx_trial_statistics.csv.gz")
# u_df <- read_csv("~/Box/SCEPTIC_fMRI/mmclock_fmri_fixed_uv_ureset_mfx_trial_statistics.csv.gz")
# use fixed-parameter U
u_df <- read_csv("/Users/andypapale/Box/SCEPTIC_fMRI/sceptic_model_fits/mmclock_fmri_fixed_uv_ureset_fixedparams_fmri_ffx_trial_statistics.csv.gz")
# Fixed params = everyone’s first-level models were estimated at fixed, group-mean parameter values.
trial_df <- trial_df %>% # trial_df trial level variables
  group_by(id, run) %>%  dplyr::mutate(rt_swing = abs(c(NA, diff(rt_csv))),
                                       rt_swing_lr = abs(log(rt_csv/lag(rt_csv))),
                                       rt_lag = lag(rt_csv) ,
                                       rt_swing_lag = lag(rt_swing),
                                       omission_lag = lag(score_csv==0),
                                       rt_vmax_lag = lag(rt_vmax),
                                       v_chosen_lag = lag(v_chosen),
                                       run_trial=1:50) %>% ungroup() #compute rt_swing within run and subject
u_df <- u_df %>% select(id, run, trial, u_chosen, u_chosen_lag, u_chosen_change,  # data frame with uncertainty variables
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
# if (unsmoothed) {
#   params_beta <- sub_df[,c("h_f1_fp", "h_f2_neg_paralimb","h_HippAntL",
#                            "DAN", "dan_r_sfg", "dan_l_sfg", "dan_parietal",
#                            "pe_f1_cort_str", "pe_f2_hipp", "pe_PH", "pe_PH_l", "pe_PH_r",
#                            "total_earnings", "LL", "alpha", "gamma", "beta")]
# } else {
#   params_beta <- sub_df[,c("h_f1_fp", "h_f2_neg_paralimb","h_HippAntL",
#                            "DAN", "dan_r_sfg", "dan_l_sfg", "dan_parietal",
#                            "pe_f1_cort_str", "pe_f2_hipp", "pe_PH",
#                            "total_earnings", "LL", "alpha", "gamma", "beta")]}
# param_cor <- corr.test(params_beta,method = 'pearson', adjust = 'none')
# 
# setwd(file.path(clock_folder, 'fmri/keuka_brain_behavior_analyses/'))

# merge into trial-level data
df <- inner_join(trial_df,sub_df, by = "id")
df$rewFunc <- relevel(as.factor(df$rewFunc),ref = "CEV")
df$rewFuncIEVsum <- df$rewFunc
contrasts(df$rewFuncIEVsum) <- contr.sum
colnames(contrasts(df$rewFuncIEVsum)) <- c('CEV','CEVR', 'DEV')

df$learning_epoch <- 'trials 1-10'
df$learning_epoch[df$run_trial>10] <- 'trials 11-50'


# obtain within-subject v_max and entropy: correlated at -.37
# isolating within subject change in entropy to address confound of related variables within subject (i.e. are all poorly performing subjects betas different)
df <- df %>% group_by(id,run) %>% mutate(v_max_wi = scale(v_max),  # grouping by id and run
                                         v_max_wi_lag = lag(v_max_wi),
                                         v_entropy_wi = scale(v_entropy), # change of entropy within run, scaled
                                         v_max_b = mean(na.omit(v_max)), # mean of subject and run
                                         v_entropy_b = mean(na.omit(v_entropy)), 
                                         rt_change = rt_csv - rt_lag,
                                         pe_max_lag = lag(pe_max), 
                                         abs_pe_max_lag = abs(pe_max_lag), 
                                         rt_vmax_change = rt_vmax - rt_vmax_lag,
                                         trial_neg_inv_sc = scale(-1/run_trial),
                                         v_chosen_change = v_chosen - lag(v_chosen)) %>% ungroup() %>% 
  mutate(rt_lag_sc = scale(rt_lag), # same scale for all subjects, scale across subjects
         rt_csv_sc = scale(rt_csv),
         rt_vmax_lag_sc = scale(rt_vmax_lag))

# correlate between-subject V and H with clusters
#b_df <- df %>% group_by(id) %>% dplyr::summarise(v_maxB = mean(v_max, na.rm = T),
#                                                v_entropyB = mean(v_entropy, na.rm = T))

# sub_df <- inner_join(sub_df, b_df, by = 'id')
# if (unsmoothed) {
#  bdf <- sub_df[,c("1","2","3","4","5","6", ...  do this in smoothed
#                   "total_earnings", "LL", "alpha", "gamma", "beta", "v_maxB", "v_entropyB")]
# } else {
#  bdf <- sub_df[,c("h_f1_fp", "h_f2_neg_paralimb",
#                   "DAN", "dan_r_sfg", "dan_l_sfg", "dan_parietal",
#                   "pe_f1_cort_str", "pe_f2_hipp", "pe_PH", "pe_PH_l", "pe_PH_r",
#                   "total_earnings", "LL", "alpha", "gamma", "beta", "v_maxB", "v_entropyB")]}
# b_cor <- corr.test(bdf,method = 'pearson', adjust = 'none')

###

###
# setwd(file.path(clock_folder, 'fmri/keuka_brain_behavior_analyses/'))
# pdf("between_subject_beta_beh_corr_fixed.pdf", width=12, height=12)
# corrplot(b_cor$r, cl.lim=c(-1,1),
#          method = "circle", tl.cex = 1.5, type = "upper", tl.col = 'black',
#          order = "hclust", diag = FALSE,
#          addCoef.col="black", addCoefasPercent = FALSE,
#          p.mat = b_cor$p, sig.level=0.05, insig = "blank")
# dev.off()

# binary splitting for plots in paper
# df$h_f1_fp_resp <- 'low'
# df$h_f1_fp_resp[df$h_f1_fp>0] <- 'high'
# df$pe_f1_cort_str_resp <- 'low'
# df$pe_f1_cort_str_resp[df$pe_f1_cort_str>0] <- 'high'
# df$pe_f2_hipp_resp <- 'low'
# df$pe_f2_hipp_resp[df$pe_f2_hipp>0] <- 'high'
# df$h_HippAntL_neg <- -df$h_HippAntL
# df$h_HippAntL_resp <- 'low'
# df$h_HippAntL_resp[df$h_HippAntL<mean(df$h_HippAntL)] <- 'high'
df$last_outcome <- NA
df$last_outcome[df$omission_lag] <- 'Omission'
df$last_outcome[!df$omission_lag] <- 'Reward'
df$last_outcome <- relevel(as.factor(df$last_outcome), ref = "Reward")

# everything below here just adapt to new mask variables MP

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
mtdf <- read_csv("/Users/andypapale/Box/SCEPTIC_fMRI/sceptic_model_fits/mmclock_meg_decay_factorize_selective_psequate_fixedparams_meg_ffx_trial_statistics.csv.gz")
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
} else {save(file = '2020-11-18-Schaeffer-vmPFC-Brain_to_Behavior_Analysis.Rdata', df, mdf)}
