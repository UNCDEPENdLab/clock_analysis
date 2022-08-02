library(tidyverse)
library(psych)
library(fmri.pipeline)
library(corrr)
library(data.table)
library(readxl)
library(lme4)
# inspect correlations of absolute PE betas to understand if interactions with RT_lag and RT_Vmax originate from different components
if (Sys.getenv("USER")=="alexdombrovski") {
  setwd("~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/betas_final/")
  # source("../get_trial_data.R")
  source("parcel_brain_behavior_functions_alex.R")
} else if (Sys.getenv("USER")=="Alex") {
  setwd("~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/betas_final/")
  # source("../get_trial_data.R")
  source("parcel_brain_behavior_functions_alex.R")
} else {
  setwd("/Users/hallquist/Data_Analysis/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/betas_final")
  source("../get_trial_data.R")
  source("/Users/hallquist/Data_Analysis/r_packages/fmri.pipeline/R/mixed_by.R")
  source("../medusa_final/plot_medusa.R")
  analysis_dir <- "~/Data_Analysis/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/betas_final"
}
analysis_dir <- "~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/betas_final"
setwd(analysis_dir)



labels_df <- setDT(read_excel("~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/MNH DAN Labels 400 Good Only 47 parcels.xlsx")) %>%
  mutate(roi_num7 = as.factor(roi7_400), 
         mask_value = as.integer(roi7_400),
         plot_label = mnh_label_400, 
         vm_gradient17 = parcel_group) %>% select(roi_num7, mask_value, plot_label, vm_gradient17, network17_400_DAN, hemi, x, y, z)
abspe_betas <- fread("/Volumes/GoogleDrive/My Drive/SCEPTIC_fMRI/wholebrain_betas/L1m-abspe_plus_rew/Schaefer_444_final_2009c_2.3mm_cope_l2.csv.gz") %>%
  filter(l2_cope_name == "overall" & !l1_cope_name  %in% c("EV_clock", "EV_feedback", "EV_rew_om")) %>% # only parametric modulators
  dplyr::select(-feat_dir, -img, -mask_name, -session, -l1_cope_number, -l2_cope_number, -l2_model, -x, -y, -z) %>%
  rename(fmri_beta = value) %>%
  # merge(label_df, by = label_join_col, all.x = TRUE)
  merge(labels_df, by = "mask_value", all = FALSE) %>% select(c(id, fmri_beta, plot_label, hemi)) %>% mutate(
  plot_label_bl  = sub("^[^_]*_", "", plot_label)
)

abspe_betas_wide <- pivot_wider(abspe_betas %>% select(c(fmri_beta, plot_label, id)), names_from = plot_label, values_from = fmri_beta) %>%
  arrange(id)



cormat <- psych::corr.test(abspe_betas_wide %>% select(-id))
setwd("./bb/plots/")
pdf("abspe_betas_network_plot.pdf", height = 20, width = 20)
network_plot(cormat$r, min_cor = 0.5, colours = c("white", "white","white", "red"))
dev.off()
pdf("MMC_abspe_betas_corrplot.pdf", height = 24, width = 24)
corrplot(cormat$r, cl.lim=c(-1,1),
         method = "circle", tl.cex = 1.5, type = "upper", tl.col = 'black',
         order = "hclust", diag = FALSE,
         addCoef.col="black", addCoefasPercent = FALSE,
         p.mat = cormat$p, sig.level=0.05, insig = "blank")
dev.off()
abspe <- nfactors(abspe_betas_wide %>% select(-id), n=5, rotate = "oblimin", diagonal = FALSE,fm = "pa", n.obs = 71, SMC = FALSE)
abspe.fa = psych::fa(cormat$r, nfactors=2, rotate = "varimax", fm = "pa")
abspe.fa = psych::fa(cormat$r, nfactors=4, rotate = "varimax", fm = "pa")
abspe_fscores <- factor.scores(abspe_betas_wide %>% select(-id), abspe.fa)$scores
abspe_betas_wide$f1_MT <- abspe_fscores[,1]
abspe_betas_wide$f4_LIP <- abspe_fscores[,2]
abspe_betas_wide$f3_SPL <- abspe_fscores[,3]
abspe_betas_wide$f2_PM <- abspe_fscores[,4]
abspe_betas_wide %>% correlate() %>% network_plot(min_cor = 0.5, colours = c("white", "white","white", "red"))
corr.test(abspe_fscores)
source("../../../get_trial_data.R")
trial_df <- get_trial_data(repo_directory = "~/code/clock_analysis/", dataset = "mmclock_fmri")

df <- inner_join(trial_df, abspe_betas_wide, by = "id")

m1 <- lmer(rt_csv_sc ~ (trial_neg_inv + rt_lag + rt_vmax_lag + v_entropy_wi + last_outcome)^2 +
           rt_lag*last_outcome*f1_MT +
           rt_lag*last_outcome*f2_PM +
           rt_lag*last_outcome*f3_SPL +
           rt_lag*last_outcome*f4_LIP +
           rt_vmax_lag:f1_MT +
           rt_vmax_lag:f2_PM +
           rt_vmax_lag:f3_SPL +
           rt_vmax_lag:f4_LIP +  
  (1 + rt_lag + rt_vmax_lag| id/run), df)
summary(m1)

# entropy_change

echange_betas <- fread("/Volumes/GoogleDrive/My Drive/SCEPTIC_fMRI/wholebrain_betas/L1m-echange/Schaefer_444_final_2009c_2.3mm_cope_l2.csv.gz") %>%
  filter(l2_cope_name == "overall" & !l1_cope_name  %in% c("EV_clock", "EV_feedback", "EV_rew_om")) %>% # only parametric modulators
  dplyr::select(-feat_dir, -img, -mask_name, -session, -l1_cope_number, -l2_cope_number, -l2_model, -x, -y, -z) %>%
  rename(fmri_beta = value) %>%
  # merge(label_df, by = label_join_col, all.x = TRUE)
  merge(labels_df, by = "mask_value", all = FALSE) %>% select(c(id, fmri_beta, plot_label, hemi)) %>% mutate(
    plot_label_bl  = sub("^[^_]*_", "", plot_label)
  )

echange_betas_wide <- pivot_wider(echange_betas %>% select(c(fmri_beta, plot_label, id)), names_from = plot_label, values_from = fmri_beta)



cormat <- psych::corr.test(abspe_betas_wide %>% select_if(is.numeric) %>% select(-id))
network_plot(cormat$r, min_cor = 0.5, colours = c("white", "white", "white",  "white", "red"))

pdf("MMC_echange_betas_corrplot.pdf", height = 24, width = 24)
corrplot(cormat$r, cl.lim=c(-1,1),
         method = "circle", tl.cex = 1.5, type = "upper", tl.col = 'black',
         order = "hclust", diag = FALSE,
         addCoef.col="black", addCoefasPercent = FALSE,
         p.mat = cormat$p, sig.level=0.05, insig = "blank")
dev.off()
echange <- nfactors(cormat$r, n=5, rotate = "oblimin", diagonal = FALSE,fm = "pa", n.obs = 71, SMC = FALSE)
echange.fa = psych::fa(cormat$r, nfactors=2, rotate = "varimax", fm = "pa")



source("../../../get_trial_data.R")
trial_df <- get_trial_data(repo_directory = "~/code/clock_analysis/", dataset = "mmclock_fmri")
