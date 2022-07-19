# plots results of parcel-wise beta-behavior survival analyses
# uses output of parcel_brain_behavior_coxme.R (this script is the "from cache" version of it)

library(data.table)
library(tidyverse)
library(afex)
library(lattice)
library(emmeans)
library(fmri.pipeline) # has mixed_By
library(coxme)
library(viridis)
library(ggnewscale)
library(RColorBrewer)
library(ggrepel)
source("~/code/Rhelpers/theme_black.R")

plots = T
inspect = F
beta_dir <- "~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/betas_final"
# source("~/code/fmri.pipeline/R/mixed_by.R")


studies = c("fmri")
# studies = c("meg", "fmri") # whether to include meg session replication

# sensitivity analyses for fMRI sample, recommend running first 2, one at a time:
censor = c(T)            # if "T", sensitivity analysis: censor first 500 ms and last second
decompose = c(T)            # if "F" sensitivity analyss: don't decompose U and V into within- vs between-trial components
omit_value = F              # do not change, needed to check once; not necessary ex-post
split_by_reward = F         # if "T", sensitivity analysis including interaction with last reward/omission; not necessary ex-post

# loop over model versions as defined above

if (omit_value) {rhs = "no_value"} else {rhs = ""}
for (censor_ends in censor) {
  for (study in studies) {
    for (decompose_within_between_trial in decompose) {
        if (!censor_ends) {
          out_dir <- file.path(beta_dir, "coxme")
        } else if (censor_ends) {
          out_dir <- file.path(beta_dir, "coxme/censored")}
        
        if (decompose_within_between_trial) {method = "decomposed"
        } else {method = "non_decomposed"}
        setwd(out_dir)
        ddf <- readRDS(file.path(paste0("beta_coxme_", method, "_", study, "_", rhs, "_trial.rds")))
      ############ Plot
      if (plots) {
        for (beta in c("abspe", "entropychange")) {
          
          pdf(paste0("value_axial_", beta, "_", method, "_", study, ".pdf"), height = 12, width = 16)
          print(ggplot(ddf %>% filter(fmri_beta == beta & term == "value_wi_t:scale(h)"),  aes(y, x, shape = network17_400_DAN)) + 
                  geom_point(size = 24, color = "white", aes(alpha = p_level_fdr, fill = Statistic)) + scale_shape_manual(values = 21:22) + geom_point(size = 24, aes(color = network17_400_DAN)) +
                  scale_fill_viridis(option = "inferno") + theme_black() + geom_text_repel(aes(label=plot_label), color = "white", force = 5, point.padding = 40, size = 8))
          dev.off()
          
          pdf(paste0("uncertainty_axial_", beta, "_", method, "_", study, ".pdf"), height = 12, width = 16)
          print(ggplot(ddf %>% filter(fmri_beta == beta & term == "scale(h):uncertainty_wi_t"),  aes(y, x, shape = network17_400_DAN)) + 
                  geom_point(size = 24, color = "white", aes(alpha = p_level_fdr, fill = Statistic)) + scale_shape_manual(values = 21:22) + geom_point(size = 24, aes(color = network17_400_DAN)) +
                  scale_fill_viridis(option = "turbo") + theme_black() + geom_text_repel(aes(label=plot_label), color = "white", force = 5, point.padding = 40, size = 8))
          dev.off()
          
          pdf(paste0("value_saggital_",beta, "_",  method, "_", study, ".pdf"), height = 12, width = 16)
          print(ggplot(ddf %>% filter(fmri_beta == beta & term == "value_wi_t:scale(h)" ),  aes(y, z, shape = network17_400_DAN)) + 
                  geom_point(size = 32, color = "white", aes(alpha = p_level_fdr, 
                                                             fill = Statistic)) + 
                  scale_shape_manual(values = 21:22) + geom_point(size = 32, color = "white") +
                  scale_fill_viridis(option = "inferno") + theme_black() + geom_text_repel(aes(label=plot_label),  color="blue"))
          dev.off()
          
          pdf(paste0("uncertainty_saggital_",beta, "_",  method, "_", study, ".pdf"), height = 12, width = 16)
          print(ggplot(ddf %>% filter(fmri_beta == beta & term == "scale(h):uncertainty_wi_t"),  aes(y, z, shape = network17_400_DAN)) + 
                  geom_point(size = 32, color = "white", aes(alpha = p_level_fdr, 
                                                             fill = Statistic)) + 
                  scale_shape_manual(values = 21:22) + geom_point(size = 32, color = "white") +
                  scale_fill_viridis(option = "turbo") + theme_black() + geom_text_repel(aes(label=plot_label),  color="blue") )
          dev.off()
          
          pdf(paste0("value_vioin_", beta, "_", method, "_", study, ".pdf"), height = 8, width = 5)
          print(ggplot(ddf %>% filter(fmri_beta == beta & term == "value_wi_t:scale(h)" )) + 
                  geom_jitter(size = 12, width = .1, height = 0,  aes(network17_400_DAN, Statistic, color = Statistic, alpha = p_level_fdr)) + 
                  geom_violin(aes(network17_400_DAN, Statistic), alpha = .2) + scale_shape_manual(values = 21:22) +
                  scale_color_viridis(option = "inferno") + theme_black() + 
                  geom_text_repel(aes(network17_400_DAN, Statistic, alpha = p_level_fdr, label=plot_label),  color="blue"))
          dev.off()
          
          
          # main figure for entropy change:
          pdf(paste0("value_violin_by_vm_gradient_", beta, "_", method, "_", study, ".pdf"), height = 8, width = 10)
          print(ggplot(ddf %>% filter(fmri_beta == beta & term == "value_wi_t:scale(h)" )) + 
                  geom_jitter(size = 12, width = .1, height = 0,  aes(vm_gradient17, Statistic, color = Statistic, alpha = p_level_fdr)) + 
                  geom_violin(aes(vm_gradient17, Statistic), alpha = .2) + scale_shape_manual(values = 21:22) +
                  scale_color_viridis(option = "inferno") + theme_black() + 
                  geom_text_repel(aes(vm_gradient17, Statistic, alpha = p_level_fdr, label=plot_label),  color="blue"))
          dev.off()
          
          
          pdf(paste0("uncertainty_violin_", beta, "_", method, "_", study, ".pdf"), height = 8, width = 5)
          print(ggplot(ddf %>% filter(fmri_beta == beta & term == "scale(h):uncertainty_wi_t")) + 
                  geom_jitter(size = 12, width = .1, height = 0,  aes(network17_400_DAN, Statistic, color = Statistic, alpha = p_level_fdr)) + 
                  geom_violin(aes(network17_400_DAN, Statistic), alpha = .2) + scale_shape_manual(values = 21:22) +
                  scale_color_viridis(option = "turbo") + theme_black() + 
                  geom_text_repel(aes(network17_400_DAN, Statistic, alpha = p_level_fdr, label=plot_label),  color="blue"))
          dev.off()
          
          pdf(paste0("uncertainty_violin_by_vm_gradient_", beta, "_", method, "_", study, ".pdf"), height = 8, width = 10)
          print(ggplot(ddf %>% filter(fmri_beta == beta & term == "scale(h):uncertainty_wi_t")) + 
                  geom_jitter(size = 12, width = .1, height = 0,  aes(vm_gradient17, Statistic, color = Statistic, alpha = p_level_fdr)) + 
                  geom_violin(aes(vm_gradient17, Statistic), alpha = .2) + scale_shape_manual(values = 21:22) +
                  scale_color_viridis(option = "turbo") + theme_black() + 
                  geom_text_repel(aes(vm_gradient17, Statistic, alpha = p_level_fdr, label=plot_label),  color="blue"))
          dev.off()
          
          pdf(paste0("uncertaintyBYtrial_violin_", beta, "_", method, "_", study, ".pdf"), height = 8, width = 5)
          print(ggplot(ddf %>% filter(fmri_beta == beta & term == "scale(h):uncertainty_wi_t:trial_neg_inv_sc")) + 
                  geom_jitter(size = 12, width = .1, height = 0,  aes(network17_400_DAN, Statistic, color = Statistic, alpha = p_level_fdr)) + 
                  geom_violin(aes(network17_400_DAN, Statistic), alpha = .2) + scale_shape_manual(values = 21:22) +
                  scale_color_viridis(option = "turbo") + theme_black() + 
                  geom_text_repel(aes(network17_400_DAN, Statistic, alpha = p_level_fdr, label=plot_label),  color="blue"))
          dev.off()
          
          pdf(paste0("uncertaintyBYtrial_violin_by_vm_gradient_", beta, "_", method, "_", study, ".pdf"), height = 8, width = 10)
          print(ggplot(ddf %>% filter(fmri_beta == beta & term == "scale(h):uncertainty_wi_t:trial_neg_inv_sc")) + 
                  geom_jitter(size = 12, width = .1, height = 0,  aes(vm_gradient17, Statistic, color = Statistic, alpha = p_level_fdr)) + 
                  geom_violin(aes(vm_gradient17, Statistic), alpha = .2) + scale_shape_manual(values = 21:22) +
                  scale_color_viridis(option = "turbo") + theme_black() + 
                  geom_text_repel(aes(vm_gradient17, Statistic, alpha = p_level_fdr, label=plot_label),  color="blue"))
          dev.off()
        }
      }
    }
  }
}


if (inspect) {
  
  summary(lm(uncertainty_wi_t ~ value_wi_t * trial_neg_inv_sc * omission_lag, trial_df))
  summary(lm(uncertainty_wi_t ~ value_wi_t * trial_neg_inv_sc * omission_lag, trial_df))
  summary(lm(uncertainty_wi ~ value_wi + trial_neg_inv_sc, trial_df))
  summary(lm(uncertainty_wi ~ value_wi, trial_df))
  
  m  <- coxme(Surv(t1,t2,response) ~ value_wi_t + uncertainty_wi_t * trial_neg_inv_sc + uncertainty_wi_t * omission_lag +
                (1|ID), surv_df)
  summary(m)
  
  m2  <- coxme(Surv(t1,t2,response) ~ value_wi_t * value_b_t + uncertainty_wi_t * uncertainty_b_t + trial_neg_inv_sc + uncertainty_wi_t * omission_lag +
                 (1|ID), surv_df)
  summary(m2)
  
  m0  <- coxme(Surv(t1,t2,response) ~  uncertainty_wi_t * trial_neg_inv_sc + uncertainty_wi_t * omission_lag +
                 (1|ID), surv_df)
  summary(m0)
  # inspect models to understand uncertainty sensitivity
  # look at the key "decomposed" models for fMRI and MEG
  ddf_fmri <- read_rds("/Volumes/GoogleDrive/My Drive/SCEPTIC_fMRI/wholebrain_betas/coxme/beta_coxme_decomposed_fmritrial.rds") %>% mutate(
    sample = "fmri",
    model_flavor = "decomposed"
  )
  ddf_meg <- read_rds("/Volumes/GoogleDrive/My Drive/SCEPTIC_fMRI/wholebrain_betas/coxme/beta_coxme_decomposed_megtrial.rds") %>% mutate(
    sample = "meg",
    model_flavor = "decomposed"
  )
  
  library(VIM)
  aggr_plot <- aggr(ddf %>% filter(term == "scale(h):uncertainty_wi_t:trial_neg_inv_sc"), col=c('navyblue','red'), numbers=TRUE, sortVars=TRUE, labels=names(data), cex.axis=.7, gap=3, ylab=c("Histogram of missing data","Pattern"))  
  
  # aggr_plot <- aggr(ddf_fmri, col=c('navyblue','red'), numbers=TRUE, sortVars=TRUE, labels=names(data), cex.axis=.7, gap=3, ylab=c("Histogram of missing data","Pattern"))
  # aggr_plot <- aggr(ddf_meg , col=c('navyblue','red'), numbers=TRUE, sortVars=TRUE, labels=names(data), cex.axis=.7, gap=3, ylab=c("Histogram of missing data","Pattern"))
  # 
  # 
  # ddf <- rbind(ddf_fmri, ddf_meg)
  # 
  toplot <- 
    ddf %>% filter(fmri_beta == "entropychange" & term == "uncertainty_wi_t")
  histogram(toplot$Statistic)
  ggplot() + 
    geom_jitter(size = 12, width = .1, height = 0,  aes(vm_gradient17, Statistic, color = Statistic, alpha = p_level_fdr)) + 
    geom_violin(aes(vm_gradient17, Statistic), alpha = .2) + scale_shape_manual(values = 21:22) +
    scale_color_viridis(option = "turbo") + theme_black() + 
    geom_text_repel(aes(vm_gradient17, Statistic, alpha = p_level_fdr, label=plot_label),  color="blue") #+ facet_wrap(~sample)
  
  ddf_meg %>% filter(fmri_beta == "abspe" & term == "h:uncertainty_wi_t") %>% View()
  
  m <- coxme(Surv(t1,t2,response) ~ value_wi_t  + uncertainty_wi_t + (1|ID), bb)
  summary(m)
  
  mt <- coxme(Surv(t1,t2,response) ~ value_wi_t  + uncertainty_wi_t * trial_neg_inv_sc + (1|ID), bb)
  summary(mt)
  
  
  mm <- coxme(Surv(t1,t2,response) ~ value_wi_t + uncertainty_wi_t + (1|ID), mbb)
  summary(mm)
  
  mmt <- coxme(Surv(t1,t2,response) ~ value_wi_t  + uncertainty_wi_t * trial_neg_inv_sc + (1|ID), mbb)
  summary(mmt)
  
}