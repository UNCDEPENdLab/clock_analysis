# plot beta-to-behavior analysis results from the 17-network Schaeffer 2018 parcellation with 47 manual DAN labels from June-July 2022
# abspe - rt_lag section makes bb parts of Fig. 3
# entropy change - rt_vmax section makes bb parts of Fig. 5

library(tidyverse)
library(lme4)
library(ggpubr)
library(car)
library(viridis)
library(ggnewscale)
library(RColorBrewer)
library(emmeans)
library(readxl)
library(readr)
library(data.table)
library(ggrepel)
library(patchwork)
source("~/code/Rhelpers/theme_black.R")
# install_github("UNCDEPENdLab/dependlab")
# library(dependlab)
# library(stringi)

signals = c("entropychange") # what L1 betas are used
all_plots = F # whether to make anatomical plots in addition to the main ones; not functional yet
# for (session in c("fmri", "meg")) {
# clock_folder <- "~/Data_Analysis/clock_analysis" #michael
clock_folder <- "~/code/clock_analysis" #alex

# # original output here:
# fmri_dir <- '/Volumes/GoogleDrive/MyDrive/SCEPTIC_fMRI/wholebrain_betas'

# just the necessary output is here:
fmri_dir <- file.path(paste0(clock_folder, "/fmri/keuka_brain_behavior_analyses/dan/betas_final/bb/"))

# DAN palette:
colors <- RColorBrewer::brewer.pal(4, "Dark2") %>% setNames(c("1" = "MT+","2" = "Premotor","3" = "Rostral PPC","4" = "Caudal PPC"))

coxme_ranef <- "rslope" # or NULL, whether to include random slopes in coxme value sensitivity models

# Manual labels from June-July 2022

labels_df <- setDT(read_excel(file.path(paste0(fmri_dir, "../../MNH DAN Labels 400 Good Only 47 parcels.xlsx")))) %>%
  mutate(roi_num7 = as.factor(roi7_400), 
         mask_value = as.integer(roi7_400),
         plot_label = mnh_label_400, 
         vm_gradient17 = parcel_group,
         vm_gradient17_names = case_when(
           vm_gradient17 == "PPCcaudal" ~ "Caudal PPC",
           vm_gradient17 == "PPCrostral" ~ "Rostral PPC",
           T ~ vm_gradient17
         ),
         vm_gradient17_names = factor(vm_gradient17_names, order = T, levels = c("MT+", "Caudal PPC", "Rostral PPC", "Premotor"))) %>% select(roi_num7, mask_value, plot_label, vm_gradient17, vm_gradient17_names, network17_400_DAN, hemi, x, y, z)
proc_df <- function(df, terms, cope) {
  # import the results of brain-behavior analyses
  df <- df$coef_df_reml %>% filter(term %in% terms) %>%
    filter(l1_cope_name == cope) %>%
    rename(roi_num7 = "mask_value") %>% mutate(roi_num7 = as.factor(roi_num7),
                                               padj_BY_term = p.adjust(p.value, method = 'BH'),
                                               p_level_fdr = as.factor(case_when(
                                                 # p_fdr > .1 ~ '0',
                                                 # p_fdr < .1 & p_fdr > .05 ~ '1',
                                                 padj_BY_term > .05 & p.value > .05 ~ '1',
                                                 padj_BY_term > .05 & p.value < .05 ~ '2',
                                                 padj_BY_term < .05 & padj_BY_term > .01 ~ '3',
                                                 padj_BY_term < .01 & padj_BY_term > .001 ~ '4',
                                                 padj_BY_term <.001 & padj_BY_term > .0001 ~ '5',
                                                 padj_BY_term <.0001 & padj_BY_term > .00001 ~ '6',
                                                 padj_BY_term <.00001 ~ '7'
                                               )),
                                               reward = case_when(
                                                 term == "rt_lag:fmri_beta" ~ "Omission",
                                                 term == "rt_lag:fmri_beta:last_outcomeReward" ~ "Reward"
                                               )) %>% inner_join(labels_df, by = "roi_num7") 
  df$p_level_fdr <- factor(df$p_level_fdr, levels = c('1', '2', '3', '4', '5', '6', '7'), labels = c("NS","p_uncorr. < .05","p < .05", "p < .01", "p < .001", "p < .0001", "p < .00001")) 
  return(df)}


if ("abspe" %in% signals) {
  setwd(file.path(fmri_dir, 'L1m-abspe'))
  # plot brain to behavior results for absolute PEs
  temp_meg <- readRDS("./parcel_maps_l2/Schaefer_400_DAN_manual_labels_47_meg__mixed_by.rds")  
  temp_fmri <- readRDS("./parcel_maps_l2/Schaefer_400_DAN_manual_labels_47_fmri_mixed_by.rds")
  abspe_terms <- c("rt_lag:fmri_beta", "fmri_beta:rt_vmax_lag", "rt_lag:fmri_beta:last_outcomeReward")
  mdf <- proc_df(temp_meg, abspe_terms, "EV_abspe") %>% mutate(session = "MEG replication")
  fdf <- proc_df(temp_fmri, abspe_terms, "EV_abspe") %>% mutate(session = "fMRI")
  df <- rbind(mdf, fdf)
  
  setwd(file.path(fmri_dir, "plots"))
  
  
  # main champagne plot, random slopes:
  # library(wesanderson)
  # pal <- wes_palette("Zissou1",  type = "continuous")
  # make insets for omission
  get_abspe_inset = function(df, study) {ggplot(df %>% filter(term %in% c("rt_lag:fmri_beta") & model_name == "slo" & session == study) %>% mutate(`RT swings` = -statistic)) + 
      geom_jitter(size = 4, width = .1, height = 0,  aes(vm_gradient17, -statistic, color = vm_gradient17_names, alpha = p_level_fdr), show.legend = F) + 
      geom_violin(aes(vm_gradient17, - statistic, color = vm_gradient17_names), alpha = .2, show.legend = F) + 
      scale_color_manual(values = colors ) +
      theme_minimal() + geom_hline(yintercept = 0, size = .3) + 
      scale_x_discrete(labels = c("MT+" = "MT+", "PPCcaudal" = "PPCc", "PPCrostral" = "PPCr", "Premotor" = "Prem.")) +
      theme(axis.text.x = element_text(size = 8, angle = 90)) + xlab(NULL) + ylab(NULL) + 
      facet_wrap(~reward)
    # facet_wrap(~reward) + guides(alpha = guide_legend(override.aes = list(color = "white"), title = expression(p[FDR])))
  }
  abspe_inset_fmri <- get_abspe_inset(df, "fMRI")
  abspe_inset_meg <- get_abspe_inset(df, "MEG replication")
  
  abspeStats <- df %>% filter(term %in% c("rt_lag:fmri_beta:last_outcomeReward")  & model_name == "slo") %>% select(plot_label, statistic, session, vm_gradient17_names) %>% 
    mutate(statistic = -statistic,
           session = sub(" .*", "", session)) %>%
    pivot_wider(names_from = session, values_from = statistic) 
  cor <- cor.test(abspeStats$fMRI, abspeStats$MEG, method = c("pearson"))
  abspe_cor_inset <- ggplot(abspeStats, aes(fMRI, MEG)) + geom_point(size = .8) + theme_light() + geom_smooth(method = "glm", alpha = .2, color = "grey") + 
    stat_cor(aes(label = ..r.label..), method = "pearson", cor.coef.name = "r", label.x = 1.5, label.y = -2, size = 2.8)
  
  
  # main panel is for reward (win-shift) with omission as inset
  pdf(paste0("abspe_rt_swings_by_vm_gradient17_", "_rslope.pdf"), height = 5, width = 10)
  print(
    ggplot(df %>% filter(term %in% c("rt_lag:fmri_beta:last_outcomeReward")  & model_name == "slo") %>% mutate(`RT swing statistic` = -statistic)) + 
      geom_jitter(size = 6, width = .1, height = 0,  aes(vm_gradient17_names, -statistic, color = vm_gradient17_names, alpha = p_level_fdr), show.legend = T) + 
      geom_violin(aes(vm_gradient17_names, color = vm_gradient17_names, -statistic), alpha = .2) + # scale_shape_manual(values = 21:25) +
      # scale_color_viridis(option = "mako") + theme_black() +  
      # scale_color_gradientn(colors = pal) + 
      scale_color_manual(values = colors ) +
      theme_minimal() + xlab(NULL) + ylab("Behavioral effect of neural response on post-reward RT swings,\n t statistic") + #labs(color = "Statistic") +
      geom_text_repel(aes(vm_gradient17_names, -statistic, alpha = p_level_fdr, label=plot_label, color = vm_gradient17_names), point.padding = 10, force = 10,  size = 2.4, show.legend = F) + 
      geom_text(aes(vm_gradient17_names, 7.2, label=vm_gradient17_names, color = vm_gradient17_names),  size = 3.5) + 
      facet_grid(~session, labeller = as_labeller(c("fMRI" = "fMRI","MEG replication" =  "MEG, out-of-session replication"))) + geom_hline(yintercept = 0, size = .2) + 
      theme(legend.key.size = unit(.4, "cm"), strip.text.x = element_text(size = 12), axis.text.x = element_blank()) +
      guides(alpha = guide_legend(title = expression(p[FDR], shape = guide_legend(override.aes = list(size = .1)))), color = "none") + 
      inset_element(abspe_inset_fmri, 0.31, -0.05, .55 ,0.24,  # left, bottom, right, top
                    align_to = "panel") + 
      inset_element(abspe_inset_meg, 0.7, -0.05, .95 ,0.24,  # left, bottom, right, top
                    align_to = "panel") + 
      inset_element(abspe_cor_inset, .97 , -0.05, 1.15, 0.24)
  )
  dev.off()
  if (all_plots){
    # anatomical detail for supplemental figures:
    # pdf(paste0("abspe_rt_swings_omission_and_reward_saggital_", session, ".pdf"), height = 12, width = 24)
    # print(ggplot(df %>% filter(term %in% c("rt_lag:fmri_beta", "rt_lag:fmri_beta:last_outcomeReward")  & model_name == "int"),  aes(y, z, shape = vm_gradient17_names, alpha = p_level_fdr)) + 
    #         geom_point(size = 20, color = "white", aes(alpha = p_level_fdr, fill = - statistic)) + scale_shape_manual(values = 21:25) + geom_point(size = 20, color = "white") +
    #         scale_fill_viridis(option = "mako") + theme_black() + 
    #         # geom_text_repel(aes(alpha = p_level_fdr, label=plot_label), color = "white", force = 5, point.padding = 40, size = 6) + 
    #         geom_text_repel(aes(alpha = p_level_fdr, label=plot_label), point.padding = 40, force = 6, color="white", size = 6) +
    #         facet_wrap(~reward) + guides(alpha = guide_legend(override.aes = list(color = "white"), title = expression(p[FDR]))))
    # dev.off()
    # 
    # pdf(paste0("abspe_rt_swings_omission_and_reward_axial_", session, ".pdf"), height = 12, width = 24)
    # print(ggplot(df %>% filter(term %in% c("rt_lag:fmri_beta", "rt_lag:fmri_beta:last_outcomeReward")  & model_name == "int"),  aes(y, -x, shape = vm_gradient17_names, alpha = p_level_fdr)) + 
    #         geom_point(size = 20, color = "white", aes(alpha = p_level_fdr, fill = - statistic)) + scale_shape_manual(values = 21:25) + geom_point(size = 20, color = "white") +
    #         scale_fill_viridis(option = "mako") + theme_black() + 
    #         # geom_text_repel(aes(alpha = p_level_fdr, label=plot_label), color = "white", force = 5, point.padding = 40, size = 6) + 
    #         geom_text_repel(aes(alpha = p_level_fdr, label=plot_label), point.padding = 40, force = 6, color="white", size = 6) +
    #         facet_wrap(~reward) + guides(alpha = guide_legend(override.aes = list(color = "white"), title = expression(p[FDR]))))
    # dev.off()
    # 
    
    
    # # test
    # summary(lm(statistic ~ vm_gradient17_names + hemi*term, df %>% filter(term %in% c("rt_lag:fmri_beta", "rt_lag:fmri_beta:last_outcomeReward"))))
    # summary(lm(statistic ~ vm_gradient17_names + hemi, df %>% filter(term %in% c("fmri_beta:rt_vmax_lag"))))
    
    # pdf(paste0("abspe_rtvmax_by_vm_gradient17_names_saggital_", session, ".pdf"), height = 12, width = 16)
    # print(ggplot(df %>% filter(term == "fmri_beta:rt_vmax_lag" & model_name == "int"),  aes(y, z, shape = vm_gradient17_names, alpha = p_level_fdr)) + 
    #         geom_point(size = 24, color = "white", aes(alpha = p_level_fdr, fill = statistic)) + scale_shape_manual(values = 21:25) + geom_point(size = 24, aes(color = vm_gradient17_names)) +
    #         scale_fill_viridis(option = "inferno") + theme_black() + geom_text_repel(aes(label=plot_label), color = "white", force = 5, point.padding = 40, size = 6) +
    #         guides(alpha = guide_legend(override.aes = list(color = "white"), title = expression(p[FDR]))))
    # dev.off()
    # pdf(paste0("abspe_rtvmax_by_vm_gradient17_axial_", session, ".pdf"), height = 12, width = 16)
    # print(ggplot(df %>% filter(term == "fmri_beta:rt_vmax_lag" & model_name == "int"),  aes(y, -x, shape = vm_gradient17, alpha = p_level_fdr)) + 
    #         geom_point(size = 24, color = "white", aes(alpha = p_level_fdr, fill = statistic)) + scale_shape_manual(values = 21:25) + geom_point(size = 24, aes(color = vm_gradient17)) +
    #         scale_fill_viridis(option = "inferno") + theme_black() + geom_text_repel(aes(label=plot_label), color = "white", force = 5, point.padding = 40, size = 6) +
    #         guides(alpha = guide_legend(override.aes = list(color = "white"), title = expression(p[FDR]))))
    # dev.off()
    # pdf(paste0("abspe_rtvmax_by_vm_gradient17_", session, ".pdf"), height = 8, width = 12)
    # 
    # print(ggplot(df %>% filter(term %in% c("fmri_beta:rt_vmax_lag")  & model_name == "int")) + 
    #         geom_jitter(size = 12, width = .1, height = 0,  aes(vm_gradient17, statistic, color = statistic, alpha = p_level_fdr)) + 
    #         geom_violin(aes(vm_gradient17, statistic), alpha = .2) + scale_shape_manual(values = 21:25) +
    #         # scale_color_viridis(option = "mako") + theme_black() + 
    #         scale_color_viridis(option = "inferno") + theme_black() + xlab("DAN region") + ylab("Convergence on RT(Vmax)") +
    #         geom_text_repel(aes(vm_gradient17, statistic, alpha = p_level_fdr, label=plot_label), point.padding = 20, force = 5, color="#4FC3F7", size = 3)  +
    #         guides(alpha = guide_legend(override.aes = list(color = "white"), title = expression(p[FDR]))))
    # # geom_text_repel(aes(vm_gradient17, - statistic, color = -statistic, alpha = p_level_fdr, label=plot_label),  color="blue")
    # dev.off()
    # 
    # # abspe exploit figure, random slopes:
    # pdf(paste0("abspe_rtvmax_by_vm_gradient17_", session, "_rslope.pdf"), height = 8, width = 12)
    # print(ggplot(df %>% filter(term %in% c("fmri_beta:rt_vmax_lag")  & model_name == "slo")) + 
    #         geom_jitter(size = 12, width = .1, height = 0,  aes(vm_gradient17, statistic, color = statistic, alpha = p_level_fdr)) + 
    #         geom_violin(aes(vm_gradient17, statistic), alpha = .2) + scale_shape_manual(values = 21:25) +
    #         # scale_color_viridis(option = "mako") + theme_black() + 
    #         scale_color_viridis(option = "inferno") + theme_black() + xlab("DAN region") + ylab("Convergence on RT(Vmax)") +
    #         geom_text_repel(aes(vm_gradient17, statistic, alpha = p_level_fdr, label=plot_label), point.padding = 20, force = 5, color="#4FC3F7", size = 3)  +
    #         guides(alpha = guide_legend(override.aes = list(color = "white"), title = expression(p[FDR]))))
    # # geom_text_repel(aes(vm_gradient17, - statistic, color = -statistic, alpha = p_level_fdr, label=plot_label),  color="blue")
    # dev.off()
  }
}

if ("entropychange" %in% signals) {
  # plot entropy change exploitation effect
  
  # import data
  setwd(file.path(fmri_dir, 'L1m-echange'))
  temp_meg <- readRDS("./parcel_maps_l2/Schaefer_400_DAN_manual_labels_47_meg__mixed_by.rds")  
  temp_fmri <- readRDS("./parcel_maps_l2/Schaefer_400_DAN_manual_labels_47_fmri_mixed_by.rds")
  terms <- c("fmri_beta:rt_vmax_lag")
  mdf <- proc_df(temp_meg, terms, "EV_entropy_change_feedback") %>% mutate(session = "MEG replication")
  fdf <- proc_df(temp_fmri, terms, "EV_entropy_change_feedback") %>% mutate(session = "fMRI")
  edf <- rbind(mdf, fdf) # different name than Abs(PE) dataframe to avoid confusion
  
  # get random slope GLM results as inset
  get_echange_inset = function(df, study, title) {
    df <- df %>% filter(term %in% c("fmri_beta:rt_vmax_lag") & model_name == "slo" & session == study)
    ggplot(df) + 
      # geom_jitter(size = 4, width = .1, height = 0,  aes(vm_gradient17, statistic, color = statistic, alpha = p_level_fdr), show.legend = F) + 
      geom_jitter(size = 4, width = .1, height = 0,  aes(vm_gradient17_names, statistic, color = vm_gradient17_names, alpha = p_level_fdr), show.legend = F) + 
      # scale_color_gradientn(colors = pal) + 
      geom_violin(aes(vm_gradient17_names, statistic, color = vm_gradient17_names), alpha = .2, show.legend = F) + ggtitle(title) +
      scale_color_manual(values = colors ) +
      theme_minimal() + geom_hline(yintercept = 0, size = .3) + 
      scale_x_discrete(labels = c("MT+" = "MT+", "Caudal PPC" = "PPCc", "Rostral PPC" = "PPCr", "Premotor" = "Prem.")) +
      theme(axis.text.x = element_text(size = 8, angle = 0), plot.title = element_text(size=9)) + xlab(NULL) + ylab(NULL) 
    # facet_wrap(~reward) + guides(alpha = guide_legend(override.aes = list(color = "white"), title = expression(p[FDR])))
  }
  echange_inset_fmri <- get_echange_inset(edf, "fMRI", "GLM") 
  echange_inset_meg <- get_echange_inset(edf, "MEG replication", "GLM")
  
  # get coxme results
  # set options as in parcel_brain_behavior_coxme.R
  censor_ends = F # whether to remove first and last second
  # if "F" sensitivity analyss: don't decompose U and V into within- vs between-trial components
  decompose_within_between_trial = F
  signals = c("entropychange")
  beta_dir <- "/Volumes/GoogleDrive/My Drive/SCEPTIC_fMRI/wholebrain_betas"
  rhs = "" # historic irrelevant analysis omitting value regressor
  if (!censor_ends) {
    out_dir <- file.path(beta_dir, "coxme")
  } else if (censor_ends) {
    out_dir <- file.path(beta_dir, "coxme/censored")}
  if (decompose_within_between_trial) {method = "decomposed"
  value_term = "value_wi_t:scale(h)"
  uncertainty_term = "scale(h):uncertainty_wi_t"
  } else {method = "non_decomposed"
  value_term = "value_wi:scale(h)"
  uncertainty_term = "scale(h):uncertainty_wi"
  }
  setwd(out_dir)
  
  meg_cox_df <- readRDS(file.path(paste0("beta_coxme_", method, "_meg", rhs, coxme_ranef, "_trial.rds"))) %>% mutate(session = "MEG replication")
  fmri_cox_df <- readRDS(file.path(paste0("beta_coxme_", method, "_fmri", rhs, coxme_ranef, "_trial.rds"))) %>% mutate(session = "fMRI")
  cox_df <- rbind(meg_cox_df, fmri_cox_df) %>% filter(fmri_beta %in% signals) %>% inner_join(labels_df %>% select(mask_value, vm_gradient17_names)) %>% group_by(session, term) %>%
    mutate(padj_BY_term = p.adjust(p, method = 'BH'),
           p_level_fdr = as.factor(case_when(
             # p_fdr > .1 ~ '0',
             # p_fdr < .1 & p_fdr > .05 ~ '1',
             padj_BY_term > .05 & p > .05 ~ '1',
             padj_BY_term > .05 & p < .05 ~ '2',
             padj_BY_term < .05 & padj_BY_term > .01 ~ '3',
             padj_BY_term < .01 & padj_BY_term > .001 ~ '4',
             padj_BY_term <.001 & padj_BY_term > .0001 ~ '5',
             padj_BY_term <.0001 & padj_BY_term > .00001 ~ '6',
             padj_BY_term <.00001 ~ '7'
           )),
           reward = case_when(
             term == "rt_lag:fmri_beta" ~ "Omission",
             term == "rt_lag:fmri_beta:last_outcomeReward" ~ "Reward"
           )) %>% ungroup() 
  cox_df$p_level_fdr <- factor(cox_df$p_level_fdr, levels = c('1', '2', '3', '4', '5', '6', '7'), labels = c("NS","p_uncorr. < .05","p < .05", "p < .01", "p < .001", "p < .0001", "p < .00001")) 
  
  if (coxme_ranef == "rslope") {yint = 4 
  cor_y_int = 3} else {yint = 18 
  cor_y_int = 13.5}  
  # sanity check: do the behavioral effects correlate across sessions anatomicall?
  
  ecStats <- cox_df %>% filter(term == value_term & fmri_beta %in% "entropychange") %>% select(plot_label, Statistic, study, vm_gradient17_names) %>% 
    pivot_wider(names_from = study, values_from = Statistic) 
  cor.test(ecStats$fmri, ecStats$meg, method = c("pearson"))
  cor_inset <- ggplot(ecStats, aes(fmri, meg)) + geom_point(size  = .8) + theme_light() + geom_smooth(method = "glm", alpha = .2, color = "grey") + 
    stat_cor(aes(label = ..r.label..), method = "pearson", cor.coef.name = "r", label.x = 0, label.y = cor_y_int, size = 2.8)
  # ggplot(ecStats, aes(plot_label, Statistic, color = session, groups = session)) + geom_point() + geom_line()
  
  setwd(file.path(fmri_dir, "plots"))
  
  
  pdf(paste0("echange_rt_vmax_by_vm_gradient17_", "_cox", coxme_ranef, "_dark.pdf"), height = 5, width = 10)
  print(ggplot(cox_df %>% filter(term == value_term) ) + 
          # geom_jitter(size = 6, width = .1, height = 0,  aes(vm_gradient17, Statistic, color = Statistic, alpha = p_level_fdr), show.legend = T) + 
          geom_jitter(size = 6, width = .1, height = 0,  aes(vm_gradient17_names, Statistic, color = vm_gradient17_names, alpha = p_level_fdr), show.legend = T) + 
          geom_violin(aes(vm_gradient17_names, color = vm_gradient17_names, Statistic), alpha = .2) + # scale_shape_manual(values = 21:25) +
          # scale_color_viridis(option = "mako") + theme_black() +  
          # scale_color_gradientn(colors = pal) + 
          scale_color_manual(values = colors ) +
          theme_minimal() + xlab(NULL) + ylab("Behavioral effect of neural response on value sensitivity, z statistic") + #labs(color = "Statistic") +
          geom_text_repel(aes(vm_gradient17_names, Statistic, alpha = p_level_fdr, label=plot_label, color = vm_gradient17_names), point.padding = 10, force = 10,  size = 2.4, show.legend = F) + 
          geom_text(aes(vm_gradient17_names, yint, label=vm_gradient17_names, color = vm_gradient17_names),  size = 3.5) + 
          facet_grid(~session, labeller = as_labeller(c("fMRI" = "fMRI","MEG replication" =  "MEG, out-of-session replication"))) + geom_hline(yintercept = 0, size = .2) + 
          theme(legend.key.size = unit(.4, "cm"), strip.text.x = element_text(size = 12), axis.text.x = element_blank()) +
          guides(alpha = guide_legend(title = expression(p[FDR], shape = guide_legend(override.aes = list(size = .1)))), color = "none") + 
          inset_element(echange_inset_fmri, .2, -.1, .45 , 0.25,  # left, bottom, right, top
                        align_to = "panel", on_top = T) + {if(coxme_ranef == "rslope")
                          inset_element(echange_inset_meg, 0.5, -0.1, 0.8 ,0.25,  # left, bottom, right, top
                                        align_to = "panel", on_top = T)} + 
          inset_element(cor_inset, 0.87 , -0.05, 1.05, 0.28)
  )
  dev.off()
  if (all_plots) {
    pdf(paste0("echange_rtvmax_by_vm_gradient17_cox_saggital_", session, ".pdf"), height = 12, width = 16)
    print(ggplot(cox_df %>% filter(term == value_term),
                 aes(y, z, alpha = p_level_fdr)) + 
            geom_point(size = 6, width = .1, height = 0,  aes(y, z, color = vm_gradient17_names, alpha = p_level_fdr), show.legend = F) + 
            scale_color_manual(values = colors ) +
            theme_minimal() +  #labs(color = "Statistic") +
            geom_text_repel(aes(y, z, alpha = p_level_fdr, label=plot_label, color = vm_gradient17_names), point.padding = 10, force = 10,  size = 2.4, show.legend = F) + 
            facet_grid(~session, labeller = as_labeller(c("fMRI" = "fMRI","MEG replication" =  "MEG, out-of-session replication"))) + geom_hline(yintercept = 0, size = .2) + 
            theme(legend.key.size = unit(.4, "cm"), strip.text.x = element_text(size = 12), axis.text.x = element_blank()) +
            guides(alpha = guide_legend(title = expression(p[FDR], shape = guide_legend(override.aes = list(size = .1))))))
    dev.off()
    pdf(paste0("echange_rtvmax_by_vm_gradient17_axial_", session, ".pdf"), height = 12, width = 16)
    print(ggplot(edf %>% filter(term == "fmri_beta:rt_vmax_lag" & model_name == "int"),  aes(y, -x, shape = vm_gradient17_names, alpha = p_level_fdr)) + 
            geom_point(size = 24, color = "white", aes(alpha = p_level_fdr, fill = statistic)) + scale_shape_manual(values = 21:25) + geom_point(size = 24, aes(color = vm_gradient17_names)) +
            scale_fill_viridis(option = "inferno") + theme_black() + geom_text_repel(aes(label=plot_label), color = "white", force = 5, point.padding = 40, size = 6) +
            guides(alpha = guide_legend(override.aes = list(color = "white"), title = expression(p[FDR]))))
    dev.off()
    pdf(paste0("echange_rtvmax_by_vm_gradient17_", session, ".pdf"), height = 8, width = 12)
    print(ggplot(edf %>% filter(term %in% c("fmri_beta:rt_vmax_lag")  & model_name == "int")) + 
            geom_jitter(size = 12, width = .1, height = 0,  aes(vm_gradient17_names, statistic, color = statistic, alpha = p_level_fdr)) + 
            geom_violin(aes(vm_gradient17_names, statistic), alpha = .2) + scale_shape_manual(values = 21:25) +
            scale_color_viridis(option = "inferno") + theme_black() + xlab("DAN region") + ylab("Convergence on RT(Vmax)") +
            geom_text_repel(aes(vm_gradient17_names, statistic, alpha = p_level_fdr, label=plot_label), point.padding = 20, force = 5, color="#4FC3F7", size = 3)  +
            guides(alpha = guide_legend(override.aes = list(color = "white"), title = expression(p[FDR]))))
    dev.off()
    
    pdf(paste0("echange_rtvmax_by_vm_gradient17_", session, "_rslope.pdf"), height = 8, width = 12)
    print(ggplot(edf %>% filter(term %in% c("fmri_beta:rt_vmax_lag")  & model_name == "slo")) + 
            geom_jitter(size = 12, width = .1, height = 0,  aes(vm_gradient17_names, statistic, color = statistic, alpha = p_level_fdr)) + 
            geom_violin(aes(vm_gradient17_names, statistic), alpha = .2) + scale_shape_manual(values = 21:25) +
            scale_color_viridis(option = "inferno") + theme_black() + xlab("DAN region") + ylab("Convergence on RT(Vmax)") +
            geom_text_repel(aes(vm_gradient17_names, statistic, alpha = p_level_fdr, label=plot_label), point.padding = 20, force = 5, color="#4FC3F7", size = 3)  +
            guides(alpha = guide_legend(override.aes = list(color = "white"), title = expression(p[FDR]))))
    dev.off()
    
  }
}
# }
