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
source("~/code/Rhelpers/theme_black.R")
# install_github("UNCDEPENdLab/dependlab")
# library(dependlab)
# library(stringi)
session = "meg"

signals = c("abspe", "echange")
# for (session in c("fmri", "meg")) {
# clock_folder <- "~/Data_Analysis/clock_analysis" #michael
clock_folder <- "~/code/clock_analysis" #alex

# # original output here:
# fmri_dir <- '/Volumes/GoogleDrive/.shortcut-targets-by-id/1ukjK6kTlaR-LXIqX6nylYOPWu1j3XGyF/SCEPTIC_fMRI/wholebrain_betas'

# just the necessary output is here:
fmri_dir <- file.path(paste0(clock_folder, "/fmri/keuka_brain_behavior_analyses/dan/betas_final/bb"))

# Manual labels from June-July 2022
# originals here
# labels_df <- setDT(read_excel("/Volumes/GoogleDrive/My Drive/SCEPTIC_fMRI/dan_medusa/schaefer_400_remap/MNH DAN Labels 400 Good Only 47 parcels.xlsx")) %>%

labels_df <- setDT(read_excel(file.path(paste0(fmri_dir, "/MNH DAN Labels 400 Good Only 47 parcels.xlsx")))) %>%
  mutate(roi_num7 = as.factor(roi7_400), 
         mask_value = as.integer(roi7_400),
         plot_label = mnh_label_400, 
         vm_gradient17 = parcel_group) %>% select(roi_num7, mask_value, plot_label, vm_gradient17, network17_400_DAN, hemi, x, y, z)

if ("abspe" %in% signals) {
  setwd(file.path(fmri_dir, 'L1m-abspe'))
  # plot brain to behavior results for absolute PEs
  if (session == "meg") {
    temp <- readRDS("./parcel_maps_l2/Schaefer_400_DAN_manual_labels_47_meg__mixed_by.rds")  
  } else if (session == "fmri") {
    ## save only the 47 DAN parcels from the 400
    # temp <- readRDS("./parcel_maps_l2/Schaefer_400_DAN_manual_labels__mixed_by.rds")
    # # library(purrr)
    # chosen47 <- unique(labels_df$mask_value)
    # temp47 <- Filter(Negate(is.null), temp) %>% purrr::map(., ~dplyr::filter(., (mask_value %in% chosen47)))  
    # saveRDS(temp47, file = "./parcel_maps_l2/Schaefer_400_DAN_manual_labels_47_fmri_mixed_by.rds")
    temp <- readRDS("./parcel_maps_l2/Schaefer_400_DAN_manual_labels_47_fmri_mixed_by.rds")
  }
  
  # import the results of brain-behavior analyses
  df <- temp$coef_df_reml %>% filter(term == "rt_lag:fmri_beta" | term == "fmri_beta:rt_vmax_lag" | term == "rt_lag:fmri_beta:last_outcomeReward") %>%
    filter(l1_cope_name == "EV_abspe") %>%
    rename(roi_num7 = "mask_value") %>% mutate(roi_num7 = as.factor(roi_num7),
                                               padj_BY_term = p.adjust(p.value, method = 'BH'),
                                               p_level_fdr = as.factor(case_when(
                                                 # p_fdr > .1 ~ '0',
                                                 # p_fdr < .1 & p_fdr > .05 ~ '1',
                                                 padj_BY_term > .05 ~ '1',
                                                 padj_BY_term < .05 & padj_BY_term > .01 ~ '2',
                                                 padj_BY_term < .01 & padj_BY_term > .001 ~ '3',
                                                 padj_BY_term <.001 & padj_BY_term > .0001 ~ '4',
                                                 padj_BY_term <.0001 & padj_BY_term > .00001 ~ '5',
                                                 padj_BY_term <.00001 ~ '6'
                                               )),
                                               reward = case_when(
                                                 term == "rt_lag:fmri_beta" ~ "Omission",
                                                 term == "rt_lag:fmri_beta:last_outcomeReward" ~ "Reward"
                                               )) %>% inner_join(labels_df, by = "roi_num7") 
  df$p_level_fdr <- factor(df$p_level_fdr, levels = c('1', '2', '3', '4', '5', '6'), labels = c("NS","p < .05", "p < .01", "p < .001", "p < .0001", "p < .00001")) 
  
  setwd(file.path(fmri_dir, "plots"))
  
  # # quick check
  # print(ggplot(df %>% filter(term == "rt_lag:fmri_beta" & model_name == "int"),  aes(y, z, shape = network17_400_DAN)) + 
  #   geom_point(size = 20, color = "white", aes(alpha = p_level_fdr, fill = - statistic)) + scale_shape_manual(values = 21:25) + geom_point(size = 20, aes(color = network17_400_DAN)) +
  #   scale_fill_viridis(option = "mako") + theme_black() + geom_text_repel(aes(label=plot_label), color = "white", force = 5, point.padding = 40, size = 8)
  # 
  # pdf(paste0("abspe_rt_swings_by_vm_gradient17_axial_", session, ".pdf"), height = 12, width = 16)
  # print(ggplot(df %>% filter(term == "rt_lag:fmri_beta" & model_name == "int"),  aes(y, -x, shape = vm_gradient17)) + 
  #   geom_point(size = 24, color = "white", aes(alpha = p_level_fdr, fill = - statistic)) + scale_shape_manual(values = 21:25) + geom_point(size = 24, aes(color = vm_gradient17)) +
  #   scale_fill_viridis(option = "mako") + theme_black() + geom_text_repel(aes(label=plot_label), color = "white", force = 5, point.padding = 40, size = 8) +
  #   guides(alpha = guide_legend(override.aes = list(color = "white"), title = expression(p[FDR]))))
  # dev.off()
  # 
  # wdf <- df %>% filter(term %in% c("rt_lag:fmri_beta", "rt_lag:fmri_beta:last_outcomeReward")) %>% 
  #   pivot_wider(names_from = term, values_from = names(df)[9:17], names_sep = "_") %>% select_all(funs(gsub(":", "_", .))) %>%
  #   group_by(model_name, plot_label) %>%
  #   mutate(statistic_mean = (statistic_rt_lag_fmri_beta + statistic_rt_lag_fmri_beta_last_outcomeReward)/2
  #   )
  # 
  # pdf(paste0("abspe_rt_swings_mean_saggital_", session, ".pdf"), height = 12, width = 16)
  # print(ggplot(wdf %>% filter(model_name == "int"),   aes(y, z, shape = vm_gradient17)) + 
  #   geom_point(size = 20, color = "white", aes(alpha = p_level_fdr, fill = - statistic)) + scale_shape_manual(values = 21:25) + geom_point(size = 20, color = "white") +
  #   scale_fill_viridis() + theme_black() + 
  #   # geom_text_repel(aes(alpha = p_level_fdr, label=plot_label), color = "white", force = 5, point.padding = 40, size = 6) + 
  #   geom_text_repel(aes(alpha = p_level_fdr, label=plot_label), point.padding = 40, force = 6, color="white", size = 7) 
  # dev.off()
  
  # main figure, champagne plot:
  pdf(paste0("abspe_rt_swings_by_vm_gradient17_", session, ".pdf"), height = 8, width = 16)
  print(ggplot(df %>% filter(term %in% c("rt_lag:fmri_beta", "rt_lag:fmri_beta:last_outcomeReward")  & model_name == "int") %>% mutate(`RT swings` = -statistic)) + 
    geom_jitter(size = 12, width = .1, height = 0,  aes(vm_gradient17, - statistic, color = `RT swings`, alpha = p_level_fdr), show.legend = T) + 
    geom_violin(aes(vm_gradient17, - statistic), alpha = .2) + scale_shape_manual(values = 21:25) +
    # scale_color_viridis(option = "mako") + theme_black() + 
    scale_color_viridis(option = "mako") + theme_black() + xlab("DAN region") + ylab("RT swing") +
    geom_text_repel(aes(vm_gradient17, -statistic, alpha = p_level_fdr, label=plot_label), point.padding = 20, force = 5, color="#4FC3F7", size = 3.5, show.legend = F) + 
    facet_wrap(~reward) + guides(alpha = guide_legend(override.aes = list(color = "white"), title = expression(p[FDR]))))
  # geom_text_repel(aes(vm_gradient17, - statistic, color = -statistic, alpha = p_level_fdr, label=plot_label),  color="blue")
  dev.off()
  
  # anatomical detail for supplemental figures:
  pdf(paste0("abspe_rt_swings_omission_and_reward_saggital_", session, ".pdf"), height = 12, width = 24)
  print(ggplot(df %>% filter(term %in% c("rt_lag:fmri_beta", "rt_lag:fmri_beta:last_outcomeReward")  & model_name == "int"),  aes(y, z, shape = vm_gradient17, alpha = p_level_fdr)) + 
    geom_point(size = 20, color = "white", aes(alpha = p_level_fdr, fill = - statistic)) + scale_shape_manual(values = 21:25) + geom_point(size = 20, color = "white") +
    scale_fill_viridis(option = "mako") + theme_black() + 
    # geom_text_repel(aes(alpha = p_level_fdr, label=plot_label), color = "white", force = 5, point.padding = 40, size = 6) + 
    geom_text_repel(aes(alpha = p_level_fdr, label=plot_label), point.padding = 40, force = 6, color="white", size = 6) +
    facet_wrap(~reward) + guides(alpha = guide_legend(override.aes = list(color = "white"), title = expression(p[FDR]))))
  dev.off()
  
  pdf(paste0("abspe_rt_swings_omission_and_reward_axial_", session, ".pdf"), height = 12, width = 24)
  print(ggplot(df %>% filter(term %in% c("rt_lag:fmri_beta", "rt_lag:fmri_beta:last_outcomeReward")  & model_name == "int"),  aes(y, -x, shape = vm_gradient17, alpha = p_level_fdr)) + 
    geom_point(size = 20, color = "white", aes(alpha = p_level_fdr, fill = - statistic)) + scale_shape_manual(values = 21:25) + geom_point(size = 20, color = "white") +
    scale_fill_viridis(option = "mako") + theme_black() + 
    # geom_text_repel(aes(alpha = p_level_fdr, label=plot_label), color = "white", force = 5, point.padding = 40, size = 6) + 
    geom_text_repel(aes(alpha = p_level_fdr, label=plot_label), point.padding = 40, force = 6, color="white", size = 6) +
    facet_wrap(~reward) + guides(alpha = guide_legend(override.aes = list(color = "white"), title = expression(p[FDR]))))
  dev.off()
  
  
  
  # # test
  # summary(lm(statistic ~ vm_gradient17 + hemi*term, df %>% filter(term %in% c("rt_lag:fmri_beta", "rt_lag:fmri_beta:last_outcomeReward"))))
  # summary(lm(statistic ~ vm_gradient17 + hemi, df %>% filter(term %in% c("fmri_beta:rt_vmax_lag"))))
  
  pdf(paste0("abspe_rtvmax_by_vm_gradient17_saggital_", session, ".pdf"), height = 12, width = 16)
  print(ggplot(df %>% filter(term == "fmri_beta:rt_vmax_lag" & model_name == "int"),  aes(y, z, shape = vm_gradient17, alpha = p_level_fdr)) + 
    geom_point(size = 24, color = "white", aes(alpha = p_level_fdr, fill = statistic)) + scale_shape_manual(values = 21:25) + geom_point(size = 24, aes(color = vm_gradient17)) +
    scale_fill_viridis(option = "inferno") + theme_black() + geom_text_repel(aes(label=plot_label), color = "white", force = 5, point.padding = 40, size = 6) +
    guides(alpha = guide_legend(override.aes = list(color = "white"), title = expression(p[FDR]))))
  dev.off()
  pdf(paste0("abspe_rtvmax_by_vm_gradient17_axial_", session, ".pdf"), height = 12, width = 16)
  print(ggplot(df %>% filter(term == "fmri_beta:rt_vmax_lag" & model_name == "int"),  aes(y, -x, shape = vm_gradient17, alpha = p_level_fdr)) + 
    geom_point(size = 24, color = "white", aes(alpha = p_level_fdr, fill = statistic)) + scale_shape_manual(values = 21:25) + geom_point(size = 24, aes(color = vm_gradient17)) +
    scale_fill_viridis(option = "inferno") + theme_black() + geom_text_repel(aes(label=plot_label), color = "white", force = 5, point.padding = 40, size = 6) +
    guides(alpha = guide_legend(override.aes = list(color = "white"), title = expression(p[FDR]))))
  dev.off()
  pdf(paste0("abspe_rtvmax_by_vm_gradient17_", session, ".pdf"), height = 8, width = 12)
  print(ggplot(df %>% filter(term %in% c("fmri_beta:rt_vmax_lag")  & model_name == "int")) + 
    geom_jitter(size = 12, width = .1, height = 0,  aes(vm_gradient17, statistic, color = statistic, alpha = p_level_fdr)) + 
    geom_violin(aes(vm_gradient17, statistic), alpha = .2) + scale_shape_manual(values = 21:25) +
    # scale_color_viridis(option = "mako") + theme_black() + 
    scale_color_viridis(option = "inferno") + theme_black() + xlab("DAN region") + ylab("Convergence on RT(Vmax)") +
    geom_text_repel(aes(vm_gradient17, statistic, alpha = p_level_fdr, label=plot_label), point.padding = 20, force = 5, color="#4FC3F7", size = 3)  +
    guides(alpha = guide_legend(override.aes = list(color = "white"), title = expression(p[FDR]))))
  # geom_text_repel(aes(vm_gradient17, - statistic, color = -statistic, alpha = p_level_fdr, label=plot_label),  color="blue")
  dev.off()
}

if ("echange" %in% signals) {
# plot entropy change exploitation effect
setwd(file.path(fmri_dir, 'L1m-echange'))
if (session == "meg") {
  temp <- readRDS("./parcel_maps_l2/Schaefer_400_DAN_manual_labels_47_meg__mixed_by.rds")  
} else if (session == "fmri") {
  # # save only the 47 DAN parcels from the 400
  # temp <- readRDS("./parcel_maps_l2/Schaefer_400_DAN_manual_labels__mixed_by.rds")
  # # library(purrr)
  # chosen47 <- unique(labels_df$mask_value)
  # temp47 <- Filter(Negate(is.null), temp) %>% purrr::map(., ~dplyr::filter(., (mask_value %in% chosen47)))
  # saveRDS(temp47, file = "./parcel_maps_l2/Schaefer_400_DAN_manual_labels_47_fmri_mixed_by.rds")
  temp <- readRDS("./parcel_maps_l2/Schaefer_400_DAN_manual_labels_47_fmri_mixed_by.rds")
}

  # name dataframe edf for Echange just in case to distinguish from abspe
edf <- temp$coef_df_reml %>% filter(term == "rt_lag:fmri_beta" | term == "fmri_beta:rt_vmax_lag" | term == "rt_lag:fmri_beta:last_outcomeReward") %>%
  filter(l1_cope_name == "EV_entropy_change_feedback") %>%
  rename(roi_num7 = "mask_value") %>% mutate(roi_num7 = as.factor(roi_num7),
                                             padj_BY_term = p.adjust(p.value, method = 'BH'),
                                             p_level_fdr = as.factor(case_when(
                                               # p_fdr > .1 ~ '0',
                                               # p_fdr < .1 & p_fdr > .05 ~ '1',
                                               padj_BY_term > .05 ~ '1',
                                               padj_BY_term < .05 & padj_BY_term > .01 ~ '2',
                                               padj_BY_term < .01 & padj_BY_term > .001 ~ '3',
                                               padj_BY_term <.001 & padj_BY_term > .0001 ~ '4',
                                               padj_BY_term <.0001 & padj_BY_term > .00001 ~ '5',
                                               padj_BY_term <.00001 ~ '6'
                                             ))) %>% inner_join(labels_df, by = "roi_num7") #%>%
setwd(file.path(fmri_dir, "plots"))

edf$p_level_fdr <- factor(edf$p_level_fdr, levels = c('1', '2', '3', '4', '5', '6'), labels = c("NS","p < .05", "p < .01", "p < .001", "p < .0001", "p < .00001")) 

pdf(paste0("echange_rtvmax_by_vm_gradient17_saggital_", session, ".pdf"), height = 12, width = 16)
print(ggplot(edf %>% filter(term == "fmri_beta:rt_vmax_lag" & model_name == "int"),  aes(y, z, shape = vm_gradient17, alpha = p_level_fdr)) + 
  geom_point(size = 24, color = "white", aes(alpha = p_level_fdr, fill = statistic)) + scale_shape_manual(values = 21:25) + geom_point(size = 24, aes(color = vm_gradient17)) +
  scale_fill_viridis(option = "inferno") + theme_black() + geom_text_repel(aes(label=plot_label), color = "white", force = 5, point.padding = 40, size = 6) +
  guides(alpha = guide_legend(override.aes = list(color = "white"), title = expression(p[FDR]))))
dev.off()
pdf(paste0("echange_rtvmax_by_vm_gradient17_axial_", session, ".pdf"), height = 12, width = 16)
print(ggplot(edf %>% filter(term == "fmri_beta:rt_vmax_lag" & model_name == "int"),  aes(y, -x, shape = vm_gradient17, alpha = p_level_fdr)) + 
  geom_point(size = 24, color = "white", aes(alpha = p_level_fdr, fill = statistic)) + scale_shape_manual(values = 21:25) + geom_point(size = 24, aes(color = vm_gradient17)) +
  scale_fill_viridis(option = "inferno") + theme_black() + geom_text_repel(aes(label=plot_label), color = "white", force = 5, point.padding = 40, size = 6) +
  guides(alpha = guide_legend(override.aes = list(color = "white"), title = expression(p[FDR]))))
dev.off()
pdf(paste0("echange_rtvmax_by_vm_gradient17_", session, ".pdf"), height = 8, width = 12)
print(ggplot(edf %>% filter(term %in% c("fmri_beta:rt_vmax_lag")  & model_name == "int")) + 
  geom_jitter(size = 12, width = .1, height = 0,  aes(vm_gradient17, statistic, color = statistic, alpha = p_level_fdr)) + 
  geom_violin(aes(vm_gradient17, statistic), alpha = .2) + scale_shape_manual(values = 21:25) +
  scale_color_viridis(option = "inferno") + theme_black() + xlab("DAN region") + ylab("Convergence on RT(Vmax)") +
  geom_text_repel(aes(vm_gradient17, statistic, alpha = p_level_fdr, label=plot_label), point.padding = 20, force = 5, color="#4FC3F7", size = 3)  +
  guides(alpha = guide_legend(override.aes = list(color = "white"), title = expression(p[FDR]))))
dev.off()
}
# }
