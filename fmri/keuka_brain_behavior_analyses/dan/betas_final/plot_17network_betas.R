# plot betas from the 17-network Schaeffer 2018 parcellation

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
source('~/code/Rhelpers/screen.lmerTest.R')
source('~/code/Rhelpers/vif.lme.R')
# library(stringi)
# session = "meg"
session = c("meg")
# clock_folder <- "~/Data_Analysis/clock_analysis" #michael
clock_folder <- "~/code/clock_analysis" #alex
# source('~/code/Rhelpers/')
fmri_dir <- '/Volumes/GoogleDrive/.shortcut-targets-by-id/1ukjK6kTlaR-LXIqX6nylYOPWu1j3XGyF/SCEPTIC_fMRI/wholebrain_betas'
# source("~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/get_trial_data.R")
labels_df <- setDT(read_excel("/Volumes/GoogleDrive/My Drive/SCEPTIC_fMRI/dan_medusa/schaefer_400_remap/MNH DAN Labels 400 Good Only 47 parcels.xlsx")) %>%
  mutate(roi_num7 = as.factor(roi7_400), 
         mask_value = as.integer(roi7_400),
         plot_label = mnh_label_400, 
         vm_gradient17 = parcel_group) %>% select(roi_num7, mask_value, plot_label, vm_gradient17, network17_400_DAN, hemi, x, y, z)

# plot PE betas across the 17-network parcellation of the 7-network DAN

setwd(file.path(fmri_dir, 'L1m-abspe'))
# brain to behavior
if (session == "meg") {
  temp <- readRDS("./parcel_maps/Schaefer_400_dan17_meg_test_mixed_by.rds")  
} else if (session == "fmri") {
  temp <- readRDS("./parcel_maps_l2/Schaefer_400_DAN_manual_labels__mixed_by.rds")
}
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
                                                                    ))) %>% inner_join(labels_df, by = "roi_num7") #%>%
  # filter( (net_num17==7 | net_num17==8))
  # filter(net_num7 ==3 & (net_num17==7 | net_num17==8))
  # filter(net_num7 ==3 ) 
df$p_level_fdr <- factor(df$p_level_fdr, levels = c('1', '2', '3', '4', '5', '6'), labels = c("NS","p < .05", "p < .01", "p < .001", "p < .0001", "p < .00001")) 

setwd(file.path(fmri_dir, "plots"))

ggplot(df %>% filter(term == "rt_lag:fmri_beta" & model_name == "int"),  aes(y, z, shape = network17_400_DAN)) + 
  geom_point(size = 20, color = "white", aes(alpha = p_level_fdr, fill = - statistic)) + scale_shape_manual(values = 21:25) + geom_point(size = 20, aes(color = network17_400_DAN)) +
  scale_fill_viridis(option = "mako") + theme_black() + geom_text_repel(aes(label=plot_label), color = "white", force = 5, point.padding = 40, size = 8)

pdf(paste0("pe_rt_swings_by_network_axial_", session, ".pdf"), height = 12, width = 16)
ggplot(df %>% filter(term == "rt_lag:fmri_beta" & model_name == "int"),  aes(y, x, shape = network17)) + 
  geom_point(size = 24, color = "white", aes(alpha = p_level_fdr, fill = - statistic)) + scale_shape_manual(values = 21:22) + geom_point(size = 24, aes(color = network17)) +
  scale_fill_viridis(option = "mako") + theme_black() + geom_text_repel(aes(label=subregion17), color = "white", force = 5, point.padding = 40, size = 8)
dev.off()

wdf <- df %>% filter(term %in% c("rt_lag:fmri_beta", "rt_lag:fmri_beta:last_outcomeReward")) %>% 
  pivot_wider(names_from = term, values_from = names(df)[9:17], names_sep = "_") %>% select_all(funs(gsub(":", "_", .))) %>%
  group_by(model_name, subregion17) %>%
  mutate(estimate_mean = (estimate_rt_lag_fmri_beta + estimate_rt_lag_fmri_beta_last_outcomeReward)/2
  )

pdf(paste0("pe_rt_swings_mean_saggital_", session, ".pdf"), height = 12, width = 16)
ggplot(wdf %>% filter(model_name == "int"),  aes(y, z, shape = network17)) + 
  geom_point(size = 32, color = "white", aes(alpha = p_level_fdr_rt_lag_fmri_beta_last_outcomeReward, 
                                             fill = - estimate_mean)) + 
  scale_shape_manual(values = 21:22) + geom_point(size = 32, color = "white") +
  scale_fill_viridis() + theme_black() + geom_text_repel(aes(label=subregion17),  color="blue") 
dev.off()

pdf(paste0("pe_rt_swings_omission_and_reward_saggital_", session, ".pdf"), height = 12, width = 24)
ggplot(df %>% filter(term %in% c("rt_lag:fmri_beta", "rt_lag:fmri_beta:last_outcomeReward")  & model_name == "int"),  aes(y, z, shape = network17)) + 
  geom_point(size = 32, color = "white", aes(alpha = p_level_fdr, fill = - estimate)) + scale_shape_manual(values = 21:22) + geom_point(size = 32, color = "white") +
  scale_fill_viridis() + theme_black() + geom_text_repel(aes(label=subregion17),  color="blue") + facet_wrap(~term)
dev.off()

pdf(paste0("pe_rt_swings_by_network_", session, ".pdf"), height = 8, width = 5)
ggplot(df %>% filter(term %in% c("rt_lag:fmri_beta")  & model_name == "int")) + 
  geom_jitter(size = 12, width = .1, height = 0,  aes(network17, - estimate, color = -estimate, alpha = p_level_fdr)) + 
  geom_violin(aes(network17, - estimate), alpha = .2) + scale_shape_manual(values = 21:22) +
  scale_color_viridis(option = "mako") + theme_black() + 
  geom_text_repel(aes(network17, - estimate, color = -estimate, alpha = p_level_fdr, label=subregion17),  color="blue")
dev.off()


# test a vs b

summary(lm(estimate ~ network17 + hemi, df %>% filter(term %in% c("rt_lag:fmri_beta"))))
summary(lm(estimate ~ network17 + hemi, df %>% filter(term %in% c("fmri_beta:rt_vmax_lag"))))

pdf(paste0("pe_rtvmax_by_network_saggital_", session, ".pdf"), height = 12, width = 16)
ggplot(df %>% filter(term == "fmri_beta:rt_vmax_lag" & model_name == "int"),  aes(y, z, shape = network17)) + 
  geom_point(size = 24, color = "white", aes(alpha = p_level_fdr, fill = statistic)) + scale_shape_manual(values = 21:22) + geom_point(size = 24, aes(color = network17)) +
  scale_fill_viridis(option = "inferno") + theme_black() + geom_text_repel(aes(label=subregion17), color = "white", force = 5, point.padding = 40, size = 8)
dev.off()
pdf(paste0("pe_rtvmax_by_network_axial_", session, ".pdf"), height = 12, width = 16)
ggplot(df %>% filter(term == "fmri_beta:rt_vmax_lag" & model_name == "int"),  aes(y, x, shape = network17)) + 
  geom_point(size = 24, color = "white", aes(alpha = p_level_fdr, fill = statistic)) + scale_shape_manual(values = 21:22) + geom_point(size = 24, aes(color = network17)) +
  scale_fill_viridis(option = "inferno") + theme_black() + geom_text_repel(aes(label=subregion17), color = "white", force = 5, point.padding = 40, size = 8)
dev.off()
# ggplot(df %>% filter(term == "rt_lag:fmri_beta" & model_name == "slo"),  aes(y, x, color = - estimate, shape = net_num17==8, alpha = abs(statistic)>1)) + 
#   geom_point(size = 32) + 
#   scale_color_viridis() + theme_black() + geom_text_repel(aes(label=subregion17), color="white")
# 
# ggplot(df %>% filter(term == "rt_lag:fmri_beta" & model_name == "int"),  aes(y, z, color = - estimate, shape = net_num17==8, alpha = abs(statistic))) + 
#   geom_point(size = 32) + 
#   scale_color_viridis() + theme_black() + geom_text_repel(aes(label=subregion17), color="white")


setwd(file.path(fmri_dir, 'L1m-echange'))
# here, mask_value = roi_num7
# brain to behavior
if (session == "meg") {
  temp <- readRDS("./parcel_maps/Schaefer_400_dan17_meg_test_mixed_by.rds")  
} else {
  temp <- readRDS("./parcel_maps/Schaefer_444_final_2009c_2.3mm_cope_l2_mixed_by.rds")
}

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
                                             ))) %>% inner_join(labels, by = "roi_num7") %>%
  # filter(net_num7 ==3 & (net_num17==7 | net_num17==8))
filter( (net_num17==7 | net_num17==8))
setwd(file.path(fmri_dir, "plots"))

# filter(net_num7 ==3 ) 
edf$p_level_fdr <- factor(edf$p_level_fdr, levels = c('1', '2', '3', '4', '5', '6'), labels = c("NS","p < .05", "p < .01", "p < .001", "p < .0001", "p < .00001")) 
pdf(paste0("echange_rtvmax_by_network_saggital_", session, ".pdf"), height = 12, width = 16)
ggplot(edf %>% filter(term == "fmri_beta:rt_vmax_lag" & model_name == "int"),  aes(y, z, shape = network17)) + 
  geom_point(size = 24, color = "white", aes(alpha = p_level_fdr, fill = statistic)) + scale_shape_manual(values = 21:22) + geom_point(size = 24, aes(color = network17)) +
  scale_fill_viridis(option = "inferno") + theme_black() + geom_text_repel(aes(label=subregion17), color = "white", force = 5, point.padding = 40, size = 8)
dev.off()
pdf(paste0("echange_rtvmax_by_network_axial_", session, ".pdf"), height = 12, width = 16)
ggplot(edf %>% filter(term == "fmri_beta:rt_vmax_lag" & model_name == "int"),  aes(y, x, shape = network17)) + 
  geom_point(size = 24, color = "white", aes(alpha = p_level_fdr, fill = statistic)) + scale_shape_manual(values = 21:22) + geom_point(size = 24, aes(color = network17)) +
  scale_fill_viridis(option = "inferno") + theme_black() + geom_text_repel(aes(label=subregion17), color = "white", force = 5, point.padding = 40, size = 8)
dev.off()
pdf(paste0("echange_rtvmax_by_network_violin_", session, ".pdf"), height = 8, width = 5)

ggplot(edf %>% filter(term %in% c("fmri_beta:rt_vmax_lag")  & model_name == "int")) + 
  geom_jitter(size = 12, width = .1, height = 0,  aes(network17, estimate, color = estimate, alpha = p_level_fdr)) + 
  geom_violin(aes(network17, estimate), alpha = .2) + scale_shape_manual(values = 21:22) +
  scale_color_viridis(option = "inferno") + theme_black() + 
  geom_text_repel(aes(network17, - estimate, color = estimate, alpha = p_level_fdr, label=subregion17),  color="blue")
dev.off()

