library(dplyr)
library(tidyr)
library(data.table)
library(ggcorrplot)
# tabulate brain-behavior effects by DAN ROI
beta_dir <- "/Volumes/GoogleDrive/My Drive/SCEPTIC_fMRI/wholebrain_betas"
setwd(beta_dir)

# MNH/AD internal labeling scheme
labels <- readxl::read_excel("~/Data_Analysis/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/MNH Dan Labels.xlsx") %>%
  dplyr::rename(mask_value=roinum) %>% select(mask_value, plot_label, Stream, Visuomotor_Gradient, lobe)

labels_200 <- read.csv(file.path(beta_dir, "schaefer_200_whereami.csv")) %>%
  dplyr::rename(mask_value = roi_num, plot_label = MNI_Glasser_HCP_v1.0) %>% dplyr::select(mask_value, plot_label) %>%
  mutate(plot_label = make.unique(sub("Focus point:\\s+", "", plot_label, perl=TRUE)))

head(labels)


# mixed_by_files <- list.files(pattern="Schaefer2018_200Parcels_7Networks_order_fonov_2.3mm_ants_cope_l2_mixed_by.rds", 
#                              full.names = TRUE, recursive=TRUE)

mixed_by_files <- list.files(pattern="Schaefer2018_200Parcels_7Networks_order_fonov_2.3mm_ants_cope_l2_mixed_by.rds", 
                             full.names = TRUE, recursive=TRUE)


effects_of_interest <- c(
  Rew="fmri_beta:last_outcomeReward",
  RTvmax="fmri_beta:rt_vmax_lag", 
  RtLag="rt_lag:fmri_beta",
  RtLagxRew="rt_lag:fmri_beta:last_outcomeReward"
)


dan_display <- lapply(mixed_by_files, function(x) {
  l1m <- sub("\\./([^/]+)/.*", "\\1", x, perl=TRUE)
  mout <- readRDS(x)$coef_df_reml %>%
    dplyr::filter(model_name == "int" & l2_cope_name == "overall") %>% # random intercept model
    dplyr::filter(!grepl("(EV_clock|EV_feedback)", l1_cope_name)) %>%
    dplyr::select(-outcome, -rhs, -effect, -l2_cope_name, -model_name) %>%
    dplyr::filter(term %in% effects_of_interest) %>%
    mutate(eff_size = case_when(
      statistic > 2 & statistic < 3 ~ "+",
      statistic > 3 & statistic < 5 ~ "++",
      statistic > 5 ~ "+++",
      statistic < -2 & statistic > -3 ~ "-",
      statistic < -3 & statistic > -5 ~ "--",
      statistic < -5 ~ "---",
      TRUE ~ ""
    )) %>%
    mutate(l1_model = !!l1m)
  
  
  mout$term_label <- factor(mout$term, levels = effects_of_interest, labels=names(effects_of_interest))
  
  dan_mout <- labels %>% left_join(mout, by="mask_value")
  #dan_mout <- labels_200 %>% left_join(mout, by="mask_value")
  
  # for display
  dan_mout_display <- dan_mout %>%
    dplyr::select(-estimate, -std.error, -df, -conf.low, -conf.high, -p.value, -padj_BY_term, -term)
  
  return(dan_mout_display)
})


#, l1_model # dump
dan_df <- rbindlist(dan_display)
dan_tab <- dan_df %>%
  dplyr::filter(l1_model %in% c("L1m-abspe", "L1m-pe", "L1m-echange", "L1m-entropy_wiz")) %>%
  dplyr::mutate(l1_cope_name = factor(
    sub("EV_", "", l1_cope_name),
    levels=c("entropy_change_feedback", "entropy_wiz_clock", "pe", "abspe"),
    labels=c("echange", "entropy", "pe", "abspe")
  )) %>%
  dplyr::select(-l1_model, -statistic) %>%
  dplyr::arrange(term_label, l1_cope_name) %>%
  pivot_wider(id_cols=c(mask_value, plot_label, Stream, Visuomotor_Gradient),
  #pivot_wider(id_cols=c(mask_value, plot_label), 
              names_from = c(term_label, l1_cope_name),
              names_glue = "{term_label}_{l1_cope_name}",
              values_from = eff_size, names_vary = "slowest")


# look at correlations among effects across DAN nodes
dan_corr <- dan_df %>%
  dplyr::filter(l1_model %in% c("L1m-abspe", "L1m-pe", "L1m-echange", "L1m-entropy_wiz")) %>%
  dplyr::mutate(l1_cope_name = factor(
    sub("EV_", "", l1_cope_name),
    levels=c("entropy_change_feedback", "entropy_wiz_clock", "pe", "abspe"),
    labels=c("echange", "entropy", "pe", "abspe")
  )) %>%
  dplyr::select(-l1_model, -eff_size) %>%
  dplyr::arrange(term_label, l1_cope_name) %>%
  pivot_wider(id_cols=c(mask_value, plot_label, Stream, Visuomotor_Gradient), 
  #pivot_wider(id_cols=c(mask_value, plot_label), 
              names_from = c(term_label, l1_cope_name),
              names_glue = "{term_label}_{l1_cope_name}",
              values_from = statistic, names_vary = "slowest")
  
cmat <- dan_corr %>% select(matches("(Rew|RTvmax|RTlag)")) %>% cor()
pdf("dan_eff_corrs.pdf", width=10, height=10)
g <- ggcorrplot(cmat, hc.order = TRUE) + ggtitle("t-stat correlations of behavioral effects across DAN nodes")
plot(g)
dev.off()


dan_corr_wide <- dan_df %>%
  dplyr::filter(l1_model %in% c("L1m-abspe", "L1m-pe", "L1m-echange", "L1m-entropy_wiz")) %>%
  dplyr::mutate(l1_cope_name = factor(
    sub("EV_", "", l1_cope_name),
    levels=c("entropy_change_feedback", "entropy_wiz_clock", "pe", "abspe"),
    labels=c("echange", "entropy", "pe", "abspe")
  )) %>%
  dplyr::select(-l1_model, -eff_size, -mask_value, -Stream, -Visuomotor_Gradient) %>%
  #dplyr::select(-l1_model, -eff_size, -mask_value) %>%
  dplyr::arrange(term_label, l1_cope_name) %>%
  pivot_wider(id_cols=c(term_label, l1_cope_name), 
              names_from = c(plot_label), # term_label, l1_cope_name),
              #names_glue = "{term_label}_{l1_cope_name}",
              values_from = statistic, names_vary = "slowest")


cmat <- dan_corr_wide %>% select(-term_label, -l1_cope_name) %>% cor()
pdf("dan_roi_corrs.pdf", width=40, height=40)
g <- ggcorrplot(cmat, hc.order = TRUE) + ggtitle("ROI t-stat correlations of behavioral effects across DAN nodes")
plot(g)
dev.off()


fwrite(dan_tab, file="dan_bb_effects_rint.csv")


## Look at effects by subnetworks within Yeo 17
schaefer_dir <- "~/Data_Analysis/schaefer_wb_parcellation"

schaefer_7 <- read.csv(file.path(schaefer_dir, "labels", "Schaefer2018_400Parcels_7Networks_order.csv")) %>%
  mutate(network=factor(network), net_num = as.numeric(network)) %>%
  rename(network7=network, net_num7=net_num, subregion7=subregion)

# this has the spatial coordinate, spatial_roi_num
schaefer_7_lookup <- read.csv(file.path(schaefer_dir, "labels", "Schaefer_400_7networks_labels.csv"))

schaefer_7 <- schaefer_7 %>% inner_join(schaefer_7_lookup, by="roi_num") %>%
  rename(roi_num7=roi_num)

schaefer_17 <- read.csv(file.path(schaefer_dir, "labels", "Schaefer2018_400Parcels_17Networks_order.csv")) %>%
  mutate(network=factor(network), net_num = as.numeric(network)) %>%
  rename(network17=network, net_num17=net_num, subregion17=subregion) %>%
  select(-hemi) # mirrored in 7

# this has the spatial coordinate, spatial_roi_num
schaefer_17_lookup <- read.csv(file.path(schaefer_dir, "labels", "Schaefer_400_17networks_labels.csv")) %>%
  select(roi_num, spatial_roi_num) # x,y,z and labels already duplicated in 7-network lookup

schaefer_17 <- schaefer_17 %>% inner_join(schaefer_17_lookup, by="roi_num") %>%
  rename(roi_num17=roi_num)

#schaefer_444 <- 

both <- schaefer_7 %>%
  select(-hemi) %>%
  inner_join(schaefer_17, by="spatial_roi_num") %>%
  rename(mask_value=roi_num7) %>%
  filter(network7 %in% c("DorsAttn", "SalVentAttn", "Cont", "Vis", "SomMot", "Default")) %>% 
  #filter(network7 == "DorsAttn") %>% 
  #filter(network17 == "DorsAttnA" | network17 == "DorsAttnB") %>%
  mutate(plot_label = make.unique(sub("Focus point:\\s+", "", MNI_Glasser_HCP_v1.0, perl=TRUE))) %>%
  mutate(lobe = if_else(x < 0, "left", "right"))

beta_dir <- "/Volumes/GoogleDrive/My Drive/SCEPTIC_fMRI/wholebrain_betas/tmp_444"
#setwd(beta_dir)

mixed_by_files <- list.files(path = beta_dir, pattern=".*Schaefer_444_final_2009c_2.3mm_cope_l2_mixed_by.rds", 
                             full.names = TRUE, recursive=TRUE)

effects_of_interest <- c(
  Rew="fmri_beta:last_outcomeReward",
  RTvmax="fmri_beta:rt_vmax_lag", 
  RtLag="rt_lag:fmri_beta",
  RtLagxRew="rt_lag:fmri_beta:last_outcomeReward"
)

dan_display <- lapply(mixed_by_files, function(x) {
  l1m <- sub("\\./([^/]+)/.*", "\\1", x, perl=TRUE)
  mout <- readRDS(x)$coef_df_reml %>%
    dplyr::filter(model_name == "int" & l2_cope_name == "overall") %>% # random intercept model
    dplyr::filter(!grepl("(EV_clock|EV_feedback)", l1_cope_name)) %>%
    dplyr::select(-outcome, -rhs, -effect, -l2_cope_name, -model_name) %>%
    dplyr::filter(term %in% effects_of_interest) %>%
    mutate(l1_model = !!l1m)
  return(mout)
})

dan_df <- rbindlist(dan_display) %>%
  #dplyr::filter(l1_model %in% c("L1m-abspe", "L1m-pe", "L1m-echange", "L1m-entropy_wiz")) %>%
  dplyr::mutate(l1_cope_name = factor(
    sub("EV_", "", l1_cope_name),
    levels=c("entropy_change_feedback", "entropy_wiz_clock", "pe", "abspe"),
    labels=c("echange", "entropy", "pe", "abspe")
  )) %>%
  dplyr::select(-l1_model) %>% #, -statistic
  inner_join(both, by="mask_value") %>%
  #inner_join(labels, by="mask_value") %>%
  mutate(sort_label = paste(lobe, sub("([LR])_(.*)", "\\2_\\1", plot_label, perl=TRUE)))

dan_e <- dan_df %>% filter(l1_cope_name %in% c("echange", "entropy"))# %>% filter(network7=="DorsAttn")
dan_pe <- dan_df %>% filter(l1_cope_name %in% c("abspe", "pe"))# %>% filter(network7=="DorsAttn")


library(patchwork)

pdf("dan_bb_17_dissociation_400.pdf", width=32, height=16)
                           
g1 <- ggplot(dan_e, aes(x=sort_label, y=estimate, ymin=estimate-std.error, ymax=estimate+std.error, color=l1_cope_name)) +
  facet_grid(term~network17, scales = "free") +
  geom_pointrange(position=position_dodge(width=0.7), size=1.5) +
  theme_bw(base_size = 18) +
  theme(axis.text.x = element_text(angle=90)) +
  scale_color_brewer(palette="Set2") +
  geom_hline(yintercept = 0)

g2 <- ggplot(dan_pe, aes(x=sort_label, y=estimate, ymin=estimate-std.error, ymax=estimate+std.error, color=l1_cope_name)) +
  facet_grid(term~network17, scales = "free") +
  geom_pointrange(position=position_dodge(width=0.7), size=1.5) +
  theme_bw(base_size = 18) +
  theme(axis.text.x = element_text(angle=90)) +
  geom_hline(yintercept = 0)

g1 + g2 + plot_layout(nrow=1)
dev.off()

pdf("dan_bb_17_dissociation_400_pe.pdf", width=32, height=16)
g2 <- ggplot(dan_pe, aes(x=sort_label, y=estimate, ymin=estimate-std.error, ymax=estimate+std.error, color=l1_cope_name)) +
  facet_grid(term~network17, scales = "free") +
  geom_pointrange(position=position_dodge(width=0.7), size=1.5) +
  theme_bw(base_size = 18) +
  theme(axis.text.x = element_text(angle=90)) +
  geom_hline(yintercept = 0)

g2
dev.off()


ggplot(dan_pe %>% filter(term=="rt_lag:fmri_beta" & l1_cope_name=="pe"), aes(x=y, y=x, color=estimate)) + geom_point(size=6) +
  scale_color_viridis_c() + brms::theme_black()

g1 <- ggplot(dan_pe %>% filter(term=="rt_lag:fmri_beta" & l1_cope_name=="pe"), 
       aes(x=y, y=z, color=-statistic, shape=network7, alpha=padj_BY_term < .01)) + geom_point(size=15) +
  scale_color_viridis_c() + brms::theme_black() +
  geom_text_repel(aes(label=subregion17), color="grey80")

g2 <- ggplot(dan_pe %>% filter(term=="rt_lag:fmri_beta" & l1_cope_name=="abspe"), 
       aes(x=y, y=z, color=-statistic, shape=network7, alpha=padj_BY_term < .01)) + geom_point(size=15) +
  scale_color_viridis_c() + brms::theme_black() +
  geom_text_repel(aes(label=subregion17), color="grey80")

g1 + g2 + plot_layout(ncol=2)

ggplot(dan_pe %>% filter(term=="rt_lag:fmri_beta" & l1_cope_name=="abspe"), 
       aes(x=y, y=z, color=-estimate, shape=lobe, alpha=padj_BY_term < .05)) + geom_point(size=15) +
  scale_color_viridis_c() + brms::theme_black() +
  geom_text_repel(aes(label=subregion17), color="grey80")


ggplot(dan_pe %>% filter(term=="rt_lag:fmri_beta:last_outcomeReward" & l1_cope_name=="abspe"), 
       aes(x=y, y=z, color=-estimate, shape=lobe, alpha=padj_BY_term < .05)) + geom_point(size=15) +
  scale_color_viridis_c() + brms::theme_black() +
  geom_text_repel(aes(label=subregion17), color="grey80")


g1 <- ggplot(dan_pe %>% filter(term=="fmri_beta:rt_vmax_lag" & l1_cope_name=="abspe"), 
       aes(x=y, y=z, color=statistic, shape=network7, alpha=abs(statistic) > 3)) + geom_point(size=15) +
  scale_color_viridis_c() + brms::theme_black() +
  geom_text(aes(label=subregion17), color="grey80") + ggtitle("abspe rtvmax")
  #scale_y_reverse()

g2 <- ggplot(dan_pe %>% filter(term=="fmri_beta:rt_vmax_lag" & l1_cope_name=="pe"), 
       aes(x=y, y=z, color=statistic, shape=network7, alpha=abs(statistic) > 3)) + geom_point(size=15) +
  scale_color_viridis_c() + brms::theme_black() +
  geom_text(aes(label=subregion17), color="grey80") + ggtitle("PE rtvmax")
  #scale_y_reverse()

g3 <- ggplot(dan_pe %>% filter(term=="rt_lag:fmri_beta" & l1_cope_name=="abspe"), 
             aes(x=y, y=z, color=statistic, shape=network17, alpha=abs(statistic) > 3)) + geom_point(size=15) +
  scale_color_viridis_c() + brms::theme_black() +
  geom_text(aes(label=subregion17), color="grey80") + ggtitle("abspe rt_lag")
  #scale_y_reverse()

g4 <- ggplot(dan_pe %>% filter(term=="rt_lag:fmri_beta" & l1_cope_name=="pe"), 
             aes(x=y, y=z, color=statistic, shape=network17, alpha=abs(statistic) > 3)) + geom_point(size=15) +
  scale_color_viridis_c() + brms::theme_black() +
  geom_text(aes(label=subregion17), color="grey80") + ggtitle("PE rt_lag")
  #scale_y_reverse()


g1 + g2 + g3 + g4 +  plot_layout(nrow=2)


###


g1 <- ggplot(dan_e %>% filter(term=="fmri_beta:rt_vmax_lag" & l1_cope_name=="echange" & network7!="SomMot"), 
             aes(x=y, y=z, fill=statistic, alpha=padj_BY_term < .05, shape=network7)) + geom_point(size=15) +
  scale_fill_viridis_c() + theme_dark() + scale_shape_manual(values=21:26) +
  geom_text(aes(label=subregion17), color="black") + ggtitle("Echange rtvmax")


g2 <- ggplot(dan_pe %>% filter(term=="fmri_beta:rt_vmax_lag" & l1_cope_name=="abspe" & network7 != "SomMot"), 
             aes(x=y, y=z, fill=statistic, alpha=padj_BY_term < .05, shape=network7)) + geom_point(size=15) +
  scale_fill_viridis_c() + theme_dark() + scale_shape_manual(values=21:26) +
  geom_text(aes(label=subregion17), color="black") + ggtitle("absPE rtvmax")

g1 + g2


g2 <- ggplot(dan_pe %>% filter(term=="rt_lag:fmri_beta" & l1_cope_name=="abspe"), 
             aes(x=y, y=x, fill=statistic, alpha=padj_BY_term < .01, shape=network7)) + geom_point(size=15) +
  scale_fill_viridis_c() + theme_dark() +
  geom_text(aes(label=subregion17), color="black") + ggtitle("abspe rtlag") +
  scale_shape_manual(values=21:26)


g2 <- ggplot(dan_pe %>% filter(term=="rt_lag:fmri_beta" & l1_cope_name=="pe"), 
             aes(x=y, y=z, color=statistic, alpha=padj_BY_term < .05)) + geom_point(size=15) +
  scale_color_viridis_c() + theme_dark() + facet_wrap(~network7) +
  geom_text(aes(label=subregion17), color="black") + ggtitle("pe rtlag")


##


g2 <- ggplot(dan_pe %>% filter(term=="fmri_beta:rt_vmax_lag" & l1_cope_name=="abspe"), 
             aes(x=y, y=x, color=statistic, alpha=abs(statistic) > 3)) + geom_point(size=15) +
  scale_color_viridis_c() + brms::theme_black() + facet_wrap(~network7) +
  geom_text(aes(label=subregion17), color="grey80") + ggtitle("absPE rtvmax")

g3 <- ggplot(dan_pe %>% filter(term=="rt_lag:fmri_beta" & l1_cope_name=="abspe"), 
             aes(x=y, y=x, color=statistic, alpha=abs(statistic) > 3)) + geom_point(size=15) +
  scale_color_viridis_c() + brms::theme_black() + facet_wrap(~network7) +
  geom_text(aes(label=subregion17), color="grey80") + ggtitle("absPE rt lag")

g1 + g2 + g3 + plot_layout(nrow=3)
  
g3 <- ggplot(dan_e %>% filter(term=="rt_lag:fmri_beta" & l1_cope_name=="echange"), 
             aes(x=y, y=x, color=statistic, shape=network7, alpha=abs(statistic) > 3)) + geom_point(size=15) +
  scale_color_viridis_c() + brms::theme_black() +
  geom_text(aes(label=subregion17), color="grey80") + ggtitle("Echange rt_lag")

g4 <- ggplot(dan_pe %>% filter(term=="rt_lag:fmri_beta" & l1_cope_name=="pe"), 
             aes(x=y, y=x, color=statistic, shape=network7, alpha=abs(statistic) > 3)) + geom_point(size=15) +
  scale_color_viridis_c() + brms::theme_black() +
  geom_text(aes(label=subregion17), color="grey80") + ggtitle("PE rt_lag")


g1 + g2 + g3 + g4 +  plot_layout(nrow=2)


###
ggplot(dan_e %>% filter(term=="rt_lag:fmri_beta" & l1_cope_name=="echange"), 
       aes(x=y, y=x, color=statistic, shape=lobe, alpha=padj_BY_term < .05)) + geom_point(size=15) +
  scale_color_viridis_c() + brms::theme_black() +
  geom_text_repel(aes(label=subregion17), color="grey80")



ggplot(dan_e %>% filter(term=="fmri_beta:rt_vmax_lag" & l1_cope_name=="entropy"), 
       aes(x=y, y=z, color=statistic, shape=lobe, alpha=padj_BY_term < .05)) + geom_point(size=15) +
  scale_color_viridis_c() + brms::theme_black() +
  geom_text_repel(aes(label=subregion17), color="grey80")



#shape=network7

x <- dan_pe %>% filter(term=="rt_lag:fmri_beta")

## meta mixed by comparing network mean B-B effects
library(afex)
m <- lm(estimate ~ network17, dan_pe %>% filter(l1_cope_name=="abspe" & term=="rt_lag:fmri_beta"))
m1 <- lm(estimate ~ network17 + hemi, dan_pe %>% filter(l1_cope_name=="abspe" & term=="rt_lag:fmri_beta"))
m <- lm(estimate ~ network17, dan_pe %>% filter(l1_cope_name=="pe" & term=="rt_lag:fmri_beta"))
summary(m1)
anova(m, m1)
library(brms)
ff2 <- brm(
  estimate | se(std.error, sigma=TRUE) ~ 1 + network17 + (1 | plot_label), 
  dan_pe %>% filter(l1_cope_name=="abspe" & term=="rt_lag:fmri_beta"),
  chains= 4, cores=4, iter = 15000,
  control = list(
    adapt_delta = 0.999,
    max_treedepth = 13
  )
)

library(emmeans)
summary(ff2)
emmeans(ff2, ~network17)
pairs(emmeans(ff2, ~network17))


ff3 <- brm(
  estimate | se(std.error, sigma=TRUE) ~ 1 + network17 + (1 | plot_label), 
  dan_pe %>% filter(l1_cope_name=="pe" & term=="rt_lag:fmri_beta"),
  chains= 4, cores=4, iter = 15000,
  control = list(
    adapt_delta = 0.999,
    max_treedepth = 13
  )
)

summary(ff3)
emmeans(ff3, ~network17)
pairs(emmeans(ff3, ~network17))





### leftovers
library(ggpubr)
library(factoextra)
to_cluster <- scale(dan_corr[, -1:-4])
fviz_nbclust(to_cluster, kmeans, nstart=100, method="wss")
res.km3 <- kmeans(scale(to_cluster), 3, nstart = 25)
res.km4 <- kmeans(scale(to_cluster), 4, nstart = 25)
res.km5 <- kmeans(scale(to_cluster), 5, nstart = 25)

rownames(to_cluster) <- dan_corr$plot_label
pdf("roi_kmeans_plots.pdf", width=12, height=10)
fviz_cluster(res.km3, data = to_cluster)
fviz_cluster(res.km4, data = to_cluster)
fviz_cluster(res.km5, data = to_cluster)
dev.off()
res.km$cluster


fviz_cluster(res.km, data = to_cluster,
             palette = c("#2E9FDF", "#00AFBB", "#E7B800"), 
             geom = "point",
             ellipse.type = "convex", 
             ggtheme = theme_bw()
)


