library(dplyr)
library(tidyr)
library(data.table)
library(ggcorrplot)
# tabulate brain-behavior effects by DAN ROI
beta_dir <- "/Volumes/GoogleDrive/My Drive/SCEPTIC_fMRI/wholebrain_betas"
setwd(beta_dir)

# MNH/AD internal labeling scheme
labels <- readxl::read_excel("~/Data_Analysis/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/MNH Dan Labels.xlsx") %>%
  dplyr::rename(mask_value=roinum) %>% select(mask_value, plot_label, Stream, Visuomotor_Gradient)

labels_200 <- read.csv(file.path(beta_dir, "schaefer_200_whereami.csv")) %>%
  dplyr::rename(mask_value = roi_num, plot_label = MNI_Glasser_HCP_v1.0) %>% dplyr::select(mask_value, plot_label) %>%
  mutate(plot_label = make.unique(sub("Focus point:\\s+", "", plot_label, perl=TRUE)))

head(labels)


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