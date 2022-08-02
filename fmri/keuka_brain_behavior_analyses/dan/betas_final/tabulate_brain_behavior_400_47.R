library(dplyr)
library(tidyr)
library(data.table)
library(ggcorrplot)
# tabulate brain-behavior effects by DAN ROI
beta_dir <- "/Volumes/GoogleDrive/My Drive/SCEPTIC_fMRI/wholebrain_betas"
setwd(beta_dir)
clock_folder <- "~/code/clock_analysis" #alex
# just the necessary output is here:
fmri_dir <- file.path(paste0(clock_folder, "/fmri/keuka_brain_behavior_analyses/dan/betas_final/bb"))
get_correlations <- F


labels_df <- setDT(read_excel(file.path(paste0(clock_folder, "/fmri/keuka_brain_behavior_analyses/dan/MNH DAN Labels 400 Good Only 47 parcels.xlsx")))) %>%
  mutate(roi_num7 = as.factor(roi7_400), 
         mask_value = as.integer(roi7_400),
         plot_label = mnh_label_400, 
         vm_gradient17 = parcel_group) %>% select(roi_num7, mask_value, plot_label, gg_label, vm_gradient17, network17_400_DAN, hemi, x, y, z)

head(labels_df)

mixed_by_files <- list.files(pattern="Schaefer_400_DAN_manual_labels_47.*mixed_by.rds", 
                             full.names = TRUE, recursive=TRUE)

effects_of_interest <- c(
  # Rew="fmri_beta:last_outcomeReward",
  RTvmax="fmri_beta:rt_vmax_lag", 
  RtLag="rt_lag:fmri_beta",
  RtLagxRew="rt_lag:fmri_beta:last_outcomeReward"
)

for (model_type in c("int", "slo")) {
  dan_display <- lapply(mixed_by_files, function(x) {
    l1m <- sub("\\./([^/]+)/.*", "\\1", x, perl=TRUE)
    mout <- readRDS(x)$coef_df_reml %>%
      dplyr::filter(model_name == model_type & l2_cope_name == "overall") %>% # random intercept model
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
        TRUE ~ "NS"
      ),
      study = case_when(
        str_detect(x, "fmri") ~ "fmri",
        str_detect(x, "meg") ~ "meg"
      )
      ) %>%
      mutate(l1_model = !!l1m)
    
    
    mout$term_label <- factor(mout$term, levels = effects_of_interest, labels=names(effects_of_interest))
    
    dan_mout <- labels_df %>% left_join(mout, by="mask_value")
    #dan_mout <- labels_200 %>% left_join(mout, by="mask_value")
    
    # for display
    dan_mout_display <- dan_mout %>%
      dplyr::select(-estimate, -std.error, -df, -conf.low, -conf.high, -p.value, -padj_BY_term, -term)
    
    return(dan_mout_display)
  })
  
  
  #, l1_model # dump
  dan_df <- rbindlist(dan_display)
  dan_tab <- dan_df %>%
    # dplyr::filter(l1_model %in% c("L1m-abspe", "L1m-pe", "L1m-echange", "L1m-entropy_wiz")) %>%
    dplyr::mutate(l1_cope_name = factor(
      sub("EV_", "", l1_cope_name),
      # levels=c("entropy_change_feedback", "entropy_wiz_clock", "pe", "abspe"),
      # labels=c("echange", "entropy", "pe", "abspe")
      levels=c("entropy_change_feedback", "abspe"),
      labels=c("echange", "abspe")
    )) %>%
    dplyr::select(term_label, l1_cope_name, gg_label, mask_value, vm_gradient17, eff_size, study, hemi) %>%
    dplyr::filter(
      (term_label == "RTvmax") | ((term_label == "RtLag" | term_label == "RtLagxRew") & l1_cope_name == "abspe")) %>%
    dplyr::arrange(term_label, l1_cope_name) %>%
    pivot_wider(
      id_cols=c(mask_value, gg_label, vm_gradient17, hemi),
      #pivot_wider(id_cols=c(mask_value, plot_label), 
      names_from = c(term_label, l1_cope_name, study),
      names_glue = "{term_label}_{l1_cope_name}_{study}",
      values_from = eff_size, names_vary = "slowest") %>%
    dplyr::arrange(vm_gradient17, gg_label, hemi)
  fwrite(dan_tab, file=paste0("./tables/dan_400_47_bb_effects_", model_type, "_.csv"))
  
  dan_tab_stat <- dan_df %>%
    # dplyr::filter(l1_model %in% c("L1m-abspe", "L1m-pe", "L1m-echange", "L1m-entropy_wiz")) %>%
    dplyr::mutate(l1_cope_name = factor(
      sub("EV_", "", l1_cope_name),
      # levels=c("entropy_change_feedback", "entropy_wiz_clock", "pe", "abspe"),
      # labels=c("echange", "entropy", "pe", "abspe")
      levels=c("entropy_change_feedback", "abspe"),
      labels=c("echange", "abspe")
    )) %>%
    dplyr::select(term_label, l1_cope_name, gg_label, mask_value, vm_gradient17, statistic, study, hemi) %>%
    dplyr::filter(
      (term_label == "RTvmax") | ((term_label == "RtLag" | term_label == "RtLagxRew") & l1_cope_name == "abspe")) %>%
    dplyr::arrange(term_label, l1_cope_name) %>%
    pivot_wider(
      id_cols=c(mask_value, gg_label, vm_gradient17, hemi),
      #pivot_wider(id_cols=c(mask_value, plot_label), 
      names_from = c(term_label, l1_cope_name, study),
      names_glue = "{term_label}_{l1_cope_name}_{study}",
      values_from = statistic, names_vary = "slowest") %>%
    dplyr::arrange(vm_gradient17, gg_label, hemi)
  fwrite(dan_tab_stat, file=paste0("./tables/dan_400_47_bb_statistics_", model_type, "_.csv"))
  
  # View(dan_tab)
  
  if (get_correlations) {
    # look at correlations among effects across DAN nodes
    dan_corr <- dan_df %>%
      # dplyr::filter(l1_model %in% c("L1m-abspe", "L1m-pe", "L1m-echange", "L1m-entropy_wiz")) %>%
      # dplyr::mutate(l1_cope_name = factor(
      #   sub("EV_", "", l1_cope_name),
      #   levels=c("entropy_change_feedback", "entropy_wiz_clock", "pe", "abspe"),
      #   labels=c("echange", "entropy", "pe", "abspe")
      # )) %>%
      dplyr::filter(
        (term_label == "RTvmax") | ((term_label == "RtLag" | term_label == "RtLagxRew") & l1_cope_name == "abspe")) %>%
      dplyr::select(-l1_model, -eff_size) %>%
      dplyr::arrange(term_label, l1_cope_name) %>%
      pivot_wider(id_cols=c(mask_value, gg_label, vm_gradient17, hemi), 
                  #pivot_wider(id_cols=c(mask_value, plot_label), 
                  names_from = c(term_label, l1_cope_name, study),
                  names_glue = "{term_label}_{l1_cope_name}_{study}",
                  values_from = statistic, names_vary = "slowest")
    
    cmat <- dan_corr %>% select(matches("(Rew|RTvmax|RTlag)")) %>% cor()
    
    # cmat <- dan_corr %>% filter(study == "fmri") select(matches("(Rew|RTvmax|RTlag)")) %>% cor()
    
    pdf(paste0("dan_400_47_eff_corrs_", model_type, "_.pdf"), width=10, height=10)
    g <- ggcorrplot(cmat, hc.order = TRUE) + ggtitle("t-stat correlations of behavioral effects across DAN nodes") 
    plot(g)
    dev.off()
    
    
    dan_corr_wide <- dan_df %>%
      dplyr::filter(
        (term_label == "RTvmax") | ((term_label == "RtLag" | term_label == "RtLagxRew") & l1_cope_name == "abspe")) %>%
      # dplyr::filter(study == "fmri") %>%
      dplyr::select(gg_label, hemi, term_label, l1_cope_name, statistic, study) %>%
      #dplyr::select(-l1_model, -eff_size, -mask_value) %>%
      dplyr::arrange(term_label, l1_cope_name) %>%
      pivot_wider(
        id_cols=c(term_label, l1_cope_name, study),
        names_from = c(hemi, gg_label), # term_label, l1_cope_name),
        #names_glue = "{term_label}_{l1_cope_name}",
        values_from = statistic, names_vary = "slowest")
    
    
    cmat <- dan_corr_wide %>% select_if(is.numeric) %>% cor()
    pdf(paste0("dan_400_47_eff_corrs_", model_type, "_.pdf"), width=40, height=40)
    g <- ggcorrplot(cmat, hc.order = TRUE, type = "upper") + ggtitle("ROI t-stat correlations of behavioral effects across DAN nodes")
    plot(g)
    dev.off()
    
    # 
    # library(corrr)
    # cormat <- dan_corr_wide %>% select_if(is.numeric) %>% correlate() %>% network_plot(min_cor = 0.7)
    # 
    # cmat <- dan_corr_wide %>% select(-term_label, -l1_cope_name) %>% cor()
    # pdf("dan_roi_corrs.pdf", width=40, height=40)
    # g <- ggcorrplot(cmat, hc.order = TRUE, type = "upper") + ggtitle("ROI t-stat correlations of behavioral effects across DAN nodes")
    # plot(g)
    # dev.off()
  } 
}
# 
# 
# dan_df <- rbindlist(dan_display) %>%
#   #dplyr::filter(l1_model %in% c("L1m-abspe", "L1m-pe", "L1m-echange", "L1m-entropy_wiz")) %>%
#   dplyr::mutate(l1_cope_name = factor(
#     sub("EV_", "", l1_cope_name),
#     levels=c("entropy_change_feedback", "entropy_wiz_clock", "pe", "abspe"),
#     labels=c("echange", "entropy", "pe", "abspe")
#   )) %>%
#   dplyr::select(-l1_model) %>% #, -statistic
#   # inner_join(labels_df %>% select(gg_label, mask_value, hemi), by="mask_value") %>%
#   #inner_join(labels, by="mask_value") %>%
#   mutate(sort_label = paste(hemi, sub("([LR])_(.*)", "\\2_\\1", gg_label, perl=TRUE)))
# 
# dan_e <- dan_df %>% filter(l1_cope_name %in% c("echange"))# %>% filter(network7=="DorsAttn")
# dan_pe <- dan_df %>% filter(l1_cope_name %in% c("abspe"))# %>% filter(network7=="DorsAttn")
# 
# 
# library(patchwork)
