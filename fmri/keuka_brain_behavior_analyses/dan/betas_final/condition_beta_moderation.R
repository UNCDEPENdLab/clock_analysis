setwd("/Users/hallquist/Data_Analysis/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/betas_final")

# rt_lag <- echange$emtrends_list$rt_lag_rewFunc %>%
#   mutate(fmri_fac = factor(fmri_beta, levels=c(-2, 0, 2), labels=c("-2SD", "Mean", "+2SD")))
labels <- readxl::read_excel("~/Data_Analysis/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/MNH Dan Labels.xlsx") %>%
  dplyr::rename(mask_value=roinum) %>% select(mask_value, plot_label, Stream)

mean_head <- function(x) { 
  if (is.numeric(x)) {
    mean(x, na.rm=T)
  } else {
    head(x, n=1)
  }
}


plot_moderation <- function(obj) {
  obj_name <- as.character(match.call())[2]
  
  for (ee in seq_along(obj$emtrends_list)) {
    ename <- names(obj$emtrends_list)[ee]
    e_df <- obj$emtrends_list[[ee]] %>%
      mutate(fmri_fac = factor(fmri_beta, levels=c(-2, 0, 2), labels=c("-2SD", "Mean", "+2SD")))
    slope_name <- grep(".trend$", names(e_df), value = TRUE)[1L]
    e_df <- e_df %>% dplyr::rename(slope = !!slope_name)
      
    if ("rewFunc" %in% names(e_df)) {
      ofile <- paste(obj_name, ename, "l1_rewFunc.pdf", sep="_")
      e_df <- e_df %>% inner_join(labels, by="mask_value") %>% dplyr::select(-rhs, -l1_cope_name, -outcome)
      if (slope_name == "rt_vmax_lag.trend") {
        gvars <- c("rewFunc", "Stream", "fmri_fac")
        e_df <- e_df %>%
          group_by(across(gvars)) %>%
          summarise_all(mean_head)

        g <- ggplot(e_df, aes(x = 1, y = slope, ymax = slope + std.error, ymin = slope - std.error, color=fmri_fac)) +
          geom_pointrange(position=position_dodge(width=1)) + facet_grid(rewFunc~Stream) + ggtitle(paste0(obj_name, " ", slope_name, "by condition effects"))  
      } else {
        gvars <- c("rewFunc", "Stream", "last_outcome", "fmri_fac")
        e_df <- e_df %>%
          group_by(across(gvars)) %>%
          summarise_all(mean_head)
        
        g <- ggplot(e_df, aes(x = last_outcome, y = slope, ymax = slope + std.error, ymin = slope - std.error, color=fmri_fac)) +
          geom_pointrange(position=position_dodge(width=1)) + facet_grid(rewFunc~Stream) + ggtitle(paste0(obj_name, " ", slope_name, "by condition effects"))  
      }
      # e_df <- e_df %>% inner_join(labels, by="mask_value") %>% dplyr::select(-rhs, -l1_cope_name, -outcome) %>%
      #   group_by(across(gvars)) %>%
      #   summarise_all(mean_head)
      # 
      
      ggsave(g, filename = ofile, width=12, height=10)
    } else {
      ofile <- paste(obj_name, ename, "l1.pdf", sep="_")
      e_df <- e_df %>% inner_join(labels, by="mask_value") %>% dplyr::select(-rhs, -l1_cope_name, -l2_cope_name, -outcome) 
      if (slope_name == "rt_vmax_lag.trend") {
        gvars <- c("Stream", "fmri_fac")
        e_df <- e_df %>%
          group_by(across(gvars)) %>%
          summarise_all(mean_head)
        
        g <- ggplot(e_df, aes(x = 1, y = slope, ymax = slope + std.error, ymin = slope - std.error, color=fmri_fac)) +
          geom_pointrange(position=position_dodge(width=1)) + facet_grid(.~Stream) + ggtitle(paste0(obj_name, " ", slope_name, "by condition effects"))  
        
      } else {
        gvars <- c("Stream", "last_outcome", "fmri_fac")
        e_df <- e_df %>%
          group_by(across(gvars)) %>%
          summarise_all(mean_head)
        
        g <- ggplot(e_df, aes(x = last_outcome, y = slope, ymax = slope + std.error, ymin = slope - std.error, color=fmri_fac)) +
          geom_pointrange(position=position_dodge(width=1)) + facet_grid(.~Stream) + ggtitle(paste0(obj_name, " ", slope_name, "by condition effects"))  
        
      }
      ggsave(g, filename = ofile, width=12, height=10)
    }
    
  }
}


echange_l1 <- readRDS("/Volumes/GoogleDrive/My Drive/SCEPTIC_fMRI/wholebrain_betas/L1m-echange/parcel_maps_l1/Schaefer2018_200Parcels_7Networks_order_fonov_2.3mm_ants_cope_l1_mixed_by.rds")
echange_l2 <- readRDS("/Volumes/GoogleDrive/My Drive/SCEPTIC_fMRI/wholebrain_betas/L1m-echange/parcel_maps_l2/Schaefer2018_200Parcels_7Networks_order_fonov_2.3mm_ants_cope_l2_mixed_by.rds")
pe_l1 <- readRDS("/Volumes/GoogleDrive/My Drive/SCEPTIC_fMRI/wholebrain_betas/L1m-pe/parcel_maps_l1/Schaefer2018_200Parcels_7Networks_order_fonov_2.3mm_ants_cope_l1_mixed_by.rds")
pe_l2 <- readRDS("/Volumes/GoogleDrive/My Drive/SCEPTIC_fMRI/wholebrain_betas/L1m-pe/parcel_maps_l2/Schaefer2018_200Parcels_7Networks_order_fonov_2.3mm_ants_cope_l2_mixed_by.rds")
entropy_l1 <- readRDS("/Volumes/GoogleDrive/My Drive/SCEPTIC_fMRI/wholebrain_betas/L1m-entropy_wiz/parcel_maps_l1/Schaefer2018_200Parcels_7Networks_order_fonov_2.3mm_ants_cope_l1_mixed_by.rds")
entropy_l2 <- readRDS("/Volumes/GoogleDrive/My Drive/SCEPTIC_fMRI/wholebrain_betas/L1m-entropy_wiz/parcel_maps_l2/Schaefer2018_200Parcels_7Networks_order_fonov_2.3mm_ants_cope_l2_mixed_by.rds")

plot_moderation(echange_l1)
plot_moderation(echange_l2)
plot_moderation(pe_l1)
plot_moderation(pe_l2)
plot_moderation(entropy_l1)
plot_moderation(entropy_l2)


# 
# rt_lag <- rt_lag %>% inner_join(labels, by="mask_value") %>% dplyr::select(-rhs, -l1_cope_name, -outcome) %>%
#   group_by(rewFunc, Stream, last_outcome, fmri_fac) %>%
#   summarise_all(mean_head)
# 
# pdf("echange_rtlag_fmri_by_rewfunc.pdf", width=12, height=10)
# g <- ggplot(rt_lag, aes(x = last_outcome, y = rt_lag.trend, ymax = rt_lag.trend + std.error, ymin = rt_lag.trend - std.error, color=fmri_fac)) +
#   geom_pointrange(position=position_dodge(width=1)) + facet_grid(rewFunc~Stream) + ggtitle("Entropy change by condition effects")
# plot(g)
# dev.off()
