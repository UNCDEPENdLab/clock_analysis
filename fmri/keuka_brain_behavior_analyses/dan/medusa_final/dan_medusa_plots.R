## parcelwise checks
ddf$coef_df_reml <- ddf$coef_df_reml %>% dplyr::filter(evt_time <= 5 & effect=="fixed") %>% 
  group_by(term, model_name) %>%
  mutate(p_FDR=p.adjust(p.value, method="fdr")) %>%
  ungroup() %>% setDT()

out_dir <- "/Volumes/GoogleDrive/My Drive/SCEPTIC_fMRI/dan_medusa"
plot_medusa(ddf, x="evt_time", y="estimate", ymin="estimate - std.error", ymax="estimate + std.error", 
            color="label", facet_by=NULL, panel_by="vm_gradient17",
            out_dir=file.path(out_dir, "rt_encode_400_final47_20Jul2022_parcelwise"), p.value="p_FDR",
            width = 15, height=15)




meg = F
if (meg) {
  ddf <- readRDS("/Volumes/GoogleDrive/My Drive/SCEPTIC_fMRI/dan_medusa/rt_encode_medusa_fmri_meg_simple_ec.rds")
} else {
  #ddf <- readRDS("/Volumes/GoogleDrive/My Drive/SCEPTIC_fMRI/dan_medusa/rt_encode_medusa_fmri_pe_posneg.rds")
  ddf <- readRDS("/Volumes/GoogleDrive/My Drive/SCEPTIC_fMRI/dan_medusa/rt_400_final47_encode_medusa_fmri.rds")
}

out_dir <- "/Volumes/GoogleDrive/My Drive/SCEPTIC_fMRI/dan_medusa/"

ddf$coef_df_reml <- ddf$coef_df_reml %>% dplyr::filter(evt_time <= 5 & effect=="fixed") %>% 
  group_by(term, model_name) %>%
  mutate(p_FDR=p.adjust(p.value, method="fdr")) %>%
  ungroup() %>% setDT()

plot_medusa(ddf, x="evt_time", y="estimate", ymin="estimate - std.error", ymax="estimate + std.error", color="vm_gradient17", facet_by="side",
            out_dir=file.path(out_dir, "rt_encode_400_final47_6Jul2022"), p.value="p_FDR")

#ddf <- readRDS("/Volumes/GoogleDrive/My Drive/SCEPTIC_fMRI/dan_medusa/clock_encode_medusa_fmri_scaled.rds")
#ddf <- readRDS("/Volumes/GoogleDrive/My Drive/SCEPTIC_fMRI/dan_medusa/clock_encode_medusa_fmri.rds")
ddf <- readRDS("/Volumes/GoogleDrive/My Drive/SCEPTIC_fMRI/dan_medusa/clock_400_final47_encode_medusa_fmri.rds")
ddf$coef_df_reml <- ddf$coef_df_reml %>% dplyr::filter(evt_time <= 5) %>% 
  filter(effect=="fixed") %>%
  group_by(term, model_name) %>%
  mutate(p_FDR=p.adjust(p.value, method="fdr")) %>%
  ungroup() %>% setDT()

plot_medusa(ddf, x="evt_time", y="estimate", ymin="estimate - std.error", ymax="estimate + std.error", color="vm_gradient17", facet_by="side",
            out_dir=file.path(out_dir, "clock_encode_400_final47_6Jul2022"), p.value="p_FDR")

# clock online
ddf <- readRDS("/Volumes/GoogleDrive/My Drive/SCEPTIC_fMRI/dan_medusa/clock_online_encode_medusa_fmri.rds")
ddf$coef_df_reml <- ddf$coef_df_reml %>% dplyr::filter(evt_time < 4) %>% 
  filter(effect=="fixed") %>%
  group_by(term) %>%
  mutate(p_FDR=p.adjust(p.value, method="fdr")) %>%
  ungroup() %>% setDT()

plot_medusa(ddf, x="evt_time", y="estimate", ymin="estimate - std.error", ymax="estimate + std.error", color="vm_gradient17", facet_by="side",
            out_dir=file.path(out_dir, "clock_encode_online_24Nov2021"), p.value="p_FDR")

# 
# 
# 
# 
# 
# 
# rdf <- readRDS("/Users/hallquist/Data_Analysis/clock_analysis/fmri/keuka_brain_behavior_analyses/rt_prediction_rt_aligned_mixed_by.rds")
# out_dir <- "/Volumes/GoogleDrive/My Drive/SCEPTIC_fMRI/dan_medusa/"
# rdf$coef_df_reml <- rdf$coef_df_reml %>% dplyr::filter(evt_time <= 5) %>% 
#   filter(effect=="fixed") %>%
#   group_by(term) %>%
#   mutate(p_FDR=p.adjust(p.value, method="fdr")) %>%
#   ungroup() %>% 
#   tidyr::separate(visuomotor_side, into=c("vm_gradient17", "side"), sep="_") %>%
#   setDT()
# 
# plot_medusa(rdf, x="evt_time", y="estimate", ymin="estimate - std.error", ymax="estimate + std.error", color="vm_gradient17", facet_by="side", 
#             out_dir=file.path(out_dir, "rt_predict"), p.value="p_FDR")
# 
# 
# ### CLOCK-ALIGNED RT PREDICTION
# cdf <- readRDS("/Users/hallquist/Data_Analysis/clock_analysis/fmri/keuka_brain_behavior_analyses/rt_prediction_clock_aligned_mixed_by.rds")
# out_dir <- "/Volumes/GoogleDrive/My Drive/SCEPTIC_fMRI/dan_medusa/"
# cdf$coef_df_reml <- cdf$coef_df_reml %>% dplyr::filter(evt_time <= 5) %>% 
#   filter(effect=="fixed") %>%
#   group_by(term) %>%
#   mutate(p_FDR=p.adjust(p.value, method="fdr")) %>%
#   ungroup() %>% 
#   tidyr::separate(visuomotor_side, into=c("vm_gradient17", "side"), sep="_") %>%
#   setDT()
# 
# plot_medusa(cdf, x="evt_time", y="estimate", ymin="estimate - std.error", ymax="estimate + std.error", color="vm_gradient17", facet_by="side", 
#             out_dir=file.path(out_dir, "clock_predict"), p.value="p_FDR")
# 
