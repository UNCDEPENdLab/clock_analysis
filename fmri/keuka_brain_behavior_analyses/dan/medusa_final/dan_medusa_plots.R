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

# manual plot of ABSPE final
get_df <- function(ddf, model_name, term) {
  require(dplyr)
  ddf$coef_df_reml %>% filter(model_name==!!model_name & effect=="fixed" & evt_time <= 5) %>%
    filter(term == !!term) %>%
    mutate(..p=p.adjust(p.value, method="fdr")) %>%
    mutate(
      p_level = case_when(
        ..p > .05 ~ '1',
        ..p < .05 & ..p > .01 ~ '2',
        ..p < .01 & ..p > .001 ~ '3',
        ..p <.001 ~ '4'),
      p_level = ordered(p_level, levels = c('1', '2', '3', '4'), labels = c("NS","p < .05", "p < .01", "p < .001"))
    ) %>%
    mutate(
      vm_gradient17 = ordered(vm_gradient17, levels=c("MT+", "PPCcaudal", "PPCrostral", "Premotor"), labels=c("MT+", "Caudal PPC", "Rostral PPC", "Premotor"))
    ) %>%
    filter(evt_time < 5) # only 0-4 post-feedback
}

ddf <- readRDS("/Volumes/GoogleDrive/My Drive/SCEPTIC_fMRI/dan_medusa/rt_400_final47_encode_medusa_fmri_rslopeonly.rds")


# final color assignments
colors <- RColorBrewer::brewer.pal(4, "Dark2") %>% setNames(c("MT+", "Premotor", "Rostral PPC", "Caudal PPC"))
colors <- colors[c("MT+", "Caudal PPC", "Rostral PPC", "Premotor")] # order along visuomotor transformation

# common elements
common <- list(
  geom_hline(yintercept = 0, size=1.5, alpha=0.6),
  geom_vline(xintercept = 0, size=1.5, alpha=0.6),
  geom_line(mapping = aes(size=NULL), size=1.0, alpha=0.6, position=position_dodge(width=0.8)),
  geom_pointrange(position=position_dodge(width=0.8)),
  scale_size_manual("FDR-corrected p", values=c(0.3, 0.7, 1.1, 1.6)),
  theme_bw(base_size=15),
  scale_color_manual("Visuomotor\ngradient", values=colors)
)

my_df <- get_df(ddf, "enc_rt_abspe_logkld_rslope", "abs_pe")


pdf("/Users/hallquist/Library/CloudStorage/OneDrive-Personal/collected_letters/papers/meg/figures/fig_3_PE_whole_brain_medusa_bb/abs_pe_medusa.pdf", width=8, height=5)
ggplot(my_df, aes(x=evt_time, y=estimate, ymin=estimate-std.error, ymax=estimate + std.error, color=vm_gradient17, group=vm_gradient17, size=p_level)) +
  common +
  ylab("abs(PE) coefficient (AU)") + xlab("Time relative to feedback (seconds)")
dev.off()



# rew > om
my_df <- get_df(ddf, "enc_rt_abspe_logkld_rslope", "outcomeReward")

pdf("/Users/hallquist/Library/CloudStorage/OneDrive-Personal/collected_letters/papers/meg/figures/fig_3_PE_whole_brain_medusa_bb/rewom_medusa.pdf", width=8, height=5)
ggplot(my_df, aes(x=evt_time, y=estimate, ymin=estimate-std.error, ymax=estimate + std.error, color=vm_gradient17, group=vm_gradient17, size=p_level)) +
  common +
  ylab("Reward > omission coefficient (AU)") + xlab("Time relative to feedback (seconds)")
dev.off()

#rew > om simple plot (no pe in model)
my_df <- get_df(ddf, "enc_rt_rewom_rslope", "outcomeReward")

pdf("/Users/hallquist/Library/CloudStorage/OneDrive-Personal/collected_letters/papers/meg/figures/fig_3_PE_whole_brain_medusa_bb/rewom_nope_medusa.pdf", width=8, height=5)
ggplot(my_df, aes(x=evt_time, y=estimate, ymin=estimate-std.error, ymax=estimate + std.error, color=vm_gradient17, group=vm_gradient17, size=p_level)) +
  common +
  ylab("Reward > omission coefficient (AU)") + xlab("Time relative to feedback (seconds)")
dev.off()


# signed pe
my_df <- get_df(ddf, "enc_rt_pe_logkld_rslope", "pe_max")

pdf("/Users/hallquist/Library/CloudStorage/OneDrive-Personal/collected_letters/papers/meg/figures/fig_3_PE_whole_brain_medusa_bb/pe_medusa.pdf", width=8, height=5)
ggplot(my_df, aes(x=evt_time, y=estimate, ymin=estimate-std.error, ymax=estimate + std.error, color=vm_gradient17, group=vm_gradient17, size=p_level)) +
  common +
  ylab("Prediction error coefficient (AU)") + xlab("Time relative to feedback (seconds)")
dev.off()


# vmax in abspe model
my_df <- get_df(ddf, "enc_rt_abspe_logkld_rslope", "v_max_wi")


pdf("/Users/hallquist/Library/CloudStorage/OneDrive-Personal/collected_letters/papers/meg/figures/fig_5_echange_whole_brain_medusa_bb_coxme/vmax_medusa.pdf", width=8, height=5)
ggplot(my_df, aes(x=evt_time, y=estimate, ymin=estimate-std.error, ymax=estimate + std.error, color=vm_gradient17, group=vm_gradient17, size=p_level)) +
  common +
  ylab("Vmax coefficient (AU)") + xlab("Time relative to feedback (seconds)")
dev.off()


my_df <- get_df(ddf, "enc_rt_abspe_logkld_rslope", "v_entropy_wi_change")


pdf("/Users/hallquist/Library/CloudStorage/OneDrive-Personal/collected_letters/papers/meg/figures/fig_5_echange_whole_brain_medusa_bb_coxme/v_entropy_wi_change_medusa.pdf", width=8, height=5)
ggplot(my_df, aes(x=evt_time, y=estimate, ymin=estimate-std.error, ymax=estimate + std.error, color=vm_gradient17, group=vm_gradient17, size=p_level)) +
  common +
  ylab("Entropy change coefficient (AU)") + xlab("Time relative to feedback (seconds)")
dev.off()










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
