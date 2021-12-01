library(dplyr)
library(ggplot2)
library(data.table)

plot_medusa <- function(coef_obj, x="evt_time", y="estimate", ymin=NULL, ymax=NULL, color=NULL, facet_by=NULL, p.value="p.value", 
                        pdf_by=c("term", "model_name"), out_dir=getwd(), plot_type="line",
                        width=9, height=7) {
  if (checkmate::test_data_frame(coef_obj$coef_df_reml)) {
    to_plot <- coef_obj$coef_df_reml
  } else if (checkmate::test_data_frame(coef_obj$coef_df_ml)) {
    to_plot <- coef_obj$coef_df_ml
  }
  
  
  if (!checkmate::test_directory_exists(out_dir)) {
    dir.create(out_dir, recursive=TRUE)
  }
  stopifnot(pdf_by %in% names(to_plot))
  
  to_plot <- to_plot %>%
    dplyr::rename(..p = !!p.value) %>%
    mutate(
      p_level = case_when(
        ..p > .05 ~ '1',
        ..p < .05 & ..p > .01 ~ '2',
        ..p < .01 & ..p > .001 ~ '3',
        ..p <.001 ~ '4'),
      p_level = ordered(p_level, levels = c('1', '2', '3', '4'), labels = c("NS","p < .05", "p < .01", "p < .001"))
    )
  
  to_plot_list <- split(to_plot, by=pdf_by)
  for (ff in seq_along(to_plot_list)) {
    this_df <- to_plot_list[[ff]]
    this_name <- make.names(names(to_plot_list[ff]))
    g <- ggplot(this_df, aes_string(x=x, y=y, color=color, ymin=ymin, ymax=ymax)) + 
      geom_line(size=1, position=position_dodge(width=0.4)) + 
      geom_pointrange(aes(size=p_level), position=position_dodge(width=0.4)) +
      scale_color_brewer(palette="Dark2", labels=c("1" = "MT+, control", "2" = "Parieto-occipital", "3" = "Post. parietal", "4" = "Frontal")) +
      geom_hline(yintercept = 0, size=1.5, alpha=0.6) +
      geom_vline(xintercept = 0, size=1.5, alpha=0.6) +
      scale_size_manual(values=c(0.5, 0.8, 1.1, 1.4)) + theme_bw(base_size=15)
    if (!is.null(facet_by)) {
      g <- g + facet_wrap({{facet_by}})
    }
    
    pdf(file.path(out_dir, paste0(this_name, ".pdf")), width=width, height=height)
    plot(g)
    dev.off()
  }
  
}

meg = F
if (meg) {
ddf <- readRDS("/Volumes/GoogleDrive/My Drive/SCEPTIC_fMRI/dan_medusa/rt_encode_medusa_fmri_meg_simple.rds")
} else {ddf <- readRDS("/Volumes/GoogleDrive/My Drive/SCEPTIC_fMRI/dan_medusa/rt_encode_medusa_fmri.rds")}
out_dir <- "/Volumes/GoogleDrive/My Drive/SCEPTIC_fMRI/dan_medusa/"
ddf$coef_df_reml <- ddf$coef_df_reml %>% dplyr::filter(evt_time <= 5 & effect=="fixed") %>% group_by(term) %>%
  mutate(p_FDR=p.adjust(p.value, method="fdr")) %>%
  ungroup() %>% setDT()

plot_medusa(ddf, x="evt_time", y="estimate", ymin="estimate - std.error", ymax="estimate + std.error", color="vm_gradient", facet_by="side", 
            out_dir=file.path(out_dir, "rt_encode"), p.value="p_FDR")

ddf <- readRDS("/Volumes/GoogleDrive/My Drive/SCEPTIC_fMRI/dan_medusa/clock_encode_medusa_fmri_scaled.rds")
ddf$coef_df_reml <- ddf$coef_df_reml %>% dplyr::filter(evt_time <= 5) %>% 
  filter(effect=="fixed") %>%
  group_by(term) %>%
  mutate(p_FDR=p.adjust(p.value, method="fdr")) %>%
  ungroup() %>% setDT()

plot_medusa(ddf, x="evt_time", y="estimate", ymin="estimate - std.error", ymax="estimate + std.error", color="vm_gradient", facet_by="side",
             out_dir=file.path(out_dir, "clock_encode_scaled"), p.value="p_FDR")




rdf <- readRDS("/Users/hallquist/Data_Analysis/clock_analysis/fmri/keuka_brain_behavior_analyses/rt_prediction_rt_aligned_mixed_by.rds")
out_dir <- "/Volumes/GoogleDrive/My Drive/SCEPTIC_fMRI/dan_medusa/"
rdf$coef_df_reml <- rdf$coef_df_reml %>% dplyr::filter(evt_time <= 5) %>% 
  filter(effect=="fixed") %>%
  group_by(term) %>%
  mutate(p_FDR=p.adjust(p.value, method="fdr")) %>%
  ungroup() %>% 
  tidyr::separate(visuomotor_side, into=c("vm_gradient", "side"), sep="_") %>%
  setDT()

plot_medusa(rdf, x="evt_time", y="estimate", ymin="estimate - std.error", ymax="estimate + std.error", color="vm_gradient", facet_by="side", 
            out_dir=file.path(out_dir, "rt_predict"), p.value="p_FDR")


### CLOCK-ALIGNED RT PREDICTION
cdf <- readRDS("/Users/hallquist/Data_Analysis/clock_analysis/fmri/keuka_brain_behavior_analyses/rt_prediction_clock_aligned_mixed_by.rds")
out_dir <- "/Volumes/GoogleDrive/My Drive/SCEPTIC_fMRI/dan_medusa/"
cdf$coef_df_reml <- cdf$coef_df_reml %>% dplyr::filter(evt_time <= 5) %>% 
  filter(effect=="fixed") %>%
  group_by(term) %>%
  mutate(p_FDR=p.adjust(p.value, method="fdr")) %>%
  ungroup() %>% 
  tidyr::separate(visuomotor_side, into=c("vm_gradient", "side"), sep="_") %>%
  setDT()

plot_medusa(cdf, x="evt_time", y="estimate", ymin="estimate - std.error", ymax="estimate + std.error", color="vm_gradient", facet_by="side", 
            out_dir=file.path(out_dir, "clock_predict"), p.value="p_FDR")

