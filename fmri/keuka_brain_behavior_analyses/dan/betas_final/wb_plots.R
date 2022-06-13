# plots for whole-brain betas
library(tidyverse)
library(ggbrain)
library(cowplot)
library(patchwork)
library(fmri.pipeline)

Sys.setenv(AFNIDIR="/Users/hallquist/abin")
setwd("~/Data_Analysis/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/betas_final")

bg_color <- "gray90"
text_color <- "black"
symmetric_legend <- FALSE
base_size <- 21 # theme
betas_dir <- "~/GDrive/SCEPTIC_fMRI/wholebrain_betas"

### OVERALL MAPS: ENTROPY, ECHANGE, PE

# entropy
# axial z = 58 IPS, FEF, Motor
# sagittal x = -7, IPS, neg vmPFC

entropy_dir <- file.path(betas_dir, "L1m-entropy_wiz")
#pdf("figures/entropy_wb_ptfce.pdf", width = 8, height = 4)
png("figures/entropy_wb_ptfce.png", width = 3.5, height = 4, units = "in", res = 600, bg = bg_color)
# ggbrain(underlay="template_brain.nii.gz", overlay = "zstat6_ptfce_fwep_0.05_1mm.nii.gz", axial_slices = list(xyz=c(58)),
#         remove_null_space = TRUE, pos_thresh = 5, neg_thresh = -5, background_color = "black", text_color = "white") # already thresholded by ptfce

entropy_gg <- ggbrain(
  underlay = file.path(entropy_dir, "template_brain.nii.gz"), 
  overlay = file.path(entropy_dir, "zstat6_ptfce_fwep_0.05_1mm.nii.gz"),
  slices = data.frame(coord = c("z = 58", "x = -7")),
  remove_null_space = TRUE, pos_thresh = 5.2, neg_thresh = -5.2,
  background_color = bg_color, text_color = text_color,
  panel_borders = FALSE, symmetric_legend = symmetric_legend, base_size = base_size,
  title = "Entropy\n(decision phase)", ncol=1
)

plot(entropy_gg)
dev.off()

######

echange_dir <- file.path(betas_dir, "L1m-echange")

# pdf("figures/entropy_wb_ptfce.pdf", width = 8, height = 4)
#png("figures/echange_wb_ptfce.png", width = 8, height = 3.5, units = "in", res = 600, bg = bg_color)
png("figures/echange_wb_ptfce.png", width = 3.5, height = 4, units = "in", res = 600, bg = bg_color)

echange_gg <- ggbrain(
  underlay = file.path(echange_dir, "template_brain.nii.gz"),
  overlay = file.path(echange_dir, "zstat6_ptfce_fwep_0.05_1mm.nii.gz"),
  slices = data.frame(coord = c("z = 54", "x = 7"), coord_labels = TRUE),
  remove_null_space = TRUE, pos_thresh = 5.2, neg_thresh = -5.2,
  background_color = bg_color, text_color = text_color,
  panel_borders = FALSE, symmetric_legend = symmetric_legend, base_size = base_size,
  underlay_contrast = "low", title = "Entropy change\n(feedback phase)", ncol=1
)

plot(echange_gg)

dev.off()

##### PE

pe_dir <- file.path(betas_dir, "L1m-pe")

# pdf("figures/entropy_wb_ptfce.pdf", width = 8, height = 4)
png("figures/pe_wb_ptfce.png", width = 5, height = 4, units = "in", res = 600, bg = bg_color)

pe_gg <- ggbrain(
  underlay = file.path(pe_dir, "template_brain.nii.gz"),
  overlay = file.path(pe_dir, "pe_overall_ptfce_fwep_0.05_1mm.nii.gz"),
  slices = data.frame(coord = c("x = -42", "y = 9", "z = 49"), coord_labels = TRUE),
  remove_null_space = TRUE, pos_thresh = 5.2, neg_thresh = -5.2,
  background_color = bg_color, text_color = text_color,
  panel_borders = FALSE, symmetric_legend = symmetric_legend, base_size = base_size,
  underlay_contrast = "low", ncol=2, nrow=2, title = "Reward prediction error\n(feedback phase)"
)

# slices = data.frame(coord = c("z = 58", "x = -7.45", "y = 10")),
# positive_colorscale = scale_fill_viridis_c(), negative_colorscale = scale_fill_distiller(palette="Blues"))
plot(pe_gg)
dev.off()

##### REW > OM


rew_dir <- file.path(betas_dir, "L1m-rew_om")

png("figures/rew_wb_ptfce.png", width = 5, height = 4, units = "in", res = 600, bg = bg_color)

rew_gg <- ggbrain(
  underlay = file.path(rew_dir, "template_brain.nii.gz"),
  overlay = file.path(rew_dir, "zstat6_ptfce_1mm.nii.gz"),
  slices = data.frame(coord = c("x = -42", "y = 9", "z = 49", "z = 55"), coord_labels = TRUE),
  remove_null_space = TRUE, pos_thresh = 5.2, neg_thresh = -5.2,
  background_color = bg_color, text_color = text_color,
  panel_borders = FALSE, symmetric_legend = symmetric_legend, base_size = base_size,
  underlay_contrast = "low", ncol=2, nrow=2, title = "Reward > Omission\n(feedback phase)"
)

plot(rew_gg)
dev.off()

#### PAIRWISE DIFFERENCES: Entropy vs. PE is not of interest
pairwise_dir <- file.path(betas_dir, "pairwise_diffs")

#### ECHANGE > PE (for figure)

# these are at 1mm, so set a higher cluster_nvox
epe_clust <- afni_3dclusterize$new(inset = file.path(pairwise_dir, "echange_m_pe_ptfce.nii.gz"), bisided = TRUE,
                                   NN = 1, clust_nvox = 400, lower_thresh = -3, upper_thresh = 3, ithr = 0)
epe_clust$run(force=TRUE)
epe_clust$get_clust_df()

png("figures/echange_gt_pe_ptfce.png", width = 6, height = 4, units = "in", res = 600, bg = bg_color)

echange_gt_pe_gg <- ggbrain(
  underlay = file.path(pairwise_dir, "template_brain.nii.gz"),
  overlay = file.path(pairwise_dir, "echange_m_pe_ptfce_cluster_data.nii.gz"),
  #slices = data.frame(coord = c("z = 52", "z = -10"), coord_labels = TRUE),
  slices = data.frame(coord = c("z = 54", "z = -10"), coord_labels = TRUE),
  remove_null_space = TRUE, pos_thresh = 3, neg_thresh = -3,
  background_color = bg_color, text_color = text_color,
  panel_borders = FALSE, symmetric_legend = symmetric_legend, base_size = base_size,
  underlay_contrast = "low", nrow=1, title = "Entropy change > RPE"
)

plot(echange_gt_pe_gg)
dev.off()

#### ECHANGE > ENTROPY

# these are at 1mm, so set a higher cluster_nvox
epe_clust <- afni_3dclusterize$new(inset = file.path(pairwise_dir, "echange_m_entropy_ptfce.nii.gz"), bisided = TRUE,
                                   NN = 1, clust_nvox = 400, lower_thresh = -3, upper_thresh = 3, ithr = 0)
epe_clust$run(force=TRUE)
epe_clust$get_clust_df()

png("figures/echange_gt_entropy_ptfce.png", width = 9, height = 4, units = "in", res = 600, bg = bg_color)

echange_gt_entropy_gg <- ggbrain(
  underlay = file.path(pairwise_dir, "template_brain.nii.gz"),
  overlay = file.path(pairwise_dir, "echange_m_entropy_ptfce_cluster_data.nii.gz"),
  #slices = data.frame(coord = c("z = 54", "x=-2.5", "z = -10"), coord_labels = TRUE),
  slices = data.frame(coord = c("z = 54", "z = -10"), coord_labels = TRUE),
  remove_null_space = TRUE, pos_thresh = 3, neg_thresh = -3,
  background_color = bg_color, text_color = text_color,
  panel_borders = FALSE, symmetric_legend = symmetric_legend, base_size = base_size,
  underlay_contrast = "low", nrow=1, title = "Entropy change > Entropy"
)

plot(echange_gt_entropy_gg)
dev.off()

#### PE > REW

# these are at 1mm, so set a higher cluster_nvox

clust_obj <- afni_3dclusterize$new(inset = file.path(pairwise_dir, "pe_m_rew_ptfce.nii.gz"), bisided = TRUE,
                                   NN = 1, clust_nvox = 500, lower_thresh = -3, upper_thresh = 3, ithr = 0)
clust_obj$run(force=TRUE)
clust_obj$get_clust_df()

# png("figures/pe_m_rew_ptfce.png", width = 4, height = 6, units = "in", res = 600, bg = "black")
# pdf("figures/pe_m_rew_ptfce.pdf", width = 4, height = 6)
png("figures/pe_m_rew_ptfce.png", width = 3, height = 2, units = "in", res = 600, bg = bg_color)

pe_gt_rew_gg <- ggbrain(
  underlay = file.path(pairwise_dir, "template_brain.nii.gz"),
  overlay = file.path(pairwise_dir, "pe_m_rew_ptfce_cluster_data.nii.gz"),
  #slices = data.frame(coord = c("z = 52", "z = -10"), coord_labels = TRUE),
  #slices = data.frame(coord = c("x = -42", "y = 9", "z = 49"), coord_labels = TRUE), # exact match to PE
  slices = data.frame(coord = c("z = 49"), coord_labels = TRUE), # single slice
  remove_null_space = TRUE, pos_thresh = 3, neg_thresh = -5,
  background_color = bg_color, text_color = text_color,
  panel_borders = FALSE, symmetric_legend = symmetric_legend, base_size = 12,
  underlay_contrast = "low", ncol=1, title="RPE > Reward"
)

plot(pe_gt_rew_gg)
dev.off()


### ECHANGE > REW

# these are at 1mm, so set a higher cluster_nvox
clust_obj <- afni_3dclusterize$new(inset = file.path(pairwise_dir, "echange_m_rew_ptfce.nii.gz"), bisided = TRUE,
                                   NN = 1, clust_nvox = 400, lower_thresh = -3, upper_thresh = 3, ithr = 0)
clust_obj$run(force=TRUE)
clust_obj$get_clust_df()

png("figures/echange_gt_rew_ptfce.png", width = 4, height = 6, units = "in", res = 600, bg = bg_color)

echange_gt_rew_gg <- ggbrain(
  underlay = file.path(pairwise_dir, "template_brain.nii.gz"),
  overlay = file.path(pairwise_dir, "echange_m_rew_ptfce_cluster_data.nii.gz"),
  # slices = data.frame(coord = c("z = 52", "z = -10"), coord_labels = TRUE),
  slices = data.frame(coord = c("z = 54", "z = -10"), coord_labels = TRUE),
  remove_null_space = TRUE, pos_thresh = 3, neg_thresh = -3,
  background_color = bg_color, text_color = text_color,
  panel_borders = FALSE, symmetric_legend = symmetric_legend, base_size = base_size,
  underlay_contrast = "low", ncol=1, title = "Entropy change - reward"
)

plot(echange_gt_rew_gg)
dev.off()

pe_inset <- pe_gt_rew_gg & guides(fill = guide_colorbar(barheight=4.5, ticks.colour = text_color))

#   fill_new = guide_colorbar(barheight = panel_height, available_aes = "fill_new"),
#   #fill_new_new = guide_colorbar(barheight = unit(0.1, "cm"), available_aes = "fill_new_new")
# )

# pe_gg & plot_annotation(
#   theme = theme(legend.justification = "top")
# ) + inset_element(pe_gt_rew_gg, left = 0.5, bottom = 0, right = 1, top = 0.5,
#                   on_top = TRUE, align_to = 'full')

# this does work, but may take too much time to be useful in code form
pdf("figures/rpe_with_inset.pdf", width = 7, height=6, bg = bg_color)
# pe_gg + inset_element(pe_gt_rew_gg, left = 0.8, bottom = 0, right = 1, top = 0.5,
#                       on_top = TRUE, align_to = 'plot')
pe_gg_with_inset <- ggdraw(
  pe_gg + plot_annotation(theme = theme(legend.justification = "top"))
) +
  draw_plot(pe_inset + plot_annotation(theme = theme(plot.background = element_rect(color="black", size = 1.1))),
            x = 0.45, y=0.01, width = 0.43, height=0.43)
plot(pe_gg_with_inset)
dev.off()

pdf("figures/echange_with_inset.pdf", width = 9, height=3, bg = bg_color)
# pe_gg + inset_element(pe_gt_rew_gg, left = 0.8, bottom = 0, right = 1, top = 0.5,
#                       on_top = TRUE, align_to = 'plot')

# shrink color bar height
p2 <- echange_gt_pe_gg & 
  guides(
    fill = guide_colorbar(barheight=3, order = 1), 
    fill_new=guide_colorbar(barheight=3, available_aes = "fill_new", order = 2)
  ) & theme_void(base_size = 11)

# need to use wrap_elements in patchwork framework to keep the overall title as plot_annotation in each panel
combo2 <- wrap_elements(plot=echange_gg, clip=F) + wrap_elements(plot=p2, clip=F) +
  plot_layout(widths=c(1, 0.6)) +
  plot_annotation(
    theme = theme(plot.background = element_rect(fill=bg_color, color = NA)),
    tag_levels = 'A'
  ) & theme(plot.background = element_rect(fill = bg_color, color = 'gray60', size=1.1))
  
      

  
#   ggdraw(
#   echange_gg + plot_annotation(theme = theme(legend.justification = "top"))
# ) +
#   draw_plot(echange_gt_pe_gg + plot_annotation(theme = theme(plot.background = element_rect(color="black", size = 1.1))),
#             x = 0.45, y=0.01, width = 0.43, height=0.43)
plot(combo2)
dev.off()

#theme(plot.background = element_rect(color="black", size = 1.1))


# redefine position of inset for final composite
pe_gg_with_inset <- ggdraw(
  pe_gg + plot_annotation(theme = theme(legend.justification = "top"))
) +
  draw_plot(pe_inset + plot_annotation(theme = theme(plot.background = element_rect(color="black", size = 1.1))),
            x = 0.45, y=0.01, width = 0.43, height=0.43)


# shrink color bars for pairwise plots to fit within panel
echange_gt_entropy_gg_small <- echange_gt_entropy_gg & 
  guides(
    fill = guide_colorbar(barheight=4.5, order = 1, ticks.colour = text_color), 
    fill_new=guide_colorbar(barheight=4.5, available_aes = "fill_new", order = 2, ticks.colour = text_color)
  )

echange_gt_pe_gg_small <- echange_gt_pe_gg & 
  guides(
    fill = guide_colorbar(barheight=4.5, order = 1, ticks.colour = text_color), 
    fill_new=guide_colorbar(barheight=4.5, available_aes = "fill_new", order = 2, ticks.colour = text_color)
  )


pdf("figures/composite.pdf", width=15, height=9)
p <- ( wrap_elements(plot=entropy_gg) + wrap_elements(plot=echange_gg) + wrap_elements(plot=pe_gg_with_inset) + plot_layout(widths=c(1,1,1.6)) ) /
  ( wrap_elements(plot=echange_gt_entropy_gg_small) + wrap_elements(plot=echange_gt_pe_gg_small) ) +
  plot_layout(heights=c(2, 1.4)) + # first composite row is 2x the second row
  plot_annotation(
    tag_levels = "A", tag_suffix = ")",
    theme = theme(plot.background = element_rect(fill="white", color = NA))
  ) &
  theme(plot.tag = element_text(face = "bold", size=21, vjust = 0),
        plot.background = element_rect(fill="white", color = "gray50", size=2) #'gray60', size=3)
        #color = 'gray60', size=1.1)
  )




# layout <- "
# AABBCC
# DDDEEE
# "
# p <- wrap_elements(plot=entropy_gg) + 
#   wrap_elements(plot=echange_gg) +
#   wrap_elements(plot=pe_gg) +
#   wrap_elements(plot=echange_gt_entropy_gg) +
#   wrap_elements(plot=echange_gt_pe_gg) +
#   plot_layout(design=layout)
# 
# # p <- plot_grid(echange_gg,
# #                entropy_gg,
# #                pe_gg,
# #                labels="AUTO"
# # )
# plot(p)

# p1 <- plot_grid(entropy_gg, echange_gg, pe_gg, nrow=1)
# p2 <- plot_grid(echange_gt_entropy_gg, echange_gt_pe_gg, nrow=1)
# p <- plot_grid(p1 + p2, ncol=1)

# p <- plot_grid(echange_gg,
#                entropy_gg,
#                pe_gg,
#                labels="AUTO"
# )
plot(p)

dev.off()
