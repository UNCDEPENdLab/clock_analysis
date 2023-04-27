# compares neural fit in the DAN across the two models
# modeled on plot_medusa.R

library(tidyverse)
library(lme4)
library(data.table)
library(readxl)
library(fmri.pipeline) # mixed_by call: use fmri.pipeline installation

out_dir <- "~/OneDrive - University of Pittsburgh/Documents/SCEPTIC_fMRI/dan_medusa"
# decon_dir <- "~/OneDrive - University of Pittsburgh/Documents/SCEPTIC_fMRI/dan_medusa"
# schaefer_dir <- "~/code/schaefer_wb_parcellation"
# # schaefer_dir <- "~/Data_Analysis/schaefer_wb_parcellation"
repo_directory <- "~/code/clock_analysis"
# repo_directory <- "~/Data_Analysis/clock_analysis"
organize <- TRUE

require(patchwork)  
require(dplyr)
require(ggplot2)
require(data.table)

colors <- RColorBrewer::brewer.pal(4, "Dark2") %>% setNames(c("1" = "MT+","2" = "Premotor","3" = "PPCrostral","4" = "PPCcaudal"))
setwd(file.path(out_dir))
if (organize) {    # organize MEDUSA results
  files <-  gsub("//", "/", list.files(pattern = ".*wm.*rds", full.names = T))
message(paste0("Found ", length(files), " files."))
csl <- lapply(files, function(x) {
  print(x)
  df <- readRDS(x) %>% Filter(Negate(is.null),.)
  coef_df <- df$coef_df_reml
  return(coef_df)
})
ddf_all_coefs <- data.table::rbindlist(csl)  %>% distinct(vm_gradient17, evt_time, model_name, term, .keep_all = T)
saveRDS(ddf_all_coefs, file=file.path(out_dir, "wm_entropy", "rt_400_final47_encode_medusa_fmri_wm_all_coefs.rds"))

csl <- lapply(files, function(x) {
  print(x)
  df <- readRDS(x) %>% Filter(Negate(is.null),.)
  fit_df <- df$fit_df
  return(fit_df)
})
ddf_all_fit <- data.table::rbindlist(csl) %>%  distinct(vm_gradient17, evt_time, model_name, .keep_all = T)
saveRDS(ddf_all_fit, file=file.path(out_dir, "wm_entropy", "rt_400_final47_encode_medusa_fmri_wm_all_fit.rds"))
}

# df <- ddf$fit_df
df <- ddf_all_fit
# g <- 
# ggplot(df, aes(evt_time, AIC, lty = model_name, color = vm_gradient17)) +
#   geom_line(size=1) + 
#   # geom_pointrange(aes(size=p_level), position=position_dodge(width=0.4)) +
#   #scale_color_brewer(palette="Dark2", labels=c("1" = "MT+, control", "2" = "Caudal post. parietal", "3" = "Rostral post. parietal", "4" = "Frontal premotor")) +
#   scale_color_manual(values = colors) +
#   geom_hline(yintercept = 0, size=1.5, alpha=0.6) +
#   geom_vline(xintercept = 0, size=1.5, alpha=0.6) +
#   scale_size_manual(values=c(0.5, 0.8, 1.1, 1.4)) + theme_bw(base_size=15)

wdf <- df %>% select(model_name, AIC, evt_time, vm_gradient17) %>%  pivot_wider(names_from = model_name, values_from = AIC) %>% mutate(
  AIC_sceptic_minus_wm = enc_rt_base - enc_rt_wm,
  AIC_sceptic_minus_wm_exp = enc_rt_base - enc_rt_wm_exp,
  AIC_sceptic_minus_wm_exp1 = enc_rt_wm_exp_sceptic - enc_rt_wm_exp,
  AIC_selective_minus_full = enc_rt_base - enc_rt_base_full,
  AIC_sceptic_full_minus_wm_exp = enc_rt_wm_exp_sceptic_full - enc_rt_wm_exp,
  AIC_selective_minus_full_wm = enc_rt_wm_exp_sceptic - enc_rt_wm_exp_sceptic_full,
  AIC_sceptic_minus_wm_exp_rslope = enc_rt_wm_exp_sceptic_rslope - enc_rt_wm_exp_rslope
)

# quick model comparison
ggplot(wdf, aes(evt_time, AIC_sceptic_minus_wm_exp_rslope, color = vm_gradient17)) +
  geom_line(size=1) + 
  # geom_pointrange(aes(size=p_level), position=position_dodge(width=0.4)) +
  #scale_color_brewer(palette="Dark2", labels=c("1" = "MT+, control", "2" = "Caudal post. parietal", "3" = "Rostral post. parietal", "4" = "Frontal premotor")) +
  scale_color_manual(values = colors) +
  geom_hline(yintercept = 0, size=1.5, alpha=0.6) +
  geom_vline(xintercept = 0, size=1.5, alpha=0.6) +
  scale_size_manual(values=c(0.5, 0.8, 1.1, 1.4)) + theme_bw(base_size=15) +  theme(legend.title=element_blank())

setwd(file.path(out_dir, "wm_entropy"))
# pdf("AIC_sceptic_wm.pdf", height = 4, width = 6)
# ggplot(wdf, aes(evt_time, AIC_sceptic_minus_wm, color = vm_gradient17)) +
#   geom_line(size=1) + 
#   # geom_pointrange(aes(size=p_level), position=position_dodge(width=0.4)) +
#   #scale_color_brewer(palette="Dark2", labels=c("1" = "MT+, control", "2" = "Caudal post. parietal", "3" = "Rostral post. parietal", "4" = "Frontal premotor")) +
#   scale_color_manual(values = colors) +
#   geom_hline(yintercept = 0, size=1.5, alpha=0.6) +
#   geom_vline(xintercept = 0, size=1.5, alpha=0.6) +
#   scale_size_manual(values=c(0.5, 0.8, 1.1, 1.4)) + theme_bw(base_size=15) +  theme(legend.title=element_blank())
# dev.off()
# pdf("AIC_sceptic_wm_exp.pdf", height = 4, width = 6)
# ggplot(wdf, aes(evt_time, AIC_sceptic_minus_wm_exp, color = vm_gradient17)) +
#   geom_line(size=1) + 
#   # geom_pointrange(aes(size=p_level), position=position_dodge(width=0.4)) +
#   #scale_color_brewer(palette="Dark2", labels=c("1" = "MT+, control", "2" = "Caudal post. parietal", "3" = "Rostral post. parietal", "4" = "Frontal premotor")) +
#   scale_color_manual(values = colors) +
#   geom_hline(yintercept = 0, size=1.5, alpha=0.6) +
#   geom_vline(xintercept = 0, size=1.5, alpha=0.6) +
#   scale_size_manual(values=c(0.5, 0.8, 1.1, 1.4)) + theme_bw(base_size=15) +  theme(legend.title=element_blank())
# dev.off()

pdf("AIC_sceptic_with_wm_exp_vs_wm_exp.pdf", height = 2.5, width = 3)
ggplot(wdf, aes(evt_time, AIC_sceptic_minus_wm_exp1, color = vm_gradient17)) +
  geom_line(size=2) + 
  # geom_pointrange(aes(size=p_level), position=position_dodge(width=0.4)) +
  #scale_color_brewer(palette="Dark2", labels=c("1" = "MT+, control", "2" = "Caudal post. parietal", "3" = "Rostral post. parietal", "4" = "Frontal premotor")) +
  scale_color_manual(values = colors) + ylab("AIC difference\nWM + value <--> only WM") + xlab("Seconds after outcome") +
  geom_hline(yintercept = 0, size=1.5, alpha=0.6) +
  geom_vline(xintercept = 0, size=1.5, alpha=0.6) +
  scale_size_manual(values=c(0.5, 0.8, 1.1, 1.4)) + theme_bw(base_size=12) +  theme(legend.title=element_blank(), legend.position="none")
dev.off()


pdf("AIC_sceptic_selective_minus_full.pdf", height = 2.9, width = 3)
ggplot(wdf, aes(evt_time, AIC_selective_minus_full, color = vm_gradient17)) +
  geom_line(size=2) + 
  # geom_pointrange(aes(size=p_level), position=position_dodge(width=0.4)) +
  #scale_color_brewer(palette="Dark2", labels=c("1" = "MT+, control", "2" = "Caudal post. parietal", "3" = "Rostral post. parietal", "4" = "Frontal premotor")) +
  scale_color_manual(values = colors) + ylab("AIC difference\ncompression <-> no compression") + xlab("Seconds after outcome") +
  geom_hline(yintercept = 0, size=1.5, alpha=0.6) +
  geom_vline(xintercept = 0, size=1.5, alpha=0.6) +
  scale_size_manual(values=c(0.5, 0.8, 1.1, 1.4)) + theme_bw(base_size=11) +  theme(legend.title=element_blank(), legend.position="none")
dev.off()
# 
# pdf("AIC_sceptic_selective_minus_full_wm.pdf", height = 4, width = 6)
# ggplot(wdf, aes(evt_time, AIC_selective_minus_full_wm, color = vm_gradient17)) +
#   geom_line(size=1) + 
#   # geom_pointrange(aes(size=p_level), position=position_dodge(width=0.4)) +
#   #scale_color_brewer(palette="Dark2", labels=c("1" = "MT+, control", "2" = "Caudal post. parietal", "3" = "Rostral post. parietal", "4" = "Frontal premotor")) +
#   scale_color_manual(values = colors) +
#   geom_hline(yintercept = 0, size=1.5, alpha=0.6) +
#   geom_vline(xintercept = 0, size=1.5, alpha=0.6) +
#   scale_size_manual(values=c(0.5, 0.8, 1.1, 1.4)) + theme_bw(base_size=15) +  theme(legend.title=element_blank())
# dev.off()
# plot_medusa <- function(coef_obj, x="evt_time", y="estimate", ymin=NULL, ymax=NULL, color=NULL, facet_by=NULL, 
#                         p.value="p.value", lty_by=NULL, pdf_by=c("term", "model_name"), panel_by=NULL, out_dir=getwd(), 
#                         plot_type="line", flip=FALSE, term_filter=NULL, width=9, height=7, include_title=TRUE)
# result=recode(result, 'Win'='1', .default=NA_character_))
to_plot <- ddf_all_coefs %>% filter(effect == "fixed" & model_name == "enc_rt_wm_exp_sceptic" & 
                                      term %in% c("v_entropy_wi_change", "wm_choice_entropy_wi_change", "wm_reward_entropy_wi_change") &
                                                    evt_time > -3 & evt_time < 7) %>% 
# to_plot <- ddf$coef_df_reml %>% filter(effect == "fixed") %>% 
  mutate(vm_gradient17 = dplyr::recode(vm_gradient17, 'PPCcaudal' = "Caudal PPC", 'PPCrostral' = "Rostral PPC"))

p <- plot_medusa(to_plot, ymin = "conf.low", ymax = "conf.high", color = "vm_gradient17", out_dir = file.path(out_dir, "wm_entropy"),
            width = 6, height = 7,  pdf_by = "model_name", facet_by = "term", include_title = F)

term_names <- c(
  `v_entropy_wi_change` = "SCEPTIC value entropy change",
  `wm_choice_entropy_wi_change` = "WM choice entropy change",
  `wm_reward_entropy_wi_change` = "WM reward entropy change"
)

p <- p + facet_wrap(~term, labeller = as_labeller(term_names), ncol = 1) + xlab("Seconds after outcome")

pdf("SCEPTIC_vs_WM_entropy_change.pdf", height = 7, width = 5.5)
print(p)
dev.off()
