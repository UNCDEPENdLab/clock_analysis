# read in L2 betas for Explore clock, examine individual differences

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
library(patchwork)
library(corrplot)
library(psych)
# devtools::install_github("UNCDEPENdLab/dependlab")

source("~/code/Rhelpers/theme_black.R")
# install_github("UNCDEPENdLab/dependlab")
# library(dependlab)
# library(stringi)

clock_folder <- "~/code/clock_analysis" #alex

# just the necessary output is here:
fmri_dir <- "~/OneDrive - University of Pittsburgh/Documents/SCEPTIC_fMRI/explore_wholebrain/"
bad_ids <- c(216845, 440311, 440336)
setwd(file.path(fmri_dir, "L2_extracted_values/L1m-v_max"))

labels <- read_csv("../../region_labels_244.csv") %>% mutate(roi_num = as.numeric(roi_num)) %>% inner_join(read_csv("../../region_lookup_244.csv"), by = "roi_num") %>% 
  inner_join(read_csv("~/code/schaefer_wb_parcellation/labels/Schaefer_200_7networks_labels.csv")) %>%  rename(roi_num7 = "roi_num")
labels17 <- read_csv("~/code/schaefer_wb_parcellation/labels/Schaefer2018_200Parcels_17Networks_order.csv") %>% 
  inner_join(read_csv("~/code/schaefer_wb_parcellation/labels/Schaefer_200_17networks_labels.csv")) %>%  
  rename(network17 = "network", roi_num17 = "roi_num") %>% select(spatial_roi_num, roi_num17, network17)

# need to merge using spatial_roi_num from _labels CSV
# 
labels <- left_join(labels, labels17, by = "spatial_roi_num")

# we are missing 2 subjects here
sub_df <- readRDS("../../../explore_medusa/data/explore_n146.rds") %>% mutate(id = registration_redcapid)

# read in V_max betas

vmax <- read_csv("transformed_Schaefer_244_final_3.125mm_cope_l2.csv.gz") %>% filter(l2_model == "intercept" & l1_cope_name == "EV_v_max" & !(id %in% bad_ids)) %>%
  mutate(roi_num7 = mask_value) %>% select(-x, -y, -z) %>% inner_join(labels, by = "roi_num7") %>% inner_join(sub_df, by = "id") %>%
  select(-c(2:11)) 
# label missing networks as subcortical
vmax$network[is.na(vmax$network)] <- "subcortical"

# inspect "functional connectivity" by examining correlation structure
# need to fix this
test <- vmax %>% 
  filter(network %in% c("Default", "Limbic", "DorsAttn", "Subcortical") | roi_num7 > 200) %>% select(c(network, subregion, value, hemi, id)) %>%  mutate(region_network_label = paste(network, subregion, hemi, sep = "_")) 
vwide <- test %>% select(c(region_network_label, value,  id)) %>% group_by(region_network_label, id) %>% mutate(row = row_number()) %>% ungroup() %>%
  pivot_wider(names_from = c(region_network_label), values_from = value) %>% filter(row == 1) #%>% select(-id)
cormat <- psych::corr.test(vwide %>% select(-id))
pdf("explore_vmax_betas_corr1.pdf", width=50, height=50)
corrplot(cormat$r, col.lim=c(-1,1), #tl.pos = 'n',
         method = "color", type = "full", #tl.col  = 'black',tl.cex = 1.5, 
         order = "alphabet", diag = FALSE,
         # addCoef.col="black", addCoefasPercent = FALSE,
         p.mat = cormat$p, sig.level=0.0001, insig = "blank"
)
dev.off()

# "functional connectivity" for only the subcortical parcels
sub_vwide <- vmax %>% filter(roi_num7 > 200) %>% select(c(network, subregion, value, hemi, id)) %>%  mutate(region_network_label = paste(network, subregion, hemi, sep = "_")) %>%
  select(c(region_network_label, value,  id)) %>% group_by(region_network_label, id) %>% mutate(row = row_number()) %>% ungroup() %>%
  pivot_wider(names_from = c(region_network_label), values_from = value) %>% filter(row == 1) %>% select(-row)
sub_just_rois <- sub_vwide %>% select(-id)
cormat <- psych::corr.test(sub_just_rois)
pdf("explore_sub_vmax_betas_corr.pdf", width=30, height=30)
corrplot(cormat$r, col.lim=c(-1,1), #tl.pos = 'n',
         method = "color", type = "full", #tl.col  = 'black',tl.cex = 1.5, 
         order = "original", diag = FALSE,
         # addCoef.col="black", addCoefasPercent = FALSE,
         p.mat = cormat$p, sig.level=0.0001, insig = "blank"
)
dev.off()


# same for DMN
dmn_vwide <- vmax %>% filter(network == "Default") %>% select(c(subregion, value, hemi, id)) %>%  mutate(region_label = paste(subregion, hemi, sep = "_")) %>%
  select(c(region_label, value,  id)) %>% group_by(region_label, id) %>% mutate(row = row_number()) %>% ungroup() %>%
  pivot_wider(names_from = c(region_label), values_from = value) %>% filter(row == 1) %>% select(-row)
dmn_just_rois <- dmn_vwide %>% select(-id)
cormat <- psych::corr.test(dmn_just_rois)
pdf("explore_dmn_vmax_betas_corr.pdf", width=30, height=30)
corrplot(cormat$r, col.lim=c(-.1,1), #tl.pos = 'n',
         method = "color", type = "full", #tl.col  = 'black',tl.cex = 1.5, 
         order = "hclust", diag = FALSE,
         # addCoef.col="black", addCoefasPercent = FALSE,
         p.mat = cormat$p, sig.level=0.0001, insig = "blank"
)
dev.off()


# plot betas by group from the DMN
dmn <- vmax %>% filter(network == "Default")
ggplot(dmn, aes(Group_a, value, color = network17)) + geom_jitter(alpha = .3) + geom_violin(draw_quantiles = c(0.5))

dmn17b <- vmax %>% filter(str_detect(network17, "DefaultB"))

# test effects of group in the DMN
# m1 <- lmer(value ~ Group + subregion + (1|id), dmn)
# Anova(m1, '3')
m2 <- lmer(value ~ Group * subregion + hemi + (1|id), dmn)
Anova(m2, '3')
summary(m2)

m2_17 <- lmer(value ~ Group * subregion + hemi + (1|id), dmn17b)
Anova(m2_17, '3')
summary(m2_17)


# plot model-estimated means
em2 <- as_tibble(emmeans(m2, specs = c("Group", "subregion")))
pdf("dmn_vmax_em_betas_by_group_and_region.pdf", width = 16, height = 8)
ggplot(em2, aes(Group, emmean, color = Group, label = subregion)) + geom_jitter(size = 6, width = .01, height = 0)  + geom_violin(draw_quantiles = 0.5, alpha = 0.2) + geom_text_repel(point.padding = 10, force = 10,  size = 2, show.legend = F)
dev.off()
pdf("dmn_vmax_betas_by_group_and_region.pdf", width = 30, height = 30)
ggplot(dmn, aes(Group_a, value, label = subregion)) + geom_jitter(size = 6, width = .01, height = 0) + geom_boxplot() + geom_hline(aes(yintercept = 0)) + geom_text_repel()
dev.off()

em2a <-  emmeans(m2, trt.vs.ctrlk ~ Group|subregion, adjust = 'none')
dmn_stats <- as_tibble(em2a$contrasts) %>% mutate(p_fdr = p.adjust(p.value, method = 'fdr'),
                                                   asterisk  = case_when(
                                                     p_fdr < .05 ~ "*",
                                                     p_fdr > .05 ~ " ",
                                                   )) %>% inner_join(labels %>% select(subregion, network17))
pdf("dmn_vmax_by_network17.pdf", height = 8, width = 24)
ggplot(dmn_stats, aes(subregion, estimate, label = asterisk, color = network17)) + geom_point() + geom_text_repel() + geom_hline(aes(yintercept = 0)) +
  geom_errorbar(aes(ymax = estimate + SE, ymin = estimate - SE)) + facet_grid(rows = vars(contrast))
dev.off()

# nothing remarkable in just the PFC
# dmn_pfc <- vmax %>% filter(network == "Default" & str_detect(subregion, "PFC"))
# m3 <- lmer(value ~ Group * subregion + (1|id), dmn_pfc)
# Anova(m3, '3')
# summary(m3)
# em3 <- as_tibble(emmeans(m3, specs = c("Group", "subregion")))
# pd = position_dodge(0.35)
# ggplot(em3, aes(Group, emmean, color = Group)) + geom_jitter() + geom_violin(draw_quantiles = 0.5, alpha = 0.2) # geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL))
# ggplot(dmn_pfc, aes(Group, value, color = Group)) + geom_jitter(alpha = 0.3) + geom_boxplot()


# same for subcortical ROIs
sub <- vmax %>% filter(roi_num > 200)
ggplot(sub, aes(Group_a, value)) + geom_jitter(alpha = .3) + geom_violin()

pdf("subcortical_vmax_betas_by_group_and_region.pdf", width = 30, height = 30)
ggplot(sub, aes(Group_a, value)) + geom_jitter(alpha = .3) + geom_boxplot() + facet_wrap(~subregion) + geom_hline(aes(yintercept = 0))
dev.off()


# test group differences in subcortical ROIs
m4 <- lmer(value ~ Group * subregion + hemi + (1|id), sub)
Anova(m4, '3')
summary(m4)
em4 <- as_tibble(emmeans(m4, specs = c("Group", "subregion", "hemi")))
pdf("explore_subcortical_vmax_by_group.pdf", height = 6, width = 16)
ggplot(em4, aes(Group, emmean, color = Group, label = subregion)) + geom_jitter(size = 6, width = .05, height = 0)  + geom_violin(draw_quantiles = 0.5, alpha = 0.2) + geom_text_repel(point.padding = 10, force = 10,  size = 2, show.legend = F)
dev.off()
ggplot(sub, aes(Group, value, color = Group)) + geom_jitter(alpha = 0.3) + geom_boxplot() 

em4p <-  emmeans(m4, trt.vs.ctrlk ~ Group|subregion, adjust = 'fdr')
em4p_group_only <-  emmeans(m4, trt.vs.ctrlk ~ Group, adjust = 'none')

# test difference between networks
m4n <- lmer(value ~ Group_a * network + hemi + (1|id), sub)
Anova(m4n, '3')
summary(m4n)
str <- sub %>% filter(!str_detect(subregion, "Cerebellum") & !str_detect(subregion, "Hippo") &
                     !str_detect(subregion, "CMN") & !str_detect(subregion, "BLA") & !str_detect(subregion, "Thalamus"))
m4_bg <- lmer(value ~ Group * subregion  * hemi + (1|id), str)
Anova(m4_bg, '3')
summary(m4_bg)
# em4_bg <-  emmeans(m4_bg, trt.vs.ctrlk ~ Group|subregion, adjust = 'none')
em4_bg <-  emmeans(m4_bg, trt.vs.ctrlk ~ Group|subregion, adjust = 'none')
bg_stats <- as_tibble(em4_bg$contrasts) %>% mutate(p_fdr = p.adjust(p.value, method = 'fdr'),
                                                   asterisk  = case_when(
                                                     p_fdr < .05 ~ "*",
                                                     p_fdr > .05 ~ " ",
                                                   ))
ggplot(bg_stats, aes(subregion, estimate, label = asterisk)) + geom_point() + geom_text_repel() +
  geom_errorbar(aes(ymax = estimate + SE, ymin = estimate - SE)) + facet_grid(rows = vars(contrast))

em4_bg_group_only <-  emmeans(m4_bg, trt.vs.ctrlk ~ Group, adjust = 'fdr')


# control for confounds
m5 <- lmer(value ~ Group * subregion + hemi + scale(education_yrs) * subregion + scale(age) * subregion + suds_dx * subregion + 
             anxiety_dx * subregion + (1|id), sub)
Anova(m5, '3')
summary(m5)

m6 <- lmer(value ~ Group * network + hemi + scale(education_yrs) * network + scale(age) * network + suds_dx * network + 
             anxiety_dx * network + (1|id), sub)
Anova(m6, '3')
summary(m6)

# individual differences in subcortical value responses: look at all correlations
# positive WTAR effects, negative effects of impulsivity in cerebellum
sub_ind <- sub_vwide %>% inner_join(sub_df, by = "id")
rois <- names(sub_ind)[2:44]
diffs <- names(sub_ind)[c(56, 68:99)]
sub_cormat <- corr.test(sub_ind %>% select(rois), sub_ind %>% select(diffs))
pdf("explore_vmax_subcortical_rois_ind_diff_corrs.pdf", height = 30, width = 30)
corrplot(sub_cormat$r, col.lim = c(-.5, .5), p.mat = sub_cormat$p, insig = "blank")
dev.off()

# same for DMN ROIs
test_dmn <- vmax %>% 
  filter(network %in% c("Default") & roi_num < 201) %>% select(c(subregion, value, hemi, id)) %>%  mutate(region_network_label = paste(subregion, hemi, sep = "_")) 
dmn_vwide <- test_dmn %>% select(c(region_network_label, value,  id)) %>% group_by(region_network_label, id) %>% mutate(row = row_number()) %>% ungroup() %>%
  pivot_wider(names_from = c(region_network_label), values_from = value) %>% filter(row == 1) 
dmn_ind <- dmn_vwide %>% inner_join(sub_df, by = "id")
rois <- names(dmn_ind)[3:48]
diffs <- names(dmn_ind)[c(56, 68:99)]

# plot correlation with age of attempt
corr.test(dmn_ind$age_at_first_attempt, dmn_ind$PFC10_L)
ggplot(dmn_ind, aes(age_at_first_attempt, PFC10_L)) + geom_point() + geom_smooth(method = "glm")


# these are characteristics of suicidal behavior
suic <- names(dmn_ind)[c(59, 60, 64)]
dmn_cormat <- corr.test(dmn_ind %>% select(rois), dmn_ind %>% select(diffs))
pdf("explore_vmax_dmn_rois_ind_diff_corrs.pdf", height = 30, width = 30)
corrplot(dmn_cormat$r, col.lim = c(-.5, .5), p.mat = dmn_cormat$p, insig = "blank", addCoef.col = "black")
dev.off()

# correlations with characteristics of suicidal behavior: effects of age at first attempt
dmn_suic_cormat <- corr.test(dmn_ind %>% select(suic), dmn_ind %>% select(rois))
pdf("explore_vmax_dmn_suic_corrs_1.pdf", height = 6, width = 30)
corrplot(dmn_suic_cormat$r, col.lim = c(-.5, .5), p.mat = dmn_suic_cormat$p, addCoef.col = "black")
dev.off()

rois <- names(sub_ind)[2:44]
sub_suic_cormat <- corr.test(sub_ind %>% select(suic), sub_ind %>% select(rois))
pdf("explore_vmax_subcortical_suic_corrs.pdf", height = 6, width = 30)
corrplot(sub_suic_cormat$r, col.lim = c(-.5, .5), p.mat = sub_suic_cormat$p, insig = "blank")
dev.off()

# follow up with formal models, DMN
ms1 <- lmer(value ~ age_at_first_attempt * subregion + hemi + (1|id), dmn)
Anova(ms1, '3')
summary(ms1)
ems1 <- as_tibble(emmeans(ms1, specs = c("age_at_first_attempt", "subregion", "hemi"), at = list(age_at_first_attempt = c(20, 60))))
pdf("explore_dmn_vmax_age_at_first_attempt.pdf", height = 8, width = 16)
ggplot(ems1, aes(as.factor(age_at_first_attempt), emmean, label = subregion)) + geom_point() + 
  geom_violin(draw_quantiles = 0.5, alpha = 0.2) + geom_text_repel(size = 1.5) + facet_wrap(~hemi)
dev.off()
ggplot(ems1, aes(age_at_first_attempt, emmean, label = subregion, groups = subregion)) + geom_point() + geom_line() +
  geom_text_repel(size = 1.5) + facet_wrap(~hemi) + xlim(15, 65) 
ggplot(sub, aes(Group, value, color = Group)) + geom_jitter(alpha = 0.3) + geom_boxplot() 


# follow up with formal models, subcortical, not as strong
ms2 <- lmer(value ~ age_at_first_attempt * subregion + hemi + (1|id), sub)
Anova(ms2, '3')
summary(ms2)
ems2 <- as_tibble(emmeans(ms2, specs = c("age_at_first_attempt", "subregion", "hemi"), at = list(age_at_first_attempt = c(20, 60))))
pdf("explore_subcortical_vmax_age_at_first_attempt.pdf", height = 8, width = 16)
ggplot(ems2, aes(age_at_first_attempt, emmean, label = subregion, groups = subregion)) + geom_point() + geom_line()+ facet_wrap(~hemi)
dev.off()

geom_text_repel(size = 1.5) + facet_wrap(~hemi) + xlim(15, 65) 
# ggplot(sub, aes(Group, value, color = Group)) + geom_jitter(alpha = 0.3) + geom_boxplot() 


# wsub_ind <- vmax %>% filter(roi_num > 200) %>% select(c(id, value, network, hemi, subregion, names(vmax)[19:72])) %>%  mutate(region_network_label = paste(network, subregion, hemi, sep = "_")) %>%
#   select(c(region_network_label, value,  id)) %>% 
#   group_by(region_network_label, id) %>% mutate(row = row_number()) %>% ungroup() %>%
#   pivot_wider(names_from = c(region_network_label), values_from = value) %>% filter(row == 1) %>% inner_join(vmax %>% select(c(id, names(vmax)[19:72])), by = "id")
# 
# s <- nfactors(sub_just_rois, n=5, rotate = "oblimin", diagonal = FALSE,fm = "pa", n.obs = 70, SMC = FALSE)
# sub.fa = fa.sort(psych::fa(sub_just_rois, nfactors=4, rotate = "varimax", fm = "pa"))
# # 2 - cerebellum, 4 - striatum, 3 - hippocampus, 1 - thalamus
# fscores <- factor.scores(sub_just_rois, sub.fa)$scores
# sub_vwide$cerebellum <- fscores[,1]
# sub_vwide$striatum <- fscores[,2]
# sub_vwide$hippocampus <- fscores[,3]
# sub_vwide$thalamus <- fscores[,4]
# subcortical <- sub_vwide %>% select(c(id, cerebellum, striatum, hippocampus, thalamus)) %>% inner_join(sub_df, by = "id")
# sub_cormat <- corr.test(subcortical %>% select (c(cerebellum, striatum, hippocampus, thalamus)), subcortical[,25:56])
# corrplot(sub_cormat$r, col.lim = c(-.3, .3))
# 
# m6 <- lm(thalamus ~ Group + scale(education_yrs)  + scale(age)  + suds_dx  + anxiety_dx , subcortical)
# Anova(m6)
# summary(m6)
# 

# random slope of region: not converging
# m4_rs <- lmer(value ~ Group * subregion + hemi + (subregion|id), sub)
# Anova(m4_rs, '3')
# summary(m4_rs)
# em4_rs <- as_tibble(emmeans(m4, specs = c("Group", "subregion", "hemi")))
# pdf("explore_subcortical_vmax_by_group_rs.pdf", height = 6, width = 16)
# ggplot(em4_rs, aes(Group, emmean, color = Group, label = subregion)) + geom_jitter() + 
#   geom_violin(draw_quantiles = 0.5, alpha = 0.2) + geom_text_repel(size = 2) + facet_wrap(~hemi)
# dev.off()
