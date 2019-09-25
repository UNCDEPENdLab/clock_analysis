###############
# Centrality and PE/entropy betas
###############
# analyze hippocampal PE and H betas vis-a-vis centrality metrics

library(tidyverse)
library(lme4)
library(emmeans)
library(readr)
library(car)

##########
# read in 

# grab voxelwise betas
pe <- read_csv("~/Box/SCEPTIC_fMRI/hippo_voxel_betas/sceptic-clock-feedback-pe_max-preconvolve_fse_groupfixed/pe_max/pe_max_atlas_betas.csv.gz")
# pe <- pe %>% mutate_at("beta", winsor,trim = .075)
pe <- pe %>% filter(beta > quantile(beta, .05) & beta < quantile(beta, .95))

h <- read_csv("~/Box/SCEPTIC_fMRI/hippo_voxel_betas/sceptic-clock-feedback-v_entropy-preconvolve_fse_groupfixed/v_entropy/v_entropy_atlas_betas.csv.gz")
# h <- h %>% mutate_at("beta", winsor,trim = .075)
h <- h %>% filter(beta > quantile(beta, .05) & beta < quantile(beta, .95))

pe$pe_beta <- pe$beta
h$h_beta <- h$beta
betas <- inner_join(pe[,c("ID", "numid", "vnum", "atlas_value", "pe_beta", "x", "y", "z")], h[,c("ID", "numid", "vnum", "atlas_value", "h_beta", "x", "y", "z")])
# ggplot(betas, aes(atlas_value, -h_beta)) + geom_smooth(method = 'gam') + facet_wrap(~x>0)

# reslice into 12 slices, summarize mean betas by slice, side
# right is x > 0
betas <- betas %>% mutate(slice = ntile(atlas_value,12), 
                          side =  case_when(
                            x > 0 ~ "right",
                            x < 0 ~ "left",
                            TRUE ~ NA_character_),
                          id = as.factor(ID))

# # sanity checks
# ggplot(betas, aes(as.factor(slice), -h_beta, color = side)) + geom_boxplot()
# ggplot(betas, aes(as.factor(slice), pe_beta, color = side)) + geom_boxplot()
# 
# # sanity check
# ggplot(betas, aes(slice, h_beta, color = side)) + geom_smooth()
# ggplot(betas, aes(slice, pe_beta, color = side)) + geom_smooth()

# load 6-slice centrality
load('~/OneDrive/collected_letters/papers/sceptic_fmri/hippo/brms_6slc_graphs.RData')
l <- gobj_lh$all_nodal
r <- gobj_rh$all_nodal
# visualize -- same as Michael's plots
# ggplot(l) + geom_boxplot(aes(node, strength_in), outlier.shape = NA) + geom_jitter(aes(node, strength_in, alpha = .1)) + 
#   scale_x_discrete(limits = c("lh2", "lh4", "lh6", "lh8", "lh10", "lh12"))
# ggplot(l) + geom_boxplot(aes(node, strength_out), outlier.shape = NA) + geom_jitter(aes(node, strength_out, alpha = .1)) + 
#   scale_x_discrete(limits = c("lh2", "lh4", "lh6", "lh8", "lh10", "lh12"))
l$id <- l$ID
l <- l %>% select(id, node, strength_in, strength_out)
# lw <- as.tibble(dcast(setDT(l), id ~ node, value.var = c("strength_in", "strength_out")))
l$side <- 'left'
r$id <- r$ID
r <- r %>% select(id, node, strength_in, strength_out)
r$side <- 'right'
# rw <- as.tibble(dcast(setDT(r), id ~ node, value.var = c("strength_in", "strength_out")))
# dc <- merge(lw,rw, by = "id")

## need to check this step!!!
c <- as_tibble(rbind(r,l))
c$slice <- as.numeric(parse_number(as.character(c$node)))
# merge betas with centrality 
cb <- as_tibble(merge(c,betas, by = c("id", "slice", "side")))
cb <- cb %>% mutate(strength_in = as.numeric(strength_in), strength_out = as.numeric(strength_out))

# # sanity checks
# ggplot(c, aes(as.factor(slice), strength_out)) + geom_boxplot(outlier.shape = NA) + geom_jitter(alpha = .2) + facet_wrap(~side)
# ggplot(c, aes(as.factor(slice), strength_in)) + geom_boxplot(outlier.shape = NA) + geom_jitter(alpha = .2) + facet_wrap(~side)
# 
# # ggplot(cb, aes(slice, scale(-h_beta), color = side)) + geom_smooth(method = 'gam') + geom_smooth(aes(slice, scale(pe_beta), color = side), method = 'gam') +
# #   scale_y_continuous(
# #     "Response to low entropy", 
# #     sec.axis = sec_axis(~ . * 1.20, name = "Response to prediction error")
# #   )
# ggplot(cb, aes(slice, pe_beta, color = side)) + geom_smooth(method = 'gam')
# 
# ggplot(cb, aes(strength_out, -h_beta)) + geom_smooth(method = 'gam')

cb$slice_f <- as.factor(cb$slice)



##########
# beta ~ centrality models
# test whether centrality predicts H beta 
mh1 <- lmer(-h_beta ~ slice_f + strength_out + side + (1|id), cb)
summary(mh1)
vif.lme(mh1)
mh2 <- lmer(-h_beta ~ slice_f * strength_out + side + (1|id), cb)
summary(mh2)
vif.lme(mh2)
h2 <- as.data.frame(emmeans::emmeans(mh2, "strength_out", by = "slice_f", at = list(strength_out = c(-.13, .39))))
ggplot(h2, aes(slice_f, emmean, group = strength_out, color = strength_out)) + geom_line() + geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL))

mh3 <- lmer(-h_beta ~ slice_f * strength_out * side + (1|id), cb)
summary(mh3)
vif.lme(mh3)
h3 <- as.data.frame(emmeans::emmeans(mh3, "strength_out", by = c("slice_f", "side"), at = list(strength_out = c(quantile(cb$strength_out, .25), cb$strength_out, .75))))
setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/hippo/figs/centrality/betas')
pdf('centrality_h_betas_no_outliers.pdf', width = 9, height = 4)
ggplot(h3, aes(slice_f, emmean, group = strength_out, color = strength_out)) + geom_line() + geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL)) + facet_wrap(~side)
dev.off()

mh3in <- lmer(-h_beta ~ slice_f * strength_in * side + (1|id), cb)
summary(mh3in)



# add REs -- that does not finish
mh4 <- lmer(-h_beta ~ slice_f * strength_out * side + (1 + slice_f + strength_out + side|id), cb)
summary(mh4)
h4 <- as.data.frame(emmeans::emmeans(mh4, "strength_out", by = c("slice_f", "side"), at = list(strength_out = c(quantile(cb$strength_out, .25), cb$strength_out, .75))))
setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/hippo/figs/centrality/betas')
pdf('centrality_h_betas.pdf', width = 12, height = 4)
ggplot(h4, aes(slice_f, emmean, group = strength_out, color = strength_out)) + geom_line() + geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL)) + facet_wrap(~side)
dev.off()

anova(mh1,mh2,mh3)

# adjust for strength_in
mh3inout <- lmer(-h_beta ~ slice_f * strength_out * side + slice_f * strength_in * side + (1|id), cb)
summary(mh3inout)
vif.lme(mh3inout)
h3inout <- as_tibble(emmeans::emmeans(mh3inout, "strength_out", by = c("slice_f", "side", "strength_in"), 
                                      at = list(strength_out = c(quantile(cb$strength_out, .25), quantile(cb$strength_out, .75)),strength_in = c(quantile(cb$strength_in, .25), quantile(cb$strength_in, .75)))))

pdf('centrality_h_betas_in_out_adjusted.pdf', width = 9, height = 9)
ggplot(h3inout, aes(slice_f, emmean, group = as.factor(strength_out), color = as.factor(strength_out))) + geom_line() + 
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL)) + facet_wrap(~ as.factor(strength_in)+ side, labeller = label_both)
dev.off()

# confusing... let's look at correlation plot
mat <- cb %>% select(h_beta, pe_beta, strength_in, strength_out, atlas_value)
cormat <- psych::corr.test(mat)

# plot the individual centrality profiles
pdf('ind_centrality_profiles.pdf', height = 10, width = 10)
ggplot(c, aes(slice,strength_out, color = side)) + geom_line() + facet_wrap(~id)
dev.off()
pdf('centrality_profiles_out.pdf', height = 6, width = 12)
ggplot(c) + geom_line(aes(slice,strength_out, color = id),alpha = .2) + geom_smooth(aes(slice,strength_out)) + theme(legend.position = "none") + facet_wrap(~side)
dev.off()
pdf('centrality_profiles_in.pdf', height = 6, width = 12)
ggplot(c) + geom_line(aes(slice,strength_in, color = id),alpha = .2) + geom_smooth(aes(slice,strength_in)) + theme(legend.position = "none") + facet_wrap(~side)
dev.off()

# pe
mpe1 <- lmer(pe_beta ~ slice_f + strength_out + side + (1|id), cb)
summary(mpe1)
vif.lme(mpe1)
mpe2 <- lmer(pe_beta ~ slice_f * strength_out + side + (1|id), cb)
summary(mpe2)
vif.lme(mpe2)
mpe3 <- lmer(pe_beta ~ slice_f * strength_out * side + (1|id), cb)
summary(mpe3)
vif.lme(mpe3)
p3 <- as.data.frame(emmeans::emmeans(mpe3, "strength_out", by = c("slice_f", "side"), at = list(strength_out = c(quantile(cb$strength_out, .25), cb$strength_out, .75))))
setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/hippo/figs/centrality')
pdf('centrality_pe_betas_no_outliers.pdf', width = 9, height = 4)
ggplot(p3, aes(slice_f, emmean, group = strength_out, color = strength_out)) + geom_line() + geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL)) + facet_wrap(~side)
dev.off()

anova(mpe1,mpe2,mpe3)



# OLD CODE, IGNORE FOR NOW
# h1 <- emmeans::emmip(mh1, axis_pos_f ~ Closeness, at = list(Closeness = c(0, .01, .02)))
# h1 <- emmeans::emmeans(mh1, "axis_pos_f", by = "Closeness", at = list(Closeness = c(0, .01, .02)))
h1 <- as.data.frame(emmeans::emmeans(mh1, "Closeness", by = "axis_pos_f", at = list(Closeness = c(0, .015, .03))))
h1$h_inv <- h1$emmean
hh <- ggplot(h1, aes(axis_pos_f, h_inv, color = Closeness)) + geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), position = position_dodge2(width =0.9))


mpe1 <- lmer(pe_beta ~ axis_pos_f + Closeness + side + (1|id), cb)
summary(mpe1)
vif.lme(mpe1)
pe1 <- as.data.frame(emmeans::emmeans(mpe1, "Closeness", by = "axis_pos_f", at = list(Closeness = c(0, .015, .03))))
pe1$pe_beta <- pe1$emmean
pp <- ggplot(pe1, aes(axis_pos_f, pe_beta, color = Closeness)) + geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), position = position_dodge2(width =0.9))

pdf('h_pe_closeness.pdf', height = 6, width = 12)
ggarrange(hh,pp)
dev.off()

# test interactions
mh2 <- lmer(-h_beta ~ axis_pos_f * scale(Closeness) + side + (1|id), cb)
summary(mh2)
vif.lme(mh2)
h2 <- as.data.frame(emmeans::emmeans(mh2, "Closeness", by = "axis_pos_f", at = list(Closeness = c(0, .015, .03))))
h2$h_inv <- h2$emmean
hh2 <- ggplot(h2, aes(axis_pos_f, h_inv, color = Closeness)) + geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), position = position_dodge2(width =0.9))
anova(mh1,mh2)

mh2c <- lmer(-h_beta ~ scale(axis_pos) * scale(Closeness) + side + (1|id), cb)
summary(mh2c)
vif.lme(mh2c)
h2c <- as.data.frame(emmeans::emmeans(mh2c, "Closeness", by = "axis_pos", at = list(Closeness = c(0, .015, .03), axis_pos = c(1,6,12))))
h2c$h_beta <- h2c$emmean
pp2c <- ggplot(h2c, aes(axis_pos, h_beta, color = Closeness)) + geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), position = position_dodge2(width =0.9))


mpe2 <- lmer(pe_beta ~ axis_pos_f * scale(Closeness) + side + (1|id), cb)
summary(mpe2)
vif.lme(mpe2)
pe2 <- as.data.frame(emmeans::emmeans(mpe2, "Closeness", by = "axis_pos_f", at = list(Closeness = c(0, .015, .03))))
pe2$pe_beta <- pe2$emmean
pp2 <- ggplot(pe2, aes(axis_pos_f, pe_beta, color = Closeness)) + geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), position = position_dodge2(width =0.9))
anova(mpe1,mpe2)

mpe2c <- lmer(pe_beta ~ scale(axis_pos)* scale(Closeness) + side + (1|id), cb)
summary(mpe2c)
vif.lme(mpe2c)
anova(mpe2,mpe2c)
pe2c <- as.data.frame(emmeans::emmeans(mpe2c, "Closeness", by = "axis_pos", at = list(Closeness = c(0, .015, .03), axis_pos = c(1,6,12))))
pe2c$pe_beta <- pe2c$emmean
pp2c <- ggplot(pe2c, aes(axis_pos, pe_beta, color = Closeness)) + geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), position = position_dodge2(width =0.9))

##########
# ICA analyses
cbi <- cb %>% group_by(id, slice, side) %>% summarise_at(vars(strength_in, strength_out, pe_beta, h_beta), mean) %>% 
  mutate(side_n = case_when(
    side == 'left' ~ 1,
    side == 'right' ~ 0
  )) %>% 
  ungroup()
library(fastICA)
library(ica)
i <- cbi %>% select(strength_in, strength_out, pe_beta, h_beta)
m <- icafast(i, 3)
plot(m$X, main = "Pre-processed data")
plot(m$X %*% m$K, main = "PCA components")
plot(m$S, main = "ICA components")
# ica_out <- cbind(signals = m$S,cbi %>% select(id,slice,side))
ica_out <- cbind(signals = m$S,cb %>% select(id,slice,side))

p1 <- ggplot(ica_out, aes(slice, side, fill = signals.1)) + geom_raster()
p2 <- ggplot(ica_out, aes(slice, side, fill = signals.2)) + geom_raster()
p3 <- ggplot(ica_out, aes(slice, side, fill = signals.3)) + geom_raster()

p4 <- ggplot(cbi, aes(slice, side, fill = strength_in)) + geom_raster()
p5 <- ggplot(cbi, aes(slice, side, fill = strength_out)) + geom_raster()
p6 <- ggplot(cbi, aes(slice, side, fill = pe_beta)) + geom_raster()
p7 <- ggplot(cbi, aes(slice, side, fill = -h_beta)) + geom_raster()

# p3 <- ggplot(ica_out, aes(slice, side, fill = signals.3)) + geom_raster()
# p4 <- ggplot(ica_out, aes(slice, side, fill = signals.4)) + geom_raster()

library(ggpubr)
pdf('ica_3_comp.pdf', height = 10, width = 5)
ggarrange(p1,p2,p3, p4, p5, p6, p7, nrow = 6, ncol = 1)
dev.off()
all <- as_tibble(cbind(m$S, cbi) %>% select_if(is.numeric))
cormat <- psych::corr.test(all)
pdf('ica_3_corrplot.pdf', height = 6, width = 6)
corrplot::corrplot(cormat$r,type = "upper")
dev.off()
cbi_scaled <- cbi %>% group_by(id) %>%  mutate_at(vars(strength_in,strength_out, pe_beta, h_beta), list(~as.vector(scale(.)))) %>% ungroup()
cbi_sc_tall <- cbi_scaled %>% gather(key = 'metric', value = 'value', strength_in,strength_out, pe_beta, h_beta)
pdf('individual_hippo_profiles.pdf', height = 30, width = 30)
ggplot(cbi_sc_tall, aes(slice, value, color = metric, lty = side)) + geom_line() + facet_wrap(~id)
dev.off()
ggplot(cbi_sc_tall, aes(slice, value, color = metric, lty = side)) + geom_smooth()
bdf <- as_tibble(cbind(m$S, cbi)) %>% mutate(ica1_hiIn_loOut = `1`, ica2_loIn_loNegH = `2`, ica3_loIn_hiNegH =`3`)
bdf_tall <- bdf %>% gather(key = 'component', value = 'value', ica1_hiIn_loOut, ica2_loIn_loNegH, ica3_loIn_hiNegH)

pdf('individual_hippo_ica_component_rasters.pdf', height = 30, width = 30)
ggplot(bdf_tall, aes(slice, id, fill = value)) + geom_tile() + facet_wrap(~component) + scale_fill_viridis_c()
dev.off()

ggplot(bdf_tall, aes(as.factor(slice), value, color = component, )) + geom_boxplot() + facet_wrap(~side)# + scale_fill_viridis_c()


pdf('ica_3_comp_slice_boxplots.pdf', height = 4, width = 8)
p2 <- ggplot(bdf_tall, aes(as.factor(slice), value, color = component, lty = side)) + geom_boxplot()
dev.off()

pdf('ica_3_comp_interpretation.pdf', height = 8, width = 8)
ggarrange(p1,p2, nrow = 1, ncol = 2)
dev.off()

# # spatial distribution: stronger gradient for 
summary(m1 <- lm(ica1_hiIn_loOut ~ as.factor(slice) * side,bdf))
Anova(m1)
summary(m2 <- lm(ica2_loIn_loNegH ~ as.factor(slice) * side,bdf))
Anova(m2)
summary(m3 <- lm(ica3_loIn_hiNegH ~ as.factor(slice) * side,bdf))
Anova(m3)
#ica on betas alone -- no, this does not make conceptual sense
# ib <- betas %>% select(pe_beta, h_beta)
# m <- icafast(ib, 2)
# bbdf <- as_tibble(cbind(m$S, betas))
# nums <- bbdf %>% select_if(is.numeric)s
# cormat <- psych::corr.test(nums)
# pdf('ica_3_corrplot.pdf', height = 6, width = 6)
# corrplot::corrplot(cormat$r,type = "upper")
# dev.off()
# 
# plot(bbdf, aes(atlas_value, `1`)) + geom_smooth()
