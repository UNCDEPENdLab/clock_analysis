###############
# Centrality and PE/entropy betas
###############

library(tidyverse)
library(lme4)
library(emmeans)
library(readr)

# read in betas
# grab voxelwise betas
pe <- read_csv("~/Box/SCEPTIC_fMRI/hippo_voxel_betas/sceptic-clock-feedback-pe_max-preconvolve_fse_groupfixed/pe_max/pe_max_atlas_betas.csv.gz")
# pe <- pe %>% mutate_at("beta", winsor,trim = .075)
pe <- pe %>% filter(beta > quantile(beta, .025) & beta < quantile(beta, .975))

h <- read_csv("~/Box/SCEPTIC_fMRI/hippo_voxel_betas/sceptic-clock-feedback-v_entropy-preconvolve_fse_groupfixed/v_entropy/v_entropy_atlas_betas.csv.gz")
# h <- h %>% mutate_at("beta", winsor,trim = .075)
h <- h %>% filter(beta > quantile(beta, .025) & beta < quantile(beta, .975))

pe$pe_beta <- pe$beta
h$h_beta <- h$beta
betas <- inner_join(pe[,c("ID", "numid", "vnum", "atlas_value", "pe_beta", "x", "y", "z")], h[,c("ID", "numid", "vnum", "atlas_value", "h_beta", "x", "y", "z")])
ggplot(betas, aes(atlas_value, -h_beta)) + geom_smooth(method = 'gam') + facet_wrap(~x>0)


# reslice into 12 slices, summarize mean betas by slice, side
# right is x > 0
betas <- betas %>% mutate(slice = ntile(atlas_value,12), 
                          side =  case_when(
                            x > 0 ~ "right",
                            x < 0 ~ "left",
                            TRUE ~ NA_character_),
                          id = as.factor(ID))

# sanity checks
ggplot(betas, aes(as.factor(slice), -h_beta, color = side)) + geom_boxplot()
ggplot(betas, aes(as.factor(slice), pe_beta, color = side)) + geom_boxplot()

# sanity check
# ggplot(betas, aes(axis_pos, h_beta, color = side)) + geom_smooth()

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

# sanity checks
ggplot(c, aes(as.factor(slice), strength_out)) + geom_boxplot(outlier.shape = NA) + geom_jitter(alpha = .1) + facet_wrap(~side)
ggplot(c, aes(as.factor(slice), strength_in)) + geom_boxplot(outlier.shape = NA) + geom_jitter(alpha = .1) + facet_wrap(~side)

ggplot(cb, aes(slice, scale(-h_beta), color = side)) + geom_smooth(method = 'gam') + geom_smooth(aes(slice, scale(pe_beta), color = side), method = 'gam') +
  scale_y_continuous(
    "Response to low entropy", 
    sec.axis = sec_axis(~ . * 1.20, name = "mpg (UK)")
  )
ggplot(cb, aes(slice, pe_beta, color = side)) + geom_smooth(method = 'gam')

ggplot(cb, aes(strength_out, -h_beta)) + geom_smooth(method = 'gam')


vif.lme <- function (fit) {
  ## adapted from rms::vif
  v <- vcov(fit)
  nam <- names(fixef(fit))
  ## exclude intercepts
  ns <- sum(1 * (nam == "Intercept" | nam == "(Intercept)"))
  if (ns > 0) {
    v <- v[-(1:ns), -(1:ns), drop = FALSE]
    nam <- nam[-(1:ns)] }
  d <- diag(v)^0.5
  v <- diag(solve(v/(d %o% d)))
  names(v) <- nam
  v }
cb$slice_f <- as.factor(cb$slice)
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



# add REs
mh4 <- lmer(-h_beta ~ slice_f * strength_out * side + (1 + slice_f * strength_out + slice_f * side + strength_out * side|id), cb)
summary(mh4)
h4 <- as.data.frame(emmeans::emmeans(mh4, "strength_out", by = c("slice_f", "side"), at = list(strength_out = c(quantile(cb$strength_out, .25), cb$strength_out, .75))))
setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/hippo/figs/centrality/betas')
pdf('centrality_h_betas_REs.pdf', width = 12, height = 4)
ggplot(h4, aes(slice_f, emmean, group = strength_out, color = strength_out)) + geom_line() + geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL)) + facet_wrap(~side)
dev.off()

anova(mh1,mh2,mh3)

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




########
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


# merge centrality metrics with voxelwise betas

