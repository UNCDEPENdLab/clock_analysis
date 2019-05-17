library(dplyr)
library(tidyverse)
library(psych)
library(ggcorrplot)
library(lme4)
library(ggpubr)
library(cowplot)
# library(sjPlot)
# library(sjmisc)
library(ggeffects)
library(gtools)
library(qgraph)
library(mlVAR)

load('~/Box Sync/SCEPTIC_fMRI/var/feedback_hipp_wide_ts.Rdata')

slices = 12

if (slices == 12) {df <- fb_wide; grp <- list(1:3,4:12); blgrp <- list(c(1:12), c(13:24))}
if (slices == 6) {df <- fb_wide6; grp <- list(1:2,3:6); blgrp <- list(c(1:6*2), c(1:6*2-1))}
rsort <- mixedsort(names(df[grep('_r', names(df))]))
lsort <- mixedsort(names(df[grep('_l', names(df))]))
blsort <- c(rsort,rev(lsort))

# sorted data
dfs <- df[,c(names(df)[1:4], rsort, lsort)]



vl1 <- mlVAR(df, vars = names(df[grep('_l', names(df))]), idvar = "id", lags = 1, dayvar = "run_trial", beepvar = "evt_time",
            estimator = "lm",
            contemporaneous = "unique", temporal = "unique",
            nCores = 1, verbose = TRUE, 
            scale = TRUE, scaleWithin = FALSE, AR = FALSE,
            iterations = "(2000)",
            chains = nCores
)
# sanity check
plot(vl1, lag = 1, layout = "circle", order = lsort, groups = grp, title = 'mean left', threshold = .02, edge.labels = T)
plot(vl1, lag = 2, subject = 3, type = "contemporaneous", order = lsort, groups = grp)
vr1 <- mlVAR(df, vars = names(df[grep('_r', names(df))]), idvar = "id", lags = 1:4, dayvar = "run_trial", beepvar = "evt_time",
             estimator = "lm",
             contemporaneous = "unique", temporal = "unique",
             nCores = 1, verbose = TRUE, 
             scale = TRUE, scaleWithin = FALSE, AR = FALSE,
             iterations = "(2000)",
             chains = nCores
)

vrl1 <- mlVAR(df, vars = names(df[grep('hipp', names(df))]), idvar = "id", lags = 1, dayvar = "run_trial", beepvar = "evt_time",
             estimator = "lm",
             contemporaneous = "unique", temporal = "unique",
             nCores = 1, verbose = TRUE, 
             scale = TRUE, scaleWithin = FALSE, AR = FALSE,
             iterations = "(2000)",
             chains = nCores
)
plot(vrl1, lag = 1, layout = "circle", subject = 1, title = 'bilateral', threshold = .4, edge.labels = T, groups = blgrp, order = blsort)
plot(vrl1, subject = 10, type = "temporal", layout = "circle", title = 'sub 1 temporal', threshold = .1)
plot(vrl1, type = "contemporaneous")

pdf(paste('fixed_effects_graphs_by_lag', as.character(slices), 'slices.pdf', sep = '_'), height = 6, width = 12)
par(mfrow = c(2,4))
for (t in 1:4) {
  plot(vr1, type = c("temporal"), lag = t, order = rsort, minimum = 0.01, layout = "circle", groups = grp, title = paste('lag', as.character(t), 'right'), edge.labels = T)
}
for (t in 1:4) {
  plot(vl1, type = c("temporal"), lag = t, order = lsort, minimum = 0.01, layout = "circle", groups = grp, title = paste('lag ', as.character(t), 'left'), edge.labels = T)
}
dev.off()

plot(vr1, type = c("contemporaneous"), order = rsort, minimum = 0.01, layout = "circle", groups = grp, title = paste( 'right'), edge.labels = T)


pdf(paste('individual_graphs_by_lag_10_subjects', as.character(slices), 'slices.pdf', sep = '_'), height = 6, width = 12)
for (sub in 1:10) {
par(mfrow = c(2,4))
for (t in 1:4) {
  plot(vr1, type = c("temporal"), lag = t, subject = sub, order = rsort, minimum = 0.1, layout = "circle", groups = grp, title = paste('subject', as.character(sub),'lag', as.character(t), 'right'))
}
for (t in 1:4) {
  plot(vl1, type = c("temporal"), lag = t, subject = sub, order = lsort, minimum = 0.1, layout = "circle", groups = grp, title = paste('subject', as.character(sub),'lag ', as.character(t), 'left'))
}
}
dev.off()

pdf(paste('individual_bl_graphs_lag1', as.character(slices), 'slices.pdf', sep = '_'), height = 6, width = 12)
par(mfrow = c(4,8))
for (sub in 1:71) {
  plot(vrl1, lag = 1, layout = "circle", subject = sub, title = paste('subject', as.character(sub)), threshold = .2, edge.labels = T, groups = blgrp, order = blsort)
}
for (sub in 1:71) {
  plot(vrl1, type = "contemporaneous", layout = "circle", subject = sub, title = paste('subject', as.character(sub)), threshold = .2, edge.labels = T, groups = blgrp, order = blsort)
}

dev.off()

plot(vr1, type = c("temporal"), lag = 1)
setwd('~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/plots/')
pdf('individual_graphs_temporal_lag1_min0.1.pdf', height = 6, width = 10)
par(mfrow = c(3,5))
for (sub in 1:71) {
plot(vr1, type = c("temporal"), lag = 1, subject = sub, order = rsort, minimum = 0.1, layout = "circle", groups = list(1:3,4:12), title = paste(as.character(sub), 'right'))
}
par(mfrow = c(3,5))
for (sub in 1:71) {
plot(vl1, type = c("temporal"), lag = 1, subject = sub, order = lsort, minimum = 0.1, layout = "circle", groups = list(1:3,4:12), title = paste(as.character(sub), 'left'))
}
dev.off()

plot(vl1, type = c("contemporaneous"), partial = T, order = lsort, layout = "circle", subject = 2, threshold = .2)

pdf("fixed_effects_group_graphs_temporal_lag1.pdf", height = 4, width = 8)
par(mfrow = c(1,2))
plot(vl1, type = c("temporal"), lag = 1, edge.labels = T, threshold = .01, layout = "circle", order = lsort, groups = list(1:3,4:12), title ='mean left')
plot(vr1, type = c("temporal"), lag = 1, edge.labels = T, threshold = .01, layout = "circle", order = rsort, groups = list(1:3,4:12), title = 'mean right')
dev.off()
########
# clock
########


load('~/Box Sync/SCEPTIC_fMRI/var/clock_hipp_wide_ts.Rdata')
cvl1 <- mlVAR(clock_wide, vars = names(clock_wide[grep('_l', names(clock_wide))]), idvar = "id", lags = 1, dayvar = "run_trial", beepvar = "evt_time",
             estimator = "lm",
             contemporaneous = "unique", temporal = "unique",
             nCores = 1, verbose = TRUE, compareToLags = 1,
             scale = TRUE, scaleWithin = FALSE, AR = FALSE,
             iterations = "(2000)",
             chains = nCores
)

cvr1 <- mlVAR(clock_wide, vars = names(clock_wide[grep('_r', names(clock_wide))]), idvar = "id", lags = 2, dayvar = "run_trial", beepvar = "evt_time",
             estimator = "lm",
             contemporaneous = "unique", temporal = "unique",
             nCores = 1, verbose = TRUE, 
             scale = TRUE, scaleWithin = FALSE, AR = FALSE,
             iterations = "(2000)",
             chains = nCores
)

rsort <- mixedsort(names(clock_wide[grep('_r', names(clock_wide))]))
lsort <- mixedsort(names(clock_wide[grep('_l', names(clock_wide))]))

setwd('~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/plots/')
pdf('individual_clock_graphs_temporal_lag1.pdf', height = 6, width = 10)
par(mfrow = c(3,5))
for (sub in 1:71) {
  plot(cvr1, type = c("temporal"), lag = 1, subject = sub, order = rsort, minimum = 0.1, layout = "circle", groups = list(1:3,4:12), title = paste(as.character(sub), 'right'))
}
par(mfrow = c(3,5))
for (sub in 1:71) {
  plot(cvl1, type = c("temporal"), lag = 1, subject = sub, order = lsort, minimum = 0.1, layout = "circle", groups = list(1:3,4:12), title = paste(as.character(sub), 'left'))
}
dev.off()

# estimate subject-wise indices
#############
# parallel to Michael's logic, kind of
# run on sorted data
vrl1_1 <- mlVAR(dfs, vars = names(dfs[grep('_', names(dfs))]), idvar = "id", lags = 1, dayvar = "run_trial", beepvar = "evt_time",
             estimator = "lm",
             contemporaneous = "unique", temporal = "unique",
             nCores = 1, verbose = TRUE, compareToLags = 1:3,
             scale = TRUE, scaleWithin = FALSE, AR = FALSE,
             iterations = "(2000)",
             chains = nCores
)
vrl1_2 <- mlVAR(dfs, vars = names(dfs[grep('_', names(dfs))]), idvar = "id", lags = 1:2, dayvar = "run_trial", beepvar = "evt_time",
                estimator = "lm",
                contemporaneous = "unique", temporal = "unique",
                nCores = 1, verbose = TRUE, compareToLags = 1:3,
                scale = TRUE, scaleWithin = FALSE, AR = FALSE,
                iterations = "(2000)",
                chains = nCores
)
vrl1_3 <- mlVAR(dfs, vars = names(dfs[grep('_', names(dfs))]), idvar = "id", lags = 1:3, dayvar = "run_trial", beepvar = "evt_time",
                estimator = "lm",
                contemporaneous = "unique", temporal = "unique",
                nCores = 1, verbose = TRUE, compareToLags = 1:3,
                scale = TRUE, scaleWithin = FALSE, AR = FALSE,
                iterations = "(2000)",
                chains = nCores
)

mlVARcompare(vrl1_1,vrl1_2,vrl1_3)
# check results
plot(vrl1, lag = 3, layout = "spring", subject = 1, title = paste('subject 1'), threshold = .2, edge.labels = T, groups = blgrp, order = blsort)

# make long df for plots
dfl <- as_tibble(matrix(nrow = 0, ncol = 0))
for (sub in 1:71) {
  gg <- getNet(vrl1, lag = lags, layout = "circle", subject = sub, title = paste('subject 1'), threshold = .2, edge.labels = T, groups = blgrp, order = blsort)
  
# gg <- getNet(vrl1, type = "temporal", lag = 1, subject = sub, layout = 'circle')
# get centrality measures for this subject
l <- qgraph::centrality(gg)
m <- as_tibble(data.frame(l))
mlong <- m[,c("OutDegree", "InDegree", "Closeness", "Betweenness", "InExpectedInfluence", "OutExpectedInfluence")]
mlong$region <- as.character(names(l$OutDegree))
mlong$id <- as.factor(unique(dfs$id)[sub])
mlong$lag <- lags
dfl <- rbind(dfl,mlong)
}
library(stringr)
numextract <- function(string){ 
  str_extract(string, "\\-*\\d+\\.*\\d*")
} 
extract.alpha <- function(x, space = ""){      
  require(stringr)
  require(purrr)
  require(magrittr)
  
  y <- strsplit(unlist(x), "[^a-zA-Z]+") 
  z <- y %>% map(~paste(., collapse = space)) %>% simplify()
  return(z)}
dfl$axis_pos <- as.numeric(numextract(dfl$region))
dfl$side <- extract.alpha(dfl$region, space = "_")
dfl$side <- recode(dfl$side, hipp_l = "left", hipp_r = "right")
i <- ggplot(na.omit(dfl), aes(axis_pos,InDegree, color = side)) + geom_smooth()
o <- ggplot(na.omit(dfl), aes(axis_pos,OutDegree, color = side)) + geom_smooth()
c <- ggplot(na.omit(dfl), aes(axis_pos,Closeness, color = side)) + geom_smooth()
b <- ggplot(na.omit(dfl), aes(axis_pos,Betweenness, color = side)) + geom_smooth()
ii <- ggplot(na.omit(dfl), aes(axis_pos,abs(InExpectedInfluence), color = side)) + geom_smooth()
oi <- ggplot(na.omit(dfl), aes(axis_pos,abs(OutExpectedInfluence), color = side)) + geom_smooth()
pdf("hipp_centrality_metrics_lag_1.pdf", height = 10, width = 16)
ggarrange(i,o,c,b,ii,oi, ncol = 3, nrow = 2)
dev.off()
# inspect individual curves
pdf("hipp_ind_closeness_ap_plots.pdf", height = 20, width = 20)
ggplot(dfl, aes(axis_pos,Closeness, color = side)) + geom_point() + facet_wrap(~id)
dev.off()
myspread <- function(df, key, value) {
  # quote key
  keyq <- rlang::enquo(key)
  # break value vector into quotes
  valueq <- rlang::enquo(value)
  s <- rlang::quos(!!valueq)
  df %>% gather(variable, value, !!!s) %>%
    unite(temp, !!keyq, variable) %>%
    spread(temp, value)
}
# spread dfl by 
dfw <- dfl %>% spread(dfl, key = )

# grab voxelwise betas
pe <- read_csv("~/Box Sync/SCEPTIC_fMRI/hippo_voxel_betas/sceptic-clock-feedback-pe_max-preconvolve_fse_groupfixed/pe_max/pe_max_atlas_betas.csv.gz")
h <- read_csv("~/Box Sync/SCEPTIC_fMRI/hippo_voxel_betas/sceptic-clock-feedback-v_entropy-preconvolve_fse_groupfixed/v_entropy/v_entropy_atlas_betas.csv.gz")
pe$pe_beta <- pe$beta
h$h_beta <- h$beta
betas <- inner_join(pe[,c("ID", "numid", "vnum", "atlas_value", "pe_beta", "x", "y", "z")], h[,c("ID", "numid", "vnum", "atlas_value", "h_beta", "x", "y", "z")])
ggplot(betas, aes(atlas_value, -h_beta)) + geom_smooth() + facet_wrap(~x>0)


# reslice into 12 slices, summarize mean betas by slice, side
# right is x > 0
betas <- betas %>% mutate(axis_pos = ntile(atlas_value,12), 
                          side =  case_when(
                            x > 0 ~ "right",
                            x < 0 ~ "left",
                            TRUE ~ NA_character_),
                          id = as.factor(ID))
# sanity check
# ggplot(betas, aes(axis_pos, h_beta, color = side)) + geom_smooth()

# merge betas with centrality 
cb <- inner_join(betas, dfl, by = c("id", "axis_pos", "side"))
# sanity checks
ggplot(cb, aes(axis_pos, Closeness, color = side)) + geom_smooth()
ggplot(cb, aes(axis_pos, -h_beta, color = side)) + geom_smooth()
ggplot(cb, aes(Closeness, pe_beta, color = side)) + geom_smooth(method = 'gam')


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
cb$axis_pos_f <- as.factor(cb$axis_pos)
# test whether centrality predicts H beta 
mh1 <- lmer(-h_beta ~ axis_pos_f + Closeness + side + (1|id), cb)
summary(mh1)
vif.lme(mh1)
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

