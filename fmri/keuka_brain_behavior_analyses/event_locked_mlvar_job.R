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
library(mlVAR)
library(gtools)
library(qgraph)

load('~/Box Sync/SCEPTIC_fMRI/var/feedback_hipp_wide_ts.Rdata')
vl1 <- mlVAR(fb_wide, vars = names(fb_wide[grep('_l', names(fb_wide))]), idvar = "id", lags = 1, dayvar = "run_trial", beepvar = "evt_time",
            estimator = "lm",
            contemporaneous = "unique", temporal = "unique",
            nCores = 1, verbose = TRUE, compareToLags = 1,
            scale = TRUE, scaleWithin = FALSE, AR = FALSE,
            iterations = "(2000)",
            chains = nCores
)

vr1 <- mlVAR(fb_wide, vars = names(fb_wide[grep('_r', names(fb_wide))]), idvar = "id", lags = 2, dayvar = "run_trial", beepvar = "evt_time",
             estimator = "lm",
             contemporaneous = "unique", temporal = "unique",
             nCores = 1, verbose = TRUE, 
             scale = TRUE, scaleWithin = FALSE, AR = FALSE,
             iterations = "(2000)",
             chains = nCores
)

rsort <- mixedsort(names(fb_wide[grep('_r', names(fb_wide))]))
lsort <- mixedsort(names(fb_wide[grep('_l', names(fb_wide))]))

plot(vr1, type = c("contemporaneous"), lag = 1, subject = 23)
setwd('~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/plots/')
pdf('individual_graphs_temporal_lag1.pdf', height = 6, width = 10)
par(mfrow = c(3,5))
for (sub in 1:71) {
plot(vr1, type = c("temporal"), lag = 1, subject = sub, order = rsort, minimum = 0.1, layout = "circle", groups = list(1:3,4:12), title = paste(as.character(sub), 'right'))
}
par(mfrow = c(3,5))
for (sub in 1:71) {
plot(vl1, type = c("temporal"), lag = 1, subject = sub, order = lsort, minimum = 0.1, layout = "circle", groups = list(1:3,4:12), title = paste(as.character(sub), 'left'))
}
dev.off()
# Sys.sleep(1.5)

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
  plot(cvr1, type = c("temporal"), lag = 1, subject = sub, order = rsort, minimum = 0.2, layout = "circle", groups = list(1:3,4:12), title = paste(as.character(sub), 'right'))
}
par(mfrow = c(3,5))
for (sub in 1:71) {
  plot(cvl1, type = c("temporal"), lag = 1, subject = sub, order = lsort, minimum = 0.2, layout = "circle", groups = list(1:3,4:12), title = paste(as.character(sub), 'left'))
}
dev.off()


# 
# vr1 <- mlVAR(fb_wide, vars = names(fb_wide[grep('_r', names(fb_wide))]), idvar = "id", lags = 1, dayvar = "run_trial", beepvar = "evt_time",
#             estimator = "lmer",
#             contemporaneous = "correlated", temporal = "fixed",
#             nCores = 8, verbose = TRUE, compareToLags = 1,
#             scale = TRUE, scaleWithin = FALSE, AR = FALSE,
#             iterations = "(2000)",
#             chains = nCores
# )
# save(vr1,"vr1.Rdata")
# layout(t(1:2))
# plot(v0, "temporal", title = "True temporal relationships", layout = "circle")
