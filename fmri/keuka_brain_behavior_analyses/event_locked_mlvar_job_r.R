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

load('~/Box Sync/SCEPTIC_fMRI/var/feedback_hipp_wide_ts.Rdata')
 
vr1 <- mlVAR(fb_wide, vars = names(fb_wide[grep('_r', names(fb_wide))]), idvar = "id", lags = 1, dayvar = "run_trial", beepvar = "evt_time",
            estimator = "lmer",
            contemporaneous = "correlated", temporal = "correlated",
            nCores = 10, verbose = TRUE, compareToLags = 1,
            scale = TRUE, scaleWithin = FALSE, AR = FALSE,
            iterations = "(2000)",
            chains = nCores
)
save(vr1,"vr1.Rdata")
# layout(t(1:2))
# plot(v0, "temporal", title = "True temporal relationships", layout = "circle")
