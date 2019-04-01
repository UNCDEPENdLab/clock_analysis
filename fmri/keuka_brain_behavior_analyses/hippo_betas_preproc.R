library(dplyr)
library(tidyverse)
library(psych)
library(corrplot)
library(lme4)

remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}

scale_within <- T

# start with entropy
setwd('~/Box Sync/skinner/projects_analyses/SCEPTIC/fMRI_paper/signals_review/MMClock_aroma_preconvolve_fse_groupfixed/hippo_voxel_betas/sceptic-clock-feedback-v_entropy-preconvolve_fse_groupfixed/v_entropy/')
hipH <- read_csv("v_entropy_atlas_betas.csv.gz")
hipH <- hipH[hipH$l2_contrast == 'overall',]
# ?winsorize
hipH$ID <- as.factor(hipH$ID)
hipH <- hipH %>% mutate_at(.vars = "beta", remove_outliers)
hipH$hemisphere <- NA
hipH$hemisphere[hipH$x>0] <- 'left'
hipH$hemisphere[hipH$x<0] <- 'right'

setwd('~/Box Sync/skinner/projects_analyses/SCEPTIC/fMRI_paper/signals_review/MMClock_aroma_preconvolve_fse_groupfixed/hippo_voxel_betas/sceptic-clock-feedback-pe_max-preconvolve_fse_groupfixed/pe_max/')
hipPE <- read_csv("pe_max_atlas_betas.csv.gz")
hipPE <- hipPE[hipPE$l2_contrast == 'overall',]
hipPE$hemisphere <- NA
hipPE$hemisphere[hipPE$x>0] <- 'left'
hipPE$hemisphere[hipPE$x<0] <- 'right'
hipPE$ID <- as.factor(hipPE$ID)

# ?winsorize
hipPE <- hipPE %>% mutate_at(.vars = "beta", remove_outliers)


setwd('~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/plots')
# sanity check: not bad
ggplot(hipH, aes(beta)) + geom_histogram()
pdf('long_axis_H_betas_RL.pdf', height = 6, width = 12)
h <- ggplot(hipH, aes(atlas_value,beta)) + geom_smooth(method = 'loess') + facet_wrap(~hemisphere)
dev.off()
ggplot(hipH, aes(x,y)) + geom_tile(aes(fill = beta)) 
ggplot(hipH, aes(x,z)) + geom_tile(aes(fill = beta)) 
# this is probably too much to ask
pdf('long_axis_H_betas_ind.pdf', height = 20, width = 20)
ggplot(hipH, aes(atlas_value,beta)) + geom_smooth(method = 'gam') + facet_wrap(~numid)
dev.off()
# sanity check: not bad
ggplot(hipPE, aes(beta)) + geom_histogram()
pdf('long_axis_PE_betas_RL.pdf', height = 6, width = 12)
ggplot(hipPE, aes(atlas_value,beta)) + geom_smooth(method = 'loess') + facet_wrap(~hemisphere)
dev.off()
ggplot(hipPE, aes(y,x)) + geom_tile(aes(fill = atlas_value)) 
ggplot(hipPE, aes(z,x)) + geom_tile(aes(fill = atlas_value)) 
# this is probably too much to ask
setwd('~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/plots')
pdf('long_axis_PE_betas_ind.pdf', height = 20, width = 20)
ggplot(hipPE, aes(atlas_value,beta, color = ID)) + geom_smooth(aes(atlas_value, beta, color = ID),method = 'gam') 
dev.off()

summary(m1 <- lmer(beta ~ atlas_value + (1|ID), hipPE))
summary(m2 <- lmer(beta ~ atlas_value + (1|ID), hipH))

ap <- 'Anterior <---> Posterior'
# plot jointly
h <- ggplot(hipH, aes(atlas_value,-beta)) + geom_smooth(method = 'loess') + 
  facet_wrap(~hemisphere) + geom_hline(yintercept=0, linetype="dashed", color = "red") + xlab(ap)
pe <- ggplot(hipPE, aes(atlas_value,beta)) + geom_smooth(method = 'loess') + 
  facet_wrap(~hemisphere) + geom_hline(yintercept=0, linetype="dashed", color = "red") + xlab(ap)
# "MAP"
mz <- ggplot(hipPE, aes(x,z)) + geom_tile(aes(fill = atlas_value)) + xlab('Right <---> Left')
my <- ggplot(hipPE, aes(-y,x)) + geom_tile(aes(fill = atlas_value)) + xlab(ap)

pdf('long_axis_PE_H_betas.pdf', height = 20, width = 10)
ggarrange(pe,h,mz,my, ncol = 1, nrow = 4, labels = c('Prediction error', 'Entropy (reversed)'))
dev.off()

ggplot(hipH, aes(atlas_value,-beta)) + geom_smooth(method = 'loess') + facet_wrap(~hemisphere) + geom_hline(yintercept=0, linetype="dashed", color = "red")
ggplot(hipH) + geom_smooth(aes(atlas_value,-beta, color = ID), method = 'gam', se = F) + geom_smooth(aes(atlas_value,-beta), method = 'gam', size = 1.5) + 
  facet_wrap(~hemisphere) + geom_hline(yintercept=0, linetype="dashed", color = "red") + xlab(ap) + theme(legend.position = "none")
ggplot(hipPE) + geom_smooth(aes(atlas_value,beta, color = ID), method = 'gam', se = F) + geom_smooth(aes(atlas_value,beta), method = 'gam', size = 1.5) + 
  facet_wrap(~hemisphere) + geom_hline(yintercept=0, linetype="dashed", color = "red") + xlab(ap) + theme(legend.position = "none")
