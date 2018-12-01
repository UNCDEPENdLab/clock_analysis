library(dplyr)
library(tidyverse)
library(psych)
library(corrplot)
library(lme4)

meta <- read_csv("~/code/clock_analysis/fmri/data/fmri_betas/v_entropy_cluster_metadata.csv")
meta_overall <- meta[meta$l2_contrast == 'overall' & meta$l3_contrast == 'Intercept',]
Hbetas <- read_csv("~/code/clock_analysis/fmri/data/fmri_betas/v_entropy_roi_betas.csv")
h <- as.tibble(Hbetas[Hbetas$l2_contrast == 'overall' & Hbetas$l3_contrast == 'Intercept',1:3]) %>% filter(cluster_number<11)
# head(merge(h,meta))
head(merge(h,rois))
h_wide <- spread(h,cluster_number,cope_value) 
# head(h_wide)
h_wide <- h_wide %>% mutate_if(is.double, winsor,trim = .075)

just_rois <- h_wide[,2:7]
# winsorize to deal with beta ouliers

# non-parametric correlations to deal with outliers
clust_cor <- cor(just_rois,method = 'spearman')
# parametric correlations on winsorised betas
clust_cor <- cor(just_rois_w,method = 'pearson')

rois <- distinct(meta_overall[,c(4,11)])
setwd('~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/')
pdf("h_cluster_corr.pdf", width=12, height=12)
corrplot(clust_cor, cl.lim=c(-1,1),
         method = "circle", tl.cex = 1.5, type = "upper", tl.col = 'black',
         order = "hclust", diag = FALSE,
         addCoef.col="orange", addCoefasPercent = FALSE,
         p.mat = 1-clust_cor, sig.level=0.5, insig = "blank")
dev.off()

h.pca = prcomp((just_rois),scale = TRUE, center = TRUE)
h_wide$h_pc1 <- h.pca$x[,1]
h_wide$h_pc2 <- h.pca$x[,2]

h.fa = psych::fa(just_rois, nfactors=2)
fscores <- factor.scores(just_rois_w, h.fa)$scores
h_wide$h_f1 <- fscores[,1]
h_wide$h_f2 <- fscores[,2]
h_wide$ParaHippL <- h_wide$`5`
h_wide$vmPFC <- h_wide$`6`

# add ids
map_df <- read_delim("~/code/clock_analysis/fmri/data/fmri_betas/map_df.txt", 
                     " ", escape_double = FALSE, trim_ws = TRUE)
pc_scores <- inner_join(h_wide[,c('subject','h_pc1','h_pc2', 'h_f1', 'h_f2', 'ParaHippL', 'vmPFC')],map_df)
# get trial_level data
trial_df <- read_csv("~/code/clock_analysis/fmri/data/mmclock_fmri_decay_factorize_selective_psequate_mfx_trial_statistics.csv.gz")
trial_df <- trial_df %>%
  group_by(id, run) %>%  dplyr::mutate(rt_swing = abs(c(NA, diff(rt_csv))), 
                                       rt_swing_lr = abs(log(rt_csv/lag(rt_csv))), 
                                       rt_lag = lag(rt_csv) , 
                                       omission_lag = lag(score_csv==0),
                                       run_trial=1:50) %>% ungroup() #compute rt_swing within run and subject

# performance
sum_df <- trial_df %>% group_by(id) %>% summarize(total_earnings = sum(score_csv)) %>% arrange(total_earnings)
beta_sum <- inner_join(pc_scores,sum_df)
plot(beta_sum$h_f1,beta_sum$total_earnings)
plot(beta_sum$h_f2,beta_sum$total_earnings)

# model parameters
params <- read_csv("~/code/clock_analysis/fmri/data/mmclock_fmri_decay_factorize_selective_psequate_mfx_sceptic_global_statistics.csv")
sub_df <- inner_join(beta_sum,params)
params_beta <- sub_df[,c("h_pc1","h_pc2","h_f1", "h_f2","ParaHippL", "vmPFC", "total_earnings", "LL", "alpha", "gamma", "beta")]
param_cor <- cor(params_beta,method = 'pearson')

setwd('~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/')
pdf("hbeta_param_corr.pdf", width=12, height=12)
corrplot(param_cor, cl.lim=c(-1,1),
         method = "circle", tl.cex = 1.5, type = "upper", tl.col = 'black',
         order = "hclust", diag = FALSE,
         addCoef.col="orange", addCoefasPercent = FALSE,
         p.mat = 1-param_cor, sig.level=0.5, insig = "blank")
dev.off()

# merge into trial-level data
df <- inner_join(trial_df,sub_df)
df$rewFunc <- relevel(as.factor(df$rewFunc),ref = "DEV")

summary(m01 <- lmer(log(rt_swing_lr) ~ scale(run_trial) * rewFunc + (1|id/run), df[df$rt_swing>0,]))
summary(m02 <- lmer(log(rt_swing_lr) ~ (scale(run_trial) + rewFunc + h_pc1 + h_pc2)^3 + (1|id/run), df[df$rt_swing>0,]))
car::Anova(m02,'3')
anova(m01,m02)


# RT swings analyses: "exploration"
summary(m1 <- lmer(log(rt_swing) ~ scale(run_trial) * scale(v_max) + rewFunc + (1|id/run), df[df$rt_swing>0,]))
# keeps flipping between log and untransformed...
summary(m2 <- lmer(log(rt_swing) ~ (scale(run_trial) + scale(v_max) + h_pc1) ^2 + rewFunc + (1|id/run), df[df$rt_swing>0,]))
car::Anova(m2, '3')
summary(m3 <- lmer(log(rt_swing) ~ (scale(run_trial) + scale(v_max) + h_pc2) ^2 + rewFunc + (1|id/run), df[df$rt_swing>0,]))
car::Anova(m3, '3')
anova(m1,m2,m3)

# EV analyses: "exploitation"
summary(m4 <- lmer(ev ~ scale(run_trial) * scale(v_max) + (1|id/run), df[df$rt_swing>0,]))
summary(m5 <- lmer(ev ~ (scale(run_trial) + scale(v_max) + h_pc1) ^2 + (1|id/run), df[df$rt_swing>0,]))
car::Anova(m5, '3')
summary(m6 <- lmer(ev ~ (scale(run_trial) + scale(v_max) + h_pc2) ^2 + (1|id/run), df[df$rt_swing>0,]))
car::Anova(m6, '3')

anova(m4,m5,m6)

# coupling between entropy and RT swings
summary(m7 <- lmer(log(rt_swing) ~ (scale(run_trial) + scale(v_max) + scale(v_entropy))^2 + rewFunc + (1|id/run), df[df$rt_swing>0,]))
summary(m8 <- lmer(log(rt_swing) ~ (scale(run_trial) + scale(v_max) + scale(v_entropy) + h_pc1) ^2 + rewFunc + (1|id/run), df[df$rt_swing>0,]))
car::Anova(m8, '3')
summary(m9 <- lmer(log(rt_swing) ~ (scale(run_trial) + scale(v_max) + scale(v_entropy) + h_pc2) ^2 + rewFunc + (1|id/run), df[df$rt_swing>0,]))
car::Anova(m9, '3')
summary(m10 <- lmer(log(rt_swing) ~ (scale(run_trial) + scale(v_max) + scale(v_entropy) +h_pc1 + h_pc2) ^2 + rewFunc + (1|id/run), df[df$rt_swing>0,]))
car::Anova(m10, '3')
summary(m10v <- lmer(log(rt_swing) ~ (scale(run_trial) + scale(v_max) + scale(v_entropy) +h_pc1 + vmPFC) ^2 + rewFunc + (1|id/run), df[df$rt_swing>0,]))
car::Anova(m10v, '3')
summary(m10p <- lmer(log(rt_swing) ~ (scale(run_trial) + scale(v_max) + scale(v_entropy) +h_pc1 + ParaHippL) ^2 + rewFunc + (1|id/run), df[df$rt_swing>0,]))
car::Anova(m10p, '3')
anova(m10,m10v,m10p)

summary(m11 <- lmer(log(rt_swing) ~ (scale(run_trial) + scale(v_max) + scale(v_entropy) +h_pc1 + h_pc2 + scale(rt_lag) + scale(total_earnings)) ^2 + rewFunc + (1|id/run), df[df$rt_swing>0,]))
car::Anova(m11, '3')

anova(m7,m8,m9, m10)

# plot out the relationships
ggplot(df, aes(v_entropy,rt_swing, color = h_pc1>0)) + geom_point() +  geom_smooth(method = "glm") + facet_wrap(~gamma>0)
ggplot(df, aes(v_max,rt_swing, color = h_pc2>0)) + geom_point() + geom_smooth(method = "glm") + facet_wrap(~gamma>0)

summary(mr1 <- lmer(scale(rt_csv) ~ (scale(run_trial) + scale(rt_vmax) + scale(rt_lag) + h_pc1 + h_pc2)^3 + rewFunc + (1|id/run), df))
car::Anova(mr1)

summary(mr2 <- lmer(scale(rt_csv) ~ (scale(run_trial) + scale(rt_vmax) + scale(rt_lag) + h_pc1 + scale(vmPFC))^3 + rewFunc + (1|id/run), df))
car::Anova(mr2)

summary(mr3 <- lmer(scale(rt_csv) ~ (scale(run_trial) + scale(rt_vmax) + scale(rt_lag) + h_pc1 + scale(ParaHippL))^3 + rewFunc + (1|id/run), df))
car::Anova(mr3)
anova(mr1,mr2,mr3)

#vif.lme(mr1) -- looks good, beginning to like this model
# plot vs. raw data
ggplot(df, aes(rt_vmax, rt_csv, color = h_pc1>0)) + geom_smooth(method = "glm") + facet_wrap(~run_trial>20)
ggplot(df, aes(rt_vmax, rt_csv, color = h_pc2>0)) + geom_smooth(method = "glm") + facet_wrap(~run_trial>20)

ggplot(df, aes(rt_lag, rt_csv, color = h_pc1>0)) + geom_smooth(method = "glm") + facet_wrap(~run_trial>20)

ggplot(df, aes(rt_lag, rt_csv, color = h_pc1>0, lty = h_pc2>0)) + geom_smooth(method = "glm")


ggplot(df, aes(run_trial, v_entropy, color = rewFunc)) + geom_smooth()
summary(m0 <- lmer(log(rt_swing) ~ scale(run_trial) * scale(v_max) * scale(v_entropy) * rewFunc + (1|id/run), df[df$rt_swing>0,]))

# plot out

# toward the goal of dissociating exploitative behavioral adjustment from H-driven exploration
# dm = dissociation model
summary(dm1 <- lmer(log(rt_swing) ~ scale(run_trial) +  scale(v_entropy) + omission_lag + (1|id/run), df[df$rt_swing>0,]))


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
