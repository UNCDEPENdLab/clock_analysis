library(dplyr)
library(tidyverse)
library(psych)
library(corrplot)
library(lme4)

# get H betas
setwd('~/Box Sync/skinner/projects_analyses/SCEPTIC/fMRI_paper/signals_review/MMClock_aroma_preconvolve_fse_groupfixed/sceptic-v_entropy-preconvolve_fse_groupfixed/v_entropy/')
meta <- read_csv("~/Box Sync/skinner/projects_analyses/SCEPTIC/fMRI_paper/signals_review/MMClock_aroma_preconvolve_fse_groupfixed/sceptic-v_entropy-preconvolve_fse_groupfixed/v_entropy/v_entropy_cluster_metadata.csv")
meta$label <- substr(meta$label,14,100)
meta_overall <- meta[meta$l2_contrast == 'overall' & meta$l3_contrast == 'Intercept' & meta$model == 'Intercept-Age',]
Hbetas <- read_csv("v_entropy_roi_betas.csv")
h <- as.tibble(Hbetas[Hbetas$l2_contrast == 'overall' & Hbetas$l3_contrast == 'Intercept' & Hbetas$model == 'Intercept',1:3]) %>% filter(cluster_number<11)
# head(merge(h,meta))
rois <- distinct(meta_overall[,c(5,12)])

# inspect distributions
hrois <- inner_join(h,meta)
hrois$labeled_cluster <- paste(hrois$cluster_number,hrois$label)
ggplot(hrois,aes(scale(cope_value))) + geom_histogram() + facet_wrap(~label)
h_wide <- spread(h,cluster_number,cope_value)
# head(h_wide)
# with group-fixed parameters we don't seem to need the winsorization step!!!
# h_wide <- h_wide %>% mutate_if(is.double, winsor,trim = .075)

just_rois <- h_wide[,2:11]
# winsorize to deal with beta ouliers

# non-parametric correlations to deal with outliers
clust_cor <- cor(just_rois,method = 'pearson')
# parametric correlations on winsorised betas
# clust_cor <- cor(just_rois_w,method = 'pearson')

setwd('~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/')
pdf("h_cluster_corr_fixed.pdf", width=12, height=12)
corrplot(clust_cor, cl.lim=c(-1,1),
         method = "circle", tl.cex = 1.5, type = "upper", tl.col = 'black',
         order = "hclust", diag = FALSE,
         addCoef.col="orange", addCoefasPercent = FALSE,
         p.mat = 1-clust_cor, sig.level=0.75, insig = "blank")
dev.off()

h.pca = prcomp((just_rois),scale = TRUE, center = TRUE)
h_wide$h_pc1 <- h.pca$x[,1]
h_wide$h_pc2 <- h.pca$x[,2]

h.fa = psych::fa(just_rois, nfactors=2)
fscores <- factor.scores(just_rois, h.fa)$scores
h_wide$h_f1 <- fscores[,1]
h_wide$h_f2 <- fscores[,2]
h_wide$ParaHippL <- h_wide$`9`
h_wide$vmPFC <- h_wide$`6`

# add value betas -- use v_max as theoretically interesting and unbiased by choice
setwd('~/Box Sync/skinner/projects_analyses/SCEPTIC/fMRI_paper/signals_review/MMClock_aroma_preconvolve_fse_groupfixed/sceptic-v_max-preconvolve_fse_groupfixed/v_max/')
vmeta <- read_csv("v_max_cluster_metadata.csv")
vmeta$label <- substr(vmeta$label,14,100)
vmeta_overall <- vmeta[vmeta$l2_contrast == 'overall' & vmeta$l3_contrast == 'Intercept' & vmeta$model == 'Intercept-Age',]
vbetas <- read_csv("v_max_roi_betas.csv")
v <- as.tibble(vbetas[vbetas$l2_contrast == 'overall' & vbetas$l3_contrast == 'Intercept' & vbetas$model == 'Intercept-Age',1:3]) %>% filter(cluster_number<15)
# head(merge(h,meta))
vrois_list <- distinct(vmeta_overall[,c(5,12)])

# inspect distributions
vrois <- inner_join(v,vmeta_overall)
vrois$labeled_cluster <- paste(vrois$cluster_number,vrois$label)
ggplot(vrois,aes(scale(cope_value))) + geom_histogram() + facet_wrap(~labeled_cluster)
v_wide <- spread(v,cluster_number,cope_value)

# few outliers, let's hold off on winsorizing for now
# h_wide <- h_wide %>% mutate_if(is.double, winsor,trim = .075)

vjust_rois <- v_wide[,2:ncol(v_wide)]
# winsorize to deal with beta ouliers

# non-parametric correlations to deal with outliers
vclust_cor <- cor(vjust_rois,method = 'pearson')
# parametric correlations on winsorised betas
# clust_cor <- cor(just_rois_w,method = 'pearson')

setwd('~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/')
pdf("v_cluster_corr_fixed.pdf", width=12, height=12)
corrplot(vclust_cor, cl.lim=c(-1,1),
         method = "circle", tl.cex = 1.5, type = "upper", tl.col = 'black',
         order = "hclust", diag = FALSE,
         addCoef.col="orange", addCoefasPercent = FALSE,
         p.mat = 1-vclust_cor, sig.level=0.75, insig = "blank")
dev.off()

v.pca = prcomp((vjust_rois),scale = TRUE, center = TRUE)
v_wide$v_pc1 <- h.pca$x[,1]
v_wide$v_pc2 <- h.pca$x[,2]

v.fa = psych::fa(vjust_rois, nfactors=3)
vfscores <- factor.scores(vjust_rois, v.fa)$scores
v_wide$v_f1 <- vfscores[,1]
v_wide$v_f2 <- vfscores[,2]

# a bit crazy, but let's put v and h into a single factor analysis
vh_wide <- inner_join(v_wide[,c(1,18:19)],h_wide[,c(1,14:15)], by = "feat_input_id")
vhjust_rois <- vh_wide[,2:ncol(vh_wide)]
vhclust_cor <- cor(vhjust_rois,method = 'pearson')
pdf("vh_cluster_corr_fixed.pdf", width=12, height=12)
corrplot(vhclust_cor, cl.lim=c(-1,1),
         method = "circle", tl.cex = 1.5, type = "upper", tl.col = 'black',
         order = "hclust", diag = FALSE,
         addCoef.col="orange", addCoefasPercent = FALSE,
         p.mat = 1-vhclust_cor, sig.level=0.75, insig = "blank")
dev.off()
# the first value (V+) factor is correlated with the entropy hippocampal factor (.55)

# OK, that went well -- let's add dAUC
setwd('~/Box Sync/skinner/projects_analyses/SCEPTIC/fMRI_paper/signals_review/MMClock_aroma_preconvolve_fse_groupfixed/sceptic-d_auc-preconvolve_fse_groupfixed/d_auc/')
dmeta <- read_csv("d_auc_cluster_metadata.csv")
dmeta$label <- substr(dmeta$label,14,100)
dmeta_overall <- dmeta[dmeta$l2_contrast == 'overall' & dmeta$l3_contrast == 'Intercept' & dmeta$model == 'Intercept-Age',]
dbetas <- read_csv("d_auc_roi_betas.csv")
d <- as.tibble(dbetas[dbetas$l2_contrast == 'overall' & dbetas$l3_contrast == 'Intercept' & dbetas$model == 'Intercept-Age',1:3]) %>% filter(cluster_number<15)
# head(merge(h,meta))
drois_list <- distinct(dmeta_overall[,c(5,12)])

# inspect distributions
# for some reason, the clusters do not seem to fit the map -- check with MNH (big clusters broken down?) !!
drois <- inner_join(d,dmeta_overall)
drois$labeled_cluster <- paste(drois$cluster_number,drois$label)
ggplot(drois,aes(scale(cope_value))) + geom_histogram() + facet_wrap(~labeled_cluster)
d_wide <- spread(d,cluster_number,cope_value)

# few outliers, let's hold off on winsorizing for now
# h_wide <- h_wide %>% mutate_if(is.double, winsor,trim = .075)

djust_rois <- d_wide[,2:ncol(d_wide)]
# winsorize to deal with beta ouliers

# non-parametric correlations to deal with outliers
dclust_cor <- cor(djust_rois,method = 'pearson')
# parametric correlations on winsorised betas
# clust_cor <- cor(just_rois_w,method = 'pearson')

setwd('~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/')
pdf("d_cluster_corr_fixed.pdf", width=12, height=12)
corrplot(dclust_cor, cl.lim=c(-1,1),
         method = "circle", tl.cex = 1.5, type = "upper", tl.col = 'black',
         order = "hclust", diag = FALSE,
         addCoef.col="orange", addCoefasPercent = FALSE,
         p.mat = 1-dclust_cor, sig.level=0.75, insig = "blank")
dev.off()
# looks like a single dimension

d.pca = prcomp((djust_rois),scale = TRUE, center = TRUE)
d_wide$d_pc1 <- d.pca$x[,1]
d_wide$d_pc2 <- d.pca$x[,2]

d.fa = psych::fa(djust_rois, nfactors=3)
dfscores <- factor.scores(djust_rois, d.fa)$scores
d_wide$d_f1 <- dfscores[,1]
d_wide$d_f2 <- dfscores[,2]
d_wide$d_f3 <- dfscores[,3]

# a bit crazy, but let's put v and h into a single factor analysis
dvh_wide <- inner_join(vh_wide,d_wide[,c(1,18:20)], by = "feat_input_id")
dvhjust_rois <- dvh_wide[,2:ncol(dvh_wide)]
dvhclust_cor <- cor(dvhjust_rois,method = 'pearson')
pdf("dvh_cluster_corr_fixed.pdf", width=12, height=12)
corrplot(dvhclust_cor, cl.lim=c(-1,1),
         method = "circle", tl.cex = 1.5, type = "upper", tl.col = 'black',
         order = "hclust", diag = FALSE,
         addCoef.col="orange", addCoefasPercent = FALSE,
         p.mat = 1-dvhclust_cor, sig.level=0.75, insig = "blank")
dev.off()

# add ids
map_df  <- as.tibble(read.csv("~/Box Sync/skinner/projects_analyses/SCEPTIC/fMRI_paper/signals_review/MMClock_aroma_preconvolve_fse_groupfixed/sceptic-v_entropy-preconvolve_fse_groupfixed/v_entropy/v_entropy-Intercept_design.txt", sep=""))

pc_scores <- inner_join(dvh_wide,map_df[,c(1:2,4:15)])
pc_scores$id <- pc_scores$ID


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
plot(beta_sum$v_f1,beta_sum$total_earnings)
plot(beta_sum$v_f2,beta_sum$total_earnings)

# model parameters
params <- read_csv("~/code/clock_analysis/fmri/data/mmclock_fmri_decay_factorize_selective_psequate_mfx_sceptic_global_statistics.csv")
sub_df <- inner_join(beta_sum,params)
params_beta <- sub_df[,c("v_f1","v_f2","h_f1", "h_f2","d_f1","d_f2","d_f3", "total_earnings", "LL", "alpha", "gamma", "beta")]
param_cor <- cor(params_beta,method = 'pearson')
param_cor <- corr.test(params_beta,method = 'pearson', adjust = 'none')

setwd('~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/')
pdf("dvhbeta_param_corr_fixed.pdf", width=12, height=12)
corrplot(param_cor$r, cl.lim=c(-1,1),
         method = "circle", tl.cex = 1.5, type = "upper", tl.col = 'black',
         order = "hclust", diag = FALSE,
         addCoef.col="orange", addCoefasPercent = FALSE,
         p.mat = param_cor$p, sig.level=0.05, insig = "blank")
dev.off()

# merge into trial-level data
df <- inner_join(trial_df,sub_df)
df$rewFunc <- relevel(as.factor(df$rewFunc),ref = "CEV")

#####
# these simple RT swing prediction models only reveal that both networks catalyze convergence regardless of condition
summary(m01 <- lmer(log(rt_swing) ~ scale(run_trial) * rewFunc + (1|id/run), df[df$rt_swing>0,]))
summary(m02 <- lmer(log(rt_swing) ~ (scale(run_trial) + rewFunc + h_f1 + h_f2)^3 + (1|id/run), df[df$rt_swing>0,]))
car::Anova(m02,'3')
anova(m01,m02)
# over-engineered, does not add much
# summary(m03 <- lmer(rt_csv ~ (scale(rt_lag) + scale(run_trial) + rewFunc + h_f1 + h_f2)^4 + (1|id/run), df[df$rt_swing>0,]))


# RT swings analyses: "exploration"
summary(m1 <- lmer(log(rt_swing) ~ scale(run_trial) * scale(v_max) + rewFunc + (1|id/run), df[df$rt_swing>0,]))
# keeps flipping between log and untransformed...
summary(m2 <- lmer(log(rt_swing) ~ (scale(run_trial) + scale(v_max) + h_f1 + h_f2) ^2 + rewFunc + (1|id/run), df[df$rt_swing>0,]))
car::Anova(m2, '3')
summary(m3 <- lmer(log(rt_swing) ~ (scale(run_trial) + scale(v_max) + v_f1 + v_f2) ^2 + rewFunc + (1|id/run), df[df$rt_swing>0,]))
car::Anova(m3, '3')
summary(m4 <- lmer(log(rt_swing) ~ (scale(run_trial) + scale(v_max) + h_f1 + h_f2 + v_f1 + v_f2) ^2 + rewFunc + (1|id/run), df[df$rt_swing>0,]))
car::Anova(m4, '3')
summary(m5 <- lmer(log(rt_swing) ~ (scale(run_trial) + scale(v_max) + h_f1 + h_f2 + v_f1 + v_f2 + d_f1 + d_f2 + d_f3) ^2 + rewFunc + (1|id/run), df[df$rt_swing>0,]))

anova(m1,m2,m3,m4,m5)
# value clusters explain more than entropy (diffAIC = 31), but each set explains unique variance; d_auc adds little

###########
# EV analyses: "exploitation" -- run removed from RE to avoid singular fit
# Differences emerge in IEV (explains why total earnings are not informative)
summary(ev1 <- lmer(ev ~ scale(run_trial) * rewFunc + (1|id), df))
summary(ev2 <- lmer(ev ~ (scale(run_trial) + rewFunc + h_f1) ^3 + (1|id), df))
summary(ev3 <- lmer(ev ~ (scale(run_trial) + rewFunc + I(-h_f2)) ^3 + (1|id), df))
summary(ev4 <- lmer(ev ~ (scale(run_trial) + rewFunc + h_f1 + I(-h_f2)) ^3 + (1|id), df))
anova(ev1,ev2,ev3,ev4)

pdf('ev_by_condition_and_brain_betas.pdf', height = 8, width = 8)
ggplot(df, aes(run_trial, ev, color = h_f2<0, lty = h_f1 >0)) + geom_smooth(method = "loess") + facet_wrap (~rewFunc)
dev.off()

# car::Anova(m5, '3')
# summary(m6 <- lmer(ev ~ (scale(run_trial) + scale(v_max) + h_f2) ^2 + (1|id/run), df[df$rt_swing>0,]))
# car::Anova(m6, '3')
#
# anova(m4,m5,m6)

# model-free look at H timecourses -- there is a bigger effect in learnable > unlearnable (esp. CEVR)
# need both networks to transition
ggplot(df, aes(run_trial, v_entropy, color = h_f1 >0, lty = h_f2 <0)) + geom_smooth(method = "loess") + facet_wrap(~rewFunc)
ggplot(df, aes(run_trial, log(rt_swing), color = h_f1 >0, lty = h_f2 <0)) + geom_smooth(method = "loess") + facet_wrap(~rewFunc)

# value clusters -- similar, but weaker
ggplot(df, aes(run_trial, log(rt_swing), color = v_f1 >0, lty = v_f2 <0)) + geom_smooth(method = "loess") + facet_wrap(~rewFunc)


# coupling between entropy and RT swings -- betas don't moderate effects of entropy.  Rather, they influence entropy as it unfolds.
# summary(m7 <- lmer(log(rt_swing) ~ (scale(run_trial) + scale(v_max) + scale(v_entropy))^2 + rewFunc + (1|id/run), df[df$rt_swing>0,]))
# summary(m8 <- lmer(log(rt_swing) ~ (scale(run_trial) + scale(v_max) + scale(v_entropy) + h_f1) ^2 + rewFunc + (1|id/run), df[df$rt_swing>0,]))
# car::Anova(m8, '3')
# summary(m9 <- lmer(log(rt_swing) ~ (scale(run_trial) + scale(v_max) + scale(v_entropy) + h_f2) ^2 + rewFunc + (1|id/run), df[df$rt_swing>0,]))
# car::Anova(m9, '3')
# summary(m10 <- lmer(log(rt_swing) ~ (scale(run_trial) + scale(v_max) + scale(v_entropy) +h_f1 + h_f2) ^2 + rewFunc + (1|id/run), df[df$rt_swing>0,]))
# car::Anova(m10, '3')
# # a bunch of relatively weak interactions, set a side for now:
# summary(m11 <- lmer(log(rt_swing) ~ (scale(run_trial) + scale(v_max) + scale(v_entropy) +h_f1 + h_f2 + v_f1 + v_f2) ^3 + rewFunc + (1|id/run), df[df$rt_swing>0,]))
# car::Anova(m11, '3')
# anova(m7,m8,m9, m10)

# plot out the relationships
ggplot(df, aes(v_entropy,rt_swing, color = h_f1>0)) + geom_point() +  geom_smooth(method = "glm") + facet_wrap(~gamma>0)
ggplot(df, aes(v_max,rt_swing, color = h_f2>0)) + geom_point() + geom_smooth(method = "glm") + facet_wrap(~gamma>0)

# just entropy
summary(rt1 <- lmer(scale(rt_csv) ~ (scale(run_trial) + scale(rt_vmax) + scale(rt_lag) + h_f1 + h_f2)^3 + rewFunc + (1|id/run), df))
car::Anova(rt1)
ggplot(df,aes(rt_vmax,rt_csv, color = h_f1>0)) + facet_wrap(~run_trial>20) + geom_smooth(method = "glm")
ggplot(df,aes(rt_vmax,rt_csv, color = h_f2<0)) + facet_wrap(~run_trial>20) + geom_smooth(method = "glm")
ggplot(df,aes(rt_vmax,rt_csv, color = h_f2<0, lty = h_f1>0)) + geom_smooth(method = "glm")

# + value
summary(rt2 <- lmer(scale(rt_csv) ~ (scale(run_trial) + scale(rt_vmax) + scale(rt_lag) + v_f1 + v_f2)^3 + rewFunc + (1|id/run), df))
car::Anova(rt2)

# + d_auc
summary(rt3 <- lmer(scale(rt_csv) ~ (scale(run_trial) + scale(rt_vmax) + scale(rt_lag) + d_f1 + d_f2 + d_f3)^3 + rewFunc + (1|id/run), df))
car::Anova(rt3)

summary(rt4 <- lmer(scale(rt_csv) ~ (scale(run_trial) + scale(rt_vmax) + scale(rt_lag) + h_f1 + h_f2 + v_f1 +v_f2 + d_f1 + d_f2 + d_f3)^3 + rewFunc + (1|id/run), df))
anova(rt2,rt1,rt3, rt4) # H predicts best, but the combination improves the prediction further (by >300 AIC points)


# vmPFC/hippocampus no longer explain more variance than factors/PCs
# summary(mr2 <- lmer(scale(rt_csv) ~ (scale(run_trial) + scale(rt_vmax) + scale(rt_lag) + h_f1 + scale(vmPFC))^3 + rewFunc + (1|id/run), df))
# car::Anova(mr2)
#
# summary(mr3 <- lmer(scale(rt_csv) ~ (scale(run_trial) + scale(rt_vmax) + scale(rt_lag) + h_f1 + scale(ParaHippL))^3 + rewFunc + (1|id/run), df))
# car::Anova(mr3)
# anova(mr1,mr2,mr3)

#vif.lme(mr1) -- looks good, beginning to like this model
# plot vs. raw data
ggplot(df, aes(rt_vmax, rt_csv, color = h_f1>0)) + geom_smooth(method = "glm") + facet_wrap(~run_trial>20)
ggplot(df, aes(rt_vmax, rt_csv, color = h_f2>0)) + geom_smooth(method = "glm") + facet_wrap(~run_trial>20)

ggplot(df, aes(rt_lag, rt_csv, color = h_f1>0)) + geom_smooth(method = "glm") + facet_wrap(~run_trial>20)

ggplot(df, aes(rt_lag, rt_csv, color = h_f1>0, lty = h_f2>0)) + geom_smooth(method = "glm")+ facet_wrap(~run_trial>20)
ggplot(df, aes(rt_vmax, rt_csv, color = h_f1>0, lty = h_f2>0)) + geom_smooth(method = "glm")+ facet_wrap(~run_trial>20)

pdf("h_timecourse_brain_fixed.pdf", width = 8, height = 8)
ggplot(df, aes(run_trial, v_entropy, color = h_f1>0, lty = h_f2>0)) + geom_smooth(method = "loess") + facet_wrap(~gamma>0)
dev.off()


ggplot(df, aes(run_trial, v_entropy, color = rewFunc)) + geom_smooth()
summary(m0 <- lmer(log(rt_swing) ~ scale(run_trial) * scale(v_max) * scale(v_entropy) * rewFunc + (1|id/run), df[df$rt_swing>0,]))

# plot out

# toward the goal of dissociating exploitative behavioral adjustment from H-driven exploration
# dm = dissociation model
# unfortunately does not work all that well due to collinearity
summary(dm1 <- lmer(log(rt_swing) ~ (scale(run_trial) + omission_lag + scale(v_max))^3 + scale(rt_lag) + (1|id/run), df[df$rt_swing>0,]))
ggplot(na.omit(df), aes(v_entropy, rt_swing, color = v_max>25, lty = omission_lag)) + geom_smooth(method = "gam")

summary(dm2 <- lmer(rt_csv ~ (scale(rt_lag) + omission_lag + scale(rt_vmax))^3 + rewFunc + (1|id/run), df[df$rt_swing>0,]))

summary(dm3 <- lmer(rt_csv ~ (scale(rt_lag) + omission_lag + scale(rt_vmax) + h_f1 + h_f2)^3 + rewFunc + (1|id/run), df[df$rt_swing>0,]))

# more focused test of dislodgment
summary(dm3a <- lmer(rt_csv ~ (scale(rt_lag) + omission_lag + scale(rt_vmax) + h_f1)^3 + rewFunc + (1|id/run), df[df$rt_swing>0,]))
summary(dm3b <- lmer(rt_csv ~ (scale(rt_lag) + omission_lag + scale(rt_vmax) + h_f2)^3 + rewFunc + (1|id/run), df[df$rt_swing>0,]))

# plot 3-way interactions
ggplot(df, aes(rt_lag, rt_csv, color = omission_lag, lty = h_f1>0)) + geom_smooth(method = "glm")
ggplot(df, aes(rt_lag, rt_csv, color = omission_lag, lty = h_f2<0)) + geom_smooth(method = "glm")
# the interesting one:
ggplot(na.omit(df), aes(rt_vmax, rt_csv, color = omission_lag, lty = h_f1>0)) + geom_smooth(method = "glm")
ggplot(na.omit(df), aes(rt_vmax, rt_csv, color = omission_lag, lty = h_f2<0)) + geom_smooth(method = "glm")

anova(dm1,dm2,dm3)

summary(dm2 <- lmer(log(rt_swing) ~ (scale(run_trial) +  scale(v_entropy) + omission_lag + scale(v_max) + h_f1 + h_f2)^3 + (1|id/run), df[df$rt_swing>0,]))
anova(dm1,dm2)


# obtain within-subject v_max and entropy 

df <- df %>% group_by(id,run) %>% mutate(v_max_wi = scale(v_max),
  v_entropy_wi = scale(v_entropy)
)
ggplot(df, aes(run_trial, v_entropy_wi, color = h_f1>0, lty = h_f2>0)) + geom_smooth(method = "loess") + facet_wrap(~gamma>0)
ggplot(df, aes( v_entropy_wi,log(rt_swing), color = h_f1>0, lty = h_f2>0)) + geom_smooth(method = "gam") + facet_wrap(~run_trial > 10)
ggplot(df, aes( v_max_wi,log(rt_swing), color = h_f1>0, lty = h_f2>0)) + geom_smooth(method = "gam") + facet_wrap(~run_trial > 10)

summary(w1 <- lmer(log(rt_swing) ~ (scale(run_trial) + v_max_wi + v_entropy_wi)^2 + rewFunc + (1|id/run), df[df$rt_swing>0,]))
summary(w2 <- lmer(log(rt_swing) ~ (scale(run_trial) + v_max_wi + v_entropy_wi + h_f1 + h_f2)^3 + rewFunc + (1|id/run), df[df$rt_swing>0,]))
car::Anova(w2)

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
