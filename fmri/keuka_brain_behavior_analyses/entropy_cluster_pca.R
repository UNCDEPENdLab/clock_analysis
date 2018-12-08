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
h <- as.tibble(Hbetas[Hbetas$l2_contrast == 'overall' & Hbetas$l3_contrast == 'Intercept' & Hbetas$model == 'Intercept-Age',1:3]) %>% filter(cluster_number<11)
# head(merge(h,meta))
rois <- distinct(meta_overall[,c(5,12)])

# inspect distributions
hrois <- inner_join(h,meta_overall)
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
         addCoef.col="black", addCoefasPercent = FALSE,
         p.mat = 1-clust_cor, sig.level=0.75, insig = "blank")
dev.off()

h.pca = prcomp((just_rois),scale = TRUE, center = TRUE)
h_wide$h_pc1 <- h.pca$x[,1]
h_wide$h_pc2 <- h.pca$x[,2]

h.fa = psych::fa(just_rois, nfactors=2)
fscores <- factor.scores(just_rois, h.fa)$scores
h_wide$h_f1_fp <- fscores[,1]
h_wide$h_f2_neg_paralimb <- fscores[,2]
h_wide$ParaHippL <- h_wide$`9`
h_wide$vmPFC <- h_wide$`6`

# add value betas -- use v_max as theoretically interesting and unbiased by choice
setwd('~/Box Sync/skinner/projects_analyses/SCEPTIC/fMRI_paper/signals_review/MMClock_aroma_preconvolve_fse_groupfixed/sceptic-v_max-preconvolve_fse_groupfixed/v_max/')
vmeta <- read_csv("v_max_cluster_metadata.csv")
vmeta$label <- substr(vmeta$label,14,100)
vmeta_overall <- vmeta[vmeta$l2_contrast == 'overall' & vmeta$l3_contrast == 'Intercept' & vmeta$model == 'Intercept-Age',]
vbetas <- read_csv("v_max_roi_betas.csv")
v <- as.tibble(vbetas[vbetas$l2_contrast == 'overall' & vbetas$l3_contrast == 'Intercept' & vbetas$model == 'Intercept-Age',1:3]) %>% filter(cluster_number<15)
vrois_list <- distinct(vmeta_overall[,c(5,12)])


# inspect distributions
vrois <- inner_join(v,vmeta_overall)
vrois$labeled_cluster <- paste(vrois$cluster_number,vrois$label)
ggplot(vrois,aes(scale(cope_value))) + geom_histogram() + facet_wrap(~labeled_cluster)
v_wide <- spread(v,cluster_number,cope_value)
# try v winsorization given the left outliers
# v_wide <- v_wide %>% mutate_if(is.double, winsor,trim = .075)

# few outliers, let's hold off on winsorizing for now
# h_wide <- h_wide %>% mutate_if(is.double, winsor,trim = .075)

# I am worried about these scores: a lot of the V+ f2 variance comes from fusiform gyri (4,11)
# let's redo using only fronto-striatal clusters (vmPFC - 10, BG - 2, vlPFC - 9,14) removing fusiform (4,11)
# and cerebellar (8) clusters
v_wide <- subset(v_wide, select = -c(`4`, `11`,`8`))
v_wide <- v_wide %>% mutate_if(is.double, winsor,trim = .075)


vjust_rois <- v_wide[,2:ncol(v_wide)]
# winsorize to deal with beta ouliers

# non-parametric correlations to deal with outliers
# parametric correlations on winsorised betas
# clust_cor <- cor(just_rois_w,method = 'pearson')
v.pca = prcomp((vjust_rois),scale = TRUE, center = TRUE)
v_wide$v_pc1 <- v.pca$x[,1]
v_wide$v_pc2 <- v.pca$x[,2]

v.fa = psych::fa(vjust_rois, nfactors=3)
vfscores <- factor.scores(vjust_rois, v.fa)$scores
v_wide$v_f1_neg_cog <- vfscores[,1]
v_wide$v_f2_paralimb <- vfscores[,2]

vclust_cor <- corr.test(v_wide %>% select_if(is.double),method = 'pearson', adjust = 'none')

setwd('~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/')
pdf("v_cluster_corr_fixed.pdf", width=12, height=12)
corrplot(vclust_cor$r, cl.lim=c(-1,1),
         method = "circle", tl.cex = 1.5, type = "upper", tl.col = 'black',
         order = "hclust", diag = FALSE,
         addCoef.col="black", addCoefasPercent = FALSE,
         p.mat = vclust_cor$p, sig.level=0.05, insig = "blank")
dev.off()

# a bit crazy, but let's put v and h into a single factor analysis
vh_wide <- inner_join(subset(v_wide, select = c("feat_input_id","v_f1_neg_cog","v_f2_paralimb")),
                      subset(h_wide, select = c("feat_input_id","h_f1_fp","h_f2_neg_paralimb")), by = "feat_input_id")
vhjust_rois <- vh_wide %>% select_if(is.double)
vhclust_cor <- corr.test(vhjust_rois,method = 'pearson', adjust = 'none')
pdf("vh_cluster_corr_fixed.pdf", width=12, height=12)
corrplot(vhclust_cor$r, cl.lim=c(-1,1),
         method = "circle", tl.cex = 1.5, type = "upper", tl.col = 'black',
         order = "hclust", diag = FALSE,
         addCoef.col="black", addCoefasPercent = FALSE,
         p.mat = vhclust_cor$p, sig.level=0.05, insig = "blank")
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
         addCoef.col="black", addCoefasPercent = FALSE,
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
dvhclust_cor <- corr.test(dvhjust_rois,method = 'pearson', adjust = 'none')
pdf("dvh_cluster_corr_fixed.pdf", width=12, height=12)
corrplot(dvhclust_cor$r, cl.lim=c(-1,1),
         method = "circle", tl.cex = 1.5, type = "upper", tl.col = 'black',
         order = "hclust", diag = FALSE,
         addCoef.col="black", addCoefasPercent = FALSE,
         p.mat = dvhclust_cor$p, sig.level=0.05, insig = "blank")
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
plot(beta_sum$h_f1_fp,beta_sum$total_earnings)
plot(beta_sum$h_f2_neg_paralimb,beta_sum$total_earnings)
plot(beta_sum$v_f1_neg_cog,beta_sum$total_earnings)
plot(beta_sum$v_f2_paralimb,beta_sum$total_earnings)

# model parameters
params <- read_csv("~/code/clock_analysis/fmri/data/mmclock_fmri_decay_factorize_selective_psequate_mfx_sceptic_global_statistics.csv")
sub_df <- inner_join(beta_sum,params)
params_beta <- sub_df[,c("v_f1_neg_cog","v_f2_paralimb","h_f1_fp", "h_f2_neg_paralimb","d_f1","d_f2","d_f3", "total_earnings", "LL", "alpha", "gamma", "beta")]
param_cor <- corr.test(params_beta,method = 'pearson', adjust = 'none')

setwd('~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/')
pdf("dvhbeta_param_corr_fixed.pdf", width=12, height=12)
corrplot(param_cor$r, cl.lim=c(-1,1),
         method = "circle", tl.cex = 1.5, type = "upper", tl.col = 'black',
         order = "hclust", diag = FALSE,
         addCoef.col="black", addCoefasPercent = FALSE,
         p.mat = param_cor$p, sig.level=0.05, insig = "blank")
dev.off()

# merge into trial-level data
df <- inner_join(trial_df,sub_df)
df$rewFunc <- relevel(as.factor(df$rewFunc),ref = "CEV")

# dichotomize betas for plotting
df$h_fp <- 'low'
df$h_fp[df$h_f1_fp>0] <- 'high'
df$low_h_paralimbic <- 'low'
df$low_h_paralimbic[df$h_f2_neg_paralimb<0] <- 'high'

df$v_paralimbic <- 'low'
df$v_paralimbic[df$v_f2_paralimb>0] <- 'high'
#NB: v_f2_paralimb is mostly vlPFC driven now...
df$low_v_fp_acc_vlpfc <- 'low'
df$low_v_fp_acc_vlpfc[df$v_f1_neg_cog<0] <- 'high'

df$learning_epoch <- 'trials 1-10'
df$learning_epoch[df$run_trial>10] <- 'trials 11-50'




# obtain within-subject v_max and entropy: correlated at -.37

df <- df %>% group_by(id,run) %>% mutate(v_max_wi = scale(v_max),
                                         v_entropy_wi = scale(v_entropy)
)
# circular, but just check to what extent each area conforms to SCEPTIC-SM: looks like there is an interaction, both need to be involved again
ggplot(df, aes(run_trial, v_entropy_wi, color = low_h_paralimbic, lty = h_fp)) + geom_smooth(method = "loess") #+ facet_wrap(~gamma>0)
# I think this means again that H estimates are more precise for high-paralimbic people, esp. late in learning
ggplot(df, aes( v_entropy_wi,log(rt_swing), color = low_h_paralimbic, lty = h_fp)) + geom_smooth(method = "gam") + facet_wrap(~learning_epoch)
ggplot(df, aes( v_max_wi,log(rt_swing), color = low_h_paralimbic, lty = h_fp)) + geom_smooth(method = "gam") + facet_wrap(~run_trial > 10)

#####
# these simple RT swing prediction models only reveal that both networks catalyze convergence regardless of condition
summary(m01 <- lmer(log(rt_swing) ~ scale(-1/run_trial) * rewFunc + (1|id/run), df[df$rt_swing>0,]))
summary(m02 <- lmer(log(rt_swing) ~ (scale(-1/run_trial) + rewFunc + h_f1_fp + h_f2_neg_paralimb)^2 + (1|id/run), df[df$rt_swing>0,]))
car::Anova(m02,'3')
anova(m01,m02)
# over-engineered, does not add much
# summary(m03 <- lmer(rt_csv ~ (scale(rt_lag) + scale(-1/run_trial) + rewFunc + h_f1_fp + h_f2_neg_paralimb)^4 + (1|id/run), df[df$rt_swing>0,]))


# RT swings analyses: "exploration"
summary(m1 <- lmer(log(rt_swing) ~ scale(-1/run_trial) * scale(v_max) + rewFunc + (1|id/run), df[df$rt_swing>0,]))
# keeps flipping between log and untransformed...
summary(m2 <- lmer(log(rt_swing) ~ (scale(-1/run_trial) + scale(v_max) + h_f1_fp + h_f2_neg_paralimb) ^3 + rewFunc + (1|id/run), df[df$rt_swing>0,]))
car::Anova(m2, '3')
summary(m3 <- lmer(log(rt_swing) ~ (scale(-1/run_trial) + scale(v_max) + v_f1_neg_cog + v_f2_paralimb) ^3 + rewFunc + (1|id/run), df[df$rt_swing>0,]))
car::Anova(m3, '3')
summary(m4 <- lmer(log(rt_swing) ~ (scale(-1/run_trial) + scale(v_max) + h_f1_fp + h_f2_neg_paralimb + v_f1_neg_cog + v_f2_paralimb) ^3 + rewFunc + (1|id/run), df[df$rt_swing>0,]))
anova(m4)
car::Anova(m4, '3')
# adding d introduces multi-collinearity
summary(m5 <- lmer(log(rt_swing) ~ (scale(-1/run_trial) + scale(v_max) + h_f1_fp + h_f2_neg_paralimb + v_f1_neg_cog + v_f2_paralimb + d_f1 + d_f2 + d_f3) ^3 + rewFunc + (1|id/run), df[df$rt_swing>0,]))

anova(m1,m2,m3,m4,m5)
# value clusters explain more than entropy (diffAIC = 31), but each set explains unique variance; d_auc adds little

###########
# EV analyses: "exploitation" -- run removed from RE to avoid singular fit
# Differences emerge in IEV (explains why total earnings are not informative)
summary(ev1 <- lmer(ev ~ scale(-1/run_trial) * rewFunc + (1|id), df))
summary(ev2 <- lmer(ev ~ (I(-scale(1/run_trial)) + rewFunc + h_f1_fp) ^3 + (1|id), df))
summary(ev3 <- lmer(ev ~ (I(-scale(1/run_trial)) + rewFunc + I(-h_f2_neg_paralimb)) ^3 + (1|id), df))
summary(ev4 <- lmer(ev ~ (I(-scale(1/run_trial)) + rewFunc + h_f1_fp + I(-h_f2_neg_paralimb)) ^4 + (1|id), df))
anova(ev1,ev2,ev3,ev4)

pdf('ev_by_condition_and_h_betas.pdf', height = 8, width = 8)
ggplot(df, aes(run_trial, ev, color = low_h_paralimbic, lty = h_fp)) + geom_smooth(method = "loess") + facet_wrap (~rewFunc)
dev.off()
pdf('ev_by_condition_and_v_betas.pdf', height = 8, width = 8)
ggplot(df, aes(run_trial, ev, color = v_paralimbic, lty = low_v_fp_acc_vlpfc)) + geom_smooth(method = "loess") + facet_wrap (~rewFunc)
dev.off()
#Hmmm...  this is very puzzling.  The v_max positive factor (vmPFC et al.) looks like a maladaptive response, even
# after removing cerebellum and fusiform!
# sanity check
# ggplot(df, aes(v_paralimbic,v_f2_paralimb)) + geom_boxplot()
# ggplot(df, aes(low_v_fp_acc_vlpfc,v_f1_neg_cog)) + geom_boxplot()
# check H timecourses
pdf('h_timecourse_by_condition_and_v_betas.pdf', height = 8, width = 8)
ggplot(df, aes(run_trial, v_entropy, color = v_paralimbic, lty = low_v_fp_acc_vlpfc)) + geom_smooth(method = "loess") + facet_wrap (~rewFunc)
dev.off()
pdf('h_timecourse_by_condition_and_h_betas.pdf', height = 8, width = 8)
ggplot(df, aes(run_trial, v_entropy, color = low_h_paralimbic, lty = h_fp)) + geom_smooth(method = "loess") + facet_wrap (~rewFunc)
dev.off()


pdf('rt_swings_by_condition_and_h_betas.pdf', height = 8, width = 8)
ggplot(df, aes(run_trial, log(rt_swing), color = low_h_paralimbic, lty = h_fp)) + geom_smooth(method = "loess") + facet_wrap(~rewFunc)
dev.off()

pdf('rt_swings_by_condition_and_v_betas.pdf', height = 8, width = 8)
ggplot(df, aes(run_trial, log(rt_swing), color = v_paralimbic, lty = low_v_fp_acc_vlpfc)) + geom_smooth(method = "loess") + facet_wrap(~rewFunc)
dev.off()


# coupling between entropy and RT swings -- betas don't moderate effects of entropy.  Rather, they influence entropy as it unfolds.
# summary(m7 <- lmer(log(rt_swing) ~ (scale(-1/run_trial) + scale(v_max) + scale(v_entropy))^2 + rewFunc + (1|id/run), df[df$rt_swing>0,]))
# summary(m8 <- lmer(log(rt_swing) ~ (scale(-1/run_trial) + scale(v_max) + scale(v_entropy) + h_f1_fp) ^2 + rewFunc + (1|id/run), df[df$rt_swing>0,]))
# car::Anova(m8, '3')
# summary(m9 <- lmer(log(rt_swing) ~ (scale(-1/run_trial) + scale(v_max) + scale(v_entropy) + h_f2_neg_paralimb) ^2 + rewFunc + (1|id/run), df[df$rt_swing>0,]))
# car::Anova(m9, '3')
# summary(m10 <- lmer(log(rt_swing) ~ (scale(-1/run_trial) + scale(v_max) + scale(v_entropy) +h_f1_fp + h_f2_neg_paralimb) ^2 + rewFunc + (1|id/run), df[df$rt_swing>0,]))
# car::Anova(m10, '3')
# # a bunch of relatively weak interactions, set a side for now:
# summary(m11 <- lmer(log(rt_swing) ~ (scale(-1/run_trial) + scale(v_max) + scale(v_entropy) +h_f1_fp + h_f2_neg_paralimb + v_f1_neg_cog + v_f2_paralimb) ^3 + rewFunc + (1|id/run), df[df$rt_swing>0,]))
# car::Anova(m11, '3')
# anova(m7,m8,m9, m10)

# plot out the relationships
ggplot(df, aes(v_entropy,rt_swing, color = h_f1_fp>0)) + geom_smooth(method = "gam")
ggplot(df, aes(v_max,rt_swing,  color = low_h_paralimbic, lty = h_fp)) + geom_smooth(method = "gam")

# just entropy -- VIFs ok
summary(rt0 <- lmer(scale(rt_csv) ~ (scale(-1/run_trial) + scale(rt_vmax) + scale(rt_lag))^3 + rewFunc + (1|id/run), df))

summary(rth <- lmer(scale(rt_csv) ~ (scale(-1/run_trial) + scale(rt_vmax) + scale(rt_lag) + h_f1_fp + h_f2_neg_paralimb)^3 + rewFunc + (1|id/run), df))
car::Anova(rt1)
ggplot(df,aes(rt_vmax,rt_csv, color = low_h_paralimbic, lty = h_fp)) + facet_wrap(~run_trial>20) + geom_smooth(method = "glm")
ggplot(df,aes(rt_lag,rt_csv, color = low_h_paralimbic, lty = h_fp)) + facet_wrap(~run_trial>20) + geom_smooth(method = "glm")

# just value -- VIFs ok
summary(rtv <- lmer(scale(rt_csv) ~ (scale(-1/run_trial) + scale(rt_vmax) + scale(rt_lag) + v_f1_neg_cog + v_f2_paralimb)^3 + rewFunc + (1|id/run), df))
# + d_auc
summary(rtd <- lmer(scale(rt_csv) ~ (scale(-1/run_trial) + scale(rt_vmax) + scale(rt_lag) + d_f1 + d_f2 + d_f3)^3 + rewFunc + (1|id/run), df))
# only H and V
summary(rtvh <- lmer(scale(rt_csv) ~ (scale(-1/run_trial) + scale(rt_vmax) + scale(rt_lag) + h_f1_fp + h_f2_neg_paralimb + v_f1_neg_cog +v_f2_paralimb)^3 + rewFunc + (1|id/run), df))
# all
summary(rtvhd <- lmer(scale(rt_csv) ~ (scale(-1/run_trial) + scale(rt_vmax) + scale(rt_lag) + h_f1_fp + h_f2_neg_paralimb  + d_f1 + d_f2 + d_f3)^3 + rewFunc + (1|id/run), df))
anova(rt0,rtv,rtd,rth,rtvh,rtvhd) # H predicts best, each addition (V, D_AUC) improves the prediction further (by >300 AIC points)
# but careful with cluster*cluster interactions -- high VIFs!


# the model to rule them all?
summary(rtwh <- lmer(scale(rt_csv) ~ (scale(run_trial) + scale(rt_vmax) + scale(rt_lag) + v_entropy_wi + v_max_wi + h_f1_fp + h_f2_neg_paralimb)^3 + rewFunc + (1|id/run), df))

# plot vs. raw data
ggplot(df, aes(rt_vmax, rt_csv, color = h_f1_fp>0)) + geom_smooth(method = "glm") + facet_wrap(~run_trial>20)
ggplot(df, aes(rt_vmax, rt_csv, color = h_f2_neg_paralimb>0)) + geom_smooth(method = "glm") + facet_wrap(~run_trial>20)

ggplot(df, aes(rt_lag, rt_csv, color = h_f1_fp>0)) + geom_smooth(method = "glm") + facet_wrap(~run_trial>20)

ggplot(df, aes(rt_lag, rt_csv, color = low_h_paralimbic, lty = h_fp)) + geom_smooth(method = "glm")+ facet_wrap(~run_trial>20)
ggplot(df, aes(rt_vmax, rt_csv, color = low_h_paralimbic, lty = h_fp)) + geom_smooth(method = "glm")+ facet_wrap(~run_trial>20)

pdf("h_timecourse_brain_fixed.pdf", width = 8, height = 8)
ggplot(df, aes(run_trial, v_entropy,color = low_h_paralimbic, lty = h_fp)) + geom_smooth(method = "loess") + facet_wrap(~gamma>0)
dev.off()


ggplot(df, aes(run_trial, v_entropy, color = rewFunc)) + geom_smooth()
summary(m0 <- lmer(log(rt_swing) ~ scale(-1/run_trial) * scale(v_max) * scale(v_entropy) * rewFunc + (1|id/run), df[df$rt_swing>0,]))

# plot out

# toward the goal of dissociating exploitative behavioral adjustment from H-driven exploration
# dm = dissociation model
# unfortunately does not work all that well due to collinearity
summary(dm1 <- lmer(log(rt_swing) ~ (scale(-1/run_trial) + omission_lag + scale(v_max))^3 + scale(rt_lag) + (1|id/run), df[df$rt_swing>0,]))
ggplot(na.omit(df), aes(v_entropy, rt_swing, color = v_max>25, lty = omission_lag)) + geom_smooth(method = "gam")

summary(dm2 <- lmer(rt_csv ~ (scale(rt_lag) + omission_lag + scale(rt_vmax))^3 + rewFunc + (1|id/run), df[df$rt_swing>0,]))

summary(dm3 <- lmer(rt_csv ~ (scale(rt_lag) + omission_lag + scale(rt_vmax) + h_f1_fp + h_f2_neg_paralimb)^3 + rewFunc + (1|id/run), df[df$rt_swing>0,]))

# more focused test of dislodgment
summary(dm3a <- lmer(rt_csv ~ (scale(rt_lag) + omission_lag + scale(rt_vmax) + h_f1_fp)^3 + rewFunc + (1|id/run), df[df$rt_swing>0,]))
summary(dm3b <- lmer(rt_csv ~ (scale(rt_lag) + omission_lag + scale(rt_vmax) + h_f2_neg_paralimb)^3 + rewFunc + (1|id/run), df[df$rt_swing>0,]))

# plot 3-way interactions
ggplot(df, aes(rt_lag, rt_csv, color = omission_lag, lty = h_f1_fp>0)) + geom_smooth(method = "glm")
ggplot(df, aes(rt_lag, rt_csv, color = omission_lag, lty = h_f2_neg_paralimb<0)) + geom_smooth(method = "glm")
# the interesting one:
ggplot(na.omit(df), aes(rt_vmax, rt_csv, color = omission_lag, lty = h_f1_fp>0)) + geom_smooth(method = "glm")
ggplot(na.omit(df), aes(rt_vmax, rt_csv, color = omission_lag, lty = h_f2_neg_paralimb<0)) + geom_smooth(method = "glm")

anova(dm1,dm2,dm3)

summary(dm2 <- lmer(log(rt_swing) ~ (scale(-1/run_trial) +  scale(v_entropy) + omission_lag + scale(v_max) + h_f1_fp + h_f2_neg_paralimb)^3 + (1|id/run), df[df$rt_swing>0,]))
anova(dm1,dm2)

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
