library(tidyverse)
library(psych)
library(corrplot)
library(lme4)

plot_dir = "~/OneDrive/collected_letters/papers/sceptic_fmri/dan/plots/rt_decode"
setwd(plot_dir)

# read in entropy change data
fname = "rt_encode_output_streams_mixed_by_entropy_change_ranef.RDS"
ddf <- as.data.frame(readRDS(fname)) %>% select(-model_name, -rhs, -outcome, -std.error, -effect, -group)
# spread by stream and time
wdf <- ddf %>% group_by(level) %>% filter(term=="v_entropy_wi_change", evt_time>0 & evt_time<5) %>% 
  select(-term) %>%
  pivot_wider(names_from = c(stream, evt_time, side), values_from = estimate) %>% ungroup()

# correlation plots

ranefs <- distinct(wdf) %>% select(-level)
ranef_cor <- corr.test(na.omit(ranefs),method = 'pearson', adjust = 'none')
# parametric correlations on winsorised betas
# clust_cor <- cor(just_rois_w,method = 'pearson')

pdf("entropy_change_ranef_corr.pdf", width=60, height=60)
corrplot(ranef_cor$r, cl.lim=c(-1,1), 
         method = "circle", tl.cex = 1.5, type = "upper", tl.col = 'black',
         order = "AOE", diag = FALSE,
         addCoef.col="black", addCoefasPercent = FALSE,
         p.mat = ranef_cor$p, sig.level=0.05, insig = "blank")
dev.off()

nr <- nfactors(ranef_cor$r, n=10, rotate = "oblimin", diagonal = FALSE,fm = "pa", n.obs = 70, SMC = FALSE)
echange.fa = psych::fa(ranefs, nfactors=2, rotate = "oblimin", fm = "pa")
echange.faba = psych::bassAckward(ranefs, nfactors=1, rotate = "oblimin", fm = "pa")
echange.fa = fa.sort(psych::fa(ranefs, nfactors=2))
echangefscores <- factor.scores(ranefs, echange.fa)$scores
wdf$echange_f1_early <- echangefscores[,1]
wdf$echange_f2_late <- echangefscores[,2]
wdf <- wdf %>% mutate(Subject = level) %>% select(Subject, echange_f1_early, echange_f2_late)

# entropy ----
# read in data
fname = "rt_encode_output_streams_mixed_by_entropy_ranef.RDS"
hddf <- as.data.frame(readRDS(fname)) %>% select(-model_name, -rhs, -outcome, -std.error, -effect, -group)
# spread by stream and time
hwdf <- hddf %>% group_by(level) %>% filter(term=="v_entropy_wi", evt_time<0) %>% 
  select(-term) %>%
  pivot_wider(names_from = c(stream, evt_time, side), values_from = estimate) %>% ungroup()

# correlation plots

ranefs <- distinct(hwdf) %>% select(-level)
ranef_cor <- corr.test(na.omit(ranefs),method = 'pearson', adjust = 'none')
# parametric correlations on winsorised betas
# clust_cor <- cor(just_rois_w,method = 'pearson')

pdf("entropy_ranef_corr.pdf", width=60, height=60)
corrplot(ranef_cor$r, cl.lim=c(-1,1), 
         method = "circle", tl.cex = 1.5, type = "upper", tl.col = 'black',
         order = "hclust", diag = FALSE,
         addCoef.col="black", addCoefasPercent = FALSE,
         p.mat = ranef_cor$p, sig.level=0.05, insig = "blank")
dev.off()

nr <- nfactors(ranef_cor$r, n=10, rotate = "oblimin", diagonal = FALSE,fm = "pa", n.obs = 70, SMC = FALSE)
e.fa = psych::fa(ranefs, nfactors=1, rotate = "oblimin", fm = "pa")
e.faba = psych::bassAckward(ranefs, nfactors=2, rotate = "oblimin", fm = "pa")
e.fa = fa.sort(psych::fa(ranefs, nfactors=1))
efscores <- factor.scores(ranefs, e.fa)$scores
wdf$e_f1 <- efscores[,1]
wdf <- wdf %>% mutate(Subject = as.integer(Subject))

# add to behavioral data file ----
behavioral_data_file = "~/code/clock_analysis/meg/MEG_n63_behavioral_data_preprocessed_trial_df.RDS"
trial_df <- as_tibble(readRDS(behavioral_data_file))  %>% inner_join(wdf)
saveRDS(trial_df, file = behavioral_data_file)
