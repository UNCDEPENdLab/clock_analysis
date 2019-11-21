# brain-to-behavior analyses with anterior and posterior hippocampal cluster betas
# first run beta_cluster_import_pca_clean.R if not run once already

library(dplyr)
library(tidyverse)
library(psych)
library(corrplot)
library(lme4)
library(ggpubr)
library(lmerTest)
library(stargazer)
library(car)
library(sjstats)
library(sjPlot)
library(emmeans)
library(cowplot)
source('~/code/Rhelpers/screen.lmerTest.R')
source('~/code/Rhelpers/vif.lme.R')
# library(stringi)

# source('~/code/Rhelpers/')
setwd('~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/')

### load data
# load('trial_df_and_vhdkfpe_clusters.Rdata')
# cleaner version with only H, PE and uncertainty trial vars
unsmoothed = F
if (unsmoothed) {
  load('trial_df_and_vh_pe_clusters_u_unsmoothed.Rdata')
} else { load('trial_df_and_vh_pe_clusters_u.Rdata') }

# vif.lme <- function (fit) {
#   ## adapted from rms::vif
#   v <- vcov(fit)
#   nam <- names(fixef(fit))
#   ## exclude intercepts
#   ns <- sum(1 * (nam == "Intercept" | nam == "(Intercept)"))
#   if (ns > 0) {
#     v <- v[-(1:ns), -(1:ns), drop = FALSE]
#     nam <- nam[-(1:ns)] }
#   d <- diag(v)^0.5
#   v <- diag(solve(v/(d %o% d)))
#   names(v) <- nam
#   v }
# 
# # check VIFs of significant effects
# screen.lmerTest <- function (mod,p=NULL) {
#   if (is.null(p)) {p <- .05}
#   c1 <- as.data.frame(coef(summary(mod))[,4:5])
#   dd <- cbind(c1[2:nrow(c1),],as.data.frame(vif.lme(mod)))
#   names(dd)[3] <- 'VIF'
#   dd$`Pr(>|t|)` <- as.numeric(dd$`Pr(>|t|)`)
#   print(dd[dd$`Pr(>|t|)`<p,c(1,3)], digits = 3)}



###############
# Main analyses including model-derived behavioral variables
# hippocampal model-based analysis
mb3hpe_hipp <-  lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + rt_vmax_lag_sc + last_outcome + 
                               v_max_wi_lag + v_entropy_wi + h_HippAntL_neg +  pe_f2_hipp)^2 + 
                     rt_lag_sc:last_outcome:h_HippAntL_neg + 
                     rt_lag_sc:last_outcome:pe_f2_hipp +
                     rt_vmax_lag_sc:trial_neg_inv_sc:h_HippAntL_neg + 
                     rt_vmax_lag_sc:trial_neg_inv_sc:pe_f2_hipp  +
                     (1|id/run), df)
screen.lmerTest(mb3hpe_hipp, .05)
summary(mb3hpe_hipp)
Anova(mb3hpe_hipp, '3')

########
# out-of-session replication with MEG behavioral data
mmb3hpe_hipp <-  lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + rt_vmax_lag_sc + last_outcome + 
                                  v_max_wi_lag + v_entropy_wi +h_HippAntL_neg +  pe_f2_hipp)^2 + 
                        rt_lag_sc:last_outcome:h_HippAntL_neg + 
                        rt_lag_sc:last_outcome:pe_f2_hipp +
                        rt_vmax_lag_sc:trial_neg_inv_sc:h_HippAntL_neg + 
                        rt_vmax_lag_sc:trial_neg_inv_sc:pe_f2_hipp  +
                        (1|id/run), mdf)
screen.lmerTest(mmb3hpe_hipp, .05)
summary(mmb3hpe_hipp)
Anova(mmb3hpe_hipp, '3')


# # NB: stargazer only works with lme4, not lmerTest
mb3hpe_hipp_lme4 <-  lme4::lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + rt_vmax_lag_sc + last_outcome + 
                                    v_max_wi_lag + v_entropy_wi + h_HippAntL_neg +  pe_f2_hipp)^2 + 
                       rt_lag_sc:last_outcome:h_HippAntL_neg + 
                       rt_lag_sc:last_outcome:pe_f2_hipp +
                       rt_vmax_lag_sc:trial_neg_inv_sc:h_HippAntL_neg + 
                       rt_vmax_lag_sc:trial_neg_inv_sc:pe_f2_hipp  +
                       (1|id/run), df)
mmb3hpe_hipp_lme4 <-  lme4::lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + rt_vmax_lag_sc + last_outcome + 
                                     v_max_wi_lag + v_entropy_wi +h_HippAntL_neg +  pe_f2_hipp)^2 + 
                        rt_lag_sc:last_outcome:h_HippAntL_neg + 
                        rt_lag_sc:last_outcome:pe_f2_hipp +
                        rt_vmax_lag_sc:trial_neg_inv_sc:h_HippAntL_neg + 
                        rt_vmax_lag_sc:trial_neg_inv_sc:pe_f2_hipp  +
                        (1|id/run), mdf)
setwd('~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/tables/')

stargazer(mb3hpe_hipp_lme4, mmb3hpe_hipp_lme4, type="html", out="hippo_mb.htm", report = "vcs*",
          digits = 1,single.row=TRUE,omit.stat = "bic",
          column.labels = c("fMRI session", "Out-of-session replication"),
          star.char = c("*", "**", "***"),
          star.cutoffs = c(0.05, 0.01, 0.001),
          notes = c("* p<0.05; ** p<0.01; *** p<0.001"),
          notes.append = F)

## AH replication forest plot
mterms <- names(fixef(mb3hpe_hipp))
setwd('../plots/')
ah <- plot_models(mb3hpe_hipp,mmb3hpe_hipp, rm.terms = mterms[c(-26, -40)], m.labels = c("fMRI", "replication"),
            show.values = T, std.est = "std2", legend.title = "Session", vline.color = "slategray3",
            wrap.labels = 20, axis.labels = c("RT(Vmax) * Ant. hippocampal low entropy response", "Trial * RT(Vmax) * Ant. hippocampal low entropy response"), 
            axis.title = "Less convergence <==> Better convergence on global max")
ah <- ah + ylim(-.01,.25) + geom_hline(yintercept = 0, color = "slategray3")
pdf("ah_beta_models_replication.pdf", height = 3, width = 5)
ah
dev.off()

## PH replication plot
setwd('~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/plots/')
ph <- plot_models(mb3hpe_hipp,mmb3hpe_hipp, rm.terms = mterms[c(-22, -39)], m.labels = c("fMRI", "replication"),
                  show.values = T,  std.est = "std2", legend.title = "Session", vline.color = "slategray3",
                  wrap.labels = 15,  axis.labels = c("Previous RT * Omission * Post. hippocampal PE response","Previous RT * Post. hippocampal PE response"),
                  axis.title = "Greater RT swing  <==>  Smaller RT swing")
ph <- ph + ggplot2::ylim(-.1,.1)
pdf("ph_beta_models_replication.pdf", height = 3, width = 5)
ph
dev.off()

# PH PE cluster suppresses the win-stay-lose-switch behaviors
vs2 <- lmer(v_chosen ~ (trial_neg_inv_sc + last_outcome + 
                          h_HippAntL_neg + pe_f2_hipp)^2 + 
              trial_neg_inv_sc*rewFunc + (1|id), df)
screen.lmerTest(vs2, .01)
vemp <- as_tibble(emmeans(vs2, ~last_outcome | pe_f2_hipp, at = list(pe_f2_hipp = c(-2,2)))) %>% mutate(`Chosen value` = emmean)
p1 <- ggplot(vemp, aes(last_outcome, `Chosen value`, color = pe_f2_hipp, group = pe_f2_hipp)) + 
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), position=position_dodge(width=0.5)) + geom_line(position=position_dodge(width=0.5)) + geom_point(position=position_dodge(width=0.5))
vema <- as_tibble(emmeans(vs2, ~last_outcome | h_HippAntL_neg, at = list(h_HippAntL_neg = c(-2,2)))) %>% mutate(`Chosen value` = emmean)
p2 <- ggplot(vema, aes(last_outcome, `Chosen value`, color = h_HippAntL_neg, group = h_HippAntL_neg)) + 
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),position=position_dodge(width=0.5)) + geom_line(position=position_dodge(width=0.5)) + geom_point(position=position_dodge(width=0.5))
setwd('/plots')
pdf('ah_ph_v_chosen_outcome.pdf', height = 3, width = 7)
ggarrange(p1,p2, ncol = 2)
dev.off()


#################
# Sensitivity analysis for the supplement

# Model-free RT analyses -- behavioral variables
mf1 <- lmer(rt_csv ~ (trial_neg_inv_sc + rt_lag_sc + last_outcome )^2  + (1|id/run), df)
summary(mf1)

##################
# best hippocampal model
mf3hpe <-  lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + last_outcome + 
                               h_HippAntL_neg + pe_f2_hipp)^2 + 
                  rt_lag_sc:last_outcome:h_HippAntL_neg + 
                  rt_lag_sc:last_outcome:pe_f2_hipp + trial_neg_inv_sc*rewFunc + (1|id/run), df)
# summary(mf3hpe)
screen.lmerTest(mf3hpe)
##
## MEG data for out-of-session replication
mmf3hpe <- lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + last_outcome + 
                               h_HippAntL_neg + pe_f2_hipp)^2 + 
                  rt_lag_sc:last_outcome:h_HippAntL_neg + 
                  rt_lag_sc:last_outcome:pe_f2_hipp + trial_neg_inv_sc*rewFunc + (1|id/run), mdf)
screen.lmerTest(mmf3hpe)
summary(mmf3hpe)
Anova(mmf3hpe,'3')
## ascertain replication -- visual check
# plot_models(mmf3hpe,mf3hpe)
# ascertain replication overall 
p1 <- plot_model(mf3hpe, show.values = T)
p2 <- plot_model(mmf3hpe, show.values = T)
pdf("model_free_beta_replication.pdf", height = 6, width = 12)
ggarrange(p1,p2,ncol = 2, labels  = c("fMRI", "MEG"))
dev.off()

# save model statistics for supplement
mf3hpe_lme4 <-  lme4::lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + last_outcome + 
                               h_HippAntL_neg + pe_f2_hipp)^2 + 
                  rt_lag_sc:last_outcome:h_HippAntL_neg + 
                  rt_lag_sc:last_outcome:pe_f2_hipp + trial_neg_inv_sc*rewFunc + (1|id/run), df)
# summary(mf3hpe)
screen.lmerTest(mf3hpe)
###############
# compare performance in fMRI and MEG sessions: better in MEG, particularly in IEV
p1 <- ggplot(df, aes(run_trial, rt_csv, color = rewFunc)) + geom_smooth(method = "gam", formula = y~splines::ns(x,4)) +  coord_cartesian(ylim=c(1300,2300))
p2 <- ggplot(mdf, aes(run_trial, rt_csv, color = rewFunc)) + geom_smooth(method = "gam", formula = y~splines::ns(x,4)) +  coord_cartesian(ylim=c(1300,2300))
ggarrange(p1,p2, labels = c("fMRI", "MEG"), align = "hv")
## MEG data for out-of-session replication
#
mmf3hpe_lme4 <-  lme4::lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + last_outcome + 
                               h_HippAntL_neg + pe_f2_hipp)^2 + 
                  rt_lag_sc:last_outcome:h_HippAntL_neg + 
                  rt_lag_sc:last_outcome:pe_f2_hipp + trial_neg_inv_sc*rewFunc + (1|id/run), mdf)

setwd('~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/tables/')
# NB: stargazer only runs with lmer, not lmerTest objects
stargazer(mf3hpe_lme4, mmf3hpe_lme4, type="html", out="hippo_mf.htm", report = "vcs*",
          digits = 1,single.row=TRUE,omit.stat = "bic",
          column.labels = c("fMRI session", "Out-of-session replication"),
          star.char = c("*", "**", "***"),
          star.cutoffs = c(0.05, 0.01, 0.001),
          notes = c("* p<0.05; ** p<0.01; *** p<0.001"),
          notes.append = F)
##############
# Sensitivity analyses (cont.):
# without covariates: effects are unchanged
summary(m0 <-  lme4::lmer(rt_csv_sc ~  rt_lag_sc*last_outcome*pe_f2_hipp +
                                  rt_vmax_lag_sc*trial_neg_inv_sc*h_HippAntL_neg + 
                                 (1|id/run), df))
summary(mm0 <-  lme4::lmer(rt_csv_sc ~  rt_lag_sc*last_outcome*pe_f2_hipp +
                      rt_vmax_lag_sc*trial_neg_inv_sc*h_HippAntL_neg + 
                      (1|id/run), mdf))
setwd('~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/tables/')
stargazer(m0, mm0, type="html", out="hippo_mb_no_covariates.htm", report = "vcs*",
          digits = 1,single.row=TRUE,omit.stat = "bic",
          column.labels = c("fMRI session", "Out-of-session replication"),
          star.char = c("*", "**", "***"),
          star.cutoffs = c(0.05, 0.01, 0.001),
          notes = c("* p<0.05; ** p<0.01; *** p<0.001"),
          notes.append = F)
#

##############
# controlling for uncertainty
mf3hpe_u <-  lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + last_outcome + scale(u_chosen_lag) +
                               h_HippAntL_neg + pe_f2_hipp)^2 + 
                  rt_lag_sc:last_outcome:h_HippAntL_neg + 
                  rt_lag_sc:last_outcome:pe_f2_hipp + trial_neg_inv_sc*rewFunc + (1|id/run), df)
while (any(grepl("failed to converge", mf3hpe_u@optinfo$conv$lme4$messages) )) {
  print(mf3hpe_u@optinfo$conv$lme4$conv)
  ss <- getME(mf3hpe_u,c("theta","fixef"))
  mf3hpe_u <- update(mf3hpe_u, start=ss)}

# summary(mf3hpe)
screen.lmerTest(mf3hpe_u)
##
## MEG data for out-of-session replication
mmf3hpe_u <-  lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + last_outcome + scale(u_chosen_lag) +
                                 h_HippAntL_neg + pe_f2_hipp)^2 + 
                    rt_lag_sc:last_outcome:h_HippAntL_neg + 
                    rt_lag_sc:last_outcome:pe_f2_hipp + trial_neg_inv_sc*rewFunc + (1|id/run), mdf)
while (any(grepl("failed to converge", mmf3hpe_u@optinfo$conv$lme4$messages) )) {
  print(mmf3hpe_u@optinfo$conv$lme4$conv)
  ss <- getME(mmf3hpe_u,c("theta","fixef"))
  mmf3hpe_u <- update(mmf3hpe_u, start=ss)}

# summary(mmf3hpe)
screen.lmerTest(mmf3hpe_u)

##############
# Sensitivity analyses (cont.):
# R vs. L PH 
#
mb3hpe_hipp_rl <-  lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + rt_vmax_lag_sc + last_outcome + 
                                      v_max_wi_lag + v_entropy_wi + h_HippAntL_neg +  pe_PH)^2 + 
                         rt_lag_sc:last_outcome:h_HippAntL_neg + 
                         rt_lag_sc:last_outcome:pe_PH +
                         rt_vmax_lag_sc:trial_neg_inv_sc:h_HippAntL_neg + 
                         rt_vmax_lag_sc:trial_neg_inv_sc:pe_PH  +
                         (1|id/run), df)
summary(mb3hpe_hipp_rl)
screen.lmerTest(mb3hpe_hipp_rl, .05)
# Anova(mmb3hpe_hipp, '3')

mmb3hpe_hipp_rl <-  lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + rt_vmax_lag_sc + last_outcome + 
                                       v_max_wi_lag + v_entropy_wi +h_HippAntL_neg +  pe_PH)^2 + 
                          rt_lag_sc:last_outcome:h_HippAntL_neg + 
                          rt_lag_sc:last_outcome:pe_PH +
                          rt_vmax_lag_sc:trial_neg_inv_sc:h_HippAntL_neg + 
                          rt_vmax_lag_sc:trial_neg_inv_sc:pe_PH  +
                          (1|id/run), mdf)
screen.lmerTest(mmb3hpe_hipp_rl, .05)
summary(mmb3hpe_hipp_rl)
Anova(mmb3hpe_hipp_rl, '3')


mb3hpe_hipp_rl <-  lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + rt_vmax_lag_sc + last_outcome + 
                                    v_max_wi_lag + v_entropy_wi + h_HippAntL_neg +  pe_PH)^2 + 
                       rt_lag_sc:last_outcome:h_HippAntL_neg + 
                       rt_lag_sc:last_outcome:pe_PH +
                       rt_vmax_lag_sc:trial_neg_inv_sc:h_HippAntL_neg + 
                       rt_vmax_lag_sc:trial_neg_inv_sc:pe_PH  +
                       (1|id/run), df)
# summary(mb3hpe_hipp)
screen.lmerTest(mb3hpe_hipp_rl, .05)
# Anova(mmb3hpe_hipp, '3')

mmb3hpe_hipp_rl <-  lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + rt_vmax_lag_sc + last_outcome + 
                                     v_max_wi_lag + v_entropy_wi +h_HippAntL_neg +  pe_PH)^2 + 
                        rt_lag_sc:last_outcome:h_HippAntL_neg + 
                        rt_lag_sc:last_outcome:pe_PH +
                        rt_vmax_lag_sc:trial_neg_inv_sc:h_HippAntL_neg + 
                        rt_vmax_lag_sc:trial_neg_inv_sc:pe_PH  +
                        (1|id/run), mdf)
while (any(grepl("failed to converge", mmb3hpe_hipp_rl@optinfo$conv$lme4$messages) )) {
  ss <- getME(mmb3hpe_hipp_rl,c("theta","fixef"))
  mmb3hpe_hipp_rl <- update(mmb3hpe_hipp_rl, start=ss)}


screen.lmerTest(mmb3hpe_hipp_rl, .05)
summary(mmb3hpe_hipp_rl)
Anova(mmb3hpe_hipp_rl, '3')

mb3hpe_hipp_r <-  lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + rt_vmax_lag_sc + last_outcome + 
                                      v_max_wi_lag + v_entropy_wi + h_HippAntL_neg +  pe_PH_r)^2 + 
                         rt_lag_sc:last_outcome:h_HippAntL_neg + 
                         rt_lag_sc:last_outcome:pe_PH_r +
                         rt_vmax_lag_sc:trial_neg_inv_sc:h_HippAntL_neg + 
                         rt_vmax_lag_sc:trial_neg_inv_sc:pe_PH_r  +
                         (1|id/run), df)
# summary(mb3hpe_hipp)
screen.lmerTest(mb3hpe_hipp_r, .05)
# Anova(mmb3hpe_hipp, '3')

mmb3hpe_hipp_r <-  lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + rt_vmax_lag_sc + last_outcome + 
                                       v_max_wi_lag + v_entropy_wi +h_HippAntL_neg +  pe_PH_r)^2 + 
                          rt_lag_sc:last_outcome:h_HippAntL_neg + 
                          rt_lag_sc:last_outcome:pe_PH_r +
                          rt_vmax_lag_sc:trial_neg_inv_sc:h_HippAntL_neg + 
                          rt_vmax_lag_sc:trial_neg_inv_sc:pe_PH_r  +
                          (1|id/run), mdf)
screen.lmerTest(mmb3hpe_hipp_r, .05)
summary(mmb3hpe_hipp_r)
Anova(mmb3hpe_hipp_r, '3')

## PH replication plot for mean of R and left
mterms <- names(fixef(mb3hpe_hipp_rl))
setwd('~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/plots/')
ph <- plot_models(mb3hpe_hipp_rl,mmb3hpe_hipp_rl, rm.terms = mterms[c(-22, -39)], m.labels = c("fMRI", "replication"),
                  show.values = T,  std.est = "std2", legend.title = "Session", vline.color = "slategray3",
                  wrap.labels = 15,  axis.labels = c("Previous RT * Omission * Post. hippocampal PE response","Previous RT * Post. hippocampal PE response"),
                  axis.title = "Greater RT swing  <==>  Smaller RT swing")
ph <- ph + ggplot2::ylim(-.1,.1)
pdf("ph_beta_models_replication_mean_LR.pdf", height = 3, width = 5)
ph
dev.off()

# visual sanity checks
# p1 <- ggplot(df, aes(rt_lag, rt_csv, color = pe_f2_hipp_resp, lty = last_outcome)) + geom_smooth(method = 'glm') #+ facet_wrap(~rewFunc)
# p2 <- ggplot(mdf, aes(rt_lag, rt_csv, color = pe_f2_hipp_resp, lty = last_outcome)) + geom_smooth(method = 'glm') #+ facet_wrap(~rewFunc)
# ggarrange(p1,p2, ncol = 1, nrow = 2, labels = c("fMRI", "MEG"))
# 
# p1 <- ggplot(df, aes(rt_vmax_lag, rt_csv, color = pe_f2_hipp_resp)) + geom_smooth(method = 'glm') + facet_wrap(~learning_epoch)
# p2 <- ggplot(mdf, aes(rt_vmax_lag, rt_csv, color = pe_f2_hipp_resp)) + geom_smooth(method = 'glm') + facet_wrap(~learning_epoch)
# ggarrange(p1,p2, ncol = 1, nrow = 2, labels = c("fMRI", "MEG"))

# p1 <- ggplot(df, aes(rt_lag, rt_csv, color = h_HippAntL_resp, lty = last_outcome)) + geom_smooth(method = 'glm') #+ facet_wrap(~learning_epoch)
# p2 <- ggplot(mdf, aes(rt_lag, rt_csv, color = h_HippAntL_resp, lty = last_outcome)) + geom_smooth(method = 'glm') #+ facet_wrap(~learning_epoch)
# ggarrange(p1,p2, ncol = 1, nrow = 2, labels = c("fMRI", "MEG"))
# 
# p1 <- ggplot(df, aes(rt_lag, rt_csv, color = h_HippAntL_resp, lty = last_outcome)) + geom_smooth(method = 'glm') #+ facet_wrap(~learning_epoch)
# p2 <- ggplot(mdf, aes(rt_lag, rt_csv, color = h_HippAntL_resp, lty = last_outcome)) + geom_smooth(method = 'glm') #+ facet_wrap(~learning_epoch)
# ggarrange(p1,p2, ncol = 1, nrow = 2, labels = c("fMRI", "MEG"))
# 
# 
# p1 <- ggplot(df, aes(rt_vmax_lag, rt_csv, color = h_HippAntL_resp)) + geom_smooth(method = 'glm') + facet_wrap(~learning_epoch)
# p2 <- ggplot(mdf, aes(rt_vmax_lag, rt_csv, color = h_HippAntL_resp)) + geom_smooth(method = 'glm') + facet_wrap(~learning_epoch)
# ggarrange(p1,p2, ncol = 1, nrow = 2, labels = c("fMRI", "MEG"))
# 

# # understand rt_vmax_change effect
# ggplot(df, aes(rt_vmax_change, rt_csv, color = pe_f2_hipp_resp)) + geom_smooth(method = "glm")

######################
# Uncertainty models #
######################

# Sanity checks
# ggplot(df, aes(run_trial, u_chosen, group = interaction(rewFunc, rt_lag>2000), color = rewFunc, lty = rt_lag>2000)) + geom_smooth()
# ggplot(df, aes(run_trial, u_chosen, group = rewFunc, color = rewFunc)) + geom_smooth()
# ggplot(df, aes(run_trial, u_chosen_change, group = rewFunc, color = rewFunc)) + geom_smooth(method = 'gam', formula = y ~ splines::ns(x, 3))

# Inspect correlations to estimate the uncertainty/value confound -- not huge
vars <- df %>% select(u_chosen, u_chosen_change, u_chosen_quantile, u_chosen_quantile_change, 
                      v_chosen, v_chosen_quantile, v_chosen_quantile_change, 
                      run_trial, magnitude, probability,rt_csv, rt_lag, rt_vmax_lag)
u_cor <- corr.test(vars,method = 'pearson', adjust = 'none')

setwd('~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/plots')
pdf("u_corr_reset.pdf", width=12, height=12)
corrplot(u_cor$r, cl.lim=c(-1,1),
         method = "circle", tl.cex = 1.5, type = "upper", tl.col = 'black',
         order = "hclust", diag = FALSE,
         addCoef.col="black", addCoefasPercent = FALSE,
         p.mat = u_cor$p, sig.level=0.05, insig = "blank")
dev.off()

# timecourses of choice uncertainty by HIPP response
# ggplot(df, aes(run_trial, rt_csv, color = rewFunc, lty = h_HippAntL_resp, group = interaction(rewFunc, h_HippAntL_resp))) + geom_smooth(method = 'gam', formula = y ~ splines::ns(x, 3)) #+ facet_wrap(~run)
p1 <- ggplot(df, aes(run_trial, u_chosen, color = rt_lag>2000, lty = pe_f2_hipp_resp, group = interaction(rt_lag>2000, pe_f2_hipp_resp))) + geom_smooth(method = 'gam', formula = y ~ splines::ns(x, 3))  + facet_wrap(last_outcome~rewFunc)
p2 <- ggplot(mdf, aes(run_trial, u_chosen, color = rt_lag>2000, lty = pe_f2_hipp_resp, group = interaction(rt_lag>2000, pe_f2_hipp_resp))) + geom_smooth(method = 'gam', formula = y ~ splines::ns(x, 3))  + facet_wrap(last_outcome~rewFunc)
ggarrange(p1,p2)

# RT timecourses by hipp response
p1 <- ggplot(df, aes(run_trial, rt_csv, color = rewFunc, lty = pe_f2_hipp_resp, group = interaction(rewFunc, pe_f2_hipp_resp))) + geom_smooth(method = 'gam', formula = y ~ splines::ns(x, 4))
p2 <- ggplot(mdf, aes(run_trial, rt_csv, color = rewFunc, lty = pe_f2_hipp_resp, group = interaction(rewFunc, pe_f2_hipp_resp))) + geom_smooth(method = 'gam', formula = y ~ splines::ns(x, 4))
ggarrange(p1,p2)

p1 <- ggplot(df, aes(run_trial, rt_csv, color = rewFunc, lty = h_HippAntL_resp, group = interaction(rewFunc, h_HippAntL_resp))) + geom_smooth(method = 'gam', formula = y ~ splines::ns(x, 4))
p2 <- ggplot(mdf, aes(run_trial, rt_csv, color = rewFunc, lty = h_HippAntL_resp, group = interaction(rewFunc, h_HippAntL_resp))) + geom_smooth(method = 'gam', formula = y ~ splines::ns(x, 4))
ggarrange(p1,p2)

# the differences in timecourse of RT swings between high and low PH are visible in fMRI, subtle in MEG
# p1 <- ggplot(df %>% filter(!is.na(last_outcome)), aes(run_trial, rt_swing, color = rewFunc, lty = pe_f2_hipp_resp, group = interaction(rewFunc, pe_f2_hipp_resp))) + geom_smooth(method = 'loess') + facet_wrap(~last_outcome)
# p2 <- ggplot(mdf %>% filter(!is.na(last_outcome)), aes(run_trial, rt_swing, color = rewFunc, lty = pe_f2_hipp_resp, group = interaction(rewFunc, pe_f2_hipp_resp))) + geom_smooth(method = 'loess')+ facet_wrap(~last_outcome)
p1 <- ggplot(df %>% filter(!is.na(last_outcome) & rewFunc!="CEV" & rewFunc!="CEVR"), aes(run_trial, rt_swing, color = pe_f2_hipp_resp)) + geom_smooth(method = 'loess') 
p2 <- ggplot(mdf %>% filter(!is.na(last_outcome) & rewFunc!="CEV" & rewFunc!="CEVR"), aes(run_trial, rt_swing, color = pe_f2_hipp_resp)) + geom_smooth(method = 'loess') 
pdf('rt_swings_by_PH_resp.pdf', height = 5, width = 10)
ggarrange(p1,p2, labels = c("fMRI", "MEG"))
dev.off()

p1 <- ggplot(df %>% filter(!is.na(last_outcome) & rewFunc!="CEV" & rewFunc!="CEVR"), aes(run_trial, rt_swing, color = h_HippAntL_resp)) + geom_smooth(method = 'loess') 
p2 <- ggplot(mdf %>% filter(!is.na(last_outcome) & rewFunc!="CEV" & rewFunc!="CEVR"), aes(run_trial, rt_swing, color = h_HippAntL_resp)) + geom_smooth(method = 'loess') 
pdf('rt_swings_by_AH_resp.pdf', height = 5, width = 10)
ggarrange(p1,p2, labels = c("fMRI", "MEG"))
dev.off()

##
# Builidng the uncertainty model for interactions with betas
# more sanity checks on quantiles (relative uncertainty) -- looks right
# ggplot(df, aes(run_trial, u_chosen)) + geom_smooth()+ facet_wrap(~rewFunc)
ggplot(df, aes(run_trial, u_chosen_quantile)) + geom_smooth()+ facet_wrap(~rewFunc)
ggplot(df, aes(run_trial, u_chosen_quantile_change)) + geom_smooth()+ facet_wrap(~rewFunc)

# Regress value out of uncertainty for plotting
m <- lmer(u_chosen_quantile ~ v_chosen + (1|ID), df)
df$u_chosen_v_partialed_out <- resid(m)
#
p1 <- ggplot(df, aes(run_trial, u_chosen_quantile, color = pe_f2_hipp_resp, lty = h_HippAntL_resp)) + geom_smooth(method = 'gam', formula = y ~ splines::ns(x,3))#+ facet_wrap(~rewFunc)
p2 <- ggplot(df, aes(run_trial, u_chosen_v_partialed_out, color = pe_f2_hipp_resp, lty = h_HippAntL_resp)) + geom_smooth(method = 'gam', formula = y ~ splines::ns(x,3))#+ facet_wrap(~rewFunc)


# most interpretable set of models
umb1 <- lmer(u_chosen_quantile_change ~ (trial_neg_inv_sc + rt_lag_sc + rt_swing + last_outcome + v_entropy_wi + h_HippAntL_neg)^2 +
              (trial_neg_inv_sc + rt_lag_sc + last_outcome + v_entropy_wi + pe_f2_hipp)^2 +
              scale(u_chosen_quantile_lag) + rt_lag_sc*rewFunc + (1|id/run), df)
screen.lmerTest(umb1, .05)

umb1v <- lmer(u_chosen_quantile_change ~ (trial_neg_inv_sc + rt_lag_sc + rt_swing + last_outcome + v_entropy_wi + h_HippAntL_neg)^2 +
              (trial_neg_inv_sc + rt_lag_sc + last_outcome + v_entropy_wi + pe_f2_hipp)^2 + scale(u_chosen_quantile_lag) + v_chosen_quantile_change + rt_lag_sc*rewFunc + (1|id/run), df)
screen.lmerTest(umb1v, .05)
# summary(umb1v)

# attempt MEG replication
mumb1 <- lmer(u_chosen_change ~ (trial_neg_inv_sc + rt_lag_sc + rt_swing + last_outcome + v_entropy_wi + h_HippAntL_neg)^2 +
               (trial_neg_inv_sc + rt_lag_sc + last_outcome + v_entropy_wi + pe_f2_hipp)^2 +
               scale(u_chosen_lag) + rt_lag_sc*rewFunc + (1|id/run), mdf)
screen.lmerTest(mumb1, .05)

mumb1v <- lmer(u_chosen_change ~ (trial_neg_inv_sc + rt_lag_sc + rt_swing + last_outcome + v_entropy_wi + h_HippAntL_neg)^2 +
                (trial_neg_inv_sc + rt_lag_sc + last_outcome + v_entropy_wi + pe_f2_hipp)^2 +
                scale(u_chosen_lag) + v_chosen_change + rt_lag_sc*rewFunc + (1|id/run), mdf)
screen.lmerTest(mumb1v, .05)

# remove entropy to minimize circularity
umb2 <- lmer(u_chosen_quantile_change ~ trial_neg_inv_sc + rt_lag_sc * rt_vmax_lag_sc * h_HippAntL_neg + 
              rt_lag_sc * rt_vmax_lag_sc * pe_f2_hipp + last_outcome + 
              scale(u_chosen_quantile_lag) +  rt_lag_sc*rewFunc + (1|id/run), df)
screen.lmerTest(umb2, .05)
# summary(umb2)
umb2v <- lmer(u_chosen_quantile_change ~ trial_neg_inv_sc + rt_lag_sc * rt_vmax_lag_sc * h_HippAntL_neg + 
                rt_lag_sc * rt_vmax_lag_sc * pe_f2_hipp + last_outcome + 
                scale(u_chosen_quantile_lag) +  v_chosen_quantile_change + rt_lag_sc*rewFunc + (1|id/run), df)
screen.lmerTest(umb2v, .05)
# summary(umb2v)

# use U instead of quantile(U) to enable MEG replication
umb3 <- lmer(u_chosen_change ~ trial_neg_inv_sc + rt_lag_sc * rt_vmax_lag_sc * h_HippAntL_neg + 
               rt_lag_sc * rt_vmax_lag_sc * pe_f2_hipp + last_outcome + 
               scale(u_chosen_lag) +  rt_lag_sc*rewFunc + (1|id/run), df)
screen.lmerTest(umb3, .05)
# summary(umb3)
umb3v <- lmer(u_chosen_change ~ trial_neg_inv_sc + rt_lag_sc * rt_vmax_lag_sc * h_HippAntL_neg + 
                rt_lag_sc * rt_vmax_lag_sc * pe_f2_hipp + last_outcome + 
                scale(u_chosen_lag) +  v_chosen_change + rt_lag_sc*rewFunc + (1|id/run), df)
screen.lmerTest(umb3v, .05)
# summary(umb3v)

# MEG replication
mumb3 <- lmer(u_chosen_change ~ trial_neg_inv_sc + rt_lag_sc * rt_vmax_lag_sc * h_HippAntL_neg + 
               rt_lag_sc * rt_vmax_lag_sc * pe_f2_hipp + last_outcome + 
               scale(u_chosen_lag) +  rt_lag_sc*rewFunc + (1|id/run), mdf)
screen.lmerTest(mumb3, .05)
# summary(mumb3)
mumb3v <- lmer(u_chosen_change ~ trial_neg_inv_sc + rt_lag_sc * rt_vmax_lag_sc * h_HippAntL_neg + 
                rt_lag_sc * rt_vmax_lag_sc * pe_f2_hipp + last_outcome + 
                scale(u_chosen_lag) +  v_chosen_change + rt_lag_sc*rewFunc + (1|id/run), mdf)
screen.lmerTest(mumb3v, .05)
# summary(mumb3v)

# plot uncertainty change following long vs short RTs
setwd('~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/plots')
p1 <- ggplot(df, aes(rt_lag, u_chosen, color = pe_f2_hipp_resp)) + geom_smooth(method = "gam", formula = y~splines::ns(x,2)) + facet_wrap(~rewFunc)
p2 <- ggplot(mdf, aes(rt_lag, u_chosen, color = pe_f2_hipp_resp)) + geom_smooth(method = "gam", formula = y~splines::ns(x,2)) + facet_wrap(~rewFunc)
# pdf('uncertainty_change_prev_RT_by_PH.pdf', height = 6, width = 12)
ggarrange(p1,p2, labels = c("fMRI", "MEG"))
# dev.off()
# summary(ub3v)

# sanity check in simple models: they both predict uncertainty-aversion
us1 <- lmer(u_chosen_quantile ~ (trial_neg_inv_sc + h_HippAntL_neg + rewFunc)^2 +
              (trial_neg_inv_sc + pe_f2_hipp + rewFunc)^2 + (1|id/run), df)
screen.lmerTest(us1, .05)

us1v <- lmer(u_chosen_quantile ~ (trial_neg_inv_sc + h_HippAntL_neg + rewFunc)^2 +
              (trial_neg_inv_sc + pe_f2_hipp + rewFunc)^2 + v_chosen_quantile + (1|id/run), df)
screen.lmerTest(us1v, .05)


mus1 <- lmer(u_chosen ~ (trial_neg_inv_sc + h_HippAntL_neg + rewFunc)^2 +
              (trial_neg_inv_sc + pe_f2_hipp + rewFunc)^2 + (1|id/run), mdf)
screen.lmerTest(mus1, .05)
mus1v <- lmer(u_chosen ~ (trial_neg_inv_sc + h_HippAntL_neg + rewFunc)^2 +
               (trial_neg_inv_sc + pe_f2_hipp + rewFunc)^2 + v_chosen + (1|id/run), mdf)
screen.lmerTest(mus1v, .05)



# demonstrate that this holds across contingencies and early/late learning
pdf('AH_entropy_uncertainty_aversion_by_cond.pdf', height = 6, width = 8)
ggplot(df %>% filter(!is.na(v_entropy_wi)), aes(run_trial, u_chosen_quantile, lty = v_entropy_wi>0, color = h_HippAntL_resp)) + 
  geom_smooth(method = 'gam', formula = y ~ splines::ns(x,3)) + facet_wrap(~rewFunc)
dev.off()


# more comprehensive models evaluated and rejected because of excessive complexity
# ub3a <- lmer(u_chosen_quantile_change ~ (trial_neg_inv_sc + rt_lag_sc + rt_vmax_lag_sc + last_outcome + v_max_wi_lag + v_entropy_wi + h_f1_fp)^3 +
#                (trial_neg_inv_sc + rt_lag_sc + rt_vmax_lag_sc + last_outcome + v_max_wi_lag + v_entropy_wi + h_HippAntL_neg)^3 +
#                (trial_neg_inv_sc + rt_lag_sc + rt_vmax_lag_sc + last_outcome + v_max_wi_lag + v_entropy_wi + pe_f1_cort_str)^3 +
#                (trial_neg_inv_sc + rt_lag_sc + rt_vmax_lag_sc + last_outcome + v_max_wi_lag + v_entropy_wi + pe_f2_hipp)^3 +
#                v_max_b + v_entropy_b + rt_lag_sc*scale(u_chosen_quantile_lag) + trial_neg_inv_sc*scale(run) + (1|id/run), df)
# screen.lmerTest(ub3a, .01)
# 
# ub3e <- lmer(u_chosen_quantile ~ (trial_neg_inv_sc + rt_lag_sc + rt_vmax_lag_sc + last_outcome + v_max_wi_lag + v_entropy_wi + h_f1_fp)^3 +
#               (trial_neg_inv_sc + rt_lag_sc + rt_vmax_lag_sc + last_outcome + v_max_wi_lag + v_entropy_wi + I(-h_f2_neg_paralimb))^3 +
#               (trial_neg_inv_sc + rt_lag_sc + rt_vmax_lag_sc + last_outcome + v_max_wi_lag + v_entropy_wi + pe_f1_cort_str)^3 +
#               (trial_neg_inv_sc + rt_lag_sc + rt_vmax_lag_sc + last_outcome + v_max_wi_lag + v_entropy_wi + pe_f2_hipp)^3 +
#               v_max_b + v_entropy_b + rt_lag_sc*scale(u_chosen_quantile_lag) + trial_neg_inv_sc*scale(run) + 
#               h_f1_fp*scale(u_chosen_quantile_lag) + I(-h_f2_neg_paralimb)*scale(u_chosen_quantile_lag) +pe_f1_cort_str*scale(u_chosen_quantile_lag) +pe_f2_hipp*scale(u_chosen_quantile_lag) + (1|id/run), edf)
# screen.lmerTest(ub3e, .01)

# plot the striking effect of HIPP on uncertainty sensitivity
p1 <- ggplot(df, aes(run_trial, v_chosen, lty = pe_f2_hipp_resp, group = pe_f2_hipp_resp)) + geom_smooth() #+ facet_wrap(~run)
p2 <- ggplot(df, aes(run_trial, u_chosen_quantile, lty = pe_f2_hipp_resp, group = (pe_f2_hipp_resp))) + geom_smooth() #+ facet_wrap(~run)
p3 <- ggplot(df, aes(run_trial, v_chosen, lty = h_HippAntL_resp, group = h_HippAntL_resp)) + geom_smooth()#+ facet_wrap(~run)
p4 <- ggplot(df, aes(run_trial, u_chosen_quantile, lty = h_HippAntL_resp, group = h_HippAntL_resp)) + geom_smooth() #+ facet_wrap(~run)
pdf("PH_AH_on_u_sensitivity.pdf", height = 8, width = 8)
ggarrange(p1,p2,p3,p4,ncol = 2, nrow = 2)
dev.off()

p1 <- ggplot(df, aes(run_trial, u_chosen_quantile, lty = pe_f1_cort_str_resp, group = interaction(pe_f1_cort_str_resp, last_outcome), color = last_outcome)) + geom_smooth() #+ facet_wrap(~run)
p2 <- ggplot(df, aes(run_trial, u_chosen_quantile, lty = pe_f2_hipp_resp, group = interaction(pe_f2_hipp_resp, last_outcome), color = last_outcome)) + geom_smooth() #+ facet_wrap(~run)
p3 <- ggplot(df, aes(run_trial, u_chosen_quantile, lty = h_f1_fp_resp, group = interaction(h_f1_fp_resp, last_outcome), color = last_outcome)) + geom_smooth() #+ facet_wrap(~run)
p4 <- ggplot(df, aes(run_trial, u_chosen_quantile, lty = h_HippAntL_resp, group = interaction(h_HippAntL_resp, last_outcome), color = last_outcome)) + geom_smooth() #+ facet_wrap(~run)

pdf("clusters_u_sensitivity_reward.pdf", height = 4, width = 8)
ggarrange(p1,p2,p3, p4, ncol = 2, nrow = 2)
dev.off()


# # simple models: do they swing in the direction of greater uncertainty?
# us1 <- lmer(u_chosen_quantile_change ~ (trial_neg_inv_sc + last_outcome)^2 + 
#               trial_neg_inv_sc*rewFunc + (1|id/run), df)
# 
# # incorporate rt swings -- hard to interpret...
# us1r <- lmer(u_chosen_quantile_change ~ (trial_neg_inv_sc + last_outcome )^2 + rt_csv_sc * rt_lag_sc +
#               trial_neg_inv_sc*rewFunc + (1|id/run), df)
# screen.lmerTest(us1r)
# us1v <- lmer(u_chosen_quantile_change ~ (trial_neg_inv_sc + last_outcome)^2 + 
#               trial_neg_inv_sc*rewFunc + scale(v_chosen) + (1|id), df)
# screen.lmerTest(us1v)
# # residualize u_chosen_quantile for v_chosen
# uv1 <- lm(u_chosen_quantile ~ 1 + v_chosen + (1|id), df)
# df$u_chosen_quantile_resid <- resid(uv1)
# uvc1 <- lm(u_chosen_quantile_change ~ 1 + v_chosen + (1|id), df)
# df$u_chosen_quantile_change_resid <- NA
# df$u_chosen_quantile_change_resid[df$trial>1] <- resid(uvc1)
# 
# ggplot(df, aes(run_trial, u_chosen_quantile_resid, color = rewFunc)) + geom_smooth()
# ggplot(df, aes(run_trial, u_chosen_quantile_change_resid, color = rewFunc)) + geom_smooth()
# 
# screen.lmerTest(us1)
# us2 <- lmer(u_chosen_quantile_change_resid ~ (trial_neg_inv_sc + last_outcome + h_f1_fp)^2 +
#               (trial_neg_inv_sc + last_outcome+ h_HippAntL_neg)^2 +
#               (trial_neg_inv_sc+ last_outcome + pe_f1_cort_str)^2 +
#             (trial_neg_inv_sc+ last_outcome + pe_f2_hipp)^2 +
#               trial_neg_inv_sc*rewFunc + u_chosen_quantile_resid + (1|id), df)
# screen.lmerTest(us2, .01)
# vs2 <- lmer(v_chosen ~ (trial_neg_inv_sc + last_outcome + h_f1_fp)^2 +
#               (trial_neg_inv_sc + last_outcome+ h_HippAntL_neg)^2 +
#               (trial_neg_inv_sc+ last_outcome + pe_f1_cort_str)^2 +
#               (trial_neg_inv_sc+ last_outcome + pe_f2_hipp)^2 +
#               trial_neg_inv_sc*rewFunc + scale(v_chosen_lag) + (1|id), df)
# screen.lmerTest(vs2, .01)
# 
# us2r <- lmer(u_chosen_quantile_resid ~ (trial_neg_inv_sc + last_outcome + h_f1_fp)^2 +
#               (trial_neg_inv_sc + last_outcome+ h_HippAntL_neg)^2 +
#               (trial_neg_inv_sc+ last_outcome + pe_f1_cort_str)^2 +
#               (trial_neg_inv_sc+ last_outcome + pe_f2_hipp)^2 +
#               trial_neg_inv_sc*rewFunc + scale(u_chosen_quantile_lag) + (1|id), df)
# screen.lmerTest(us2r, .01)
# summary(us2r)
# Anova(us2r, '3')
# understand the time course of uncertainty residualized for value
ggplot(df, aes(run_trial, u_chosen_quantile_resid, color = pe_f2_hipp_resp, lty = last_outcome)) + geom_smooth(method = 'gam', formula = y~splines::ns(x,3))
ggplot(df, aes(run_trial, u_chosen_quantile, color = rewFunc, lty = last_outcome)) + geom_smooth(method = 'gam', formula = y~splines::ns(x,3))
ggplot(df, aes(run_trial, v_chosen, color = rewFunc, lty = last_outcome)) + geom_smooth(method = 'gam', formula = y~splines::ns(x,3))

p1 <- ggplot(df %>% filter(!is.na(last_outcome)), aes(pe_f2_hipp, u_chosen_quantile_change, color = last_outcome)) + geom_smooth(method = 'gam') #+ 
  geom_hline(yintercept = 1437.59)#+ facet_wrap(~run)
p2 <- ggplot(df %>% filter(!is.na(last_outcome)), aes(pe_f1_cort_str, u_chosen_quantile_change, color = last_outcome)) + geom_smooth(method = 'gam') #+ 
  geom_hline(yintercept = 1437.59)#+ facet_wrap(~run)
pdf("pe_clusters_u_reward.pdf", height = 4, width = 8)
ggarrange(p1,p2, ncol = 2, nrow = 1)
dev.off()
# # control for value: this changes results
# us3 <- lmer(u_chosen_quantile ~ (trial_neg_inv_sc + last_outcome + h_f1_fp)^2 +
#               (trial_neg_inv_sc + last_outcome + h_HippAntL_neg)^2 +
#               (trial_neg_inv_sc + last_outcome + pe_f1_cort_str)^2 +
#               (trial_neg_inv_sc + last_outcome + pe_f2_hipp)^2 +
#               trial_neg_inv_sc*rewFunc + scale(v_chosen) + (1|id), df)
# screen.lmerTest(us3, .05)
# Anova(us3, '3')
# anova(us1,us2,us3)

p1 <- ggplot(df %>% filter(!is.na(last_outcome)), aes(pe_f2_hipp, v_chosen, color = last_outcome)) + geom_smooth(method = 'gam') + 
  geom_hline(yintercept = 27.62)#+ facet_wrap(~run)
p2 <- ggplot(df %>% filter(!is.na(last_outcome)), aes(pe_f1_cort_str, v_chosen, color = last_outcome)) + geom_smooth(method = 'gam') + 
  geom_hline(yintercept = 27.62)#+ facet_wrap(~run)
pdf("pe_clusters_v_reward.pdf", height = 4, width = 8)
ggarrange(p1,p2, ncol = 2, nrow = 1)
dev.off()

p1 <- ggplot(df %>% filter(!is.na(last_outcome)), aes(pe_f2_hipp, rt_csv, color = last_outcome)) + geom_smooth(method = 'gam') 
p2 <- ggplot(df %>% filter(!is.na(last_outcome)), aes(pe_f1_cort_str, rt_csv, color = last_outcome)) + geom_smooth(method = 'gam')
p3 <- ggplot(df %>% filter(!is.na(last_outcome)), aes(-h_HippAntL, rt_csv, color = last_outcome)) + geom_smooth(method = 'gam') 
pdf("pe_ah_clusters_rt_reward.pdf", height = 8, width = 20)
ggarrange(p1,p2, p3, ncol = 3, nrow = 1)
dev.off()


# do the PH PE people RT-swing more post rewards?!!!

# predict RT with u_chosen_quantile and HIPP
urs1 <- lmer(rt_csv ~ (trial_neg_inv_sc + rt_lag_sc + last_outcome + u_chosen_quantile)^2 + 
               trial_neg_inv_sc * rewFunc +  (1|id/run), df)
screen.lmerTest(urs1, .01)
vif(urs1)
Anova(urs1)
urs2 <- lmer(rt_csv ~ (trial_neg_inv_sc + rt_lag_sc + last_outcome + u_chosen_quantile + h_f1_fp)^3 +
               (trial_neg_inv_sc + rt_lag_sc + last_outcome + u_chosen_quantile + h_HippAntL_neg)^3 + 
               (trial_neg_inv_sc + rt_lag_sc + last_outcome + u_chosen_quantile + pe_f1_cort_str)^3 + 
               (trial_neg_inv_sc + rt_lag_sc + last_outcome + u_chosen_quantile + pe_f2_hipp)^3 + 
               trial_neg_inv_sc * rewFunc +
               (1|id/run), df)
screen.lmerTest(urs2, .01)
Anova(urs2)
# unpack PH*u_chosen_quantile
ggplot(df, aes(rt_lag, rt_csv, lty = last_outcome, color = pe_f2_hipp_resp)) + geom_smooth(method = "gam",
                                                                                            formula = y ~ splines::ns(x,2))
ggplot(df, aes(rt_lag, rt_csv, lty = last_outcome, color = h_f1_fp>0)) + geom_smooth(method = "gam",
                                                                                           formula = y ~ splines::ns(x,2))

# ideas for improving uncertainty analyses:
# try ML Cox survival with time-varying within-trial U
# look at relative rather than absolute uncertainty of the choice

### plot single subject
sdf <- df[df$id==10811 & df$emotion == "scram" & (df$rewFunc=="IEV" | df$rewFunc=="DEV"),] %>% filter(run_trial>1)
pdf("single_subject_swings_10811.pdf", height = 3, width = 5)
ggplot(sdf, aes(run_trial, rt_csv/1000, color = ev, size = score_csv)) + geom_point() + facet_wrap(~rewFunc) + 
  xlab("Trial") + ylab("Response time, s") + scale_color_viridis_c(option = "plasma", name = "Expected value") + theme_dark()
dev.off()
#
save(file = 'vhd_u_meg_models.Rdata', list = ls(all.names = TRUE))

# # lm to elucidate interactions -- effects very similar to lmer
# lmf3hpe <-  lm(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + last_outcome + 
#                                h_HippAntL_neg + pe_f2_hipp)^2 + 
#                   rt_lag_sc:last_outcome:h_HippAntL_neg + 
#                   rt_lag_sc:last_outcome:pe_f2_hipp + trial_neg_inv_sc*rewFunc, df)
# lmfa <- as_tibble(emtrends(lmf3hpe, var = "rt_lag_sc", specs = c("last_outcome", "h_HippAntL_neg"), at = list(h_HippAntL_neg = c(-2,2), rt_lag_sc = c(-1,2))))
# lmfp <- as_tibble(emtrends(lmf3hpe,var = "rt_lag_sc", specs = c("last_outcome", "pe_f2_hipp"), at = list(pe_f2_hipp = c(-2,2))))
# lmmf3hpe <- lm(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + last_outcome + 
#                                h_HippAntL_neg + pe_f2_hipp)^2 + 
#                   rt_lag_sc:last_outcome:h_HippAntL_neg + 
#                   rt_lag_sc:last_outcome:pe_f2_hipp + trial_neg_inv_sc*rewFunc, mdf)
# lmmfa <- as_tibble(emtrends(lmmf3hpe,var = "rt_lag_sc", specs = c("last_outcome", "h_HippAntL_neg"), at = list(h_HippAntL_neg = c(-2,2))))
# lmmfp <- as_tibble(emtrends(lmmf3hpe,var = "rt_lag_sc", specs = c("last_outcome", "pe_f2_hipp"), at = list(pe_f2_hipp = c(-2,2))))
# 
# # PH increases exploration post-reward, AH increases stickiness post-omission
# p1 <- ggplot(lmfa, aes(last_outcome, rt_lag_sc.trend, group = h_HippAntL_neg, color = h_HippAntL_neg)) + geom_line() + geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL))
# p2 <- ggplot(lmmfa, aes(last_outcome, rt_lag_sc.trend, group = h_HippAntL_neg, color = h_HippAntL_neg)) + geom_line() + geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL))
# p3 <- ggplot(lmfp, aes(last_outcome, rt_lag_sc.trend, group = pe_f2_hipp, color = pe_f2_hipp)) + geom_line() + geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL))
# p4 <- ggplot(lmmfp, aes(last_outcome, rt_lag_sc.trend, group = pe_f2_hipp, color = pe_f2_hipp)) + geom_line() + geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL))
