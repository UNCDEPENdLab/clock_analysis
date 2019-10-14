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
library(stringi)

# source('~/code/Rhelpers/')
setwd('~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/')
# load('trial_df_and_vhdkfpe_clusters.Rdata')
# cleaner version with only H, PE and uncertainty trial vars
unsmoothed = F
if (unsmoothed) {
  load('trial_df_and_vh_pe_clusters_u_unsmoothed.Rdata')
} else { load('trial_df_and_vh_pe_clusters_u.Rdata') }

######
# end of preprocessing

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

# check VIFs of significant effects
screen.lmerTest <- function (mod,p=NULL) {
  if (is.null(p)) {p <- .05}
  c1 <- as.data.frame(coef(summary(mod))[,4:5])
  dd <- cbind(c1[2:nrow(c1),],as.data.frame(vif.lme(mod)))
  names(dd)[3] <- 'VIF'
  dd$`Pr(>|t|)` <- as.numeric(dd$`Pr(>|t|)`)
  print(dd[dd$`Pr(>|t|)`<p,c(1,3)], digits = 3)}
#####
## "model-free analyses"


## let's set decay aside for now

# MF RT analyses
mf1 <- lmer(rt_csv ~ (trial_neg_inv_sc + rt_lag_sc + last_outcome )^2  + rt_lag_sc + (1|id/run), df)
summary(mf1)
#################

## earlier code
# best model
# mf3hpedh3 <-  lmer(rt_csv ~ (scale(-1/run_trial) + scale(rt_lag) +   omission_lag + 
#                                h_f1_fp + I(-h_f2_neg_paralimb) + I(-dh_f_neg_vmpfc_precun) + pe_f1_cort_str + pe_f2_hipp)^2 + 
#                      scale(rt_lag):omission_lag:h_f1_fp + scale(rt_lag):omission_lag:I(-h_f2_neg_paralimb) + 
#                      scale(rt_lag):omission_lag:I(-dh_f_neg_vmpfc_precun) +
#                      scale(rt_lag):omission_lag:pe_f1_cort_str + scale(rt_lag):omission_lag:pe_f2_hipp +
#                      scale(rt_lag):omission_lag:I(-dh_f_neg_vmpfc_precun) +
#                      (1|id/run), df)
# #screen.lmerTest(mf3hpedh3)
# 
# # swap in hippocampal low-H cluster for vmpfcHipp factor -- fit is worse
# mf3hpedh3hipp <-  lmer(rt_csv ~ (scale(-1/run_trial) + scale(rt_lag) +   omission_lag + 
#                                    I(-h_HippAntL) +  pe_f2_hipp)^2 + 
#                          scale(rt_lag):omission_lag:h_f1_fp + scale(rt_lag):omission_lag:I(-h_HippAntL) + 
#                          scale(rt_lag):omission_lag:pe_f2_hipp + (1|id/run), df)
# screen.lmerTest(mf3hpedh3hipp)
# mmf3hpedh3hipp <-  lmer(rt_csv ~ (scale(-1/run_trial) + scale(rt_lag) +   omission_lag + 
#                                    I(-h_HippAntL) +  pe_f2_hipp)^2 + 
#                          scale(rt_lag):omission_lag:I(-h_HippAntL) + 
#                          scale(rt_lag):omission_lag:pe_f2_hipp + scale(-1/run_trial)*rewFunc + (1|id/run), mdf)
# screen.lmerTest(mmf3hpedh3hipp)
# 
# stargazer(mf1, mf2h, mf2v, mf2pe,mf2dh,mf3hpedh3, type="html", out="mf.htm", report = "vcs*",
#           digits = 1,single.row=TRUE,omit.stat = "bic",
#           column.labels = c("Null", "Value", "Prediction error", "Entropy", "Entropy change", "Best model"),
#           star.char = c("*", "**", "***"),
#           star.cutoffs = c(0.05, 0.01, 0.001),
#           notes = c("* p<0.05; ** p<0.01; *** p<0.001"), 
#           notes.append = F)
# m <- mf3hpedh3hipp
# library(emmeans)
# 
# anova(mf1,mf2h,mf2v,mf2d,mf2pe, mf2pef2only,mf2kf, mf2dh,mf2dhp, mf3hpedh3, mf3hpedh3hipp)
# 
##################
# best hippocampal model
mf3hpe <-  lmer(rt_csv ~ (trial_neg_inv_sc + rt_lag_sc + last_outcome + 
                            h_HippAntL_neg + pe_f2_hipp)^2 + 
                     rt_lag_sc:last_outcome:h_HippAntL_neg + 
                     rt_lag_sc:last_outcome:pe_f2_hipp + trial_neg_inv_sc*rewFunc + (1|id/run), df[df$run_trial>2,])
summary(mf3hpe)
screen.lmerTest(mf3hpe)
##
## MEG data
mmf3hpe <-  lmer(rt_csv ~ (trial_neg_inv_sc + rt_lag_sc + last_outcome + 
                             h_HippAntL_neg + pe_f2_hipp)^2 + 
                  rt_lag_sc:last_outcome:h_HippAntL_neg + 
                  rt_lag_sc:last_outcome:pe_f2_hipp + trial_neg_inv_sc*rewFunc + (1|id/run), mdf)
screen.lmerTest(mmf3hpe)

mfterms <- names(fixef(mf3hpe))
# mmterms <- attributes(terms(mmf3hpe))$term.labels
# mterms <- stri_replace_all(mterms, "last_outcomeReward", fixed="last_outcome")

pdf("ph_beta_models_replication.pdf", height = 3, width = 5)
plot_models(mf3hpe,mmf3hpe, rm.terms = mfterms[c(-16, -24)], m.labels = c("fMRI", "replication"),
            show.values = T, axis.lim = c(-55,10), std.est = "std2", legend.title = "Session", vline.color = "slategray3",
            wrap.labels = 10,  axis.labels = c("Previous RT * Reward * PH","Previous RT * PH"),
            axis.title = "Greater RT swing  <==>  Smaller RT swing")
dev.off()

# 
# p1 <- plot_model(mf3hpe, show.values = T)
# p2 <- plot_model(mmf3hpe, show.values = T)
# pdf("model_free_beta_replication.pdf", height = 6, width = 12)
# ggarrange(p1,p2,ncol = 2, labels  = c("fMRI", "MEG"))
# dev.off()
###############
# purely hippocampal model-based
mb3hpe_hipp <-  lmer(rt_csv ~ (trial_neg_inv_sc + rt_lag_sc + rt_vmax_lag_sc + last_outcome + 
                               v_max_wi_lag + v_entropy_wi + h_HippAntL_neg +  pe_f2_hipp)^2 + 
                     rt_lag_sc:last_outcome:h_HippAntL_neg + 
                     rt_lag_sc:last_outcome:pe_f2_hipp +
                     rt_vmax_lag_sc:trial_neg_inv_sc:h_HippAntL_neg + 
                     rt_vmax_lag_sc:trial_neg_inv_sc:pe_f2_hipp  +
                     (1|id/run), df)
summary(mb3hpe_hipp)
screen.lmerTest(mb3hpe_hipp, .05)
Anova(mmb3hpe_hipp, '3')
anova(mf3hpe, mb3hpe_hipp)
########
# out-of-session replication with MEG behavioral data
mmb3hpe_hipp <-  lmer(rt_csv ~ (trial_neg_inv_sc + rt_lag_sc + rt_vmax_lag_sc + last_outcome + 
                                  v_max_wi_lag + v_entropy_wi +h_HippAntL_neg +  pe_f2_hipp)^2 + 
                        rt_lag_sc:last_outcome:h_HippAntL_neg + 
                        rt_lag_sc:last_outcome:pe_f2_hipp +
                        rt_vmax_lag_sc:trial_neg_inv_sc:h_HippAntL_neg + 
                        rt_vmax_lag_sc:trial_neg_inv_sc:pe_f2_hipp  +
                        (1|id/run), mdf)
screen.lmerTest(mmb3hpe_hipp, .05)
summary(mmb3hpe_hipp)
Anova(mmb3hpe_hipp, '3')
mterms <- attributes(terms(mb3hpe_hipp))$term.labels
mmterms <- attributes(terms(mb3hpe_hipp))$term.labels
library(stringi)
mterms <- stri_replace_all(mterms, "last_outcomeReward", fixed="last_outcome")

pdf("ah_beta_models_replication.pdf", height = 3, width = 5)
plot_models(mb3hpe_hipp,mmb3hpe_hipp, rm.terms = mterms[c(-25, -39)], m.labels = c("fMRI", "replication"),
            show.values = T, axis.lim = c(-50,175), std.est = "std2", legend.title = "Session", vline.color = "slategray3",
            wrap.labels = 50, axis.labels = c("RT(Vmax) * AH", "Trial * RT(Vmax) * AH"), 
            axis.title = "Less convergence <==> More convergence on global max")
dev.off()

# visual checks
p1 <- ggplot(df, aes(rt_lag, rt_csv, color = pe_f2_hipp_resp, lty = last_outcome)) + geom_smooth(method = 'glm') #+ facet_wrap(~rewFunc)
p2 <- ggplot(mdf, aes(rt_lag, rt_csv, color = pe_f2_hipp_resp, lty = last_outcome)) + geom_smooth(method = 'glm') #+ facet_wrap(~rewFunc)
ggarrange(p1,p2, ncol = 1, nrow = 2, labels = c("fMRI", "MEG"))

p1 <- ggplot(df, aes(rt_vmax_lag, rt_csv, color = pe_f2_hipp_resp)) + geom_smooth(method = 'glm') + facet_wrap(~learning_epoch)
p2 <- ggplot(mdf, aes(rt_vmax_lag, rt_csv, color = pe_f2_hipp_resp)) + geom_smooth(method = 'glm') + facet_wrap(~learning_epoch)
ggarrange(p1,p2, ncol = 1, nrow = 2, labels = c("fMRI", "MEG"))

## initial replication code
# mmb3hpe_hipp <-  lmer(rt_csv ~ (scale(-1/run_trial) + scale(rt_lag) + scale(rt_vmax_lag) + omission_lag + 
#                                  v_max_wi_lag + v_entropy_wi +rt_vmax_change +  I(-h_HippAntL) +  pe_f2_hipp)^2 + 
#                        scale(rt_lag):omission_lag:I(-h_HippAntL) + 
#                        scale(rt_lag):omission_lag:pe_f2_hipp +
#                        scale(rt_vmax_lag):scale(-1/run_trial):I(-h_HippAntL) + 
#                        scale(rt_vmax_lag):scale(-1/run_trial):pe_f2_hipp  +
#                        v_max_b + v_entropy_b + (1|id/run), mdf)
# screen.lmerTest(mmb3hpe_hipp, .05)

p1 <- ggplot(df, aes(rt_lag, rt_csv, color = h_HippAntL_resp, lty = last_outcome)) + geom_smooth(method = 'glm') #+ facet_wrap(~learning_epoch)
p2 <- ggplot(mdf, aes(rt_lag, rt_csv, color = h_HippAntL_resp, lty = last_outcome)) + geom_smooth(method = 'glm') #+ facet_wrap(~learning_epoch)
ggarrange(p1,p2, ncol = 1, nrow = 2, labels = c("fMRI", "MEG"))

p1 <- ggplot(df, aes(rt_vmax_lag, rt_csv, color = h_HippAntL_resp)) + geom_smooth(method = 'glm') + facet_wrap(~learning_epoch)
p2 <- ggplot(mdf, aes(rt_vmax_lag, rt_csv, color = h_HippAntL_resp)) + geom_smooth(method = 'glm') + facet_wrap(~learning_epoch)
ggarrange(p1,p2, ncol = 1, nrow = 2, labels = c("fMRI", "MEG"))



# understand rt_vmax_change effect
ggplot(df, aes(rt_vmax_change, rt_csv, color = pe_f2_hipp_resp)) + geom_smooth(method = "glm")

########
# uncertainty models
# ggplot(df, aes(run_trial, u_chosen, group = run, color = run)) + geom_smooth(method = 'gam', formula = y ~ splines::ns(x, 2))
# ggplot(df, aes(run_trial, u_chosen, group = interaction(run, rt_lag>2000), color = run, lty = rt_lag>2000)) + geom_smooth()
ggplot(df, aes(run_trial, u_chosen, group = interaction(rewFunc, rt_lag>2000), color = rewFunc, lty = rt_lag>2000)) + geom_smooth()
ggplot(df, aes(run_trial, u_chosen, group = rewFunc, color = rewFunc)) + geom_smooth()


ggplot(df, aes(run_trial, u_chosen_change, group = rewFunc, color = rewFunc)) + geom_smooth(method = 'gam', formula = y ~ splines::ns(x, 3))
vars <- df %>% select(u_chosen, u_chosen_change, run_trial, magnitude, probability, rt_lag, rt_vmax_lag)
u_cor <- corr.test(vars,method = 'pearson', adjust = 'none')

setwd('~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/')
pdf("u_corr_reset.pdf", width=12, height=12)
corrplot(u_cor$r, cl.lim=c(-1,1),
         method = "circle", tl.cex = 1.5, type = "upper", tl.col = 'black',
         order = "hclust", diag = FALSE,
         addCoef.col="black", addCoefasPercent = FALSE,
         p.mat = u_cor$p, sig.level=0.05, insig = "blank")
dev.off()

uf1 <- lmer(u_chosen ~ (scale(run_trial) + rt_lag_sc + rewFunc)^2 + 
              scale(run_trial)*scale(run) + (1|id/run), df)
vif(uf1)

# check: does entropy decay here? -- need to import the right entropy, value, etc.
ggplot(df, aes(run_trial, v_entropy_wi)) + geom_smooth(method = 'gam', formula = y ~ splines::ns(x, 4)) #+ facet_wrap(~run)
ggplot(df, aes(run_trial, v_entropy_wi)) + geom_smooth(method = 'loess')

# select learnable conditions, which are sampled across various runs
ldf <- df %>% filter(rewFunc== 'IEV' | rewFunc=='DEV')
fdf <- df %>% filter(run_trial>5)
edf <- df %>% filter(run_trial<6)
ggplot(df, aes(run_trial, u_chosen, color = rewFunc, lty = h_HippAntL_resp, group = interaction(rewFunc, h_HippAntL_resp))) + geom_smooth() #+ facet_wrap(~run)

ggplot(df %>% filter(run_trial>5), aes(run_trial, u_chosen, color = rewFunc, lty = h_HippAntL_resp, group = interaction(rewFunc, h_HippAntL_resp))) + geom_smooth() #+ facet_wrap(~run)

ggplot(df, aes(run_trial, rt_csv, color = rewFunc, lty = h_HippAntL_resp, group = interaction(rewFunc, h_HippAntL_resp))) + geom_smooth(method = 'gam', formula = y ~ splines::ns(x, 3)) #+ facet_wrap(~run)
ggplot(df, aes(run_trial, u_chosen, color = rt_lag>2000, lty = pe_f2_hipp_resp, group = interaction(rt_lag>2000, pe_f2_hipp_resp))) + geom_smooth(method = 'gam', formula = y ~ splines::ns(x, 3))  + facet_wrap(~rewFunc)
ggplot(df, aes(run_trial, rt_csv, color = rewFunc, lty = pe_f2_hipp_resp, group = interaction(rewFunc, pe_f2_hipp_resp))) + geom_smooth(method = 'gam', formula = y ~ splines::ns(x, 3)) + facet_wrap(~rewFunc)
ggplot(df, aes(run_trial, rt_csv, color = rewFunc, lty = h_HippAntL_resp, group = interaction(rewFunc, h_HippAntL_resp))) + geom_smooth(method = 'gam', formula = y ~ splines::ns(x, 3)) + facet_wrap(~rewFunc)
ggplot(df, aes(run_trial, rt_csv, color = rewFunc, lty = h_f1_fp_resp, group = interaction(rewFunc, h_f1_fp_resp))) + geom_smooth(method = 'gam', formula = y ~ splines::ns(x, 3)) + facet_wrap(~rewFunc)

# U vs. V
ggplot(df %>% filter(v_chosen>0), aes(v_chosen, u_chosen)) + geom_smooth() #+ facet_wrap(~run)
ggplot(df %>% filter(run_trial>5), aes(y_chosen, u_chosen)) + geom_smooth() #+ facet_wrap(~run)


# reasonable basis for builidng the model for interactions with betas
# NB: use run_trial, not -1/run_trial in U analyses
ub1 <- lmer(u_chosen ~ (scale(run_trial) + rt_lag_sc + rt_vmax_lag_sc + last_outcome + v_max_wi_lag + v_entropy_wi)^2 + 
              scale(u_chosen_lag) + 
              v_max_b + v_entropy_b + scale(run_trial)*scale(run) + (1|id/run), df)
screen.lmerTest(ub1,.01)
# # excluding the first 5 trials -- same results, no need to do that
# ub1f <- lmer(u_chosen ~ (scale(run_trial) + rt_lag_sc + rt_vmax_lag_sc + last_outcome + v_max_wi_lag + v_entropy_wi)^2 + 
#                last_outcome * scale(u_chosen_lag) + rt_vmax_lag_sc * scale(u_chosen_lag) + rt_lag_sc*scale(u_chosen_lag) + 
#                v_max_b + v_entropy_b + scale(run_trial)*scale(run) + (1|id/run), fdf)
# screen.lmerTest(ub1f,.01)
# ub1e <- lmer(u_chosen ~ (scale(run_trial) + rt_lag_sc + rt_vmax_lag_sc + last_outcome + v_max_wi_lag + v_entropy_wi)^2 + 
#                last_outcome * scale(u_chosen_lag) + rt_vmax_lag_sc * scale(u_chosen_lag) + rt_lag_sc*scale(u_chosen_lag) + 
#                v_max_b + v_entropy_b + scale(run_trial)*scale(run) + (1|id/run), edf)
# screen.lmerTest(ub1e,.001)


# cannot allow u_chosen_lag to interact with trial because of multi-collinearity
ub2 <- lmer(u_chosen_change ~ (scale(run_trial) + rt_lag_sc + rt_vmax_lag_sc + last_outcome + v_max_wi_lag + v_entropy_wi)^2 + 
              v_max_b + v_entropy_b + scale(run_trial)*scale(run) + scale(u_chosen_lag) + (1|id/run), df)
screen.lmerTest(ub2,.01)

# brute force approach to betas
# u_chosen looks more interpretable than u_chosen_change
ub3 <- lmer(u_chosen ~ (scale(run_trial) + rt_lag_sc + rt_vmax_lag_sc + last_outcome + v_entropy_wi + h_f1_fp)^3 +
              (scale(run_trial) + rt_lag_sc + rt_vmax_lag_sc + last_outcome + v_entropy_wi + h_HippAntL_neg)^3 +
              (scale(run_trial) + rt_lag_sc + rt_vmax_lag_sc + last_outcome + v_entropy_wi + pe_f1_cort_str)^3 +
              (scale(run_trial) + rt_lag_sc + rt_vmax_lag_sc + last_outcome + v_entropy_wi + pe_f2_hipp)^3 +
              scale(u_chosen_lag) + scale(run_trial)*scale(run)  +  (1|id/run), df)
screen.lmerTest(ub3, .01)

# unpack entropy* trial
ggplot(df %>% filter(!is.na(v_entropy_wi)), aes(run_trial, u_chosen, lty = v_entropy_wi>0, color = v_chosen>27)) + 
  geom_smooth(method = 'gam', formula = y ~ splines::ns(x,2)) #+ facet_wrap(~rewFunc)

ggplot(df %>% filter(!is.na(v_entropy_wi)), aes(run_trial, v_chosen, lty = v_entropy_wi>0, color = u_chosen>1341)) + 
  geom_smooth(method = 'gam', formula = y ~ splines::ns(x,3)) + facet_wrap(~rewFunc)


pdf('AH_entropy_uncertainty_aversion_by_cond.pdf', height = 6, width = 8)
ggplot(df %>% filter(!is.na(v_entropy_wi)), aes(run_trial, u_chosen, lty = v_entropy_wi>0, color = h_HippAntL_resp)) + 
  geom_smooth(method = 'gam', formula = y ~ splines::ns(x,2)) + facet_wrap(~rewFunc)
dev.off()

# 
# ub3a <- lmer(u_chosen_change ~ (scale(run_trial) + rt_lag_sc + rt_vmax_lag_sc + last_outcome + v_max_wi_lag + v_entropy_wi + h_f1_fp)^3 +
#                (scale(run_trial) + rt_lag_sc + rt_vmax_lag_sc + last_outcome + v_max_wi_lag + v_entropy_wi + h_HippAntL_neg)^3 +
#                (scale(run_trial) + rt_lag_sc + rt_vmax_lag_sc + last_outcome + v_max_wi_lag + v_entropy_wi + pe_f1_cort_str)^3 +
#                (scale(run_trial) + rt_lag_sc + rt_vmax_lag_sc + last_outcome + v_max_wi_lag + v_entropy_wi + pe_f2_hipp)^3 +
#                v_max_b + v_entropy_b + rt_lag_sc*scale(u_chosen_lag) + scale(run_trial)*scale(run) + (1|id/run), df)
# screen.lmerTest(ub3a, .01)
# 

# ub3e <- lmer(u_chosen ~ (scale(run_trial) + rt_lag_sc + rt_vmax_lag_sc + last_outcome + v_max_wi_lag + v_entropy_wi + h_f1_fp)^3 +
#               (scale(run_trial) + rt_lag_sc + rt_vmax_lag_sc + last_outcome + v_max_wi_lag + v_entropy_wi + I(-h_f2_neg_paralimb))^3 +
#               (scale(run_trial) + rt_lag_sc + rt_vmax_lag_sc + last_outcome + v_max_wi_lag + v_entropy_wi + pe_f1_cort_str)^3 +
#               (scale(run_trial) + rt_lag_sc + rt_vmax_lag_sc + last_outcome + v_max_wi_lag + v_entropy_wi + pe_f2_hipp)^3 +
#               v_max_b + v_entropy_b + rt_lag_sc*scale(u_chosen_lag) + scale(run_trial)*scale(run) + 
#               h_f1_fp*scale(u_chosen_lag) + I(-h_f2_neg_paralimb)*scale(u_chosen_lag) +pe_f1_cort_str*scale(u_chosen_lag) +pe_f2_hipp*scale(u_chosen_lag) + (1|id/run), edf)
# screen.lmerTest(ub3e, .01)

# plot the striking effect of HIPP on uncertainty sensitivity
p1 <- ggplot(df, aes(run_trial, u_chosen, lty = pe_f2_hipp_resp, group = pe_f2_hipp_resp)) + geom_smooth() #+ facet_wrap(~run)
p2 <- ggplot(df, aes(run_trial, u_chosen, lty = pe_f2_hipp_resp, group = interaction(run,pe_f2_hipp_resp) , color = run)) + geom_smooth() #+ facet_wrap(~run)
p3 <- ggplot(df, aes(run_trial, u_chosen, lty = h_HippAntL_resp, group = h_HippAntL_resp)) + geom_smooth()#+ facet_wrap(~run)
p4 <- ggplot(df, aes(run_trial, u_chosen, lty = h_HippAntL_resp, group = interaction(run, h_HippAntL_resp), color = run)) + geom_smooth() #+ facet_wrap(~run)
pdf("PH_AH_on_u_sensitivity.pdf", height = 8, width = 8)
ggarrange(p1,p2,p3,p4,ncol = 2, nrow = 2)
dev.off()

p1 <- ggplot(df, aes(run_trial, u_chosen, lty = pe_f1_cort_str_resp, group = pe_f1_cort_str_resp)) + geom_smooth() #+ facet_wrap(~run)
p2 <- ggplot(df, aes(run_trial, u_chosen, lty = pe_f2_hipp_resp, group = pe_f2_hipp_resp)) + geom_smooth() #+ facet_wrap(~run)
p3 <- ggplot(df, aes(run_trial, u_chosen, lty = h_f1_fp_resp, group = h_f1_fp_resp)) + geom_smooth() #+ facet_wrap(~run)
p4 <- ggplot(df, aes(run_trial, u_chosen, lty = h_HippAntL_resp, group = h_HippAntL_resp)) + geom_smooth() #+ facet_wrap(~run)

pdf("clusters_u_sensitivity.pdf", height = 4, width = 8)
ggarrange(p1,p2,p3, p4, ncol = 2, nrow = 2)
dev.off()

# 
# ub4 <- lmer(u_chosen_change ~ (scale(run_trial) + rt_lag_sc + rt_vmax_lag_sc + last_outcome + 
#                                  v_max_wi_lag + v_entropy_wi + scale(rt_vmax_change) +  
#                                  h_f1_fp + I(-h_f2_neg_paralimb) + pe_f1_cort_str + pe_f2_hipp)^2 + 
#               rt_lag_sc:last_outcome:h_f1_fp + rt_lag_sc:last_outcome:I(-h_f2_neg_paralimb) + 
#               rt_lag_sc:last_outcome:pe_f1_cort_str + rt_lag_sc:last_outcome:pe_f2_hipp +
#               rt_vmax_lag_sc:scale(run_trial):h_f1_fp + rt_vmax_lag_sc:scale(run_trial):I(-h_f2_neg_paralimb) + 
#               rt_vmax_lag_sc:scale(run_trial):pe_f1_cort_str + rt_vmax_lag_sc:scale(run_trial):pe_f2_hipp  +
#               v_max_b + v_entropy_b + scale(run_trial)*scale(run) + (1|id/run), df)
# summary(ub4)
# 

# simple models: do they swing in the direction of greater uncertainty?
us1 <- lmer(u_chosen ~ (scale(run_trial) + last_outcome)^2 + 
              scale(run_trial)*rewFunc + (1|id), df)
us1v <- lmer(u_chosen ~ (scale(run_trial) + last_outcome)^2 + 
              scale(run_trial)*rewFunc + scale(v_chosen) + (1|id), df)
# residualize u_chosen for v_chosen
uv1 <- lm(u_chosen ~ 1 + v_chosen + (1|id), df)
df$u_chosen_resid <- resid(uv1)
screen.lmerTest(us1)
us2 <- lmer(u_chosen ~ (scale(run_trial) + last_outcome + h_f1_fp)^2 +
              (scale(run_trial) + last_outcome+ h_HippAntL_neg)^2 +
              (scale(run_trial)+ last_outcome + pe_f1_cort_str)^2 +
            (scale(run_trial)+ last_outcome + pe_f2_hipp)^2 +
              scale(run_trial)*rewFunc + scale(u_chosen_lag) + (1|id), df)
screen.lmerTest(us2, .01)
us2r <- lmer(u_chosen_resid ~ (scale(run_trial) + last_outcome + h_f1_fp)^2 +
              (scale(run_trial) + last_outcome+ h_HippAntL_neg)^2 +
              (scale(run_trial)+ last_outcome + pe_f1_cort_str)^2 +
              (scale(run_trial)+ last_outcome + pe_f2_hipp)^2 +
              scale(run_trial)*rewFunc + scale(u_chosen_lag) + (1|id), df)
screen.lmerTest(us2r, .01)
# understand the time course of uncertainty residualized for value
ggplot(df, aes(run_trial, u_chosen_resid, color = rewFunc, lty = last_outcome)) + geom_smooth(method = 'gam', formula = y~splines::ns(x,3))
ggplot(df, aes(run_trial, u_chosen, color = rewFunc, lty = last_outcome)) + geom_smooth(method = 'gam', formula = y~splines::ns(x,3))
ggplot(df, aes(run_trial, v_chosen, color = rewFunc, lty = last_outcome)) + geom_smooth(method = 'gam', formula = y~splines::ns(x,3))

Anova(us2, '3')
p1 <- ggplot(df %>% filter(!is.na(last_outcome)), aes(pe_f2_hipp, u_chosen_change, color = last_outcome)) + geom_smooth(method = 'gam') #+ 
  geom_hline(yintercept = 1437.59)#+ facet_wrap(~run)
p2 <- ggplot(df %>% filter(!is.na(last_outcome)), aes(pe_f1_cort_str, u_chosen_change, color = last_outcome)) + geom_smooth(method = 'gam') #+ 
  geom_hline(yintercept = 1437.59)#+ facet_wrap(~run)
pdf("pe_clusters_u_reward.pdf", height = 4, width = 8)
ggarrange(p1,p2, ncol = 2, nrow = 1)
dev.off()
# control for value: this changes results
us3 <- lmer(u_chosen ~ (scale(run_trial) + last_outcome + h_f1_fp)^2 +
              (scale(run_trial) + last_outcome + h_HippAntL_neg)^2 +
              (scale(run_trial) + last_outcome + pe_f1_cort_str)^2 +
              (scale(run_trial) + last_outcome + pe_f2_hipp)^2 +
              scale(run_trial)*rewFunc + scale(v_chosen) + (1|id), df)
screen.lmerTest(us3, .05)
Anova(us3, '3')
anova(us1,us2,us3)

# is there a similar effect with v_chosen?
vs2 <- lmer(v_chosen ~ (scale(run_trial) + last_outcome + 
                          h_f1_fp + h_HippAntL_neg + pe_f1_cort_str + pe_f2_hipp)^2 + 
              scale(run_trial)*rewFunc + (1|id), df)
screen.lmerTest(vs2, .01)
p1 <- ggplot(df %>% filter(!is.na(last_outcome)), aes(pe_f2_hipp, v_chosen, color = last_outcome)) + geom_smooth(method = 'gam') + 
  geom_hline(yintercept = 27.62)#+ facet_wrap(~run)
p2 <- ggplot(df %>% filter(!is.na(last_outcome)), aes(pe_f1_cort_str, v_chosen, color = last_outcome)) + geom_smooth(method = 'gam') + 
  geom_hline(yintercept = 27.62)#+ facet_wrap(~run)
pdf("pe_clusters_v_reward.pdf", height = 4, width = 8)
ggarrange(p1,p2, ncol = 2, nrow = 1)
dev.off()

p1 <- ggplot(df %>% filter(!is.na(last_outcome)), aes(pe_f2_hipp, rt_csv-rt_lag, color = last_outcome)) + geom_smooth(method = 'gam') + 
  geom_hline(yintercept = 0)#+ facet_wrap(~rewFunc)
p2 <- ggplot(df %>% filter(!is.na(last_outcome)), aes(pe_f1_cort_str, rt_csv-rt_lag, color = last_outcome)) + geom_smooth(method = 'gam') + 
  geom_hline(yintercept = 0)#+ facet_wrap(~rewFunc)
p3 <- ggplot(df %>% filter(!is.na(last_outcome)), aes(-h_HippAntL, rt_csv-rt_lag, color = last_outcome)) + geom_smooth(method = 'gam') + 
  geom_hline(yintercept = 0)#+ facet_wrap(~rewFunc)
pdf("pe_ah_clusters_rt_reward.pdf", height = 8, width = 20)
ggarrange(p1,p2, p3, ncol = 3, nrow = 1)
dev.off()


# do the PH PE people RT-swing more post rewards?!!!

# predict RT with u_chosen and HIPP
urs1 <- lmer(rt_csv ~ (scale(run_trial) + rt_lag_sc + last_outcome)^2 + 
               scale(run_trial) * rewFunc +  (1|id/run), df)
screen.lmerTest(urs1, .01)
vif(urs1)
Anova(urs1)
urs2 <- lmer(rt_csv ~ (scale(run_trial) + rt_lag_sc + last_outcome + h_f1_fp)^3 +
               (scale(run_trial) + rt_lag_sc + last_outcome + h_HippAntL_neg)^3 + 
               (scale(run_trial) + rt_lag_sc + last_outcome + pe_f1_cort_str)^3 + 
               (scale(run_trial) + rt_lag_sc + last_outcome + pe_f2_hipp)^3 + 
               scale(run_trial) * rewFunc +
               (1|id/run), df)
screen.lmerTest(urs2, .01)
Anova(urs2)
# unpack PH*u_chosen
ggplot(df, aes(rt_lag, rt_csv, lty = last_outcome, color = pe_f2_hipp_resp)) + geom_smooth(method = "gam",
                                                                                            formula = y ~ splines::ns(x,2))
ggplot(df, aes(rt_lag, rt_csv, lty = last_outcome, color = h_f1_fp>0)) + geom_smooth(method = "gam",
                                                                                           formula = y ~ splines::ns(x,2))

# ideas for improving uncertainty analyses:
# try ML Cox survival with time-varying within-trial U
# look at relative rather than absolute uncertainty of the choice

### plot single subject
sdf <- df[df$id==10811 & df$emotion == "scram" & (df$rewFunc=="IEV" | df$rewFunc=="DEV"),]
pdf("single_subject_swings_10811.pdf", height = 3, width = 5)
ggplot(sdf, aes(run_trial, rt_csv/1000, color = ev)) + geom_point(size = 4) + facet_wrap(~rewFunc) + 
  xlab("Trial") + ylab("Response time, s") + scale_color_viridis_c(option = "plasma", name = "Expected value") + theme_dark()
dev.off()
#
save(file = 'vhd_u_meg_models.Rdata', list = ls(all.names = TRUE))
# load('vhd_models.Rdata')
