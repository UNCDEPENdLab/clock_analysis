library(dplyr)
library(tidyverse)
library(psych)
library(corrplot)
library(lme4)
# library(lmerTest)
library(stargazer)
# source('~/code/Rhelpers/')
setwd('~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/')
load('trial_df_and_vhdkfpe_clusters.Rdata')
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
mf1 <- lmer(rt_csv ~ (scale(-1/run_trial) + scale(rt_lag) + omission_lag )^2  + scale(rt_lag) + (1|id/run), df)
#screen.lmerTest(mf1)
mf2h <- lmer(rt_csv ~ (scale(-1/run_trial) + scale(rt_lag) +   omission_lag  + h_f1_fp + I(-h_f2_neg_paralimb))^2 + 
               scale(rt_lag):omission_lag:h_f1_fp + scale(rt_lag):omission_lag:I(-h_f2_neg_paralimb) +
               (1|id/run), df)
#screen.lmerTest(mf2h)
mf2v <- lmer(rt_csv ~ (scale(-1/run_trial) + scale(rt_lag) +   omission_lag  + I(-v_f1_neg_cog) + v_f2_paralimb)^2 + 
               scale(rt_lag):omission_lag:I(-v_f1_neg_cog) + scale(rt_lag):omission_lag:v_f2_paralimb +
               (1|id/run), df)
#screen.lmerTest(mf2v)
# # add KLD -- prefronto-striatal makes them a bit faster, but does not interact with anything
# mf2k <- lmer(rt_csv ~ (scale(-1/run_trial) + scale(rt_lag) +   omission_lag  + k_f1_IPL_ventr_stream + k_f2_prefrontal_bg)^2 +   (1|id/run), df)
# #screen.lmerTest(mf2k)

# KLD at feedback: kf_f2_vmpfc_precun favors exploitation, but fit is comparatively poor
mf2kf <- lmer(rt_csv ~ (scale(-1/run_trial) + scale(rt_lag) + omission_lag  + kf_f1_pos + I(-kf_f2_vmpfc_precun))^2 +
                scale(rt_lag):omission_lag:kf_f1_pos + scale(rt_lag):omission_lag:I(-kf_f2_vmpfc_precun) +
              (1|id/run), df)
#screen.lmerTest(mf2kf)

# entropy change -- nothing spectacular, but vmpfc-precuneus promotes exploitation
# the inter-correlation between dh_f1_co_bg and dh_f2_dan is likely a problem (r=.73)
mf2dh <- lmer(rt_csv ~ (scale(-1/run_trial) + scale(rt_lag) +   omission_lag  + dh_f1_co_bg + dh_f2_dan + I(-dh_f_neg_vmpfc_precun))^2 +
                scale(rt_lag):omission_lag:dh_f1_co_bg + scale(rt_lag):omission_lag:dh_f2_dan + scale(rt_lag):omission_lag:I(-dh_f_neg_vmpfc_precun) +
                (1|id/run), df)
#screen.lmerTest(mf2dh)

# set DHP aside for now, too many factors
mf2dhp <- lmer(rt_csv ~ (scale(-1/run_trial) + scale(rt_lag) +   omission_lag  + dhp_f1_all)^2 +  
                 scale(rt_lag):omission_lag:dhp_f1_all + (1|id/run), df)
#screen.lmerTest(mf2dhp)

# PE
mf2pe <- lmer(rt_csv ~ (scale(-1/run_trial) + scale(rt_lag) +   omission_lag  + pe_f1_cort_str + pe_f2_hipp)^2 + 
                scale(rt_lag):omission_lag:pe_f1_cort_str + scale(rt_lag):omission_lag:pe_f2_hipp +
                (1|id/run), df)
#screen.lmerTest(mf2pe)

# but the two PE factors are inter-correlated at r=.52
# post. hippocampal effects are even stronger after taking the first factor out!
mf2pef2only <- lmer(rt_csv ~ (scale(-1/run_trial) + scale(rt_lag) +   omission_lag  + pe_f2_hipp)^2 + 
                scale(rt_lag):omission_lag:pe_f2_hipp +
                (1|id/run), df)
#screen.lmerTest(mf2pef2only)

ggplot(df,aes(rt_csv,rt_lag, lty = pe_f2_hipp_resp, color = omission_lag)) + geom_smooth(method = 'glm')

# many effects, poor fit (too many clusters?)
mf2d <- lmer(rt_csv ~ (scale(-1/run_trial) + scale(rt_lag) +   omission_lag + d_f1_FP_SMA + d_f2_VS + d_f3_ACC_ins)^2 + 
               scale(rt_lag):omission_lag:d_f1_FP_SMA + scale(rt_lag):omission_lag:d_f2_VS + scale(rt_lag):omission_lag:d_f3_ACC_ins + 
               (1|id/run), df)
#screen.lmerTest(mf2d)

#################
# best model
mf3hpedh3 <-  lmer(rt_csv ~ (scale(-1/run_trial) + scale(rt_lag) +   omission_lag + 
                    h_f1_fp + I(-h_f2_neg_paralimb) + I(-dh_f_neg_vmpfc_precun) + pe_f1_cort_str + pe_f2_hipp)^2 + 
                     scale(rt_lag):omission_lag:h_f1_fp + scale(rt_lag):omission_lag:I(-h_f2_neg_paralimb) + 
                     scale(rt_lag):omission_lag:I(-dh_f_neg_vmpfc_precun) +
                     scale(rt_lag):omission_lag:pe_f1_cort_str + scale(rt_lag):omission_lag:pe_f2_hipp +
                     scale(rt_lag):omission_lag:I(-dh_f_neg_vmpfc_precun) +
                    (1|id/run), df)
#screen.lmerTest(mf3hpedh3)

# swap in hippocampal low-H cluster for vmpfcHipp factor -- fit is worse
mf3hpedh3hipp <-  lmer(rt_csv ~ (scale(-1/run_trial) + scale(rt_lag) +   omission_lag + 
                               h_f1_fp + I(-h_HippAntL) + I(-dh_f_neg_vmpfc_precun) + pe_f1_cort_str + pe_f2_hipp)^2 + 
                     scale(rt_lag):omission_lag:h_f1_fp + scale(rt_lag):omission_lag:I(-h_HippAntL) + 
                     scale(rt_lag):omission_lag:I(-dh_f_neg_vmpfc_precun) +
                     scale(rt_lag):omission_lag:pe_f1_cort_str + scale(rt_lag):omission_lag:pe_f2_hipp +
                     scale(rt_lag):omission_lag:I(-dh_f_neg_vmpfc_precun) +
                     (1|id/run), df)
#screen.lmerTest(mf3hpedh3hipp)

stargazer(mf1, mf2h, mf2v, mf2pe,mf2dh,mf3hpedh3, type="html", out="mf.htm", report = "vcs*",
          digits = 1,single.row=TRUE,omit.stat = "bic",
          column.labels = c("Null", "Value", "Prediction error", "Entropy", "Entropy change", "Best model"),
          star.char = c("*", "**", "***"),
          star.cutoffs = c(0.05, 0.01, 0.001),
          notes = c("* p<0.05; ** p<0.01; *** p<0.001"), 
          notes.append = F)


anova(mf1,mf2h,mf2v,mf2d,mf2pe, mf2pef2only,mf2kf, mf2dh,mf2dhp, mf3hpedh3, mf3hpedh3hipp)

##################

# Model-based RT models -- an improvement of 4000 AIC points over MF
vv <- df[,c("rt_csv","rt_lag", "rt_vmax", "rt_vmax_lag", "rt_vmax_change", "run_trial", "omission_lag", "v_max_wi", "v_max_wi_lag", "abs_pe_max_lag", "v_entropy_wi")]
w_cor <- corr.test(vv)
pdf("within_sub_predictor_corr.pdf", width=12, height=12)
corrplot(w_cor$r, cl.lim=c(-1,1),
         method = "circle", tl.cex = 1.5, type = "upper", tl.col = 'black',
         order = "hclust", diag = FALSE,
         addCoef.col="black", addCoefasPercent = FALSE,
         p.mat = w_cor$p, sig.level=0.05, insig = "blank")
dev.off()
# NB: omissions have two effects: RT shortening and non-directional RT swing
# collinearity remediated by taking lags of rt_vmax and v_max, conceptually the same variables, just discard the effect of last outcome
# staying away from brain*brain interactions, too complex
df <- df[!is.na(df$rt_vmax_change),]
mb1 <- lmer(rt_csv ~ (scale(-1/run_trial) + scale(rt_lag) + scale(rt_vmax_lag) + omission_lag + v_max_wi_lag + v_entropy_wi)^2 + 
              v_max_b + v_entropy_b + (1|id/run), df)
#screen.lmerTest(mb1)

# mb1pe <- lmer(rt_csv ~ (scale(-1/run_trial) + scale(rt_lag) + scale(rt_vmax_lag) + omission_lag + v_max_wi_lag + abs(pe_max_lag) + v_entropy_wi)^2 + v_max_b + v_entropy_b + scale(rt_lag) + (1|id/run), df)
# #screen.lmerTest(mb1pe)

mb1a <- lmer(rt_csv ~ (scale(-1/run_trial) + scale(rt_lag) + scale(rt_vmax_lag) + v_max_wi_lag + rt_vmax_change + v_entropy_wi)^2 + 
               v_max_b + v_entropy_b + (1|id/run), df)
#screen.lmerTest(mb1a)

mb2h <- lmer(rt_csv ~ (scale(-1/run_trial) + scale(rt_lag) + scale(rt_vmax_lag) + omission_lag + v_max_wi_lag + v_entropy_wi +rt_vmax_change +  h_f1_fp + I(-h_f2_neg_paralimb))^2 + 
                      scale(rt_lag):omission_lag:h_f1_fp + scale(rt_lag):omission_lag:I(-h_f2_neg_paralimb) +
                      scale(rt_vmax_lag):v_max_wi_lag:h_f1_fp + scale(rt_vmax_lag):v_max_wi_lag:I(-h_f2_neg_paralimb) +
                      scale(-1/run_trial):scale(rt_vmax_lag):h_f1_fp +  scale(-1/run_trial):scale(rt_vmax_lag):I(-h_f2_neg_paralimb) +
                      v_max_b + v_entropy_b +  (1|id/run), df)
#screen.lmerTest(mb2h)
mb2v <- lmer(rt_csv ~ (scale(-1/run_trial) + scale(rt_lag) + scale(rt_vmax_lag) + omission_lag + v_max_wi_lag + v_entropy_wi + rt_vmax_change + I(-v_f1_neg_cog) + v_f2_paralimb)^2 + 
               scale(rt_lag):omission_lag:I(-v_f1_neg_cog) + scale(rt_lag):omission_lag:v_f2_paralimb +
               scale(rt_vmax_lag):v_max_wi_lag:I(-v_f1_neg_cog) + scale(rt_vmax_lag):v_max_wi_lag:v_f2_paralimb +
               scale(-1/run_trial):scale(rt_vmax_lag):I(-v_f1_neg_cog) +  scale(-1/run_trial):scale(rt_vmax_lag):v_f2_paralimb +
               v_max_b + v_entropy_b +   (1|id/run), df)
#screen.lmerTest(mb2v)
# # add KLD -- prefronto-striatal makes them a bit faster, but does not interact with anything
# mb2k <- lmer(rt_csv ~ (scale(-1/run_trial) + scale(rt_lag) + scale(rt_vmax_lag) + omission_lag + v_max_wi_lag + v_entropy_wi + k_f1_IPL_ventr_stream + k_f2_prefrontal_bg)^2 + v_max_b + v_entropy_b +   (1|id/run), df)
# #screen.lmerTest(mb2k)

# KLD at feedback: kf_f2_vmpfc_precun favors exploitation, but fit is comparatively poor
# nothing remarkable for kf_f1_pos
mb2kf <- lmer(rt_csv ~ (scale(-1/run_trial) + scale(rt_lag) + scale(rt_vmax_lag) + omission_lag + v_max_wi_lag + v_entropy_wi + rt_vmax_change + kf_f1_pos + I(-kf_f2_vmpfc_precun))^2 +
  scale(rt_lag):omission_lag:kf_f1_pos + scale(rt_lag):omission_lag:I(-kf_f2_vmpfc_precun) +
  scale(rt_vmax_lag):v_max_wi_lag:kf_f1_pos + scale(rt_vmax_lag):v_max_wi_lag:I(-kf_f2_vmpfc_precun) +
  scale(-1/run_trial):scale(rt_vmax_lag):kf_f1_pos + scale(-1/run_trial):scale(rt_vmax_lag):I(-kf_f2_vmpfc_precun) +
  v_max_b + v_entropy_b +  (1|id/run), df)
#screen.lmerTest(mb2kf)

# entropy change -- nothing spectacular, but vmpfc-precuneus promotes exploitation 
mb2dh <- lmer(rt_csv ~ (scale(-1/run_trial) + scale(rt_lag) + scale(rt_vmax_lag) + omission_lag + v_max_wi_lag + v_entropy_wi + rt_vmax_change + dh_f1_co_bg + dh_f2_dan + I(-dh_f_neg_vmpfc_precun))^2 +
                scale(rt_lag):omission_lag:dh_f1_co_bg + scale(rt_lag):omission_lag:dh_f2_dan + scale(rt_lag):omission_lag:I(-dh_f_neg_vmpfc_precun) +
                scale(rt_vmax_lag):v_max_wi_lag:dh_f1_co_bg + scale(rt_vmax_lag):v_max_wi_lag:scale(rt_vmax_lag):v_max_wi_lag:dh_f2_dan + scale(rt_vmax_lag):v_max_wi_lag:I(-dh_f_neg_vmpfc_precun) +
                scale(-1/run_trial):scale(rt_vmax_lag):dh_f1_co_bg + scale(-1/run_trial):scale(rt_vmax_lag):dh_f2_dan + scale(-1/run_trial):scale(rt_vmax_lag):I(-dh_f_neg_vmpfc_precun) +
                v_max_b + v_entropy_b +   (1|id/run), df)
#screen.lmerTest(mb2dh)
# intercorrelated factors, test separately:
mb2dh1 <- lmer(rt_csv ~ (scale(-1/run_trial) + scale(rt_lag) + scale(rt_vmax_lag) + omission_lag + v_max_wi_lag + v_entropy_wi + rt_vmax_change + dh_f1_co_bg + I(-dh_f_neg_vmpfc_precun))^2 +
                scale(rt_lag):omission_lag:dh_f1_co_bg + scale(rt_lag):omission_lag:I(-dh_f_neg_vmpfc_precun) +
                scale(rt_vmax_lag):v_max_wi_lag:dh_f1_co_bg + scale(rt_vmax_lag):v_max_wi_lag:I(-dh_f_neg_vmpfc_precun) +
                scale(-1/run_trial):scale(rt_vmax_lag):dh_f1_co_bg + scale(-1/run_trial):scale(rt_vmax_lag):I(-dh_f_neg_vmpfc_precun) +
                v_max_b + v_entropy_b +   (1|id/run), df)
#screen.lmerTest(mb2dh1)
# dh_f2_dan dampens the effect of rt_vmax_change (there was also a small effect of v2 above)!
mb2dh2 <- lmer(rt_csv ~ (scale(-1/run_trial) + scale(rt_lag) + scale(rt_vmax_lag) + omission_lag + v_max_wi_lag + v_entropy_wi + rt_vmax_change + dh_f2_dan + I(-dh_f_neg_vmpfc_precun))^2 +
                scale(rt_lag):omission_lag:dh_f2_dan + scale(rt_lag):omission_lag:I(-dh_f_neg_vmpfc_precun) +
                scale(rt_vmax_lag):v_max_wi_lag:scale(rt_vmax_lag):v_max_wi_lag:dh_f2_dan + scale(rt_vmax_lag):v_max_wi_lag:I(-dh_f_neg_vmpfc_precun) +
                scale(-1/run_trial):scale(rt_vmax_lag):dh_f2_dan + scale(-1/run_trial):scale(rt_vmax_lag):I(-dh_f_neg_vmpfc_precun) +
                v_max_b + v_entropy_b +   (1|id/run), df)
#screen.lmerTest(mb2dh2)


# dhp_f1_all
mb2dhp <- lmer(rt_csv ~ (scale(-1/run_trial) + scale(rt_lag) + scale(rt_vmax_lag) + omission_lag + v_max_wi_lag + v_entropy_wi + rt_vmax_change + dhp_f1_all + I(-dh_f_neg_vmpfc_precun))^2 +
                scale(rt_lag):omission_lag:dhp_f1_all + scale(rt_lag):omission_lag:I(-dh_f_neg_vmpfc_precun) +
                scale(rt_vmax_lag):v_max_wi_lag:dhp_f1_all + scale(rt_vmax_lag):v_max_wi_lag:I(-dh_f_neg_vmpfc_precun) +
                scale(-1/run_trial):scale(rt_vmax_lag):dhp_f1_all + scale(-1/run_trial):scale(rt_vmax_lag):I(-dh_f_neg_vmpfc_precun) +
                v_max_b + v_entropy_b +   (1|id/run), df)
#screen.lmerTest(mb2dhp)

# set DHP aside for now, too many factors
# mb2dhp <- lmer(rt_csv ~ (scale(-1/run_trial) + scale(rt_lag) + scale(rt_vmax_lag) + omission_lag + v_max_wi_lag + v_entropy_wi + 
#                           dhp_f1_dlpfc_r + dhp_f2_str + dhp_f3_ventr_stream_cerebell + dhp_f4_prefront_l)^2 + v_max_b + v_entropy_b +   (1|id/run), df)
# #screen.lmerTest(mb2dhp)

# PE
mb2pe <- lmer(rt_csv ~ (scale(-1/run_trial) + scale(rt_lag) + scale(rt_vmax_lag) + omission_lag + v_max_wi_lag + v_entropy_wi + rt_vmax_change + pe_f1_cort_str + pe_f2_hipp)^2 + 
                scale(rt_lag):omission_lag:pe_f1_cort_str + scale(rt_lag):omission_lag:pe_f2_hipp +
                scale(rt_vmax_lag):v_max_wi_lag:pe_f1_cort_str + scale(rt_vmax_lag):v_max_wi_lag:pe_f2_hipp +
                scale(-1/run_trial):scale(rt_vmax_lag):pe_f1_cort_str + scale(-1/run_trial):scale(rt_vmax_lag):pe_f2_hipp +
                v_max_b + v_entropy_b +   (1|id/run), df)
#screen.lmerTest(mb2pe)

mb2pe_hipp <- lmer(rt_csv ~ (scale(-1/run_trial) + scale(rt_lag) + scale(rt_vmax_lag) + omission_lag + v_max_wi_lag + v_entropy_wi + rt_vmax_change +  pe_f2_hipp)^2 + 
                scale(rt_lag):omission_lag:pe_f2_hipp +
                scale(rt_vmax_lag):v_max_wi_lag:pe_f2_hipp +
                scale(-1/run_trial):scale(rt_vmax_lag):pe_f2_hipp +
                v_max_b + v_entropy_b + (1|id/run), df)
screen.lmerTest(mb2pe_hipp)

# add h_HippAntL
mb2hipp_pa <- lmer(rt_csv ~ (scale(-1/run_trial) + scale(rt_lag) + scale(rt_vmax_lag) + omission_lag + v_max_wi_lag + v_entropy_wi + rt_vmax_change + h_HippAntL + pe_f2_hipp)^2 + 
                     scale(rt_lag):omission_lag:h_HippAntL +
                     scale(rt_lag):omission_lag:pe_f2_hipp +
                     scale(rt_vmax_lag):v_max_wi_lag:h_HippAntL +
                     scale(rt_vmax_lag):v_max_wi_lag:pe_f2_hipp +
                     scale(-1/run_trial):scale(rt_vmax_lag):h_HippAntL +
                     scale(-1/run_trial):scale(rt_vmax_lag):pe_f2_hipp +
                     v_max_b + v_entropy_b + (1|id/run), df)
screen.lmerTest(mb2hipp_pa)


# # add trial-wise PEmax -- not feasible
# mb2pe2 <- lmer(rt_csv ~ (scale(-1/run_trial) + scale(rt_lag) + scale(rt_vmax_lag) + omission_lag + v_max_wi_lag + v_entropy_wi + abs_pe_max_lag +  pe_f1_cort_str + pe_f2_hipp)^3 + 
#                 scale(rt_lag):omission_lag:pe_f1_cort_str + scale(rt_lag):omission_lag:pe_f2_hipp +
#                 scale(rt_vmax_lag):v_max_wi_lag:pe_f1_cort_str + scale(rt_vmax_lag):v_max_wi_lag:pe_f2_hipp +
#                 scale(-1/run_trial):scale(rt_vmax_lag):pe_f1_cort_str + scale(-1/run_trial):scale(rt_vmax_lag):pe_f2_hipp +
#                 v_max_b + v_entropy_b +   (1|id/run), df)
# #screen.lmerTest(mb2pe2)


# many effects, poor fit (too many clusters?)
mb2d <- lmer(rt_csv ~ (scale(-1/run_trial) + scale(rt_lag) + scale(rt_vmax_lag) + omission_lag + v_max_wi_lag  + v_entropy_wi+ rt_vmax_change + d_f1_FP_SMA + d_f2_VS + d_f3_ACC_ins)^2 + 
               scale(rt_lag):omission_lag:d_f1_FP_SMA + scale(rt_lag):omission_lag:d_f2_VS + scale(rt_lag):omission_lag:d_f3_ACC_ins + 
               scale(rt_vmax_lag):v_max_wi_lag:d_f1_FP_SMA + scale(rt_vmax_lag):v_max_wi_lag:d_f2_VS +scale(rt_vmax_lag):v_max_wi_lag:d_f3_ACC_ins + 
               scale(-1/run_trial):scale(rt_vmax_lag):d_f1_FP_SMA + scale(-1/run_trial):scale(rt_vmax_lag):d_f2_VS + scale(-1/run_trial):scale(rt_vmax_lag):d_f3_ACC_ins +
               v_max_b + v_entropy_b +  (1|id/run), df)
#screen.lmerTest(mb2d)


#################
# second-best model
mb3hpedh3 <-  lmer(rt_csv ~ (scale(-1/run_trial) + scale(rt_lag) + scale(rt_vmax_lag) + omission_lag + 
                     v_max_wi_lag + v_entropy_wi +rt_vmax_change +  h_f1_fp + I(-h_f2_neg_paralimb) + I(-dh_f_neg_vmpfc_precun) + pe_f1_cort_str + pe_f2_hipp)^2 + 
                     scale(rt_lag):omission_lag:h_f1_fp + scale(rt_lag):omission_lag:I(-h_f2_neg_paralimb) + 
                     scale(rt_lag):omission_lag:I(-dh_f_neg_vmpfc_precun) +
                     scale(rt_lag):omission_lag:pe_f1_cort_str + scale(rt_lag):omission_lag:pe_f2_hipp +
                     scale(rt_lag):omission_lag:I(-dh_f_neg_vmpfc_precun) +
                     scale(rt_vmax_lag):scale(-1/run_trial):h_f1_fp + scale(rt_vmax_lag):scale(-1/run_trial):I(-h_f2_neg_paralimb) + 
                     scale(rt_vmax_lag):scale(-1/run_trial):pe_f1_cort_str + scale(rt_vmax_lag):scale(-1/run_trial):pe_f2_hipp  +
                     scale(rt_vmax_lag):scale(-1/run_trial):I(-dh_f_neg_vmpfc_precun) +
                     v_max_b + v_entropy_b + (1|id/run), df)
#screen.lmerTest(mb3hpedh3, .05)

# add in DHP -- best model
mb3hpedh_p4 <-  lmer(rt_csv ~ (scale(-1/run_trial) + scale(rt_lag) + scale(rt_vmax_lag) + omission_lag + 
                               v_max_wi_lag + v_entropy_wi + scale(rt_vmax_change) +  
                                 h_f1_fp + I(-h_f2_neg_paralimb) + I(-dh_f_neg_vmpfc_precun) + pe_f1_cort_str + pe_f2_hipp + dhp_f1_all)^2 + 
                     scale(rt_lag):omission_lag:h_f1_fp + scale(rt_lag):omission_lag:I(-h_f2_neg_paralimb) + scale(rt_lag):omission_lag:dhp_f1_all +
                     scale(rt_lag):omission_lag:I(-dh_f_neg_vmpfc_precun) +
                     scale(rt_lag):omission_lag:pe_f1_cort_str + scale(rt_lag):omission_lag:pe_f2_hipp +
                     scale(rt_lag):omission_lag:I(-dh_f_neg_vmpfc_precun) +
                     scale(rt_vmax_lag):scale(-1/run_trial):h_f1_fp + scale(rt_vmax_lag):scale(-1/run_trial):I(-h_f2_neg_paralimb) + 
                     scale(rt_vmax_lag):scale(-1/run_trial):pe_f1_cort_str + scale(rt_vmax_lag):scale(-1/run_trial):pe_f2_hipp  +
                     scale(rt_vmax_lag):scale(-1/run_trial):I(-dh_f_neg_vmpfc_precun) + scale(rt_vmax_lag):scale(-1/run_trial):dhp_f1_all +
                     v_max_b + v_entropy_b + (1|id/run), df)
#screen.lmerTest(mb3hpedh_p4, .05)

# remove pe_f1 to make sure pe_f2_hipp effects are not driven by multicollinearity
mb3hpedh_p5 <-  lmer(rt_csv ~ (scale(-1/run_trial) + scale(rt_lag) + scale(rt_vmax_lag) + omission_lag + 
                                 v_max_wi_lag + v_entropy_wi +rt_vmax_change +  
                                 h_f1_fp + I(-h_f2_neg_paralimb) + I(-dh_f_neg_vmpfc_precun) +  pe_f2_hipp + dhp_f1_all)^2 + 
                       scale(rt_lag):omission_lag:h_f1_fp + scale(rt_lag):omission_lag:I(-h_f2_neg_paralimb) + scale(rt_lag):omission_lag:dhp_f1_all +
                       scale(rt_lag):omission_lag:I(-dh_f_neg_vmpfc_precun) +
                       scale(rt_lag):omission_lag:pe_f2_hipp +
                       scale(rt_lag):omission_lag:I(-dh_f_neg_vmpfc_precun) +
                       scale(rt_vmax_lag):scale(-1/run_trial):h_f1_fp + scale(rt_vmax_lag):scale(-1/run_trial):I(-h_f2_neg_paralimb) + 
                       scale(rt_vmax_lag):scale(-1/run_trial):pe_f2_hipp  +
                       scale(rt_vmax_lag):scale(-1/run_trial):I(-dh_f_neg_vmpfc_precun) + scale(rt_vmax_lag):scale(-1/run_trial):dhp_f1_all +
                       v_max_b + v_entropy_b + (1|id/run), df)
#screen.lmerTest(mb3hpedh_p5, .05)

# test with current rt_vmax instead of change
mb3hpedh_p_cur <-  lmer(rt_csv ~ (scale(-1/run_trial) + scale(rt_lag) + scale(rt_vmax) + omission_lag + 
                                 v_max_wi_lag + v_entropy_wi +
                                 h_f1_fp + I(-h_f2_neg_paralimb) + I(-dh_f_neg_vmpfc_precun) + pe_f1_cort_str + pe_f2_hipp + dhp_f1_all)^2 + 
                       scale(rt_lag):omission_lag:h_f1_fp + scale(rt_lag):omission_lag:I(-h_f2_neg_paralimb) + scale(rt_lag):omission_lag:dhp_f1_all +
                       scale(rt_lag):omission_lag:I(-dh_f_neg_vmpfc_precun) +
                       scale(rt_lag):omission_lag:pe_f1_cort_str + scale(rt_lag):omission_lag:pe_f2_hipp +
                       scale(rt_lag):omission_lag:I(-dh_f_neg_vmpfc_precun) +
                       scale(rt_vmax):scale(-1/run_trial):h_f1_fp + scale(rt_vmax_lag):scale(-1/run_trial):I(-h_f2_neg_paralimb) + 
                       scale(rt_vmax):scale(-1/run_trial):pe_f1_cort_str + scale(rt_vmax_lag):scale(-1/run_trial):pe_f2_hipp  +
                       scale(rt_vmax):scale(-1/run_trial):I(-dh_f_neg_vmpfc_precun) + scale(rt_vmax_lag):scale(-1/run_trial):dhp_f1_all +
                       v_max_b + v_entropy_b + (1|id/run), df)
#screen.lmerTest(mb3hpedh_p_cur, .05)

# remove pe_f1 to make sure pe_f2_hipp effects are not driven by multicollinearity
mb3hpedh_p_cur2 <-  lmer(rt_csv ~ (scale(-1/run_trial) + scale(rt_lag) + scale(rt_vmax) + omission_lag + 
                                 v_max_wi_lag + v_entropy_wi +
                                 h_f1_fp + I(-h_f2_neg_paralimb) + I(-dh_f_neg_vmpfc_precun) +  pe_f2_hipp + dhp_f1_all)^2 + 
                       scale(rt_lag):omission_lag:h_f1_fp + scale(rt_lag):omission_lag:I(-h_f2_neg_paralimb) + scale(rt_lag):omission_lag:dhp_f1_all +
                       scale(rt_lag):omission_lag:I(-dh_f_neg_vmpfc_precun) +
                       scale(rt_lag):omission_lag:pe_f2_hipp +
                       scale(rt_lag):omission_lag:I(-dh_f_neg_vmpfc_precun) +
                       scale(rt_vmax_lag):scale(-1/run_trial):h_f1_fp + scale(rt_vmax_lag):scale(-1/run_trial):I(-h_f2_neg_paralimb) + 
                       scale(rt_vmax_lag):scale(-1/run_trial):pe_f2_hipp  +
                       scale(rt_vmax_lag):scale(-1/run_trial):I(-dh_f_neg_vmpfc_precun) + scale(rt_vmax_lag):scale(-1/run_trial):dhp_f1_all +
                       v_max_b + v_entropy_b + (1|id/run), df)
#screen.lmerTest(mb3hpedh_p_cur2, .05)

# and test the rt_vmax*omission*ROI interaction to see if they get more dislodged
mb3hpedh_p_cur3 <-  lmer(rt_csv ~ (scale(-1/run_trial) + scale(rt_lag) + scale(rt_vmax) + omission_lag + 
                           v_max_wi_lag + v_entropy_wi +
                           h_f1_fp + I(-h_f2_neg_paralimb) + I(-dh_f_neg_vmpfc_precun) + pe_f1_cort_str + pe_f2_hipp + dhp_f1_all)^2 + 
                           scale(rt_lag):omission_lag:h_f1_fp + scale(rt_lag):omission_lag:I(-h_f2_neg_paralimb) + scale(rt_lag):omission_lag:dhp_f1_all +
                           scale(rt_lag):omission_lag:I(-dh_f_neg_vmpfc_precun) +
                           scale(rt_lag):omission_lag:pe_f1_cort_str + scale(rt_lag):omission_lag:pe_f2_hipp +
                           scale(rt_lag):omission_lag:I(-dh_f_neg_vmpfc_precun) +
                           scale(rt_vmax):scale(-1/run_trial):h_f1_fp + scale(rt_vmax_lag):scale(-1/run_trial):I(-h_f2_neg_paralimb) + 
                           scale(rt_vmax):scale(-1/run_trial):pe_f1_cort_str  + scale(rt_vmax):scale(-1/run_trial):pe_f2_hipp  +
                           scale(rt_vmax):scale(-1/run_trial):I(-dh_f_neg_vmpfc_precun) + scale(rt_vmax):scale(-1/run_trial):dhp_f1_all +
                           scale(rt_vmax):omission_lag:h_f1_fp + scale(rt_vmax):omission_lag:I(-h_f2_neg_paralimb) + scale(rt_vmax):omission_lag:dhp_f1_all +
                           scale(rt_vmax):omission_lag:I(-dh_f_neg_vmpfc_precun) +
                           scale(rt_vmax):omission_lag:pe_f1_cort_str + scale(rt_vmax):omission_lag:pe_f2_hipp +
                           scale(rt_vmax):omission_lag:I(-dh_f_neg_vmpfc_precun) +
                           v_max_b + v_entropy_b + (1|id/run), df)
#screen.lmerTest(mb3hpedh_p_cur3, .05)


##################
anova(mb1a,mb2h,mb2v,mb2d,mb2pe,mb2kf, mb2dh, mb2dhp, mb3hpedh3, mb3hpedh_p4, mb3hpedh_p_cur2, mb3hpedh_p_cur3)
# at the end of the day, the h clusters explain the most
# pe is better than KLD or entropy change

stargazer(mb1a,mb2v,mb2pe,mb2h, mb2dh, mb3hpedh3, type="html", out="mb.htm", report = "vcs*",
          digits = 1,single.row=TRUE,omit.stat = "bic",
          column.labels = c("Null", "Value", "Prediction error", "Entropy", "Entropy change", "Best model"),
          star.char = c("*", "**", "***"),
          star.cutoffs = c(0.05, 0.01, 0.001),
          notes = c("* p<0.05; ** p<0.01; *** p<0.001"), 
          notes.append = F)


# understand rt_vmax_change effect
ggplot(df, aes(rt_vmax_change, rt_csv, color = pe_f2_hipp_resp)) + geom_smooth(method = "glm")

ggplot(df, aes(rt_csv, color = pe_f2_hipp_resp)) + geom_smooth(method = "glm")


###########
## PLOTS ##
#####
## "model-based analyses": plots

# behavioral exegisis on Dauc betas
# crescendo contributions to RT swings

# analogously, for KLD:
ggplot(df,aes(run_trial,rt_csv, color = k_f1_IPL_ventr_stream_resp, lty = k_f2_prefrontal_bg_resp)) + geom_smooth() + facet_wrap(~rewFunc)
ggplot(df,aes(run_trial,rt_csv, color = k_f1_IPL_ventr_stream_resp)) + geom_smooth() + facet_wrap(~rewFunc)
ggplot(df,aes(run_trial,rt_csv, lty = k_f2_prefrontal_bg_resp)) + geom_smooth() + facet_wrap(~rewFunc)

# analogously, for KLD at feedback
# not seeing much
ggplot(df,aes(run_trial,rt_csv, color = kf_f1_pos_resp, lty = kf_f3_str_front_ins_resp, size = kf_f2_vmpfc_precun_resp)) + geom_smooth() + facet_wrap(~rewFunc)
ggplot(df,aes(run_trial,rt_csv, color = kf_f1_pos_resp)) + geom_smooth() + facet_wrap(~rewFunc)
ggplot(df,aes(run_trial,rt_csv, lty = kf_f3_str_front_ins_resp)) + geom_smooth() + facet_wrap(~rewFunc)
ggplot(df,aes(run_trial,rt_csv, color = kf_f2_vmpfc_precun_resp)) + geom_smooth() + facet_wrap(~rewFunc)

ggplot(df,aes(run_trial,log(rt_swing), color = kf_f1_pos_resp)) + geom_smooth() + facet_wrap(~rewFunc)
# well, maybe with striato-fronto_insular and RT swings
ggplot(df,aes(run_trial,log(rt_swing), lty = kf_f3_str_front_ins_resp)) + geom_smooth() + facet_wrap(~rewFunc)
# and perhaps 
ggplot(df,aes(run_trial,log(rt_swing), size = kf_f2_vmpfc_precun_resp)) + geom_smooth() + facet_wrap(~rewFunc)


# and for PE
ggplot(df,aes(run_trial,rt_csv, color = pe_f1_cort_str_resp, lty = pe_f2_hipp_resp)) + geom_smooth() + facet_wrap(~rewFunc)
ggplot(df,aes(run_trial,rt_csv, color = pe_f1_cort_str_resp)) + geom_smooth() + facet_wrap(~rewFunc)
ggplot(df,aes(run_trial,rt_csv, lty = pe_f2_hipp_resp)) + geom_smooth() + facet_wrap(~rewFunc)

#hippocampal regions -- big effect of posterior Hipp/PE on IEV RTs in low antHippH subjects...
ggplot(df,aes(run_trial,rt_csv, color = h_HippAntL_resp)) + geom_smooth() + facet_wrap(~rewFunc)
ggplot(df,aes(run_trial,rt_csv, lty = pe_f2_hipp_resp)) + geom_smooth() + facet_wrap(~rewFunc)

ggplot(df,aes(rt_lag,rt_csv, color = h_HippAntL_resp)) + geom_smooth(method = 'loess') #+ facet_wrap(~rewFunc)
ggplot(df,aes(rt_lag,rt_csv, lty = pe_f2_hipp_resp)) + geom_smooth(method = 'loess') 

ggplot(df,aes(rt_vmax,rt_csv, color = h_HippAntL_resp)) + geom_smooth(method = 'loess') 
ggplot(df,aes(rt_vmax,rt_csv, lty = pe_f2_hipp_resp)) + geom_smooth(method = 'loess') 


# hippo RT swings
ggplot(df,aes(run_trial,log(rt_swing)/rt_lag, color = h_HippAntL_resp)) + geom_smooth(method = 'loess') + facet_wrap(~rewFunc)
ggplot(df,aes(run_trial,log(rt_swing)/rt_lag, lty = pe_f2_hipp_resp)) + geom_smooth(method = 'loess') + facet_wrap(~rewFunc)

ggplot(df[!is.na(df$rt_vmax),],aes(rt_lag,log(rt_change), color = h_HippAntL_resp)) + geom_smooth(method = 'glm') + 
  facet_wrap(~rt_vmax>15)
ggplot(df[!is.na(df$rt_vmax),],aes(rt_lag,log(rt_change), lty = pe_f2_hipp_resp)) + geom_smooth(method = 'glm') + 
  facet_wrap(~rt_vmax>15)



ggplot(df,aes(rt_lag,rt_csv, color = h_HippAntL_resp)) + geom_smooth() + facet_wrap(~rewFunc)
ggplot(df,aes(rt_lag,rt_csv, color = omission_lag, lty = pe_f2_hipp_resp)) + geom_smooth() + facet_wrap(~rewFunc)

ggplot(df,aes(run_trial,log(rt_swing), color = h_HippAntL_resp)) + geom_smooth() #+ facet_wrap(~rewFunc)
ggplot(df,aes(run_trial,log(rt_swing), color = omission_lag, lty = pe_f2_hipp_resp)) + geom_smooth(method = 'loess') #+ facet_wrap(~rewFunc)

ggplot(df,aes(rt_lag,log(rt_swing), color = pe_f1_cort_str_resp, lty = pe_f2_hipp_resp)) + geom_smooth(method = 'loess') #+ facet_wrap(~rewFunc)

ggplot(df,aes(pe_f2_hipp,log(rt_swing), color = pe_f1_cort_str_resp)) + geom_smooth(method = 'loess') + facet_wrap(~rewFunc)
ggplot(df,aes(-h_HippAntL,log(rt_swing), color = h_fp)) + geom_smooth(method = 'loess') + facet_wrap(~rewFunc)


# not a huge impact of KLD betas on RT swings...
# hold off on further plotting for now (2/19/19)
ggplot(df,aes(v_max_wi,log(rt_swing), lty = k_f1_IPL_crbl_MFG_ITG_resp)) + geom_smooth() + facet_wrap(~rewFunc)
ggplot(df,aes(v_max_wi,log(rt_swing), size = k_f2_visual_resp)) + geom_smooth() + facet_wrap(~rewFunc)
ggplot(df,aes(v_max_wi,log(rt_swing), color = k_f3_OFC_resp)) + geom_smooth() + facet_wrap(~rewFunc)
# does KL control H?
p1 <- ggplot(df[df$run>1,],aes(run_trial,v_entropy, color = d_f1_FP_SMAresp)) + geom_smooth() + facet_wrap(~rewFunc) + scale_color_brewer(palette = 'Set3', direction = -1)+ scale_x_continuous(breaks = c(1,50)) + guides(color=FALSE) + theme(axis.title.x=element_blank())
p2 <- ggplot(df[df$run>1,],aes(run_trial,v_entropy, color = d_f2_VSresp)) + geom_smooth() + facet_wrap(~rewFunc) + scale_color_brewer(palette = 'Dark2', direction = 1)+ scale_x_continuous(breaks = c(1,50)) + guides(color=FALSE)+ theme(axis.title.x=element_blank())
p3 <- ggplot(df[df$run>1,],aes(run_trial,v_entropy, color = d_f3_ACC_ins_resp)) + geom_smooth() + facet_wrap(~rewFunc) + scale_color_brewer(palette = 'Set1', direction = 1)+ scale_x_continuous(breaks = c(1,50)) + guides(color=FALSE)
p4 <- ggplot(df[df$run>1,],aes(run_trial,v_entropy, color = low_v_fp_acc_vlpfc)) + geom_smooth() + facet_wrap(~rewFunc) + scale_color_brewer(palette = 'Paired', direction = -1)+ scale_x_continuous(breaks = c(1,50)) + guides(color=FALSE)+ theme(axis.title.x=element_blank())
p5 <- ggplot(df[df$run>1,],aes(run_trial,v_entropy, color = v_paralimbic)) + geom_smooth() + facet_wrap(~rewFunc) + scale_color_brewer(palette = 'Pastel1', direction = 1)+ scale_x_continuous(breaks = c(1,50)) + guides(color=FALSE)+ theme(axis.title.x=element_blank())
p6 <- ggplot(df[df$run>1,],aes(run_trial,v_entropy, color = h_fp)) + geom_smooth() + facet_wrap(~rewFunc) + scale_color_brewer(palette = 'Pastel2', direction = -1)+ scale_x_continuous(breaks = c(1,50)) + guides(color=FALSE)+ theme(axis.title.x=element_blank())
p7 <- ggplot(df[df$run>1,],aes(run_trial,v_entropy, color = low_h_paralimbic)) + geom_smooth() + facet_wrap(~rewFunc) + scale_color_brewer(palette = 'Accent', direction = -1)+ scale_x_continuous(breaks = c(1,50)) + guides(color=FALSE)+ theme(axis.title.x=element_blank())
pdf("all_betas_H_timecourse.pdf", width = 12, height = 18)
ggarrange(p4,p5,p6,p7,p1,p2,p3, ncol = 2,nrow = 7, common.legend = F, align = c('v'))
dev.off()

t1 <- ggplot(df[df$run>1,],aes(run_trial,log(rt_swing), color = d_f1_FP_SMAresp)) + geom_smooth() + facet_wrap(~rewFunc) + scale_color_brewer(palette = 'Set3', direction = -1) + scale_x_continuous(breaks = c(1,50)) + guides(color=FALSE)+ theme(axis.title.x=element_blank())
t2 <- ggplot(df[df$run>1,],aes(run_trial,log(rt_swing), color = d_f2_VSresp)) + geom_smooth() + facet_wrap(~rewFunc) + scale_color_brewer(palette = 'Dark2', direction = 1)+ scale_x_continuous(breaks = c(1,50)) + guides(color=FALSE)+ theme(axis.title.x=element_blank())
t3 <- ggplot(df[df$run>1,],aes(run_trial,log(rt_swing), color = d_f3_ACC_ins_resp)) + geom_smooth() + facet_wrap(~rewFunc) + scale_color_brewer(palette = 'Set1', direction = 1)+ scale_x_continuous(breaks = c(1,50)) + guides(color=FALSE)
# just double-check the d3 effect
# ggplot(df[df$run>1 & df$rt_swing>1,],aes(run_trial,(rt_swing), color = d_f3_ACC_ins_resp)) + geom_smooth() + facet_wrap(~rewFunc) + scale_color_brewer(palette = 'Set1', direction = 1)+ scale_x_continuous(breaks = c(1,50)) + guides(color=FALSE)
t4 <- ggplot(df[df$run>1,],aes(run_trial,log(rt_swing), color = low_v_fp_acc_vlpfc)) + geom_smooth() + facet_wrap(~rewFunc) + scale_color_brewer(palette = 'Paired', direction = -1)+ scale_x_continuous(breaks = c(1,50)) + guides(color=FALSE)+ theme(axis.title.x=element_blank())
t5 <- ggplot(df[df$run>1,],aes(run_trial,log(rt_swing), color = v_paralimbic)) + geom_smooth() + facet_wrap(~rewFunc) + scale_color_brewer(palette = 'Pastel1', direction = 1)+ scale_x_continuous(breaks = c(1,50)) + guides(color=FALSE)+ theme(axis.title.x=element_blank())
t6 <- ggplot(df[df$run>1,],aes(run_trial,log(rt_swing), color = h_fp)) + geom_smooth() + facet_wrap(~rewFunc) + scale_color_brewer(palette = 'Pastel2', direction = -1)+ scale_x_continuous(breaks = c(1,50)) + guides(color=FALSE)+ theme(axis.title.x=element_blank())
t7 <- ggplot(df[df$run>1,],aes(run_trial,log(rt_swing), color = low_h_paralimbic)) + geom_smooth() + facet_wrap(~rewFunc) + scale_color_brewer(palette = 'Accent', direction = -1)+ scale_x_continuous(breaks = c(1,50)) + guides(color=FALSE)+ theme(axis.title.x=element_blank())
pdf("all_betas_RT_swing_timecourse.pdf", width = 12, height = 18)
ggarrange(t4,t5,t6,t7,t1,t2,t3, ncol = 2,nrow = 7, common.legend = F, align = c('v'))
dev.off()

r1 <- ggplot(df[df$run>1,],aes(run_trial,rt_csv, color = d_f1_FP_SMAresp)) + geom_smooth() + facet_wrap(~rewFunc) + scale_color_brewer(palette = 'Set3', direction = -1)+ scale_x_continuous(breaks = c(1,50))+ theme(axis.title.x=element_blank())
r2 <- ggplot(df[df$run>1,],aes(run_trial,rt_csv, color = d_f2_VSresp)) + geom_smooth() + facet_wrap(~rewFunc) + scale_color_brewer(palette = 'Dark2', direction = 1)+ scale_x_continuous(breaks = c(1,50))+ theme(axis.title.x=element_blank())
r3 <- ggplot(df[df$run>1,],aes(run_trial,rt_csv, color = d_f3_ACC_ins_resp)) + geom_smooth() + facet_wrap(~rewFunc) + scale_color_brewer(palette = 'Set1', direction = 1)+ scale_x_continuous(breaks = c(1,50))
r4 <- ggplot(df[df$run>1,],aes(run_trial,rt_csv, color = low_v_fp_acc_vlpfc)) + geom_smooth() + facet_wrap(~rewFunc) + scale_color_brewer(palette = 'Paired', direction = -1)+ scale_x_continuous(breaks = c(1,50))+ theme(axis.title.x=element_blank())
r5 <- ggplot(df[df$run>1,],aes(run_trial,rt_csv, color = v_paralimbic)) + geom_smooth() + facet_wrap(~rewFunc) + scale_color_brewer(palette = 'Pastel1', direction = 1)+ scale_x_continuous(breaks = c(1,50))+ theme(axis.title.x=element_blank())
r6 <- ggplot(df[df$run>1,],aes(run_trial,rt_csv, color = h_fp)) + geom_smooth() + facet_wrap(~rewFunc) + scale_color_brewer(palette = 'Pastel2', direction = -1)+ scale_x_continuous(breaks = c(1,50))+ theme(axis.title.x=element_blank())
r7 <- ggplot(df[df$run>1,],aes(run_trial,rt_csv, color = low_h_paralimbic)) + geom_smooth() + facet_wrap(~rewFunc) + scale_color_brewer(palette = 'Accent', direction = -1)+ scale_x_continuous(breaks = c(1,50))+ theme(axis.title.x=element_blank())
pdf("all_betas_RT_timecourse.pdf", width = 12, height = 18)
ggarrange(r4,r5,r6,r7,r1,r2,r3, ncol = 2,nrow = 7, common.legend = F, align = c('v'))
dev.off()

allplots <- ggarrange(p4, t4, r4,
          p5, t5, r5,
          p6, t6, r6,
          p7, t7, r7,
          p1, t1, r1,
          p2, t2, r2,
          p3, t3, r3, ncol = 3,nrow = 7, common.legend = F)


pdf("all_betas_H_rt_swing_RT_timecourse.pdf", width = 12, height = 16)
ggarrange(p4, t4, r4,
          p5, t5, r5,
          p6, t6, r6,
          p7, t7, r7,
          p1, t1, r1,
          p2, t2, r2,
          p3, t3, r3, ncol = 3,nrow = 7, common.legend = F)
dev.off()


ggplot(df,aes(run_trial,v_entropy_wi, lty = d_f1_FP_SMAresp)) + geom_smooth() + facet_wrap(~rewFunc)
ggplot(df,aes(run_trial,v_entropy_wi, size = d_f2_VSresp)) + geom_smooth() + facet_wrap(~rewFunc)
ggplot(df,aes(run_trial,v_entropy_wi, color = d_f3_ACC_ins_resp)) + geom_smooth() + facet_wrap(~rewFunc) + scale_color_brewer(palette = 'BrBG', direction = 1)


# reward omission
ggplot(na.omit(df),aes(d_f1_FP_SMAresp,log(rt_swing), lty = last_outcome)) + geom_boxplot()
ggplot(na.omit(df),aes(d_f2_VSresp,log(rt_swing), lty = last_outcome)) + geom_boxplot()
ggplot(na.omit(df),aes(d_f3_ACC_ins_resp,log(rt_swing), lty = last_outcome)) + geom_boxplot()

ggplot(na.omit(df),aes(last_outcome,log(rt_swing), lty = d_f1_FP_SMAresp, size = d_f2_VSresp, color = d_f3_ACC_ins_resp)) + geom_boxplot()
ggplot(na.omit(df),aes(d_f2_VSresp,log(rt_swing), lty = last_outcome)) + geom_boxplot()
ggplot(na.omit(df),aes(d_f3_ACC_ins_resp,log(rt_swing), lty = last_outcome)) + geom_boxplot()


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

pdf('rt_by_condition_and_h_betas.pdf', height = 8, width = 8)
ggplot(df, aes(run_trial, rt_csv, color = low_h_paralimbic, lty = h_fp)) + geom_smooth(method = "loess") + facet_wrap(~rewFunc)
dev.off()

pdf('rt_by_condition_and_v_betas.pdf', height = 8, width = 8)
ggplot(df, aes(run_trial, rt_csv, color = v_paralimbic, lty = low_v_fp_acc_vlpfc)) + geom_smooth(method = "loess") + facet_wrap(~rewFunc)
dev.off()

ggplot(df, aes(v_max,rt_swing,  color = low_h_paralimbic, lty = h_fp)) + geom_smooth(method = "gam")

# plot vs. raw data
ggplot(df, aes(rt_lag, rt_csv, color = low_h_paralimbic, lty = h_fp)) + geom_smooth(method = "glm")+ facet_wrap(~run_trial>20)
ggplot(df, aes(rt_vmax, rt_csv, color = low_h_paralimbic, lty = h_fp)) + geom_smooth(method = "glm")+ facet_wrap(~run_trial>20)
pdf("h_timecourse_brain_fixed.pdf", width = 8, height = 8)
ggplot(df, aes(run_trial, v_entropy,color = low_h_paralimbic, lty = h_fp)) + geom_smooth(method = "loess") + facet_wrap(~gamma>0)
dev.off()


ggplot(df, aes(run_trial, v_entropy, color = rewFunc)) + geom_smooth()
summary(m0 <- lmer(log(rt_swing) ~ scale(-1/run_trial) * scale(v_max) * scale(v_entropy) * rewFunc + (1|id/run), df[df$rt_swing>0,]))

# plot out


save(file = 'vhd_models.Rdata', list = ls(all.names = TRUE))
 # load('vhd_models.Rdata')
