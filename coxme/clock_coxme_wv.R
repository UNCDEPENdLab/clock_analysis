# runs mixed-effects Cox models on clock data
# when running the first time, first run compute_sceptic_fmri_statistics.R
# basedir <- "~/Data_Analysis"
basedir <- "~/code"
setwd(file.path(basedir, "clock_analysis/coxme"))
library(readr)
library(ggplot2)
library(tidyr)
library(tibble)
library(reshape2)
library(emmeans)
library(dplyr)
library(lme4)
library(survival)
library(coxme)
library(survminer)
# library(OIsurv)
library(ggpubr)
# devtools::install_github('junkka/ehahelper') # requires gfortran, $ brew cask install gfortran
library(ehahelper)
library(broom)
library(broom.mixed)
library(car)
load ("fMRI_MEG_coxme_objects_Nov15_2020")


# coxme on wv smiles
# summary(cox_wv1 <- coxme(Surv(t1,t2,response) ~ rtlag_sc + wv3b0a1 + trial_neg_inv_sc*rewFunc + 
#                           value_wi + uncertainty_wi + 
#                            (1|ID), bb))

# summary(cox_wv2 <- coxme(Surv(t1,t2,response) ~ rtlag_sc + wv3b1a1 + trial_neg_inv_sc*rewFunc + 
#                            value_wi + uncertainty_wi + 
#                            (1|ID), bb))
# summary(cox_wv3 <- coxme(Surv(t1,t2,response) ~ rtlag_sc + wv3b1a2 + trial_neg_inv_sc*rewFunc + 
#                            value_wi + uncertainty_wi + 
#                            (1|ID), bb))

# MEG replication

###########
# simple models 
# with single factors
# _f suffix denotes model run on full interval (including the suspect ends)
###########

 ###################
 # Main analysis (interval ends filtered out)
 ###################
 summary(wv_ge <- coxme(Surv(t1,t2,response) ~ wvs1b1a1*general_entropy + wvs2b1a1*general_entropy + wvs3b1a1*general_entropy + 
                                 value_wi*general_entropy + uncertainty_wi*general_entropy + 
                                 (1|ID), fbb))
 Anova(wv_ge, '3')
 
 # reduced model w/o value or uncertainty
 summary(wv_only_ge <- coxme(Surv(t1,t2,response) ~ wvs1b1a1*general_entropy + wvs2b1a1*general_entropy + wvs3b1a1*general_entropy + 
                          (1|ID), fbb))
 
 summary(ge_rtlag <- coxme(Surv(t1,t2,response) ~ rtlag*general_entropy +
                          value_wi*general_entropy + uncertainty_wi*general_entropy + 
                          (1|ID), fbb))
 Anova(ge_rtlag, '3')
 
 
 summary(wv_vlpfc <- coxme(Surv(t1,t2,response) ~ wvs1b1a1*entropy_vlPFC + wvs2b1a1*entropy_vlPFC + wvs3b1a1*entropy_vlPFC + 
                                    value_wi*entropy_vlPFC + uncertainty_wi*entropy_vlPFC + 
                                    (1|ID), fbb))
 summary(wv_vlpfc)
 Anova(wv_vlpfc, '3')

 # MEG replication
 summary(mwv_ge <- coxme(Surv(t1,t2,response) ~ wvs1b1a1*general_entropy + wvs2b1a1*general_entropy + wvs3b1a1*general_entropy + 
                          value_wi*general_entropy + uncertainty_wi*general_entropy + 
                          (1|ID), mfbb))
 Anova(mwv_ge, '3')
 
 # reduced model w/o value or uncertainty
 summary(mwv_only_ge <- coxme(Surv(t1,t2,response) ~ wvs1b1a1*general_entropy + wvs2b1a1*general_entropy + wvs3b1a1*general_entropy + 
                               (1|ID), mfbb))
 
 
 summary(mge_rtlag <- coxme(Surv(t1,t2,response) ~ rt_lag*general_entropy +
                             value_wi*general_entropy + uncertainty_wi*general_entropy + 
                             (1|ID), mfbb))
 Anova(mge_rtlag, '3')
 
 
 summary(mwv_vlpfc <- coxme(Surv(t1,t2,response) ~ wvs1b1a1*entropy_vlPFC + wvs2b1a1*entropy_vlPFC + wvs3b1a1*entropy_vlPFC + 
                             value_wi*entropy_vlPFC + uncertainty_wi*entropy_vlPFC + 
                             (1|ID), mfbb))
 Anova(mwv_vlpfc, '3')
 
  
# add the fef and med_par factors with general entropy
 
 summary(wv_ge_fef_par <- coxme(Surv(t1,t2,response) ~ wvs1b1a1*general_entropy + wvs2b1a1*general_entropy + wvs3b1a1*general_entropy +
                                  wvs1b1a1*fef + wvs2b1a1*fef + wvs3b1a1*fef +
                                  wvs1b1a1*med_par + wvs2b1a1*med_par + wvs3b1a1*med_par +
                                  value_wi*general_entropy + uncertainty_wi*general_entropy + 
                                  value_wi*fef + uncertainty_wi*fef + 
                                  value_wi*med_par + uncertainty_wi*med_par + 
                                  (1|ID), fbb))
 Anova(wv_ge_fef_par, '3')
 
 # interactions with trial
 summary(wv_ge_trial <- coxme(Surv(t1,t2,response) ~ wvs1b1a1*general_entropy*trial + wvs2b1a1*general_entropy*trial + wvs3b1a1*general_entropy*trial + 
                                       value_wi*general_entropy + uncertainty_wi*general_entropy*trial + 
                                       (1|ID), fbb))
 Anova(wv_ge_trial, '3')
 
 # add specific factors
 
 summary(wv_alldan_middle_trial <- coxme(Surv(t1,t2,response) ~ wvs1b1a1*general_entropy*trial + 
                                           wvs2b1a1*general_entropy*trial + 
                                           wvs3b1a1*general_entropy*trial + 
                                           value_wi*general_entropy + 
                                           uncertainty_wi*general_entropy*trial + 
                                           wvs1b1a1*fef*trial + 
                                           wvs2b1a1*fef*trial + 
                                           wvs3b1a1*fef*trial + 
                                           value_wi*fef + 
                                           uncertainty_wi*fef*trial + 
                                           wvs1b1a1*med_par*trial + 
                                           wvs2b1a1*med_par*trial + 
                                           wvs3b1a1*med_par*trial + 
                                           value_wi*med_par + 
                                           uncertainty_wi*med_par*trial + 
                                           (1|ID), fbb))
 Anova(wv_alldan_middle_trial, '3')
 # curious effect of med_par on uncertainty-seeking 
 # this also holds in a simple model above (wv_med_par, z=6.65)
 # check w/o other factors:
 summary(wv_medpar_middle_trial <- coxme(Surv(t1,t2,response) ~ 
                                           wvs1b1a1*med_par*trial + 
                                           wvs2b1a1*med_par*trial + 
                                           wvs3b1a1*med_par*trial + 
                                           value_wi*med_par + 
                                           uncertainty_wi*med_par*trial + 
                                           (1|ID), fbb))
 Anova(wv_medpar_middle_trial, '3')
 
 ########################
 # Effects of PE clusters
 ########################
 # only IPS PE
 summary(wv_pe_ips <- coxme(Surv(t1,t2,response) ~ wvs1b1a1*pe_ips + wvs2b1a1*pe_ips + wvs3b1a1*pe_ips +
                              (1|ID), fbb))
 Anova(wv_pe_ips, '3')
 
 
 summary(wv_pe_ips <- coxme(Surv(t1,t2,response) ~ wvs1b1a1*pe_ips + wvs2b1a1*pe_ips + wvs3b1a1*pe_ips + 
                          value_wi*pe_ips + uncertainty_wi*pe_ips + 
                          (1|ID), fbb))
 Anova(wv_pe_ips, '3')

 summary(wv_pe_ips_v <- coxme(Surv(t1,t2,response) ~ value_wi*wvs1b1a1*pe_ips + value_wi*wvs2b1a1*pe_ips + value_wi*wvs3b1a1*pe_ips + 
                               uncertainty_wi*pe_ips + 
                              (1|ID), fbb))
 Anova(wv_pe_ips_v, '3')
 
 summary(wv_h_dan_v <- coxme(Surv(t1,t2,response) ~ value_wi*wvs1b1a1*general_entropy + value_wi*wvs2b1a1*general_entropy + value_wi*wvs3b1a1*general_entropy + 
                                uncertainty_wi*general_entropy + 
                                (1|ID), fbb))
 Anova(wv_h_dan_v, '3')
 
  
 summary(rtlag_pe_ips <- coxme(Surv(t1,t2,response) ~ rtlag*pe_ips + 
                              value_wi*pe_ips + uncertainty_wi*pe_ips + 
                              (1|ID), fbb))
 Anova(rtlag_pe_ips, '3')
 
  
 # all cortical PE
 summary(wv_pe_cort <- coxme(Surv(t1,t2,response) ~ wvs1b1a1*pe_f1_cort_hipp + wvs2b1a1*pe_f1_cort_hipp + wvs3b1a1*pe_f1_cort_hipp + 
                              value_wi*pe_f1_cort_hipp + uncertainty_wi*pe_f1_cort_hipp + 
                              (1|ID), fbb))
 Anova(wv_pe_cort, '3')
 
 
 # striatal PE
 summary(wv_pe_str <- coxme(Surv(t1,t2,response) ~ wvs1b1a1*pe_f3_str + wvs2b1a1*pe_f3_str + wvs3b1a1*pe_f3_str + 
                              value_wi*pe_f3_str + uncertainty_wi*pe_f3_str + 
                              (1|ID), fbb))
 Anova(wv_pe_str, '3')
 
 # cerebellar PE
 summary(wv_pe_cerebell <- coxme(Surv(t1,t2,response) ~ wvs1b1a1*pe_f2_cerebell + wvs2b1a1*pe_f2_cerebell + wvs3b1a1*pe_f2_cerebell + 
                              value_wi*pe_f2_cerebell + uncertainty_wi*pe_f2_cerebell + 
                              (1|ID), fbb))
 Anova(wv_pe_cerebell, '3')
 
 # hippocampal PE
 summary(wv_pe_ph <- coxme(Surv(t1,t2,response) ~ wvs1b1a1*pe_PH_r + wvs2b1a1*pe_PH_r + wvs3b1a1*pe_PH_r + 
                                   value_wi*pe_PH_r + uncertainty_wi*pe_PH_r + 
                                   (1|ID), fbb))
 Anova(wv_pe_ph, '3')
 
 ######################
 ## PE: MEG replication
 ######################
 summary(mwv_pe_ips <- coxme(Surv(t1,t2,response) ~ wvs1b1a1*pe_ips + wvs2b1a1*pe_ips + wvs3b1a1*pe_ips + 
                              value_wi*pe_ips + uncertainty_wi*pe_ips + 
                              (1|ID), mfbb))
 Anova(mwv_pe_ips, '3')
 
 summary(mrtlag_pe_ips <- coxme(Surv(t1,t2,response) ~ rt_lag*pe_ips + 
                                 value_wi*pe_ips + uncertainty_wi*pe_ips + 
                                 (1|ID), mfbb))
 Anova(mrtlag_pe_ips, '3')
 
 # all cortical PE
 summary(mwv_pe_cort <- coxme(Surv(t1,t2,response) ~ wvs1b1a1*pe_f1_cort_hipp + wvs2b1a1*pe_f1_cort_hipp + wvs3b1a1*pe_f1_cort_hipp + 
                               value_wi*pe_f1_cort_hipp + uncertainty_wi*pe_f1_cort_hipp + 
                               (1|ID), mfbb))
 Anova(mwv_pe_cort, '3')
 
 
 # striatal PE
 summary(mwv_pe_str <- coxme(Surv(t1,t2,response) ~ wvs1b1a1*pe_f3_str + wvs2b1a1*pe_f3_str + wvs3b1a1*pe_f3_str + 
                              value_wi*pe_f3_str + uncertainty_wi*pe_f3_str + 
                              (1|ID), mfbb))
 Anova(mwv_pe_str, '3')
 
 # cerebellar PE
 summary(mwv_pe_cerebell <- coxme(Surv(t1,t2,response) ~ wvs1b1a1*pe_f2_cerebell + wvs2b1a1*pe_f2_cerebell + wvs3b1a1*pe_f2_cerebell + 
                                   value_wi*pe_f2_cerebell + uncertainty_wi*pe_f2_cerebell + 
                                   (1|ID), mfbb))
 Anova(mwv_pe_cerebell, '3')
 
 summary(mwv_pe_ph <- coxme(Surv(t1,t2,response) ~ wvs1b1a1*pe_PH_r + wvs2b1a1*pe_PH_r + wvs3b1a1*pe_PH_r + 
                             value_wi*pe_PH_r + uncertainty_wi*pe_PH_r + 
                             (1|ID), mfbb))
 Anova(mwv_pe_ph, '3')
 
 
 # just the SFG blobs
 summary(wv_sfgl_middle <- coxme(Surv(t1,t2,response) ~ wvs1b1a1*dan_l_sfg + wvs2b1a1*dan_l_sfg + wvs3b1a1*dan_l_sfg + 
                                   value_wi*dan_l_sfg + uncertainty_wi*dan_l_sfg + 
                                   (1|ID), fbb))
 Anova(wv_sfgl_middle, '3')
 
 summary(wv_sfgr_middle <- coxme(Surv(t1,t2,response) ~ wvs1b1a1*dan_r_sfg + wvs2b1a1*dan_r_sfg + wvs3b1a1*dan_r_sfg + 
                                   value_wi*dan_r_sfg + uncertainty_wi*dan_r_sfg + 
                                   (1|ID), fbb))
 Anova(wv_sfgr_middle, '3')
 # 
 
# DAN mean beta
wv_dan1_f <- coxme(Surv(t1,t2,response) ~ wvs1b1a1*DAN + wvs2b1a1*DAN + wvs3b1a1*DAN +
                           value_wi*DAN + uncertainty_wi*DAN + 
                           (1|ID), bb)
summary(wv_dan1_f)
Anova(wv_dan1_f, '3')

# fef from bifactor model
summary(wv_fef_f <- coxme(Surv(t1,t2,response) ~ wvs1b1a1*fef + wvs2b1a1*fef + wvs3b1a1*fef + 
                           value_wi*fef + uncertainty_wi*fef + 
                           (1|ID), bb))
Anova(wv_fef_f, '3')

# parietal from bifactor model
summary(wv_med_par_f <- coxme(Surv(t1,t2,response) ~ wvs1b1a1*med_par + wvs2b1a1*med_par + wvs3b1a1*med_par + 
                          value_wi*med_par + uncertainty_wi*med_par + 
                          (1|ID), bb))
Anova(wv_med_par_f, '3')

# general DAN from bifactor model
summary(wv_ge_f <- coxme(Surv(t1,t2,response) ~ wvs1b1a1*general_entropy + wvs2b1a1*general_entropy + wvs3b1a1*general_entropy + 
                          value_wi*general_entropy + uncertainty_wi*general_entropy + 
                          (1|ID), bb))
Anova(wv_ge, '3')

# vlPFC mean
summary(wv_vlpfc_f <- coxme(Surv(t1,t2,response) ~ wvs1b1a1*entropy_vlPFC + wvs2b1a1*entropy_vlPFC + wvs3b1a1*entropy_vlPFC + 
                         value_wi*entropy_vlPFC + trial_neg_inv_sc*uncertainty_wi*entropy_vlPFC + 
                         (1|ID), bb))
Anova(wv_vlpfc_f, '3')

 # try interaction with trials

# summary(cox_wv2 <- coxme(Surv(t1,t2,response) ~ rtlag_sc + wv3b1a1 + trial_neg_inv_sc*rewFunc + 
#                            value_wi + uncertainty_wi + 
#                            (1|ID), bb))
# summary(cox_wv3 <- coxme(Surv(t1,t2,response) ~ rtlag_sc + wv3b1a2 + trial_neg_inv_sc*rewFunc + 
#                            value_wi + uncertainty_wi + 
#                            (1|ID), bb))
# 
# # main analysis
# summary(cox_hipp1a <- coxme(Surv(t1,t2,response) ~ rtlag_sc*pe_f2_hipp + rtlag_sc*AH_sc + trial_neg_inv_sc*rewFunc + trial_neg_inv_sc*uncertainty_wi + trial_neg_inv_sc*value_wi +
#                               value_wi*pe_f2_hipp + value_wi*AH_sc + uncertainty_wi*AH_sc +uncertainty_wi*pe_f2_hipp + (1|ID), sdf))
# out <- capture.output(summary(cox_hipp1a))
# setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/hippo/supp/')
# write.csv(out,file = "coxme_hipp_summary.csv", sep = ":\t")
# 
# 
# # sensitivity analysis: remove no-go zones in the first 1000 ms and last 500 ms
# summary(cox_hipp1 <- coxme(Surv(t1,t2,response) ~ rtlag_sc*pe_f2_hipp + rtlag_sc*AH_sc + trial_neg_inv_sc*rewFunc + trial_neg_inv_sc*uncertainty_wi + trial_neg_inv_sc*value_wi +
#                             value_wi*pe_f2_hipp + value_wi*AH_sc + uncertainty_wi*AH_sc +uncertainty_wi*pe_f2_hipp + (1|ID), fdf))
# 
# # add linear and quadratic time
# summary(cox_hipp1b <- coxme(Surv(t1,t2,response) ~ scale(bin) + I(scale(bin)^2) + rtlag_sc*pe_f2_hipp + rtlag_sc*AH_sc + trial_neg_inv_sc*rewFunc + trial_neg_inv_sc*uncertainty_wi + trial_neg_inv_sc*value_wi +
#                               value_wi*pe_f2_hipp + value_wi*AH_sc + uncertainty_wi*AH_sc +uncertainty_wi*pe_f2_hipp + (1|ID), sdf))
# 
# 
# devtools::install_github('junkka/ehahelper')
# library(ehahelper)
# library(broom)
# library(broom.mixed)
# tidy_cox_hipp <- tidy(cox_hipp1a, exponentiate = F)[,1:5]
# stargazer(tidy_cox_hipp, type = 'html', out = 'tidy_cox_hipp.html', summary = F, digits = 3, digits.extra = 10)
# write_csv(tidy_cox_hipp, 'tidy_cox_hipp.csv')
# p1 <- plot_model(cox_hipp1a,transform = 'exp', terms = c("pe_f2_hipp:value_wi", "AH_sc:value_wi", "AH_sc:uncertainty_wi", "pe_f2_hipp:uncertainty_wi"), show.values = T, show.p = T, value.offset = .3 )
# pdf('AH_PH_uncertainty_value_coxme.pdf', height = 3, width = 4.5)
# p1 + ylim(c(.95, 1.05)) + ylab("Effect on hazard of response, A.U.") + scale_x_discrete(labels = c("PH * value", "PH * uncertainty", "AH * value", "AH * uncertainty")) + labs(title = "") +
#   geom_hline(yintercept = 1)
# dev.off()
# # interactions with trial_neg_inv_sc: inferior model
# summary(cox_hipp3 <- coxme(Surv(t1,t2,response) ~ rtlag_sc*pe_f2_hipp + rtlag_sc*AH_sc + trial_neg_inv_sc*rewFunc + trial_neg_inv_sc*value_wi +
#                              value_wi*pe_f2_hipp + value_wi*AH_sc + uncertainty_wi*AH_sc + trial_neg_inv_sc*uncertainty_wi*pe_f2_hipp + (1|ID), fdf))
# summary(cox_hipp3)
# Anova(cox_hipp3, '3')
# 
