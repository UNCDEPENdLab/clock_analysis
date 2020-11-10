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
load(file="clock_for_coxme_value_only_070518.RData")

# sanity checks on within-trial matrices

# add hippocampal betas
setwd(file.path(basedir, 'clock_analysis/fmri/keuka_brain_behavior_analyses/'))
load('trial_df_and_vh_pe_clusters_u.Rdata')
hdf <- df %>% select (ID, DAN, dan_parietal, dan_l_sfg, dan_r_sfg, general_entropy,
med_par, fef, entropy_vlPFC) %>% unique()
sdf <- sdf %>% inner_join(hdf, by = "ID")

# mark events (responses)
sdf$response <- round(sdf$rt/1000, digits = 1)==sdf$t2

# # KS across conditions
# 
# ks.test(df$rt_csv[df$rewFunc=='CEV'], "punif", 1, 4000)
# ks.test(df$rt_csv[df$rewFunc=='CEVR'], "punif", 1, 4000)
# ks.test(df$rt_csv[df$rewFunc=='IEV'], "punif", 1, 4000)
# ks.test(df$rt_csv[df$rewFunc=='DEV'], "punif", 1, 4000)
# 
# # Why is this happening?
# hist(df$rt_csv[df$rewFunc=='CEV'])
# hist(df$rt_csv[df$rewFunc=='CEVR'])
# hist(df$rt_csv[df$rewFunc=='IEV'])
# hist(df$rt_csv[df$rewFunc=='DEV'])
# library(texmex)
# attach(mtcars)
# par(mfrow=c(2,2))
# plot(df$rt_csv[df$rewFunc=='DEV'], edf(df$rt_csv[df$rewFunc=='DEV']))
# plot(df$rt_csv[df$rewFunc=='IEV'], edf(df$rt_csv[df$rewFunc=='IEV']))
# plot(df$rt_csv[df$rewFunc=='CEV'], edf(df$rt_csv[df$rewFunc=='CEV']))
# plot(df$rt_csv[df$rewFunc=='CEVR'], edf(df$rt_csv[df$rewFunc=='CEVR']))
# uniform <- runif(1000, 1,4000)
# because IEV has the most uniform EDF of RTs overall

# # inspect piecewize hazard functions
# library(muhaz)
# # piecewise <- pehaz(ddf$latency, delta=ddf$quit, width=NA, min.time=0, max.time=20.1)
# # plot(piecewise, xlab="Time, seconds", ylab="Quit Rate")
# pwI <- pehaz(na.omit(bdf$rt[bdf$rewFunc=='IEV']), width=200, min.time=0, max.time=4000)
# pwD <- pehaz(na.omit(bdf$rt[bdf$rewFunc=='DEV']), width=200, min.time=0, max.time=4000)
# pwR <- pehaz(na.omit(bdf$rt[bdf$rewFunc=='CEVR']), width=200, min.time=0, max.time=4000)
# pwC <- pehaz(na.omit(bdf$rt[bdf$rewFunc=='CEV']), width=200, min.time=0, max.time=4000)
# 
# h <- c(pwI$Hazard,pwD$Hazard,pwC$Hazard,pwR$Hazard)
# times <- c(pwI$Cuts[2:length(pwI$Cuts)],pwD$Cuts[2:length(pwD$Cuts)],pwC$Cuts[2:length(pwC$Cuts)],pwR$Cuts[2:length(pwR$Cuts)])
# condition <- c(rep('IEV',length(pwI$Hazard)), rep('DEV',length(pwD$Hazard)),rep('CEV',length(pwC$Hazard)),rep('CEVR',length(pwR$Hazard)))
# H <- as.tibble(cbind(h,times,condition))
# H$h <- as.numeric(H$h)
# H$logh <- log(H$h)
# H$times <- as.numeric(H$times)
# # pdf("piecewise_hazard_fx_all_groups.pdf", width = 12, height = 4)
# ggplot(H,aes(times, h, color = condition)) + geom_line()
# ggplot(H,aes(times, logh, color = condition)) + geom_line()
# # not what I expected: shared underlying hazard across contingencies
# 
# badfit <- survfit(Surv(t2) ~ rewFunc, type = 'right', origin = .1, data=sdf)
# badfit1 <- survfit(Surv(rt) ~ rewFunc,  data=bdf)
# 
# 
# plot(badfit1, mark.time=FALSE, lty=1:4,
#      xlab="ms", ylab="Proportion still waiting")
# legend(3000, .85, c("CEV", "CEVR", "DEV", "IEV"),
#        lty=1:4, bty='n')

# number bins
sdf <- sdf %>% group_by(ID, run, trial) %>% dplyr::mutate(bin = 1:n(), time = bin/10) %>% ungroup() %>% 
  group_by(ID) %>% mutate(value_wi = scale(value),
                          uncertainty_wi = scale(uncertainty),
                          value_b = mean(value),
                          uncertainty_b = mean(uncertainty), 
                          trial_neg_inv_sc = scale(-1/trial)) %>% ungroup() 
# filter out no-go zones at the edges
fdf <- sdf %>% filter(bin >10 & bin <35)

# diagnostics on uncertainty and value distributions

# ggplot(sdf %>% filter(bin >10 & bin <35), aes(time, value_wi, color = rewFunc)) + geom_smooth(method = 'gam', formula = y~splines::ns(x,2)) 
# ggplot(sdf %>% filter(bin >10 & bin <35), aes(time, uncertainty_wi, color = rewFunc)) + geom_smooth(method = 'gam', formula = y~splines::ns(x,2))
# corr.test(fdf$value_wi, fdf$uncertainty_wi)
# corr.test(fdf$value, fdf$uncertainty)
# 


##Build 0/1 WV smiles (selection history with different buffer sizes)
sdf <- sdf %>% select(ID, run, trial, bin, everything())

#use trial-wise data frame for getting lagged timesteps
trialwise <- sdf %>% filter(bin==1) %>% 
  group_by(ID, run) %>%
  mutate(timesteplag1=dplyr::lag(timestep, n=1, order_by=trial),
         timesteplag2=dplyr::lag(timestep, n=2, order_by=trial),
         timesteplag3=dplyr::lag(timestep, n=3, order_by=trial),
         timesteplag4=dplyr::lag(timestep, n=4, order_by=trial)) %>%
  ungroup() %>% select(ID, run, trial, timesteplag1, timesteplag2, timesteplag3, timesteplag4)

sdf <- sdf %>% select(-timesteplag) %>% left_join(trialwise, by=c("ID", "run", "trial"))


get_wv_smile <- function(microdf, nlags=3, nbefore=0, nafter=1, spec=FALSE) {
  smile <- rep(0, nrow(microdf))
  lcols <- ifelse(isTRUE(spec), paste0("timesteplag", nlags), paste0("timesteplag", 1:nlags))
  wvprefix <- ifelse(isTRUE(spec), "wvs", "wv")
  vv <- microdf[1,lcols] #just first row (bin) of this trial is needed
  if (nrow(vv) > 0L) {
    allpos <- lapply(vv, function(x) {
      if (is.na(x)) {
        return(NULL)
      } else {
        return(seq.int(x-nbefore, x+nafter))
      }
    })
    
    topopulate <- unname(unlist(allpos))
    smile[topopulate[topopulate <= nrow(microdf)]] <- 1
  }
  microdf[[paste0(wvprefix, nlags, "b", nbefore, "a", nafter)]] <- smile
  return(microdf)
}

sdf$splitbasis <- with(sdf, paste(ID, run, trial, sep="."))

splitdf <- split(sdf, sdf$splitbasis)
splitdf <- lapply(splitdf, function(microdf) {
  microdf %>% get_wv_smile(n=3, nbefore=0, nafter=1) %>%
    get_wv_smile(n=3, nbefore=1, nafter=1) %>%
    get_wv_smile(n=3, nbefore=1, nafter=2) %>%
    get_wv_smile(n=4, nbefore=0, nafter=1) %>%
    get_wv_smile(n=4, nbefore=1, nafter=1) %>%
    get_wv_smile(n=4, nbefore=1, nafter=2) %>%
    get_wv_smile(n=1, nbefore=1, nafter=1, spec=TRUE) %>%
    get_wv_smile(n=2, nbefore=1, nafter=1, spec=TRUE) %>%
    get_wv_smile(n=3, nbefore=1, nafter=1, spec=TRUE) %>%
    get_wv_smile(n=4, nbefore=1, nafter=1, spec=TRUE)
})

bb <- bind_rows(splitdf)

fbb <- bb %>% filter(bin >10 & bin <35)

#spot check
bb %>% filter(trial > 5) %>% select(bin, timestep, timesteplag1, timesteplag2, wv3b0a1)
# table(bb$wv3b0a1)
# table(bb$wv4b0a1)

### end wv smiles


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

###########
# simple models 
# with single factors
###########

# DAN mean beta
wv_dan1 <- coxme(Surv(t1,t2,response) ~ wvs1b1a1*DAN + wvs2b1a1*DAN + wvs3b1a1*DAN +
                           value_wi*DAN + uncertainty_wi*DAN + 
                           (1|ID), bb)
summary(wv_dan1)
Anova(wv_dan1, '3')

# fef from bifactor model
summary(wv_fef <- coxme(Surv(t1,t2,response) ~ wvs1b1a1*fef + wvs2b1a1*fef + wvs3b1a1*fef + 
                           value_wi*fef + uncertainty_wi*fef + 
                           (1|ID), bb))
Anova(wv_fef, '3')

# parietal from bifactor model
summary(wv_med_par <- coxme(Surv(t1,t2,response) ~ wvs1b1a1*med_par + wvs2b1a1*med_par + wvs3b1a1*med_par + 
                          value_wi*med_par + uncertainty_wi*med_par + 
                          (1|ID), bb))
Anova(wv_med_par, '3')

# general DAN from bifactor model
summary(wv_ge <- coxme(Surv(t1,t2,response) ~ wvs1b1a1*general_entropy + wvs2b1a1*general_entropy + wvs3b1a1*general_entropy + 
                          value_wi*general_entropy + uncertainty_wi*general_entropy + 
                          (1|ID), bb))
Anova(wv_ge, '3')

# vlPFC mean
summary(wv_vlpfc <- coxme(Surv(t1,t2,response) ~ wvs1b1a1*entropy_vlPFC + wvs2b1a1*entropy_vlPFC + wvs3b1a1*entropy_vlPFC + 
                         value_wi*entropy_vlPFC + trial_neg_inv_sc*uncertainty_wi*entropy_vlPFC + 
                         (1|ID), bb))
Anova(wv_vlpfc, '3')

###################
# 
# Main analysis: filter the ends
summary(wv_ge_middle <- coxme(Surv(t1,t2,response) ~ wvs1b1a1*general_entropy + wvs2b1a1*general_entropy + wvs3b1a1*general_entropy + 
                         value_wi*general_entropy + uncertainty_wi*general_entropy + 
                         (1|ID), fbb))
Anova(wv_ge_middle, '3')

summary(wv_vlpfc_middle <- coxme(Surv(t1,t2,response) ~ wvs1b1a1*entropy_vlPFC + wvs2b1a1*entropy_vlPFC + wvs3b1a1*entropy_vlPFC + 
                            value_wi*entropy_vlPFC + uncertainty_wi*entropy_vlPFC + 
                            (1|ID), fbb))
summary(wv_vlpfc_middle)
Anova(wv_vlpfc_middle, '3')

summary(wv_med_par_middle <- coxme(Surv(t1,t2,response) ~ wvs1b1a1*med_par + wvs2b1a1*med_par + wvs3b1a1*med_par + 
                              value_wi*med_par + uncertainty_wi*med_par + 
                              (1|ID), fbb))
Anova(wv_med_par_middle, '3')


# add the fef and med_par factors with general entropy, God bless us!
# _c will stand for circumcised
summary(wv_ge_fef_par_c <- coxme(Surv(t1,t2,response) ~ wvs1b1a1*general_entropy + wvs2b1a1*general_entropy + wvs3b1a1*general_entropy +
                                 wvs1b1a1*fef + wvs2b1a1*fef + wvs3b1a1*fef +
                                 wvs1b1a1*med_par + wvs2b1a1*med_par + wvs3b1a1*med_par +
                                value_wi*general_entropy + uncertainty_wi*general_entropy + 
                                 value_wi*fef + uncertainty_wi*fef + 
                                 value_wi*med_par + uncertainty_wi*med_par + 
                                (1|ID), fbb))
Anova(wv_ge_fef_par_c, '3')

# interactions with trial
summary(wv_ge_middle_trial <- coxme(Surv(t1,t2,response) ~ wvs1b1a1*general_entropy*trial + wvs2b1a1*general_entropy*trial + wvs3b1a1*general_entropy*trial + 
                                value_wi*general_entropy + uncertainty_wi*general_entropy*trial + 
                                (1|ID), fbb))
Anova(wv_ge_middle_trial, '3')

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


# just the SFG blobs
summary(wv_sfgl_middle <- coxme(Surv(t1,t2,response) ~ wvs1b1a1*dan_l_sfg + wvs2b1a1*dan_l_sfg + wvs3b1a1*dan_l_sfg + 
                                   value_wi*dan_l_sfg + uncertainty_wi*dan_l_sfg + 
                                   (1|ID), fbb))
Anova(wv_sfgl_middle, '3')

summary(wv_sfgr_middle <- coxme(Surv(t1,t2,response) ~ wvs1b1a1*dan_r_sfg + wvs2b1a1*dan_r_sfg + wvs3b1a1*dan_r_sfg + 
                                  value_wi*dan_r_sfg + uncertainty_wi*dan_r_sfg + 
                                  (1|ID), fbb))
Anova(wv_sfgr_middle, '3')
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
