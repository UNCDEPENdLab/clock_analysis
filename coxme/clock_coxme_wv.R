# runs mixed-effects Cox models on clock data
# when running the first time, first run compute_sceptic_fmri_statistics.R

setwd("~/code/clock_analysis/coxme")
library(readr)
library(ggplot2)
library(tidyr)
library(tibble)
library(reshape2)
library(emmeans)
# library(factoextra)
# library(ggfortify)
# library(RColorBrewer)
# library(MASS)
# library(readr)
# library(VIM)
# library(mice)
# library(multcompView)
# library(stargazer)
library(dplyr)
library(lme4)
library(survival)
library(coxme)
library(survminer)
# library(OIsurv)
library(ggpubr)

load(file="clock_for_coxme_value_only_070518.RData")

# sanity checks on within-trial matrices

# add hippocampal betas
setwd('~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/')
load('trial_df_and_vh_pe_clusters_u.Rdata')
hdf <- df %>% select (ID, pe_f2_hipp, h_HippAntL_neg) %>% unique()
sdf <- sdf %>% inner_join(hdf, by = "ID")

# mark events (responses)
sdf$response <- round(sdf$rt/1000, digits = 1)==sdf$t2

# KS across conditions

ks.test(df$rt_csv[df$rewFunc=='CEV'], "punif", 1, 4000)
ks.test(df$rt_csv[df$rewFunc=='CEVR'], "punif", 1, 4000)
ks.test(df$rt_csv[df$rewFunc=='IEV'], "punif", 1, 4000)
ks.test(df$rt_csv[df$rewFunc=='DEV'], "punif", 1, 4000)

hist(df$rt_csv[df$rewFunc=='CEV'])
hist(df$rt_csv[df$rewFunc=='CEVR'])
hist(df$rt_csv[df$rewFunc=='IEV'])
hist(df$rt_csv[df$rewFunc=='DEV'])


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
sdf <- sdf %>% group_by(ID, run, trial) %>% mutate(bin = 1:n(), time = bin/10) %>% ungroup() %>% 
  group_by(ID) %>% mutate(value_wi = scale(value),
                          uncertainty_wi = scale(uncertainty),
                          value_b = mean(value),
                          uncertainty_b = mean(uncertainty), 
                          trial_neg_inv_sc = scale(-1/trial)) %>% ungroup() %>% 
  mutate(rtlag_sc = scale(rtlag),
         AH_sc = scale(h_HippAntL_neg))
# filter out no-go zones at the edges
fdf <- sdf %>% filter(bin >10 & bin <35)

# diagnostics on uncertainty and value distributions

ggplot(sdf %>% filter(bin >10 & bin <35), aes(time, value_wi, color = rewFunc)) + geom_smooth(method = 'gam', formula = y~splines::ns(x,2)) 
ggplot(sdf %>% filter(bin >10 & bin <35), aes(time, uncertainty_wi, color = rewFunc)) + geom_smooth(method = 'gam', formula = y~splines::ns(x,2))
corr.test(fdf$value_wi, fdf$uncertainty_wi)
corr.test(fdf$value, fdf$uncertainty)


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
# # devtools::install_github('junkka/ehahelper')
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
