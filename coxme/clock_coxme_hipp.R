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

# add hippocampal betas
setwd('~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/')
load('trial_df_and_vh_pe_clusters_u.Rdata')
hdf <- df %>% select (ID, pe_f2_hipp, h_HippAntL_neg) %>% unique()
sdf <- sdf %>% inner_join(hdf, by = "ID")

# mark events (responses)
sdf$response <- round(sdf$rt/1000, digits = 1)==sdf$t2


# inspect piecewize hazard functions
library(muhaz)
# piecewise <- pehaz(ddf$latency, delta=ddf$quit, width=NA, min.time=0, max.time=20.1)
# plot(piecewise, xlab="Time, seconds", ylab="Quit Rate")
pwI <- pehaz(na.omit(bdf$rt[bdf$rewFunc=='IEV']), width=200, min.time=0, max.time=4000)
pwD <- pehaz(na.omit(bdf$rt[bdf$rewFunc=='DEV']), width=200, min.time=0, max.time=4000)
pwR <- pehaz(na.omit(bdf$rt[bdf$rewFunc=='CEVR']), width=200, min.time=0, max.time=4000)
pwC <- pehaz(na.omit(bdf$rt[bdf$rewFunc=='CEV']), width=200, min.time=0, max.time=4000)

h <- c(pwI$Hazard,pwD$Hazard,pwC$Hazard,pwR$Hazard)
times <- c(pwI$Cuts[2:length(pwI$Cuts)],pwD$Cuts[2:length(pwD$Cuts)],pwC$Cuts[2:length(pwC$Cuts)],pwR$Cuts[2:length(pwR$Cuts)])
condition <- c(rep('IEV',length(pwI$Hazard)), rep('DEV',length(pwD$Hazard)),rep('CEV',length(pwC$Hazard)),rep('CEVR',length(pwR$Hazard)))
H <- as.tibble(cbind(h,times,condition))
H$h <- as.numeric(H$h)
H$logh <- log(H$h)
H$times <- as.numeric(H$times)
# pdf("piecewise_hazard_fx_all_groups.pdf", width = 12, height = 4)
ggplot(H,aes(times, h, color = condition)) + geom_line()
ggplot(H,aes(times, logh, color = condition)) + geom_line()
# not what I expected: shared underlying hazard across contingencies

badfit <- survfit(Surv(t2) ~ rewFunc, type = 'right', origin = .1, data=sdf)
badfit1 <- survfit(Surv(rt) ~ rewFunc,  data=bdf)


plot(badfit1, mark.time=FALSE, lty=1:4,
     xlab="ms", ylab="Proportion still waiting")
legend(3000, .85, c("CEV", "CEVR", "DEV", "IEV"),
       lty=1:4, bty='n')

# add AH and PH, first without run
summary(cox_hipp <- coxme(Surv(t1,t2,response) ~ rtlag*pe_f2_hipp + rtlag*h_HippAntL_neg + trial*rewFunc + trial*uncertainty +
                      value*pe_f2_hipp + value*h_HippAntL_neg + uncertainty*h_HippAntL_neg +uncertainty*pe_f2_hipp + (1|ID), sdf))
summary(cox_hipp_sc <- coxme(Surv(t1,t2,response) ~ rtlag*pe_f2_hipp + rtlag*h_HippAntL_neg + trial*rewFunc + trial*scale(uncertainty) +
                            scale(value)*pe_f2_hipp + scale(value)*h_HippAntL_neg + scale(uncertainty)*h_HippAntL_neg + scale(uncertainty)*pe_f2_hipp + (1|ID), sdf))

out <- capture.output(summary(cox_hipp))
setwd('~/OneDrive/collected_letters/papers/sceptic_fmri/hippo/supp/')
write.csv(out,file = "coxme_hipp_summary.csv", sep = ":\t")

# devtools::install_github('junkka/ehahelper')
library(ehahelper)
library(broom)
library(broom.mixed)
tidy_cox_hipp <- tidy(cox_hipp, exponentiate = F)[,1:5]
stargazer(tidy_cox_hipp, type = 'html', out = 'tidy_cox_hipp.html', summary = F, digits = 3, digits.extra = 10)
write_csv(tidy_cox_hipp, 'tidy_cox_hipp.csv')
p1 <- plot_model(cox_hipp,transform = 'exp', terms = c("pe_f2_hipp:value", "h_HippAntL_neg:value", "h_HippAntL_neg:uncertainty", "pe_f2_hipp:uncertainty"), axis.lim = c(.99, 1.01))
p1 + ylim(c(.99, 1.025))

# no interesting interactions with trial, inferior model
summary(cox_hipp1 <- coxme(Surv(t1,t2,response) ~ rtlag*pe_f2_hipp + rtlag*h_HippAntL_neg + trial*rewFunc + 
                             trial*value*pe_f2_hipp + trial*value*h_HippAntL_neg + trial*uncertainty*h_HippAntL_neg + trial*uncertainty*pe_f2_hipp + (1|ID), sdf))


# interaction with trial
summary(c5a <- coxme(Surv(t1,t2,response) ~ scale(rtlag) + scale(value)*scale(trial) + scale(value)*scale(uncertainty) + scale(entropylag) +  (1|ID/run), sdf))


# limit analysis to middle 3s to eliminate speed constraints and avoidance of interval end: U-aversion holds
# remove distance from the edge
summary(c6 <- coxme(Surv(t1,t2,response) ~ rtlag + trial + value + uncertainty  + (1|ID/run), sdf[sdf$t1>.5 & sdf$t1<3.5,]))

# limit analysis to middle 2s to eliminate speed constraints and avoidance of interval end: U-aversion holds
summary(c6a <- coxme(Surv(t1,t2,response) ~ rtlag + trial + value + uncertainty + (1|ID/run), sdf[sdf$t1>1 & sdf$t1<3,]))

# limit to IEV to r/o uncertainty/value tradeoff explanation: U-aversion holds
summary(c7 <- coxme(Surv(t1,t2,response) ~ rtlag + trial + value + uncertainty +  (1|ID/run), sdf[sdf$rewFunc=='IEV' & sdf$t1>1 & sdf$t1<3,]))

# late in learning

summary(c7 <- coxme(Surv(t1,t2,response) ~ rtlag + value + uncertainty + entropylag + distfromedgelag +  (1|ID/run), sdf[sdf$rewFunc=='IEV',]))
library(s)

summary(test <- lmer(scale(uncertainty) ~ scale(value)*as.factor(trial) + (1|ID/run), sdf[sdf$rewFunc=='IEV',]))
summary(test <- lmer(scale(uncertainty) ~ scale(value)*as.factor(trial) + (1|ID/run), sdf[sdf$rewFunc=='DEV',]))
summary(test <- lmer(scale(uncertainty) ~ scale(value)*as.factor(trial) + (1|ID/run), sdf[sdf$rewFunc=='CEVR',]))

# limit to CEVR to r/o probability/magnitude tradeoff explanation: U-aversion holds
summary(c8 <- coxme(Surv(t1,t2,response) ~ rtlag + trial + value + uncertainty +  (1|ID/run), sdf[sdf$rewFunc=='CEVR',]))

summary(check <- lmer(uncertainty ~ value + (1|ID/run/trial), sdf))

summary(check <- lmer(uncertainty ~ t2 + scale(value) + (1|ID/run/trial), sdf[sdf$rewFunc=='DEV' & sdf$t1>.5,]))
summary(check <- lmer(uncertainty ~ t2 + scale(value) + (1|ID/run/trial), sdf[sdf$rewFunc=='IEV' & sdf$t1<1,]))
summary(check <- lmer(scale(uncertainty) ~ scale(value) + (1|ID/run/trial), adf[adf$rewFunc=='IEV',]))
summary(check <- lmer(scale(uncertainty) ~ scale(value) + (1|ID/run/trial), adf[adf$rewFunc=='DEV' & adf$t1>1 & adf$t1<3,]))
summary(check <- lmer(scale(uncertainty) ~ scale(value) + (1|ID/run/trial), adf[adf$rewFunc=='IEV' & adf$t1>3 & adf$t1<4,]))

# they are least correlated in the very end of IEV...
summary(c9 <- coxme(Surv(t1,t2,response) ~ rtlag + value + uncertainty + entropylag + distfromedgelag +  (1|ID/run), sdf[sdf$rewFunc=='IEV' & sdf$t1>3 & sdf$t1<4,]))

d <- adf[adf$rewFunc=='IEV' & adf$t1>3.75 & adf$t1<4,]
cor.test(d$uncertainty,d$value)

anova(c2,c3,c4, c5)

# apa library did not work
# library(apa)
out <- capture.output(summary(c5))

write.csv(out,file = "coxme_c5_summary.csv", sep = ":\t")
cat("coxme_c5", out, file="coxme_c5_summary.txt", fill = TRUE, sep=":\t", append=FALSE)

stargazer(out, type="html", out="coxme_c5_summary.htm", digits = 2,single.row=FALSE,omit.stat = "bic",
          star.char = c("+", "*", "**", "***"),
          star.cutoffs = c(0.1, 0.05, 0.01, 0.001),
          notes = c("+ p<0.1; * p<0.05; ** p<0.01; *** p<0.001"), 
          notes.append = F)


summary(c5)
#########
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
