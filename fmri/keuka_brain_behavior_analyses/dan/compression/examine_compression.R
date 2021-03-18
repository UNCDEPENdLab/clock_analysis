library(readr)
library(tidyverse)
library(data.table)
library(lme4)
library(broom.mixed)
library(ggplot2)
library(emmeans)
library(lmerTest)
source("~/Downloads/screen.lmerTest.R")
source("/Users/michael/Data_Analysis/GKY_analysis/functions/vif.lme.R")
#setwd("/proj/mnhallqlab/users/michael/fmri.pipeline/R")
setwd("/Users/michael/Data_Analysis/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/compression")

#source("get_mmy3_trial_df.R")
#source("/Users/hallquist/Data_Analysis/r_packages/fmri.pipeline/R/get_mmy3_trial_df.R")
#setwd(file.path(getMainDir(), "users", "michael", "sceptic_decon"))

#source("/Users/michael/Data_Analysis/r_packages/fmri.pipeline/R/rle_dt.R")

#compress_results <- read_csv("Schaefer_DorsAttn_2.3mm_clock_long_compression.csv.gz") %>% as.data.table()
#compress_results <- read_csv("Schaefer_DorsAttn_2.3mm_whole_trial_compression.csv.gz") %>% as.data.table()
compress_results <- read_csv("Schaefer_DorsAttn_2.3mm_rt8_compression.csv.gz") %>% filter(decon_compress_nt > 8) %>% as.data.table()
compress_results[, atlas := NULL]

#swiss cheese problem for now
xtabs(~id + run, compress_results) #not too whole-y

z_score <- function(x) { (x - mean(x, na.rm=TRUE))/sd(x, na.rm=TRUE) }

trial_df <- readRDS("/Users/michael/Data_Analysis/clock_analysis/fmri/mmy3_trial_df_selective_groupfixed.rds") %>%
#trial_df <- get_mmy3_trial_df(model="selective", groupfixed=TRUE) %>%
  mutate(rt_time=clock_onset + rt_csv/1000, #should be pretty redundant with isi_onset, but still
    rt_vmax=rt_vmax/10, #to put into seconds
    rt_vmax_cum=clock_onset + rt_vmax,
    rt_cur = rt_csv/1000,
    inv_run_trial = -100/run_trial,
    outcome=factor(score_vba > 0, levels=c(FALSE, TRUE), labels=c("Omission", "Reward"))) %>%
  group_by(id, run) %>%
  mutate(
    rt_lag=lag(rt_cur, order_by=run_trial),
    rt_next=lead(rt_cur, 1L, order_by=run_trial),
    rt_swing=abs(rt_cur - rt_lag),
    iti_lag=lag(iti_ideal, order_by=run_trial),
    outcome_lag=lag(outcome, 1L, order_by=run_trial),
    v_entropy_wi = v_entropy - mean(v_entropy, na.rm=TRUE),
    v_entropy_bw = mean(v_entropy, na.rm=TRUE),
    v_entropy_wi_z = z_score(v_entropy_wi),
    v_entropy_wi_z_lag = lag(v_entropy_wi_z, 1L, order_by=run_trial),
    v_max_wi_z = z_score(v_max),
    v_max_wi_z_lag = lag(v_max_wi_z, 1L, order_by=run_trial)
    ) %>% ungroup() %>%
  select(id, run, trial, run_trial, inv_run_trial, rewFunc, emotion, outcome_lag,
         rt_cur, rt_lag, rt_next, ev, rt_vmax, rt_vmax_lag, rt_swing,
         v_entropy, v_entropy_change, v_entropy_bw, v_entropy_wi, v_entropy_wi_z, v_entropy_wi_z_lag,
         v_max_wi_z, v_max_wi_z_lag,
         pe_max, iti_ideal, iti_lag)

combined <- left_join(compress_results, trial_df, by=c("id", "run", "trial")) %>%
  rename(acom=decon_acompress_0.9, com=decon_compress_0.9) %>%
  group_by(id, run) %>%
  mutate(
    acom_bw=mean(acom, na.rm=TRUE),
    acom_wi_z = z_score(acom),
    acom_wi_z_lag=lag(acom_wi_z, 1, order_by=run_trial),
    com_lag=lag(com, 1, order_by=run_trial)
  ) %>% ungroup()


#scale things at once
combined <- combined %>% mutate(
  across(.cols=c(rt_cur, rt_lag, rt_next, rt_vmax, rt_vmax_lag, inv_run_trial, run_trial, acom,
                 v_entropy, v_entropy_bw, v_entropy_change, 
                 pe_max, iti_ideal, iti_lag, acom_bw), .fns=list(z=z_score))
)

#contents of acom
m1 <- lmer(acom ~ iti_lag + iti_ideal + rt_cur + rt_lag + (1|id/run) + (1|label), combined %>% filter(decon_compress_nt==16))
summary(m1)

#combined_rle <- rle_dt$new(data=combined, keys=c("id", "run", "trial", "atlas_value", "rewFunc", "emotion"), optimize_order = FALSE)

avs <- sort(unique(combined$atlas_value))

# roi_labels <- read.table("/Users/michael/Data_Analysis/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/Schaefer2018_200Parcels_7Networks_order_manual.txt", sep="\t")
# names(roi_labels)[1:3] <- c("atlas_value", "schaefer_label", "manual_label")
# roi_labels <- roi_label

labels <- as_tibble(read_delim("~/Data_Analysis/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/Schaefer2018_200Parcels_DAN_order_manual.txt", 
                               "\t", escape_double = FALSE, col_names = FALSE, 
                               trim_ws = TRUE)) %>% select(1:4)
names(labels) <- c("atlas_value", "label_long", "label_short", "entropy_signal")
labels <- labels %>% filter(grepl("DorsAtt", label_long))  %>% mutate(
  side  = case_when(
    grepl("LH", label_long) ~ "L",
    grepl("RH", label_long) ~ "R"
  ),
  label_short_side = paste(label_short, side, sep ="_"),
  label_long1 = substr(label_long, 23, 100),
  label = case_when(
    label_short!="0" ~ label_short_side,
    label_short=="0" ~ label_long1
  )
) %>% select(c(label, side, atlas_value))

combined <- combined %>% left_join(labels, by="atlas_value")


#force cleanup
#rm(compress_results, combined)
#gc()



# toplot <- combined_rle$get() %>% filter(atlas_value %in% c(41, 145, 33, 141, 31, 135)) %>%
#   mutate(atlas_value=factor(atlas_value, levels=c(33, 41, 141, 145, 31, 135), labels=c("L_LIPd", "L_FEF", "R_LIPd", "R_FEF", "L_MT", "R_MT")))
# 
# 
# rfef <- combined_rle$get() %>% filter(atlas_value == 145)

# toplot <- combined %>% filter(atlas_value %in% c(41, 145, 33, 141, 31, 135)) %>%
#   mutate(atlas_value=factor(atlas_value, levels=c(33, 41, 141, 145, 31, 135), labels=c("L_LIPd", "L_FEF", "R_LIPd", "R_FEF", "L_MT", "R_MT")))

rfef <- combined %>% filter(atlas_value == 145)



m_rfef <- lmer(acom ~ v_entropy + (1 | id/run), data=rfef)

m_rfef2 <- lmer(acom ~ v_entropy_change*rewFunc + rt_cur + 
                  run_trial + decon_compress_nt + (1 | id/run), data=rfef)

#include prior and current ITI as covariates to control for how much of each period
#is factored into the number of timepoints... seems like this may be too collinear!
m_rfef3 <- lmer(acom ~ v_entropy_change*rewFunc + rt_cur + 
                 run_trial + decon_compress_nt + iti_ideal + iti_lag + (1 | id/run), data=rfef)

car::vif(m_rfef2)

summary(m_rfef2)
summary(m_rfef3)
car::Anova(m_rfef3, type=3)
emtrends(m_rfef3, ~rewFunc, var="v_entropy_change")

#look at bw/wi
#rewFunc is confounded with v_entropy_bw
m_rfef4 <- lmer(acom ~ v_entropy_wi*rewFunc + rt_cur + v_entropy_bw + #*rewFunc
                  v_entropy_change + run_trial + rt_cur*iti_ideal + rt_cur*iti_lag + (1 | id/run), data=rfef)

#decon_compress_nt

#entropy bw effect -- people who have higher average entropy in a run tend to have more
#compressed representations in FEF...
summary(m_rfef4)
car::Anova(m_rfef4, type=3)
car::vif(m_rfef4, type=3)

#for whole-trial, confirmed that iti_lag is not influential in nt, whereas iti and rt account for 95% of variance
summary(lm(decon_compress_nt ~ iti_ideal + rt_cur, rfef)) #iti_lag + 

#1) ITI-only windowing -- maybe 3+ seconds? (~36%)
#2) decompose NT into online/offline, or pre, online, post
#3) window something like +3 from last feedback (or clock, whichever is first) through + 10 (next clock as collider)
#4) windows of equal length aligned to different events (mid-clock, rt+1, mid-iti)
#5) compare against rest?
#6) carve out and concatenate ITI periods, windowed SVD, ignoring trial

#continuity of compression? yes, big-time
m_rfef5 <- lmer(acom ~ acom_lag*v_entropy_change + v_entropy_wi*rewFunc + rt_cur + v_entropy_bw + #*rewFunc
                  v_entropy_change + run_trial + decon_compress_nt + iti_ideal + iti_lag + (1 | id/run), data=rfef)

#should not include nt and iti in same model under whole-trial
m_rfef5 <- lmer(acom ~ acom_lag*v_entropy_change + v_entropy_wi*rewFunc + rt_cur + v_entropy_bw + #*rewFunc
                  v_entropy_change + run_trial + iti_ideal + iti_lag + (1 | id/run), data=rfef)


summary(m_rfef5)
car::Anova(m_rfef5, type=3)
car::vif(m_rfef5, type=3)


m_rfef6 <- lmer(acom ~ acom_lag + v_entropy_wi*rewFunc + rt_cur + v_entropy_bw + #*rewFunc
                  v_entropy_change + run_trial + decon_compress_nt + iti_ideal + iti_lag + (1 | id/run), data=rfef)

car::Anova(m_rfef6, type=3)
car::vif(m_rfef4, type=3)

#rt modeling
m_rt1 <- lmer(rt_cur ~ acom_lag + v_entropy_wi + v_entropy_bw + rewFunc + rt_vmax_lag +
                  v_entropy_change + run_trial + (1 | id/run), data=rfef) # iti_ideal + iti_lag +

summary(m_rt1)
car::Anova(m_rt1, type=3)
car::vif(m_rt1, type=3)

#what is relationship of entropy_wi with acom_lag?
# m_e1 <- lmer(v_entropy_wi ~ acom_lag + run_trial + (1 | id/run), data=rfef) # iti_ideal + iti_lag +
# summary(m_e1) #huh, not significant... and in RT model where both are main effects, we are okay

#bring in AR(1) RT
#huh, so more compressed representations in R FEF go with longer RTs
m_rt2 <- lmer(rt_cur ~ rt_lag * acom_lag + rewFunc + rt_vmax_lag * acom_bw + rt_vmax_lag * acom_lag +  #v_entropy_wi + v_entropy_bw + v_entropy_change +
                run_trial + decon_compress_nt + (1 | id/run), data=rfef) # iti_ideal + iti_lag +
#run*acom_bw +

m_rt2 <- lmer(rt_cur ~ iti_ideal + iti_lag + rewFunc + run_trial + (1 | id/run), data=rfef) # 


summary(m_rt2)
car::Anova(m_rt2, type=3)
car::vif(m_rt2, type=3)

#1) one complication of decon_compress_nt: it contains information about the current RT!
#     - this yields criterion contamination


m_rt2 <- lmer(rt_next ~ rt_lag * acom_lag + rewFunc + rt_vmax_lag * acom_bw + rt_vmax_lag * acom_lag +  #v_entropy_wi + v_entropy_bw + v_entropy_change +
                run_trial + decon_compress_nt + (1 | id/run), data=rfef) # iti_ideal + iti_lag +


###
#main beta-behavior analysis from NComms paper
library(lmerTest)
mb3hpe_hipp <-  lmer(rt_cur ~ (inv_run_trial_z + rt_lag_z + rt_vmax_lag_z + outcome_lag + 
                                    v_max_wi_z_lag + v_entropy_wi_z + iti_ideal_z + acom_lag_z + acom_bw_z)^2 + (1|id/run), rfef)
screen.lmerTest(mb3hpe_hipp, .05)
summary(mb3hpe_hipp)
Anova(mb3hpe_hipp, '3')
car::vif(mb3hpe_hipp)

tt <- tidy(mb3hpe_hipp) %>% filter(abs(statistic) > 2)
###


pdf("compress_preview_approx.pdf", width=15, height=15)
ggplot(toplot, aes(x=run_trial, y=acom, color=rewFunc)) + geom_smooth(method="loess") +
  facet_wrap(~atlas_value)
dev.off()


mlist <- list()
coef_list <- list()
avs <- sort(unique(combined$label))
#ANOMALY 1: IF iti_ideal_z is included, we see a big interaction with ACOM_LAG

for (aa in avs) {
  df <- combined %>% filter(label==!!aa & !is.na(v_entropy))
  mm <- lmer(rt_cur ~ (inv_run_trial_z + rt_lag_z + rt_vmax_lag_z + outcome_lag + 
                         v_max_wi_z_lag + v_entropy_wi_z + acom_wi_z_lag + acom_bw_z)^2 + (1|id/run), df) #iti_lag_z + 
  mlist[[as.character(aa)]] <- mm
  coef_list[[as.character(aa)]] <- tidy(mm) 
}

#plot matrix of coefficients
coef_df <- do.call(rbind, lapply(1:length(coef_list), function(x) { coef_list[[x]]$region <- names(coef_list)[x]; return(coef_list[[x]]) }))

coef_df <- coef_df %>% filter(grepl(".*acom.*", term)) %>% filter(p.value < .1)

ggplot(coef_df, aes(x=term, y=region, fill=statistic)) + geom_tile() + scale_fill_viridis_c() + theme(axis.text.x = element_text(angle=90))
ggsave(file="acom_preview_no_itilag.pdf", width=15, height=12)

screen.lmerTest(mlist$`1_f_6a_R`)

cor(df$acom_wi_z, df$rt_swing, use="pairwise")


#CCF examination
df <- combined %>% filter(label=="1_f_6a_R")

#g1 <- ggplot(df, aes(x=run_trial, y=v_entropy_wi_z)) + geom_smooth(method="loess") + facet_wrap(~rewFunc)
g1 <- ggplot(df, aes(x=run_trial, y=rt_swing)) + geom_smooth(method="loess") + facet_wrap(~rewFunc)
g2 <- ggplot(df, aes(x=run_trial, y=v_entropy_wi_z)) + geom_smooth(method="loess") + facet_wrap(~rewFunc)
g3 <- ggplot(df, aes(x=run_trial, y=acom_wi_z)) + geom_smooth(method="loess") + facet_wrap(~rewFunc)
g4 <- ggplot(df, aes(x=run_trial, y=acom)) + geom_smooth(method="loess") + facet_wrap(~rewFunc)

library(cowplot)
pdf("trial compression comparisons with acomraw 2.pdf", width=15, height=9)
cowplot::plot_grid(g1, g2, g3, g4, nrow=1)
dev.off()

#library(lattice)
pnl <- function(x, y = x, maxlag=10) { 
  par(new = TRUE)
  ccf(x, y, lag.max=maxlag, na.action=na.pass)
}
# loess.pnl <- function(x, y, ...) { 
#   #panel.grid(h = -1, v = -1) 
#   par(new=TRUE)
#   panel.xyplot(x, y, ...) 
#   panel.loess(x, y, ..., col = "black")
#   panel.rug(x = x[is.na(y)], y = y[is.na(x)])
# }
loess.pnl <- function(x, y, span=1, lwd=2, ...) { 
  panel.smooth(x, y, span=span, lwd=lwd, iter=10, ...) 
}

#run variance is already partialed out by run centering (avoid singular fit)
acom_resid <- lmer(acom_wi_z ~ 1 + run_trial + rt_cur + (1|id), df, na.action=na.exclude)
summary(acom_resid)

acom_resid_iti <- lmer(acom_wi_z ~ 1 + run_trial + iti_ideal + rt_cur + (1|id), df, na.action=na.exclude)
summary(acom_resid_iti)

toplot <- df %>% select(id, run, run_trial, rewFunc, acom_wi_z, v_entropy_change_z, v_entropy_wi_z, rt_swing, rt_lag_z, pe_max) %>% # %>% na.omit()
  mutate(acom_resid=resid(acom_resid), acom_resid_iti=resid(acom_resid_iti))

plot_split <- split(toplot, list(toplot$id, toplot$run))
#pdf("ccfs/acom_ccfs.pdf", width=12, height=12)
acom_entropy_ccfs <- list()
acom_entropy_ccfs_resid <- list()
for (ss in 1:length(plot_split)) {
  subdf <- plot_split[[ss]]
  if (nrow(subdf) == 0L) { next }
  #pdf(paste0("ccfs/id", ss$id, "_run", ss$run, "ccf.pdf"), width=12, height=12)
  #pairs(subdf %>% select(-id, -run, -run_trial, -rewFunc), upper.panel = pnl, lower.panel=loess.pnl, diag.panel = pnl, cex.labels = 1, gap=3, 
  #      main=paste0("id ", subdf$id, ", run ", subdf$run, ", cont ", subdf$rewFunc))
  #dev.off()
  ccres <- ccf(subdf$acom_wi_z, subdf$v_entropy_wi_z, lag.max=8, na.action=na.pass, plot = FALSE)
  acom_entropy_ccfs[[ss]] <- data.frame(id=subdf$id[1], run=subdf$run[1], rewFunc=subdf$rewFunc[1], lag=as.vector(ccres$lag), ccf=as.vector(ccres$acf))
  
  ccres <- ccf(subdf$acom_resid_iti, subdf$v_entropy_wi_z, lag.max=8, na.action=na.pass, plot = FALSE)
  acom_entropy_ccfs_resid[[ss]] <- data.frame(id=subdf$id[1], run=subdf$run[1], rewFunc=subdf$rewFunc[1], lag=as.vector(ccres$lag), ccf=as.vector(ccres$acf))
}
#dev.off()

ccdf <- bind_rows(acom_entropy_ccfs)
ccdf_resid <- bind_rows(acom_entropy_ccfs_resid)

pdf("acom entropy ccf avg.pdf", width=8, height=8)
ggplot(ccdf, aes(x=lag, y=ccf)) + stat_summary(fun.data = mean_cl_boot, geom = "pointrange", fun.args=list(conf.int=.9)) +
  xlab("ACom_wi - Entropy_wi (future entropy [neg] to past entropy [pos]") + facet_wrap(~rewFunc) + geom_hline(yintercept=0)
dev.off()

pdf("acom entropy ccf resid avg.pdf", width=8, height=8)
ggplot(ccdf_resid, aes(x=lag, y=ccf)) + stat_summary(fun.data = mean_cl_boot, geom = "pointrange", fun.args=list(conf.int=.9)) +
  xlab("ACom_wi - Entropy_wi (future entropy [neg] to past entropy [pos]") + facet_wrap(~rewFunc) + geom_hline(yintercept=0)
dev.off()


ccres <- ccf(subdf$acom_wi_z, subdf$v_entropy_wi_z, lag.max=8, na.action=na.pass, plot = FALSE)
ccres
cvec <- c()
for (ll in -8:8) {
  #cvec <- c(cvec, cor(subdf$v_entropy_wi_z, Hmisc::Lag(subdf$acom_wi_z, ll), use="pairwise"))
  cvec <- c(cvec, cor(Hmisc::Lag(subdf$v_entropy_wi_z, ll), subdf$acom_wi_z, use="pairwise"))
}

#-1 is entropy on next trial, +1 is entropy on prior trial
cbind(lag_m1=Hmisc::Lag(subdf$v_entropy_wi_z, -1), lag_0=subdf$v_entropy_wi_z, lag_p1=Hmisc::Lag(subdf$v_entropy_wi_z, 1))

#conclusion: for CCF, negative lags mean current acom correlated with future entropy; positive lags mean current acom correlated with past entropy
cor(cvec, as.vector(ccres$acf))
cbind(cvec, ccres$acf, lag=-8:8)

cor(lag(subdf$v_entropy_wi_z, 1), subdf$acom_wi_z, use="pairwise")

#DECON model-free analysis
m1 <- lmer(acom ~ outcome_lag + (1|id/run), df) #+ rt_lag + rt_cur +
summary(m1)
m1 <- lmer(acom ~ iti_ideal + rt_cur + decon_nt + (1|id/run), df) #+ rt_lag + rt_cur +
summary(m1)
m1 <- lmer(acom ~ iti_ideal + rt_cur + v_entropy_wi_z + (1|id/run), df) #+ rt_lag + rt_cur +
summary(m1)

m1 <- lmer(ev ~ acom_z*run_trial_z + iti_ideal + rt_cur + v_entropy_wi_z + (1|id/run), df) #+ rt_lag + rt_cur +
vif.lme(m1)
summary(m1)

avs <- sort(unique(combined$label))

library(sjPlot)
coef_list <- list()
for (aa in avs) {
  df <- combined %>% filter(label==!!aa)
  #decon_compress_nt is essentially equal to iti_ideal + rt_cur, so don't include it
  #v_entropy_wi_z
  #m1 <- lmer(ev ~ acom_z*run_trial_z + iti_ideal + rt_cur*rewFunc + acom_z*rewFunc + (1|id/run), df) #+ rt_lag + rt_cur +
  m1 <- lmer(ev ~ acom_wi_z*run_trial_z + acom_bw_z*run_trial_z + iti_ideal + rt_cur + (1|id/run), df) #+ rt_lag + rt_cur +
  summary(m1)
  #vif.lme(m1)
  #plot_model(m1, type="pred", terms = c("run_trial_z", "acom_z"))
  coef_list[[as.character(aa)]] <- tidy(m1) %>% mutate(roi=!!aa)
}

coef_df <- bind_rows(coef_list)
coef_df <- coef_df %>% filter(grepl(".*acom.*", term)) %>% filter(p.value < .1)

ggplot(coef_df, aes(x=term, y=roi, fill=statistic)) + geom_tile() + scale_fill_viridis_c() + theme(axis.text.x = element_text(angle=90))
ggsave(file="acom_ev_effects_decomposed.pdf", width=15, height=12)

plot_model(m1, type="pred", terms = c("run_trial_z", "acom_wi_z"))


## decon as outcome models
avs <- sort(unique(combined$label))
coef_list <- list()
for (aa in avs) {
  df <- combined %>% filter(label==!!aa)
  #decon_compress_nt is essentially equal to iti_ideal + rt_cur, so don't include it
  #v_entropy_wi_z
  m1 <- lmer(acom_z ~ rt_lag_z + run_trial_z + iti_ideal + rt_cur + acom_wi_z_lag*outcome_lag*run_trial_z + rewFunc*acom_wi_z_lag*v_entropy_wi_z_lag  + (1|id/run), df) #+ rt_lag + rt_cur +
  #m1 <- lmer(ev ~ acom_wi_z*run_trial_z + acom_bw_z*run_trial_z + iti_ideal + rt_cur + (1|id/run), df) #+ rt_lag + rt_cur +
  summary(m1)
  car::Anova(m1, 3)
  emtrends(m1, ~rewFunc, var="acom_wi_z_lag")
  vif.lme(m1)
  #plot_model(m1, type="pred", terms = c("run_trial_z", "acom_z"))
  coef_list[[as.character(aa)]] <- tidy(m1) %>% mutate(roi=!!aa)
}

#m1 <- lmer(acom_z ~ rt_lag_z + run_trial_z + iti_ideal + rt_cur + acom_wi_z_lag*outcome_lag*run_trial_z + acom_wi_z_lag*v_entropy_wi_z_lag + rewFunc + (1|id/run), df) #+ rt_lag + rt_cur +

coef_df <- bind_rows(coef_list)
coef_df <- coef_df %>% filter(grepl(".*(acom|entropy).*", term)) %>% filter(term != "acom_wi_z_lag") #%>% filter(p.value < .05)

ggplot(coef_df, aes(x=term, y=roi, fill=statistic)) + geom_tile() + scale_fill_viridis_c() + theme(axis.text.x = element_text(angle=90))
ggsave(file="acom_outcome.pdf", width=15, height=12)

##hand-code mlVAR-type approach
m1 <-  bf(lh1 ~  lh1_l + lh3_l + lh5_l + lh7_l + lh9_l + lh11_l + (1|int|id) + (0 + lh1_l + lh3_l + lh5_l + lh7_l + lh9_l + lh11_l |slo|id))
m2 <-  bf(lh3 ~  lh1_l + lh3_l + lh5_l + lh7_l + lh9_l + lh11_l + (1|int|id) + (0 + lh1_l + lh3_l + lh5_l + lh7_l + lh9_l + lh11_l |slo|id))

lm3_lh_6slc_odd <- brm(b1 + b3 + b5 + b7 + b9 + b11 + set_rescor(TRUE), data = fb_l_lags,
                       chains = 4, cores = 4, autocor = cor_arma(~run_trial|id, 1), iter=2500, init_r=0.9)

summary(lm3_lh_6slc_odd)


library(mlVAR)
mlv_list <- list()
for (aa in avs) {
  df <- combined %>% filter(label==!!aa)
  #decon_compress_nt is essentially equal to iti_ideal + rt_cur, so don't include it
  #v_entropy_wi_z
  fit1 <- mlVAR(df, vars = c("acom", "v_entropy"), idvar = "id", lags = 3, dayvar="run",
                temporal = "correlated", contemporaneous="correlated")
  
  mlv_list[[as.character(aa)]] <- fit1
}


fit1 <- mlVAR(combined, vars = c("acom", "v_entropy"), idvar = "id", lags = 3, 
              temporal = "correlated", contemporaneous="fixed")

df <- combined %>% filter(label=="1_f_6a_R")

acom_resid <- lmer(acom_wi_z ~ 1 + run_trial + rt_cur + (1|id), df, na.action=na.exclude)
summary(acom_resid)

df <- df %>% mutate(acom_resid=resid(acom_resid)) %>%
  group_by(id, run) %>%
  mutate(
    acom_bw=mean(acom, na.rm=TRUE),
    acom_wi_z = z_score(acom),
    acom_wi_z_lag=lag(acom_wi_z, 1, order_by=run_trial),
    
  )

get_acom_resid <- function(df) { 
  mm <- lmer(acom_wi_z ~ 1 + rt_cur + (1|id), df, na.action=na.exclude)
  resid(mm)
}

run_compress <- combined %>% 
  mutate(acom_residrt=get_acom_resid(.)) %>%group_by(id, run, label) %>% 
  summarize(acom_median=median(acom, na.rm=TRUE), acom_mean=mean(acom, na.rm=TRUE),
         acom_residrt_median=median(acom_residrt, na.rm=TRUE), acom_residrt_mean=mean(acom_residrt, na.rm=TRUE)) %>%
  ungroup()

write_csv(run_compress, file="run_level_compression.csv.gz")

