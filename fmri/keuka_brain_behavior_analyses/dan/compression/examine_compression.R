library(readr)
library(tidyverse)
library(data.table)
library(lme4)
library(broom.mixed)
library(ggplot2)
library(emmeans)

#setwd("/proj/mnhallqlab/users/michael/fmri.pipeline/R")
setwd("/Users/michael/Data_Analysis/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/compression")

#source("get_mmy3_trial_df.R")
#source("/Users/hallquist/Data_Analysis/r_packages/fmri.pipeline/R/get_mmy3_trial_df.R")
#setwd(file.path(getMainDir(), "users", "michael", "sceptic_decon"))

#source("/Users/michael/Data_Analysis/r_packages/fmri.pipeline/R/rle_dt.R")

compress_results <- read_csv("Schaefer_DorsAttn_2.3mm_clock_long_compression.csv.gz") %>% as.data.table()
compress_results[, atlas := NULL]

trial_df <- readRDS("/Users/michael/Data_Analysis/clock_analysis/fmri/mmy3_trial_df_selective_groupfixed.rds") %>%
#trial_df <- get_mmy3_trial_df(model="selective", groupfixed=TRUE) %>%
  mutate(rt_time=clock_onset + rt_csv/1000, #should be pretty redundant with isi_onset, but still
    rt_vmax=rt_vmax/10, #to put into seconds
    rt_vmax_cum=clock_onset + rt_vmax,
    rt_sec = rt_csv/1000) %>%
  group_by(id, run) %>%
  mutate(
    rt_lag=as.vector(scale(lag(rt_sec, order_by=run_trial))),
    iti_lag=lag(iti_ideal, order_by=run_trial),
    v_entropy_wi = v_entropy - mean(v_entropy, na.rm=TRUE),
    v_entropy_bw = mean(v_entropy, na.rm=TRUE)
    ) %>% ungroup() %>%
  select(id, run, trial, run_trial, rewFunc, emotion, rt_sec, rt_lag, rt_next, ev, rt_vmax_lag, 
         v_entropy, v_entropy_change, pe_max, iti_ideal, iti_lag,
         v_entropy_bw, v_entropy_wi)

combined <- left_join(compress_results, trial_df, by=c("id", "run", "trial")) %>%
  group_by(id, run) %>%
  mutate(
    acom_bw=mean(decon_acompress_0.9, na.rm=TRUE),
    acom_lag=as.vector(scale(lag(decon_acompress_0.9, 1, order_by=run_trial))),
    com_lag=as.vector(scale(lag(decon_compress_0.9, 1, order_by=run_trial)))
  ) %>% ungroup() %>%
  mutate(
    acom_bw=scale(acom_bw),
    rt_vmax_lag=scale(rt_vmax_lag)
  )

#combined_rle <- rle_dt$new(data=combined, keys=c("id", "run", "trial", "atlas_value", "rewFunc", "emotion"), optimize_order = FALSE)

avs <- sort(unique(combined$atlas_value))

#force cleanup
#rm(compress_results, combined)
#gc()



# toplot <- combined_rle$get() %>% filter(atlas_value %in% c(41, 145, 33, 141, 31, 135)) %>%
#   mutate(atlas_value=factor(atlas_value, levels=c(33, 41, 141, 145, 31, 135), labels=c("L_LIPd", "L_FEF", "R_LIPd", "R_FEF", "L_MT", "R_MT")))
# 
# 
# rfef <- combined_rle$get() %>% filter(atlas_value == 145)

toplot <- combined %>% filter(atlas_value %in% c(41, 145, 33, 141, 31, 135)) %>%
  mutate(atlas_value=factor(atlas_value, levels=c(33, 41, 141, 145, 31, 135), labels=c("L_LIPd", "L_FEF", "R_LIPd", "R_FEF", "L_MT", "R_MT")))


rfef <- combined %>% filter(atlas_value == 145)



m_rfef <- lmer(decon_acompress_0.9 ~ v_entropy + (1 | id/run), data=rfef)

m_rfef2 <- lmer(decon_acompress_0.9 ~ v_entropy_change*rewFunc + rt_sec + 
                  run_trial + decon_compress_nt + (1 | id/run), data=rfef)

#include prior and current ITI as covariates to control for how much of each period
#is factored into the number of timepoints... seems like this may be too collinear!
m_rfef3 <- lmer(decon_acompress_0.9 ~ v_entropy_change*rewFunc + rt_sec + 
                 run_trial + decon_compress_nt + iti_ideal + iti_lag + (1 | id/run), data=rfef)

car::vif(m_rfef2)

summary(m_rfef2)
summary(m_rfef3)
car::Anova(m_rfef3, type=3)
emtrends(m_rfef3, ~rewFunc, var="v_entropy_change")

#look at bw/wi
#rewFunc is confounded with v_entropy_bw
m_rfef4 <- lmer(decon_acompress_0.9 ~ v_entropy_wi*rewFunc + rt_sec + v_entropy_bw + #*rewFunc
                  v_entropy_change + run_trial + rt_sec*iti_ideal + rt_sec*iti_lag + (1 | id/run), data=rfef)

#decon_compress_nt

#entropy bw effect -- people who have higher average entropy in a run tend to have more
#compressed representations in FEF...
summary(m_rfef4)
car::Anova(m_rfef4, type=3)
car::vif(m_rfef4, type=3)

summary(lm(decon_compress_nt ~ iti_ideal + iti_lag + rt_sec, rfef))

#1) ITI-only windowing -- maybe 3+ seconds? (~36%)
#2) decompose NT into online/offline, or pre, online, post
#3) window something like +3 from last feedback (or clock, whichever is first) through + 10 (next clock as collider)
#4) windows of equal length aligned to different events (mid-clock, rt+1, mid-iti)
#5) compare against rest?
#6) carve out and concatenate ITI periods, windowed SVD, ignoring trial

#continuity of compression? yes, big-time
m_rfef5 <- lmer(decon_acompress_0.9 ~ acom_lag*v_entropy_change + v_entropy_wi*rewFunc + rt_sec + v_entropy_bw + #*rewFunc
                  v_entropy_change + run_trial + decon_compress_nt + iti_ideal + iti_lag + (1 | id/run), data=rfef)

summary(m_rfef5)
car::Anova(m_rfef5, type=3)
car::vif(m_rfef4, type=3)


m_rfef6 <- lmer(decon_acompress_0.9 ~ acom_lag + v_entropy_wi*rewFunc + rt_sec + v_entropy_bw + #*rewFunc
                  v_entropy_change + run_trial + decon_compress_nt + iti_ideal + iti_lag + (1 | id/run), data=rfef)

car::Anova(m_rfef6, type=3)
car::vif(m_rfef4, type=3)

#rt modeling
m_rt1 <- lmer(rt_sec ~ acom_lag + v_entropy_wi + v_entropy_bw + rewFunc + rt_vmax_lag +
                  v_entropy_change + run_trial + (1 | id/run), data=rfef) # iti_ideal + iti_lag +

summary(m_rt1)
car::Anova(m_rt1, type=3)
car::vif(m_rt1, type=3)

#what is relationship of entropy_wi with acom_lag?
# m_e1 <- lmer(v_entropy_wi ~ acom_lag + run_trial + (1 | id/run), data=rfef) # iti_ideal + iti_lag +
# summary(m_e1) #huh, not significant... and in RT model where both are main effects, we are okay

#bring in AR(1) RT
#huh, so more compressed representations in R FEF go with longer RTs
m_rt2 <- lmer(rt_sec ~ rt_lag * acom_lag + rewFunc + rt_vmax_lag * acom_bw + rt_vmax_lag * acom_lag +  #v_entropy_wi + v_entropy_bw + v_entropy_change +
                run_trial + decon_compress_nt + (1 | id/run), data=rfef) # iti_ideal + iti_lag +
#run*acom_bw +

m_rt2 <- lmer(rt_sec ~ iti_ideal + iti_lag + rewFunc + run_trial + (1 | id/run), data=rfef) # 


summary(m_rt2)
car::Anova(m_rt2, type=3)
car::vif(m_rt2, type=3)

m_rt2 <- lmer(rt_next ~ rt_lag * acom_lag + rewFunc + rt_vmax_lag * acom_bw + rt_vmax_lag * acom_lag +  #v_entropy_wi + v_entropy_bw + v_entropy_change +
                run_trial + decon_compress_nt + (1 | id/run), data=rfef) # iti_ideal + iti_lag +



pdf("compress_preview_approx.pdf", width=15, height=15)
ggplot(toplot, aes(x=run_trial, y=decon_acompress_0.9, color=rewFunc)) + geom_smooth(method="loess") +
  facet_wrap(~atlas_value)
dev.off()


mlist <- list()
coef_list <- list()
for (aa in avs) {
  df <- combined %>% filter(atlas_value==!!aa & !is.na(v_entropy))
  mm <- lmer(decon_acompress_0.9 ~ v_entropy + (1 | id/run), data=df)
  mlist[[as.character(aa)]] <- mm
  coef_list[[as.character(aa)]] <- tidy(mm) 
}

