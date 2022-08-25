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

# install_github("UNCDEPENdLab/dependlab")
# library(dependlab)
source('~/code/Rhelpers/screen.lmerTest.R')
source('~/code/Rhelpers/vif.lme.R')
# library(stringi)

# clock_folder <- "~/Data_Analysis/clock_analysis" #michael
clock_folder <- "~/code/clock_analysis" #alex
# source('~/code/Rhelpers/')
setwd(file.path(clock_folder, 'fmri/keuka_brain_behavior_analyses/dan'))
tab_dir <- "~/OneDrive/collected_letters/papers/meg/tables/table_s1/"

### load data
source("~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/get_trial_data.R")
# get_trial_data <- function(repo_directory=NULL, dataset="mmclock_fmri", groupfixed=TRUE) 
mdf <- get_trial_data(repo_directory = clock_folder, dataset = "mmclock_meg", groupfixed = T)
fdf <- get_trial_data(repo_directory = clock_folder, dataset = "mmclock_fmri", groupfixed = T)


wm_fmri1 <- lme4::lmer(rt_csv_sc ~
             rt_lag * omission_lag  +
               rt_lag2 * omission_lag2  +
               rt_lag3 * omission_lag3  +
               rt_lag4 * omission_lag4  +
               rt_lag5 * omission_lag5  +
             (1|id/run), fdf %>% filter(rt_csv<4000 & !is.na(rt_vmax_lag_sc) & !is.na(trial_neg_inv_sc)))
summary(wm_fmri1)
Anova(wm_fmri1)
vif(wm_fmri1)

wm_rtvmax_fmri1 <- lme4::lmer(rt_csv_sc ~
                         rt_lag * omission_lag  +
                         rt_lag2 * omission_lag2  +
                         rt_lag3 * omission_lag3  +
                         rt_lag4 * omission_lag4  +
                         rt_lag5 * omission_lag5  +
                         rt_vmax_lag_sc * trial_neg_inv_sc  +
                           (1|id/run), fdf %>% filter(rt_csv<4000 & !is.na(rt_vmax_lag_sc) & !is.na(trial_neg_inv_sc)))

summary(wm_rtvmax_fmri1)

wm_rtvmax_fmri1_rs <- lme4::lmer(rt_csv_sc ~
                                rt_lag * omission_lag  +
                                rt_lag2 * omission_lag2  +
                                rt_lag3 * omission_lag3  +
                                rt_lag4 * omission_lag4  +
                                rt_lag5 * omission_lag5  +
                                rt_vmax_lag_sc * trial_neg_inv_sc  +
                                (rt_lag * omission_lag  +
                                   rt_lag2 * omission_lag2  +
                                   rt_lag3 * omission_lag3  +
                                   rt_lag4 * omission_lag4  +
                                   rt_lag5 * omission_lag5  +
                                   rt_vmax_lag_sc * trial_neg_inv_sc |id/run), fdf %>% filter(rt_csv<4000 & !is.na(rt_vmax_lag_sc) & !is.na(trial_neg_inv_sc)))

summary(wm_rtvmax_fmri1_rs)


anova(wm_fmri1, wm_rtvmax_fmri1)

wm_meg1 <- lme4::lmer(rt_csv_sc ~
                         rt_lag * omission_lag  +
                         rt_lag2 * omission_lag2  +
                         rt_lag3 * omission_lag3  +
                         rt_lag4 * omission_lag4  +
                         rt_lag5 * omission_lag5  +
                         (1|id/run), mdf %>% filter(rt_csv<4000 & !is.na(rt_vmax_lag_sc) & !is.na(trial_neg_inv_sc)))
summary(wm_meg1)
Anova(wm_meg1)
vif(wm_meg1)

wm_rtvmax_meg1 <- lme4::lmer(rt_csv_sc ~
                                rt_lag * omission_lag  +
                                rt_lag2 * omission_lag2  +
                                rt_lag3 * omission_lag3  +
                                rt_lag4 * omission_lag4  +
                                rt_lag5 * omission_lag5  +
                                rt_vmax_lag_sc * trial_neg_inv_sc  +
                                (1|id/run), mdf %>% filter(rt_csv<4000 & !is.na(rt_vmax_lag_sc) & !is.na(trial_neg_inv_sc)))

summary(wm_rtvmax_meg1)
Anova(wm_rtvmax_meg1)
anova(wm_meg1, wm_rtvmax_meg1)
tab_model(wm_rtvmax_fmri1, wm_rtvmax_meg1, order.terms = c(1,13,2,4,6,8,10,3,5,7,9,11,14:18,12,19), 
          show.stat = T, show.ci = F, show.est = F,  show.re.var = F, dv.labels = c("fMRI session", "MEG session"), show.icc = F,  show.ngroups = F,
          pred.labels = 
            c("Intercept", 
              "RT_lag1", "Omission_lag1", 
              "RT_lag2", "Omission_lag2", 
              "RT_lag3", "Omission_lag3", 
              "RT_lag4", "Omission_lag4", 
              "RT_lag5", "Omission_lag5",
              "RT_Vmax", "Trial",
              "RT_lag1 * Omission_lag1", 
              "RT_lag2 * Omission_lag2", 
              "RT_lag3 * Omission_lag3", 
              "RT_lag4 * Omission_lag4", 
              "RT_lag5 * Omission_lag5", 
              "RT_Vmax * Trial"
              ), file = "dan_table_s1_wm.html", p.style = "numeric", p.threshold = c(.05, .01, .001))
library(rio)
setwd(tab_dir)
saveXML(tab, file = "dan_table_s1_wm.html", format = "html")
export2html(tab, file = "dan_table_s1_wm.html")
# save_html(paste(as.character(tab), collapse = "\n"), file = "dan_table_s1_wm.html", background = "white", lang = "en")


# same with Coxme


