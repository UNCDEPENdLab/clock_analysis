library(dplyr)
library(tidyverse)
library(psych)
library(ggcorrplot)
library(lme4)
library(ggpubr)
library(cowplot)
# library(sjPlot)
# library(sjmisc)
#library(mlVAR)

library(ggeffects)
library(viridis)
library(car)
library(data.table)
library(emmeans)

# read in, process; go with "long" [-1:10] clock windows for now, will censor later
#####################
#plots = F
reprocess = FALSE #load from cache
#analyze = F

#repo_dir <- "~/code/clock_analysis"
repo_dir <- "~/Data_Analysis/clock_analysis"

source(file.path(repo_dir, "fmri/keuka_brain_behavior_analyses/load_medusa_data.R"))

fb_comb$bin_num_f <- as.factor(fb_comb$bin_num)

load(file.path(repo_dir, "/fmri/keuka_brain_behavior_analyses/trial_df_and_vhdkfpe_clusters.Rdata"))

attr(df, "labels") <- NULL #somehow this is holding a 560-row data.frame
df <- df %>% dplyr::ungroup()

# read in behavioral data
# select relevant columns for compactness
df <- df %>% select(id, run, run_trial, rewFunc,emotion, rt_csv, score_csv, rt_next, rt_vmax, rt_vmax_lag,rt_vmax_change, v_max_wi, v_entropy_wi, v_entropy_b, v_entropy, v_max_b, Age, Female)
# add deconvolved hippocampal timeseries
d <- merge(df, fb_comb, by = c("id", "run", "run_trial")) %>% mutate(v_entropy_wi_z = as.vector(scale(v_entropy_wi)))

library(permuco)
tt <- aovperm(decon_interp ~ bin_num_f*evt_time_f*v_entropy_wi_z + rewFunc.x + rt_csv.y + Error(id/(bin_num_f*evt_time_f)),
    data = d, method = "Rd_kheradPajouh_renaud")

# replicate lm decoding analyses
# scale(-1/run_trial)*rewFunc + reward + scale(rt_csv) + scale(rt_vmax_lag) + scale(rt_vmax_change) + v_entropy_wi


#analyze_df <- fb_comb %>% filter(iti_prev > 1 & iti_ideal > 8 & evt_time < 9)
analyze_df <- fb_comb #fit all data

dm1 <- lmer(decon_interp ~ bin_num_f*evt_time_f*scale(-1/run_trial)*rewFunc + 
        bin_num_f*evt_time_f*scale(rt_csv) + 
        bin_num_f*evt_time_f*scale(rt_vmax_lag) +
        bin_num_f*evt_time_f*scale(rt_vmax_change) + 
        bin_num_f*evt_time_f*scale(v_entropy_wi) +
        bin_num_f*evt_time_f*scale(v_entropy_wi_change) +
        (1 | id/run) + (1 | side), analyze_df)

summary(dm1)
car::Anova(dm1, '3')
vif(dm1)
library(emmeans)
r1 <- emtrends(dm1, var = 'rt_vmax_change', specs = c('bin_num_f','evt_time_f'), data = analyze_df)
r1 <- as.data.frame(r1)
ggplot(r1, aes(evt_time_f, bin_num_f, color = rt_vmax_change.trend)) + geom_tile()

# reduce this monstrosity to just one effect of interest
dm2 <- lmer(decon_interp ~ 
        bin_num_f*evt_time_f*scale(rt_vmax_lag) +
        (1 | id/run) + (1 | side), fb_comb %>% filter (iti_prev>1 & iti_ideal>8 & evt_time < 9))
summary(dm2)
car::Anova(dm2, '3')
r2 <- emtrends(dm2, var = 'rt_vmax_lag', specs = c('bin_num_f','evt_time_f'), data = fb_comb %>% filter (iti_prev>1 & iti_ideal>8 & evt_time < 9))
r2 <- as.data.frame(r2)
ggplot(r2, aes(evt_time_f, bin_num_f, color = rt_vmax_lag.trend)) + geom_tile()

dm3 <- lmer(decon_interp ~ 
        bin_num_f*evt_time_f*scale(v_entropy_wi_change) +
        (1 | id/run) + (1 | side), fb_comb %>% filter (iti_prev>1 & iti_ideal>8 & evt_time < 9))
summary(dm3)
car::Anova(dm3, '3')
r3 <- emtrends(dm3, var = 'rt_vmax_lag', specs = c('bin_num_f','evt_time_f'), data = fb_comb %>% filter (iti_prev>1 & iti_ideal>8 & evt_time < 9))
r3 <- as.data.frame(r3)
ggplot(r2, aes(evt_time_f, bin_num_f, color = rt_vmax_lag.trend)) + geom_tile()

# not even reward??
dm4 <- lmer(decon_interp ~ 
        bin_num_f*evt_time_f*reward + side +
        (1 | id/run), analyze_df)
summary(dm4)
car::Anova(dm4, '3')
em <- emmeans(dm4, ~bin_num_f*evt_time_f | reward)
print(em)

vv <- vcov(dm4)

#no evidence of correlation in the fixed effects
summary(vv[lower.tri(vv)])

