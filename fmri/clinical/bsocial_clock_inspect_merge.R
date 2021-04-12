# inspects bSocial behavioral data
# merges with demographic/clinical subject characteristics

library(dplyr)
library(tidyverse)
library(psych)
library(corrplot)
library(lme4)
library(stargazer)

clock_folder <- "~/code/clock_analysis"
data_folder <- "~/Box/SCEPTIC_fMRI/bsocial/behav/"

# load data
setwd(data_folder)
load("behav_concat.rdata")
sub_df <- read_csv("bsoc_clock_demo.csv") %>% select(!X1) %>% mutate(ID = as.factor(ID)
                                                                     )
df <- as_tibble(df) %>% mutate(ID = as.factor(redcapID),
                               rt = as.numeric(rt),
                               run = ceiling(trial/50)) %>% merge(sub_df) %>%
  # get lags
  group_by(ID, run) %>% mutate(rt_lag = lag(rt),
                               reward = score>0,
                               reward_lag = lag(reward),
                               run_trial = trial - (run-1)*50,
                               rt_swing = abs(rt - rt_lag))

# check design
table(df$run, df$rewFunc, df$emotion)

# diagnostic plots
setwd("~/code/clock_analysis/fmri/clinical/bsocial_behav/")
pdf("bsocial_behav_diags.pdf", height = 30, width = 30)
ggplot(df, aes(trial, rt, color = rewFunc)) + facet_wrap(~ID) + geom_line()
dev.off()

# learning: looks pretty good
ggplot(df, aes(run_trial, rt, color = rewFunc)) + geom_smooth()
ggplot(df, aes(run_trial, rt_swing, color = rewFunc)) + geom_smooth()

# simple behavioral model for validation 
m1 <- lmer(rt ~ rt_lag*reward_lag + rewFunc*run_trial + (1|ID), df)


