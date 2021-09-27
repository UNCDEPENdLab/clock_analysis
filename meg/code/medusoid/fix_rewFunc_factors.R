# correct rewFunc factor coding error in the MEG dataset
library(tidyverse)
library(emmeans)
library(lme4)
behavioral_data_file <- "~/code/clock_analysis/meg/MEG_n63_behavioral_data_preprocessed_trial_df.RDS"
df <- readRDS(behavioral_data_file)

# check behavior in MEG session
df$rewFunc_char <- as.character(df$rewFunc)

ggplot(df, aes(rt_csv, ev)) + geom_smooth() + facet_wrap(~rewFunc)

ggplot(df, aes(run_trial, rt_csv, color = rewFunc_char)) + geom_smooth()
# CRAZY, DEV = CEVR, IEV = CEV, CEV = DEV, CEVR = IEV

