# adds full-maintenance ("fixed") entropy stats to the behavioral trial df
# calculates entropy change, pos/neg entropy change

library(tidyverse)
plots = F
# read in existing behavioral trial-level data
behavioral_data_file <- "~/code/clock_analysis/meg/MEG_n63_behavioral_data_preprocessed_trial_df.RDS"
trial_df <- readRDS(behavioral_data_file) 
# 
# # read in full entropy stats
# full_df <- read.csv("~/code/clock_analysis/meg/data/mmclock_meg_fixed_fixedparams_meg_ffx_trial_statistics.csv.gz")
# 
# # calculate within- vs. between-run entropy, lags, change, positive, negative
# full_df <- full_df %>% 
#   rename(v_entropy_full = v_entropy, pe_max_full = pe_max, rt_vmax_full = rt_vmax) %>% 
#   group_by(id,run) %>% mutate(v_entropy_wi_full = as.vector(scale(v_entropy_full)),
#                               v_entropy_b_full = mean(na.omit(v_entropy_full)),
#                               v_entropy_wi_change_full = lead(v_entropy_wi_full) - v_entropy_wi_full) %>%
#   ungroup() %>% mutate(entropy_change_pos_wi_full = case_when(
#     v_entropy_wi_change_full > 0 ~ v_entropy_wi_change_full,
#     TRUE ~ 0),
#     entropy_change_neg_wi_full = case_when(
#       v_entropy_wi_change_full < 0 ~ -v_entropy_wi_change_full,
#       TRUE ~ 0
#     ),
#     # recode for compatibility
#     id  = as.integer(substr(id,1,5)),
#     run = as.numeric(run),
#     trial = as.numeric(trial),
#     score_csv = as.numeric(score_csv),
#     rt_csv = as.numeric(rt_csv)/1000) %>% select(id, run, trial, asc_trial, rewFunc, emotion, rt_csv, score_csv, pe_max_full,
#                   v_entropy_full, v_entropy_wi_full, v_entropy_b_full, v_entropy_wi_change_full,
#                   entropy_change_pos_wi_full, entropy_change_neg_wi_full)
# # str(full_df)
# 
# # now we can merge; remove condition coded as a character string instead of factor
# trial_df <- merge(trial_df, full_df, by = c("id", "run", "trial", "rt_csv", "score_csv")) %>% 
#   rename(rewFunc = rewFunc.x) %>% select(!rewFunc.y)
str(trial_df)

# plot for sanity checks
if (plots) {
library(ggpubr)
p1 <- ggplot(trial_df %>% filter(run>1), aes(run_trial, v_entropy_wi)) + geom_smooth()
p2 <- ggplot(trial_df %>% filter(run>1), aes(run_trial, v_entropy_wi_full)) + geom_smooth()
p3 <- ggplot(trial_df %>% filter(run>1), aes(run_trial, v_entropy_wi_change)) + geom_smooth()
p4 <- ggplot(trial_df %>% filter(run>1), aes(run_trial, v_entropy_wi_change_full)) + geom_smooth()
setwd("~/OneDrive/collected_letters/papers/meg/plots/wholebrain")
pdf("full_vs_selective_entropy_sanity_chesk.pdf", height = 14, width = 20)
ggarrange(p1,p2,p3,p4)
dev.off()
}
psych::corr.test(trial_df$v_entropy_wi, trial_df$v_entropy_wi_full)
psych::corr.test(trial_df %>% select(v_entropy_wi, v_entropy_wi_full, v_entropy_wi_change, v_entropy_wi_change_full))

# OK to save updated df
saveRDS(object = trial_df, file = behavioral_data_file)
