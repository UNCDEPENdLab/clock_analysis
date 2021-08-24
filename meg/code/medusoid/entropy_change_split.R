# splits entropy change into positive and negative deflections
# this is necessary to understand whether mu-suppression takes place to increases vs. mu-synchronization for 
# decreases of entropy


behavioral_data_file <- "~/code/clock_analysis/meg/MEG_n63_behavioral_data_preprocessed_trial_df.RDS"
trial_df <- readRDS(behavioral_data_file) %>% as.data.frame(lapply(trial_df, function(x) {
  if (inherits(x, "matrix")) { x <- as.vector(x) }
  return(x)
})) %>% mutate(entropy_change_pos_wi = case_when(
  v_entropy_wi_change > 0 ~ v_entropy_wi_change,
  TRUE ~ 0
),
entropy_change_neg_wi = case_when(
  v_entropy_wi_change < 0 ~ -v_entropy_wi_change,
  TRUE ~ 0
))

# examine correlations
psych::corr.test(trial_df$entropy_change_neg_wi, trial_df$entropy_change_pos_wi)
# correlation is only -0.23
psych::corr.test(trial_df$v_entropy_wi, trial_df$entropy_change_pos_wi)
# -0.23
psych::corr.test(trial_df$v_entropy_wi, trial_df$entropy_change_neg_wi)
# 0.27
psych::corr.test(trial_df$abs_pe, trial_df$entropy_change_neg_wi)
# 0.34
psych::corr.test(trial_df$abs_pe, trial_df$entropy_change_pos_wi)
# 0.25
psych::corr.test(trial_df$trial_neg_inv_sc, trial_df$entropy_change_neg_wi)
# -0.03
psych::corr.test(trial_df$trial_neg_inv_sc, trial_df$entropy_change_pos_wi)
# -0.08
psych::corr.test(as.numeric(trial_df$outcome), trial_df$entropy_change_neg_wi, method = "spearman")
# 0.39
psych::corr.test(as.numeric(trial_df$outcome), trial_df$entropy_change_pos_wi, method = "spearman")
# -0.05

# No major multi-colinearity problems 

# check plots: 
ggplot(trial_df, aes(run_trial, entropy_change_pos_wi, color = rewFunc)) + geom_smooth()
ggplot(trial_df, aes(run_trial, entropy_change_neg_wi, color = rewFunc)) + geom_smooth()