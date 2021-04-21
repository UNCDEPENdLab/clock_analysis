library(dplyr)
library(data.table)
library(lme4)
source("~/Data_Analysis/r_packages/fmri.pipeline/R/mixed_by.R")

repo_directory <- "~/Data_Analysis/clock_analysis"
load(file.path(repo_directory, "fmri/keuka_brain_behavior_analyses/dan/trial_df_and_vh_pe_hd_clusters_u.Rdata"))


brain_df <- readRDS(file.path(repo_directory, "fmri/keuka_brain_behavior_analyses/dan/encode_ranefs.rds"))
#now read in the ranefs and generate MEDuSA results for the behavioral model

#lacks a few things like v_max_wi_lag
#trial_df <- readRDS(file.path(repo_directory, "fmri/keuka_brain_behavior_analyses/dan/fmri_trial_df_medusa.RDS"))

df <- mdf #switch to meg as behavior

#for simplicity, change all matrix columns (wi-centering and scaling) to numeric
df <- as.data.frame(lapply(df, function(x) {
  if (inherits(x, "matrix")) { x <- as.vector(x) }
  return(x)
}))

#narrow trial_df to retain only crucial columns
df <- df %>%  dplyr::select(
  id, run, run_trial, rt_csv_sc, trial_neg_inv_sc,
  rt_lag_sc, rt_vmax_lag_sc, last_outcome, v_max_wi_lag, v_entropy_wi)

setDT(df)

#reshape to have ranefs as columns? No, not useful
# test <- brain_df %>% select(stream, side, evt_time, term, level, estimate) %>%
#   pivot_wider(names_from="term", values_from = "estimate") %>%
#   rename(id=level) %>% as_tibble(.name_repair = "universal")
#
# library(corrr)
# xx <- brain_df %>% group_by(stream, evt_time, term) %>%
#   dlookr::correlate("estimate")
#
# ggplot(xx, aes(x=var1, y=var2, fill=coef_corr)) + geom_tile() + facet_grid(stream ~ evt_time)
#

#cleanup
brain_df <- brain_df %>%
  select(stream, side, evt_time, term, level, estimate) %>%
  rename(id=level, brain_blup=estimate, brain_ranef=term) %>%
  mutate(brain_ranef=make.names(brain_ranef), id=as.integer(id)) %>%
  filter(brain_ranef %in% c("scale.abs_pe.", "v_entropy_wi")) #just pe and entropy slopes for now

#combine ranefs with behavior
combined_df <- merge(brain_df, df, by="id", allow.cartesian = TRUE)


behav_formula <- formula( ~ (trial_neg_inv_sc + rt_lag_sc + rt_vmax_lag_sc + last_outcome +
                               v_max_wi_lag + v_entropy_wi + brain_blup)^2 +
                            rt_lag_sc:last_outcome:brain_blup +
                            rt_vmax_lag_sc:trial_neg_inv_sc:brain_blup + (1|id/run))

#combine ranefs with trial_df

splits <- c("stream", "side", "evt_time", "brain_ranef")

# test_df <- combined_df %>% filter(
#   stream=="dorso-dorsal" & side=="L" & evt_time=="-4" & brain_ranef=="scale.abs_pe."
# )
# mtest <- lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + rt_vmax_lag_sc + last_outcome +
#      v_max_wi_lag + v_entropy_wi + brain_blup)^2 +
#   rt_lag_sc:last_outcome:brain_blup +
#   rt_vmax_lag_sc:trial_neg_inv_sc:brain_blup + (1|id/run), data=test_df)

behav_brain_eblupeffect <- mixed_by(
  combined_df, outcomes = "rt_csv_sc", rhs_model_formulae = behav_formula,
  split_on = splits, ncores = 16, refit_on_nonconvergence = 5, padjust_by = "term")

#saveRDS(behav_brain_eblupeffect, file=file.path(repo_directory, "mixed_by_entropy_pe_ranef_effects.rds"))
saveRDS(behav_brain_eblupeffect, file=file.path(repo_directory, "mixed_by_entropy_pe_ranef_effects_MEG.rds"))
