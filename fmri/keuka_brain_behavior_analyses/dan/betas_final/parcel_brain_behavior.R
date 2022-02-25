# final parcel-wise analyses of DAN brain-to-behavior
library(data.table)
library(tidyverse)
library(afex)
library(lattice)
library(emmeans)
library(fmri.pipeline) # has mixed_By
library(dplyr)
library(readr)
library(tidyr)
library(purrr)

#analysis_dir <- "~/Data_Analysis/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/betas_final"
analysis_dir <- "/proj/mnhallqlab/projects/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/betas_final"
setwd(analysis_dir)

source("../get_trial_data.R")
source("parcel_brain_behavior_functions.R")

# location of whole brain betas for analysis
# beta_dir <- "/Volumes/GoogleDrive/My Drive/SCEPTIC_fMRI/wholebrain_betas"
beta_dir <- file.path(analysis_dir, "wholebrain_betas")

# MNH/AD internal labeling scheme
labels <- readxl::read_excel(file.path(analysis_dir, "..", "MNH Dan Labels.xlsx")) %>%
  dplyr::rename(mask_value=roinum) %>% select(mask_value, plot_label)

# BALSA parcel labels from whereami
labels_200 <- read.csv(file.path(beta_dir, "schaefer_200_whereami.csv")) %>%
  dplyr::rename(mask_value = roi_num, plot_label = MNI_Glasser_HCP_v1.0) %>% dplyr::select(mask_value, plot_label) %>%
  mutate(plot_label = sub("Focus point:\\s+", "", plot_label, perl=TRUE))

#trial_df <- get_trial_data(repo_directory = "~/Data_Analysis/clock_analysis") %>%
trial_df <- get_trial_data(repo_directory = "/proj/mnhallqlab/projects/clock_analysis") %>%
  group_by(id, run) %>%
  mutate(v_max_wi_lag = lag(v_max_wi, order_by=run_trial)) %>%
  ungroup() %>%
  dplyr::select(id, run, run_trial, trial_neg_inv, rt_csv, rt_lag, v_entropy_wi, v_max_wi_lag, 
                rt_vmax_lag, last_outcome)

echange_l2_copes <- file.path(beta_dir, "L1m-echange/Schaefer2018_200Parcels_7Networks_order_fonov_2.3mm_ants_cope_l2.csv.gz")
int <- formula(~ (trial_neg_inv + rt_lag + v_max_wi_lag + v_entropy_wi + fmri_beta + last_outcome)^2 +
                        rt_lag:last_outcome:fmri_beta +
                        rt_vmax_lag*trial_neg_inv*fmri_beta +
                        (1 | id/run)
)


slo <- formula(~ (trial_neg_inv + rt_lag + v_max_wi_lag + v_entropy_wi + fmri_beta + last_outcome)^2 +
                 rt_lag:last_outcome:fmri_beta +
                 rt_vmax_lag*trial_neg_inv*fmri_beta +
                 (1 + rt_lag + rt_vmax_lag | id/run)
)

efiles <- list.files(beta_dir, pattern = "Schaefer2018_200Parcels_7Networks_order_fonov_2.3mm_ants_cope_l2.csv.gz", 
                    recursive = TRUE, full.names = TRUE)

save.image(file="parcel_input_snapshot.RData")

#efiles <- efiles[3:6]
for (ee in efiles) {
  job <- R_batch_job$new(
    job_name = "parcel_bb", batch_directory = getwd(), scheduler = "slurm",
    input_rdata_file = "parcel_input_snapshot.RData",
    n_nodes = 1, n_cpus = 16, wall_time = "12:00:00",
    mem_total = "64G",
    r_code = glue("to_plot <- mixed_by_betas('{ee}', labels_200, trial_df, mask_file = 'Schaefer2018_200Parcels_7Networks_order_fonov_1mm_ants.nii.gz',
                            rhs_form = fmri.pipeline:::named_list(int, slo), ncores = 16, afni_dir = '/proj/mnhallqlab/sw/afni',
                            split_on = c('l1_cope_name', 'l2_cope_name', 'mask_value'))"),
    r_packages = c("fmri.pipeline", "tidyverse", "data.table", "sfsmisc"),
    batch_code = c("module use /proj/mnhallqlab/sw/modules", "module load r/4.0.3_depend")
  )
  
  job$submit()


  # local execution
  # to_plot <- mixed_by_betas(ee, labels_200, trial_df, mask_file = "Schaefer2018_200Parcels_7Networks_order_fonov_1mm_ants.nii.gz",
  #                           rhs_form = fmri.pipeline:::named_list(int, slo), ncores = 16,
  #                           split_on = c("l1_cope_name", "l2_cope_name", "mask_value"))
}


# these are betas I manually extracted by using colMeans, readNifti, and mask indices in the old/traditional voxelwise approach.
# see entropy_betas.R. These are specifically the entropy change overall betas from the echange model

manual_betas <- read.csv("/Volumes/GoogleDrive/My Drive/SCEPTIC_fMRI/wholebrain_betas/beta_debugging/schaefer200_echange_overall_l2_means.csv") %>%
  mutate(id = 1:n()) %>%
  pivot_longer(names_to = "roi", cols = -id, values_to = "fmri_beta") %>%
  mutate(mask_value = as.integer(sub("^V", "", roi))) %>%
  dplyr::select(-roi)

onesamp_betas(manual_betas, mask_file="Schaefer2018_200Parcels_7Networks_order_fonov_1mm_ants.nii.gz",
              out_dir = "/Volumes/GoogleDrive/My Drive/SCEPTIC_fMRI/wholebrain_betas/beta_debugging")

man_means <- manual_betas %>% group_by(mask_value) %>% summarise(mean = mean(fmri_beta))

# these are betas for the overall echange contrast extracted at the time of the extract_glm_betas_in_mask
echange_intermediate <- readRDS("/Volumes/GoogleDrive/My Drive/SCEPTIC_fMRI/wholebrain_betas/beta_debugging/overall_echange_stat_result.rds") %>%
  unnest(img_stats) %>% select(-img, -img_exists, -x, -y, -z)

int_means <- echange_intermediate %>% group_by(mask_value) %>% summarise(mean = mean(value))

onesamp_betas(echange_intermediate, mask_file="Schaefer2018_200Parcels_7Networks_order_fonov_1mm_ants.nii.gz", dv = "value",
              out_dir = "/Volumes/GoogleDrive/My Drive/SCEPTIC_fMRI/wholebrain_betas/beta_debugging", img_prefix = "echange_intermediate")



#setwd("/Users/michael/Data_Analysis/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/betas_final")

######



###


#beta_dir <- "/Volumes/GoogleDrive/My Drive/SCEPTIC_fMRI/final_betas"



echange_l2_copes <- file.path(beta_dir, "L1m-echange/Schaefer2018_200Parcels_7Networks_order_fonov_2.3mm_ants_cope_l2.csv.gz")

mixed_by_betas(echange_l2_copes, labels_200, trial_df)

# subject-level
#echange_l2_copes <- fread(file.path(beta_dir, "L1m-echange/Schaefer_DorsAttn_2.3mm_cope_l2.csv.gz")) %>%
#echange_l2_copes <- fread(file.path(beta_dir, "L1m-echange/zstat6_ptfce_Schaefer2018_200Parcels_7Networks_order_fonov_2.3mm_ants_overlap_cope_l2.csv.gz")) %>%
#echange_l2_copes <- fread(file.path(beta_dir, "L1m-echange/Schaefer2018_200Parcels_7Networks_order_fonov_2.3mm_ants_cope_l2.csv.gz")) %>%
entropy_l2_copes <- fread(file.path(beta_dir, "L1m-entropy_wiz/Schaefer2018_200Parcels_7Networks_order_fonov_2.3mm_ants_cope_l2.csv.gz")) %>%
  #filter(l1_cope_name == "EV_entropy_change_feedback" & l2_cope_name == "overall") %>%
  filter(l1_cope_name == "EV_entropy_wiz_clock" & l2_cope_name == "overall") %>%
  select(-feat_dir, -img) %>%
  rename(fmri_beta=value) %>%
  left_join(labels_200, by="mask_value")

# meg effects of interest
# meg_wide <- readRDS("MEG_betas_wide_echange_vmax_reward_Nov15_2021.RDS")

# keep as long to have meg beta as split
meg_ranefs <- readRDS("MEG_betas_echange_vmax_reward_Nov15_2021.RDS") %>%
  mutate(id=as.numeric(id)) %>% rename(meg_beta = avg)

# run-level
echange_run_copes <- fread(file.path(beta_dir, "L1m-echange/Schaefer_DorsAttn_2.3mm_cope_l1.csv.gz")) %>%
  dplyr::select(-feat_dir, -img) %>%
  dplyr::filter(l1_cope_name == "EV_entropy_change_feedback") %>%
  dplyr::rename(run=run_number)

# very mild winsorization -- .25 %
echange_run_copes <- echange_run_copes %>%
  mutate(fmri_beta = DescTools::Winsorize(value, probs=c(.0025, .9975))) %>%
  left_join(labels, by="mask_value")

echange_sub_copes <- echange_run_copes %>%
  group_by(id, mask_value, plot_label) %>%
  summarize(fmri_beta=mean(fmri_beta)) %>%
  ungroup()

ee <- echange_sub_copes %>% rename(avg_l1=fmri_beta) %>% select(id, mask_value, avg_l1) %>%
  left_join(echange_l2_copes %>% select(id, mask_value, fmri_beta))

# correlation of l2 overall contrast with avg l1 beta
tt <- ee %>% group_by(id) %>%
  do({
    df <- .
    data.frame(r=cor(df$fmri_beta, df$avg_l1))
  }) %>% ungroup()


# sub_sep <- merge(echange_l2_copes, meg_ranefs, by="id") %>%
#   left_join(labels, by="mask_value")

sub_sep <- left_join(echange_l2_copes, meg_ranefs, by="id")


tt <- sub_sep %>% group_by(plot_label, reg_region) %>%
  do({
    df <- .
    data.frame(r=cor(df$fmri_beta, df$meg_beta))
  }) %>% ungroup()

ggplot(tt, aes(x=reg_region, y=plot_label, fill=r)) + geom_tile() + scale_fill_viridis_c()


histogram(echange_run_copes$value)
histogram(echange_run_copes$fmri_beta)
histogram(echange_l2_copes$fmri_beta)

histogram(~fmri_beta | id, echange_run_copes)


#combo <- echange_run_copes %>% left_join(trial_set, by=c("id", "run"))
#combo <- echange_sub_copes %>% left_join(trial_set, by=c("id"))
#combo <- echange_l2_copes %>% left_join(trial_set, by=c("id"))
combo <- entropy_l2_copes %>% left_join(trial_set, by=c("id"))


model_base <- formula(~ (trial_neg_inv + rt_lag + v_max_wi_lag + v_entropy_wi + fmri_beta + last_outcome)^2 +
                        rt_lag:last_outcome:fmri_beta +
                        rt_vmax_lag*trial_neg_inv*fmri_beta +
                        (1 | id/run)
                      )

#test_df <- combo %>% dplyr::filter(plot_label=="R_pVIP") %>%
test_df <- combo %>% dplyr::filter(plot_label==" R_Orbital_Frontal_Complex") %>%
  mutate_if(is.numeric, scale)

model_base_nobrain <- formula(~ (trial_neg_inv + rt_lag + v_max_wi_lag + v_entropy_wi + last_outcome)^2 + (1 | id/run))
ff <- update.formula(model_base_nobrain, "rt_csv ~ .")
test_model <- lmer(ff, data=test_df)

summary(test_model)


ff <- update.formula(model_base, "rt_csv ~ .")
test_model <- lmer(ff, data=test_df)
summary(test_model)


splits <- "mask_value"
ddf <- mixed_by(combo, outcomes = "rt_csv", rhs_model_formulae = list(main=model_base),
                split_on = splits, scale_predictors = c("trial_neg_inv", "rt_lag", "rt_vmax_lag", "v_max_wi_lag", "fmri_beta"),
                tidy_args = list(effects = c("fixed", "ran_vals", "ran_pars", "ran_coefs"), conf.int = TRUE), 
                calculate = c("parameter_estimates_reml"), ncores = 16, refit_on_nonconvergence = 5, padjust_by = "term",
                emtrends_spec = list(
                  rt_lag = list(outcome = "rt_csv", model_name = "main", var = "rt_lag", specs = c("last_outcome", "fmri_beta"), 
                                at=list(fmri_beta = c(-2, 0, 2))) # z scores
                ))


to_plot <- ddf$coef_df_reml %>%
  filter(effect=="fixed") %>%
  group_by(term) %>%
  mutate(p_FDR=p.adjust(p.value, method="fdr")) %>%
  ungroup() %>% 
  #left_join(labels, by="mask_value") %>%
  left_join(labels_200, by="mask_value") %>%
  setDT()
  

#write.csv(to_plot %>% select(-rhs, -effect), file="fmri_brainbehavior_parcel_betas.csv", row.names=FALSE)
#write.csv(to_plot %>% select(-rhs, -effect), file="fmri_brainbehavior_parcel_betas_ptfce.csv", row.names=FALSE)
#write.csv(to_plot %>% select(-rhs, -effect), file="fmri_brainbehavior_parcel_betas_200.csv", row.names=FALSE)
write.csv(to_plot %>% select(-rhs, -effect), file="fmri_brainbehavior_parcel_entropy_betas_200.csv", row.names=FALSE)


# plot_medusa(ddf, x="plot_label", y=1, color="estimate",
#             out_dir="parcel_betas", p.value="p_FDR", plot_type="heat")

plot_medusa(to_plot, x="plot_label", y="estimate", ymin="estimate - std.error", ymax="estimate + std.error",
            out_dir="parcel_betas_wb", p.value="p_FDR", plot_type="line", flip = TRUE, term_filter="fmri_beta")

###

emt <- ddf$emtrends_list$rt_lag
emt <- emt[,-2:-3]

emt <- emt %>% left_join(labels, by="mask_value") %>%
  mutate(fmri_beta=factor(fmri_beta, levels=c(-2, 0, 2), labels=c("low", "mean", "high")))
ggplot(emt, aes(x=plot_label, y=rt_lag.trend, color=fmri_beta, ymin=rt_lag.trend - std.error, ymax = rt_lag.trend + std.error)) +
  geom_pointrange() + coord_flip() + facet_wrap(~last_outcome)

## MEG -> fMRI


with_meg <- combo %>% inner_join(meg_ranefs, by="id")


splits <- c("mask_value", "reg_region")
model_meg <- formula(~ (trial_neg_inv + rt_lag + rt_vmax_lag + v_max_wi_lag + v_entropy_wi + fmri_beta + meg_beta + last_outcome)^2 +
                       rt_lag:last_outcome:fmri_beta +
                       rt_vmax_lag:trial_neg_inv:fmri_beta +
                       rt_lag:last_outcome:meg_beta +
                       rt_vmax_lag:trial_neg_inv:meg_beta +
                       (1 | id/run)
)

model_nofmri <- formula(~ (trial_neg_inv + rt_lag + rt_vmax_lag + v_max_wi_lag + v_entropy_wi + meg_beta + last_outcome)^2 +
                         rt_lag:last_outcome:meg_beta +
                         rt_vmax_lag:trial_neg_inv:meg_beta +
                         (1 | id/run)
)


model_nomeg <- formula(~ (trial_neg_inv + rt_lag + rt_vmax_lag + v_max_wi_lag + v_entropy_wi + fmri_beta + last_outcome)^2 +
                       rt_lag:last_outcome:fmri_beta +
                       rt_vmax_lag:trial_neg_inv:fmri_beta +
                       (1 | id/run)
)


meg_df <- mixed_by(with_meg, outcomes = "rt_csv", rhs_model_formulae = list(main=model_meg, nomeg=model_nomeg),
                split_on = splits, scale_predictors = c("trial_neg_inv", "rt_lag", "rt_vmax_lag", "v_max_wi_lag", "fmri_beta", "meg_beta"),
                tidy_args = list(effects = c("fixed", "ran_vals", "ran_pars", "ran_coefs"), conf.int = TRUE), 
                calculate = c("parameter_estimates_reml"), ncores = 16, refit_on_nonconvergence = 5, padjust_by = "term")


test_df <- with_meg %>% dplyr::filter(reg_region == "entropy_change_late_beta" & plot_label=="R_pVIP") %>%
  mutate_if(is.numeric, scale)
ff <- update.formula(model_meg, "rt_csv ~ .")
#ff <- update.formula(model_nofmri, "rt_csv ~ .")
test_model <- lmer(ff, data=test_df)
summary(test_model)

meg_ee <- emtrends(test_model, specs = ~ last_outcome * meg_beta, var="rt_lag", at=list(meg_beta=c(-2, 0, 2)))
ggplot(as.data.frame(meg_ee), aes(x=last_outcome, y=rt_lag.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=as.factor(meg_beta))) + 
         geom_pointrange(position=position_dodge(width=1))

fmri_ee <- emtrends(test_model, specs = ~ last_outcome * fmri_beta, var="rt_lag", at=list(fmri_beta=c(-2, 0, 2)))
ggplot(as.data.frame(fmri_ee), aes(x=last_outcome, y=rt_lag.trend, ymin=asymp.LCL, ymax=asymp.UCL, color=as.factor(fmri_beta))) + 
  geom_pointrange(position=position_dodge(width=1))


plot_df <- meg_df$coef_df_reml %>%
  filter(effect=="fixed") %>%
  group_by(term) %>%
  mutate(p_FDR=p.adjust(p.value, method="fdr")) %>%
  ungroup() %>% 
  left_join(labels, by="mask_value") %>%
  mutate(statistic = if_else(abs(statistic) < 1, NA_real_, statistic)) %>%
  #filter(abs(statistic) > 1) %>%
  setDT()

plot_medusa(plot_df, x="plot_label", y="reg_region", color="statistic", term_filter="(meg_beta|fmri_beta)",
            out_dir="parcel_betas_meg", p.value="p_FDR", plot_type="heat")




# add MEG late beta .6-.8 suppression, early theta synchronization in separate models.

# 
#   mb_dan1 <-  lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + rt_vmax_lag_sc + last_outcome + 
#                                   v_max_wi_lag + v_entropy_wi + dan_parietal + pe_ips)^2 + 
#                      rt_lag_sc:last_outcome:dan_parietal + 
#                      rt_lag_sc:last_outcome:pe_ips +
#                      rt_vmax_lag_sc:trial_neg_inv_sc:dan_parietal + 
#                      rt_vmax_lag_sc:trial_neg_inv_sc:pe_ips  +
#                      (1|id/run), df %>% filter(rt_csv<4000))





# 
# mb_dan1 <-  lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + rt_vmax_lag_sc + last_outcome + 
#                                 v_max_wi_lag + v_entropy_wi + dan_parietal + pe_ips)^2 + 
#                    rt_lag_sc:last_outcome:dan_parietal + 
#                    rt_lag_sc:last_outcome:pe_ips +
#                    rt_vmax_lag_sc:trial_neg_inv_sc:dan_parietal + 
#                    rt_vmax_lag_sc:trial_neg_inv_sc:pe_ips  +
#                    (1|id/run), df %>% filter(rt_csv<4000))
# screen.lmerTest(mb_dan1, .05)
# summary(mb_dan1)
# Anova(mb_dan1, '3')
# 
# # basic model with only the general entropy factor:
# mb_dan1a <-  lmer(rt_csv_sc ~ (trial_neg_inv_sc + rt_lag_sc + rt_vmax_lag_sc + last_outcome + 
#                                  v_max_wi_lag + v_entropy_wi + general_entropy)^2 + 
#                     rt_lag_sc:last_outcome:general_entropy + 
#                     rt_lag_sc:trial_neg_inv_sc: v_entropy_wi:general_entropy + 
#                     rt_vmax_lag_sc:trial_neg_inv_sc:general_entropy + 
#                     (1|id/run), df %>% filter(rt_csv<4000))
# screen.lmerTest(mb_dan1a, .05)
# summary(mb_dan1a)
# 
