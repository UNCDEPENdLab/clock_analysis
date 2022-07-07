# final parcel-wise analyses of DAN brain-to-behavior
library(data.table)
library(tidyverse)
library(afex)
library(lattice)
session  = "fmri"
if (Sys.getenv("USER")=="alexdombrovski") {
  setwd("~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/betas_final")
  
  source("../get_trial_data.R")
  # source("../medusa_final/plot_medusa.R")
  source("~/code/fmri.pipeline/R/mixed_by.R")
  # labels <- readxl::read_excel("~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/MNH Dan Labels.xlsx") %>%
  #   dplyr::rename(mask_value=roinum) %>% select(mask_value, plot_label)
  
  beta_dir <- "/Volumes/GoogleDrive/My Drive/SCEPTIC_fMRI/final_betas"
  # trial_df <- get_trial_data(repo_directory = "~/code/clock_analysis")
  # trial_df <- get_trial_data(repo_directory = "~/code/clock_analysis", dataset = "mmclock_fmri") %>% mutate(id = as.integer(id))# for MEG
} else {
  setwd("/Users/hallquist/Data_Analysis/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/betas_final")
  
  source("../get_trial_data.R")
  #source("/Users/hallquist/Data_Analysis/r_packages/fmri.pipeline/R/mixed_by.R")
  source("../medusa_final/plot_medusa.R")
  # labels <- readxl::read_excel("/Users/hallquist/Data_Analysis/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/MNH Dan Labels.xlsx") %>%
  #   dplyr::rename(mask_value=roinum) %>% select(mask_value, plot_label)
  # trial_df <- get_trial_data(repo_directory = "/Users/hallquist/Data_Analysis/clock_analysis")
}
beta_dir <- "/Volumes/GoogleDrive/My Drive/SCEPTIC_fMRI/final_betas"

library(emmeans)
library(fmri.pipeline) # has mixed_By
library(dplyr)
library(readr)
library(tidyr)
library(purrr)
library(glue)

#analysis_dir <- "~/Data_Analysis/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/betas_final"
analysis_dir <- "~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/betas_final"
setwd(analysis_dir)

source("../get_trial_data.R")
source("parcel_brain_behavior_functions_alex.R")

# location of whole brain betas for analysis
beta_dir <- "/Volumes/GoogleDrive/My Drive/SCEPTIC_fMRI/wholebrain_betas"
# beta_dir <- file.path(analysis_dir, "wholebrain_betas")

labels_df <- setDT(read_excel("/Volumes/GoogleDrive/My Drive/SCEPTIC_fMRI/dan_medusa/schaefer_400_remap/MNH DAN Labels 400 Good Only 47 parcels.xlsx")) %>%
  mutate(roi_num7 = as.factor(roi7_400), 
         mask_value = as.integer(roi7_400),
         plot_label = mnh_label_400, 
         vm_gradient17 = parcel_group) %>% select(roi_num7, mask_value, plot_label, vm_gradient17, network17_400_DAN, hemi, x, y, z)

#trial_df <- get_trial_data(repo_directory = "~/Data_Analysis/clock_analysis") %>%
if (session == "meg") {study = "mmclock_meg"} else if (session == "fmri") {study = "mmclock_fmri"}
trial_df <- setDT(get_trial_data(repo_directory = "~/code/clock_analysis", dataset = study)) %>%
# trial_df <- trial_df %>%
  group_by(id, run) %>%
  mutate(v_max_wi_lag = lag(v_max_wi, order_by = run_trial),
         id = as.integer(id)) %>%
  ungroup() %>%
  dplyr::rename(run_number = run) %>%
  dplyr::select(id, run_number, run_trial, trial_neg_inv, rt_csv, rt_lag, v_entropy_wi, v_max_wi_lag, 
                rt_vmax_lag, last_outcome, rewFunc) 
# take a single subject
# trial_df <- trial_df %>% filter(id == 10637)

#echange_l2_copes <- file.path(beta_dir, "L1m-echange/Schaefer2018_400Parcels_7Networks_order_fonov_2.3mm_ants_cope_l2.csv.gz")
# echange_l2_copes <- file.path(beta_dir, "L1m-echange/Schaefer_444_final_2009c_2.3mm_cope_l2.csv.gz")


int <- formula(~ (trial_neg_inv + rt_lag + v_max_wi_lag + v_entropy_wi + fmri_beta + last_outcome)^2 +
                 rt_lag:last_outcome:fmri_beta +
                 rt_vmax_lag*trial_neg_inv*fmri_beta +
                 (1 | id/run_number)
)


slo <- formula(~ (trial_neg_inv + rt_lag + v_max_wi_lag + v_entropy_wi + fmri_beta + last_outcome)^2 +
                 rt_lag:last_outcome:fmri_beta +
                 rt_vmax_lag*trial_neg_inv*fmri_beta +
                 (1 + rt_lag + rt_vmax_lag | id/run_number)
)


# # int only with rewFunc moderation
# rewFunc <- formula(~ (trial_neg_inv + rt_lag + v_max_wi_lag + v_entropy_wi + fmri_beta + last_outcome + rewFunc)^2 +
#                      rt_lag:last_outcome:fmri_beta +
#                      rt_vmax_lag*trial_neg_inv*fmri_beta +
#                      rt_lag*last_outcome*fmri_beta*rewFunc +
#                      rt_vmax_lag*fmri_beta*rewFunc +
#                      (1 | id/run_number)
# )


emtrends_spec <- list(
  rt_lag_int = list(
    outcome = "rt_csv", model_name = "int", var = "rt_lag", specs = c("last_outcome", "fmri_beta"),
    at = list(fmri_beta = c(-2, 0, 2))
  ), # z scores
  rt_lag_slo = list(
    outcome = "rt_csv", model_name = "slo", var = "rt_lag", specs = c("last_outcome", "fmri_beta"),
    at = list(fmri_beta = c(-2, 0, 2))
  ), # z scores
  rt_lag_rewFunc = list(
    outcome = "rt_csv", model_name = "rewFunc", var = "rt_lag", specs = c("last_outcome", "fmri_beta", "rewFunc"),
    at = list(fmri_beta = c(-2, 0, 2)) # z scores
  ),
  rt_vmax_int = list(
    outcome = "rt_csv", model_name = "int", var = "rt_vmax_lag", specs = c("last_outcome", "fmri_beta"),
    at = list(fmri_beta = c(-2, 0, 2)) # z scores
  ),
  rt_vmax_slo = list(
    outcome = "rt_csv", model_name = "slo", var = "rt_vmax_lag", specs = c("last_outcome", "fmri_beta"),
    at = list(fmri_beta = c(-2, 0, 2)) # z scores
  ),
  rt_vmax_rewFunc = list(
    outcome = "rt_csv", model_name = "rewFunc", var = "rt_vmax_lag", specs = c("fmri_beta", "rewFunc"),
    at = list(fmri_beta = c(-2, 0, 2)) # z scores
  )
)

# save snapshot of environment to use inside estimation
save.image(file="parcel_input_snapshot.RData")


efiles_l2 <- list.files(beta_dir,
                        pattern = "Schaefer_444_final_2009c_2.3mm_cope_l2.csv.gz",
                        recursive = TRUE, full.names = TRUE
)

# manually filter just echange
efiles_l2 <- efiles_l2[3]

# # test on one parcel
# labels_df <- labels_df[1,]
# efiles_l2 <- efiles_l2[1]

for (ee in efiles_l2) {
  job <- R_batch_job$new(
    job_name = "parcel_bb", batch_directory = getwd(), scheduler = "local",
    input_rdata_file = "parcel_input_snapshot.RData",
    n_nodes = 1, n_cpus = 18, wall_time = "23:00:00",
    mem_total = "64G",
    #r_code = glue("to_plot <- mixed_by_betas('{ee}', labels_df, trial_df, mask_file = 'Schaefer2018_400Parcels_7Networks_order_fonov_1mm_ants.nii.gz',
    r_code = glue("to_plot <- mixed_by_betas('{ee}', labels_df, trial_df, mask_file = 'Schaefer_444_final_2.3mm.nii.gz',
                            rhs_model_formulae = fmri.pipeline:::named_list(int, slo), ncores = 18, afni_dir = '~/abin',
                            out_prefix = 'Schaefer_400_dan_parcels_fmri_',
    
                                               #                         rhs_model_formulae = fmri.pipeline:::named_list(int, slo), ncores = 16, afni_dir = '/proj/mnhallqlab/sw/afni',
                                               #                         calculate = c('parameter_estimates_reml', 'fit_statistics'),
                                               #                         beta_level = 2L, focal_contrast = 'overall', emtrends_spec = emtrends_spec,
                                               split_on = c('l1_cope_name', 'l2_cope_name', 'mask_value'))"),
    r_packages = c("fmri.pipeline", "tidyverse", "data.table", "sfsmisc"),
    batch_code = c("module use /proj/mnhallqlab/sw/modules", "module load r/4.1.2_depend")
  )
  
  job$submit()
  # original glue:
  # r_code = glue("to_plot <- mixed_by_betas('{ee}', labels_200, trial_df, mask_file = 'Schaefer2018_200Parcels_7Networks_order_fonov_1mm_ants.nii.gz',
  
  
  # local execution
  # to_plot <- mixed_by_betas(ee, labels_df, trial_df, mask_file = "Schaefer_444_final_2.3mm.nii.gz",
  #                           rhs_model_formulae = fmri.pipeline:::named_list(int, slo), ncores = 16,
  #                           split_on = c("l1_cope_name", "l2_cope_name", "mask_value"))

}



efiles_l1 <- list.files(beta_dir,
  pattern = "Schaefer_444_final_2009c_2.3mm_cope_l1.csv.gz",
  recursive = TRUE, full.names = TRUE
)

efiles_l1 <- efiles_l1[c(6, 9, 11)] # entropy, pe, echange
l1_contrasts <- c("EV_entropy_change_feedback", "EV_entropy_wiz_clock", "EV_pe")

# l1 brain-behavior
for (ee in seq_along(efiles_l1)) {
  job <- R_batch_job$new(
    job_name = "parcel_bb", batch_directory = getwd(), scheduler = "slurm",
    input_rdata_file = "parcel_input_snapshot.RData",
    n_nodes = 1, n_cpus = 16, wall_time = "12:00:00",
    mem_total = "96G",
    r_code = glue("to_plot <- mixed_by_betas('{efiles_l1[ee]}', labels_200, trial_df, mask_file = 'Schaefer2018_200Parcels_7Networks_order_fonov_1mm_ants.nii.gz',
                                             rhs_model_formulae = fmri.pipeline:::named_list(rewFunc), ncores = 16, afni_dir = '/proj/mnhallqlab/sw/afni',
                                             calculate = c('parameter_estimates_reml', 'fit_statistics'),
                                             trial_join_col = c('id', 'run_number'), beta_level = 1L, focal_contrast = '{l1_contrasts[ee]}', emtrends_spec = emtrends_spec,
                                             split_on = c('l1_cope_name', 'mask_value'))"),
    r_packages = c("fmri.pipeline", "tidyverse", "data.table", "sfsmisc"),
    batch_code = c("module use /proj/mnhallqlab/sw/modules", "module load r/4.1.2_depend")
  )

  #eval(parse(text = job$r_code))
  job$submit()

  # local execution
  # to_plot <- mixed_by_betas(ee, labels_df, trial_df, mask_file = "Schaefer2018_200Parcels_7Networks_order_fonov_1mm_ants.nii.gz",
  #                           rhs_model_formulae = fmri.pipeline:::named_list(int, slo), ncores = 16,
  #                           split_on = c("l1_cope_name", "l2_cope_name", "mask_value"))
}




####




# these are betas I manually extracted by using colMeans, readNifti, and mask indices in the old/traditional voxelwise approach.
# see entropy_betas.R. These are specifically the entropy change overall betas from the echange model

manual_betas <- read.csv("/Volumes/GoogleDrive/My Drive/SCEPTIC_fMRI/wholebrain_betas/beta_debugging/schaefer200_echange_overall_l2_means.csv") %>%
  mutate(id = 1:n()) %>%
  pivot_longer(names_to = "roi", cols = -id, values_to = "fmri_beta") %>%
  mutate(mask_value = as.integer(sub("^V", "", roi))) %>%
  dplyr::select(-roi)

onesamp_betas(manual_betas, mask_file="Schaefer2018_400Parcels_7Networks_order_fonov_1mm_ants.nii.gz",
              out_dir = "/Volumes/GoogleDrive/My Drive/SCEPTIC_fMRI/wholebrain_betas/beta_debugging")

man_means <- manual_betas %>% group_by(mask_value) %>% summarise(mean = mean(fmri_beta))

# these are betas for the overall echange contrast extracted at the time of the extract_glm_betas_in_mask
echange_intermediate <- readRDS("/Volumes/GoogleDrive/My Drive/SCEPTIC_fMRI/wholebrain_betas/beta_debugging/overall_echange_stat_result.rds") %>%
  unnest(img_stats) %>% select(-img, -img_exists, -x, -y, -z)

# int_means <- echange_intermediate %>% group_by(mask_value) %>% summarise(mean = mean(value))
# 
# onesamp_betas(echange_intermediate, mask_file="Schaefer2018_400Parcels_7Networks_order_fonov_1mm_ants.nii.gz", dv = "value",
#               out_dir = "/Volumes/GoogleDrive/My Drive/SCEPTIC_fMRI/wholebrain_betas/beta_debugging", img_prefix = "echange_intermediate")
# 


#setwd("/Users/michael/Data_Analysis/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/betas_final")

######



###


#beta_dir <- "/Volumes/GoogleDrive/My Drive/SCEPTIC_fMRI/final_betas"



echange_l2_copes <- file.path(beta_dir, "L1m-echange/Schaefer_444_final_2009c_2.3mm_cope_l2.csv.gz")

mixed_by_betas(echange_l2_copes, labels_df, trial_df)

# subject-level
#echange_l2_copes <- fread(file.path(beta_dir, "L1m-echange/Schaefer_DorsAttn_2.3mm_cope_l2.csv.gz")) %>%
#echange_l2_copes <- fread(file.path(beta_dir, "L1m-echange/zstat6_ptfce_Schaefer2018_400Parcels_7Networks_order_fonov_2.3mm_ants_overlap_cope_l2.csv.gz")) %>%
#echange_l2_copes <- fread(file.path(beta_dir, "L1m-echange/Schaefer_444_final_2009c_2.3mm_cope_l2.csv.gz")) %>%
entropy_l2_copes <- fread(file.path(beta_dir, "L1m-entropy_wiz/Schaefer_444_final_2009c_2.3mm_cope_l2.csv.gz")) %>%
  #filter(l1_cope_name == "EV_entropy_change_feedback" & l2_cope_name == "overall") %>%
  filter(l1_cope_name == "EV_entropy_wiz_clock" & l2_cope_name == "overall") %>%
  select(-feat_dir, -img) %>%
  rename(fmri_beta=value) %>%
  left_join(labels_df, by="mask_value")

# meg effects of interest
# meg_wide <- readRDS("MEG_betas_wide_echange_vmax_reward_Nov15_2021.RDS")

# keep as long to have meg beta as split

meg_ranefs <- readRDS("MEG_betas_echange_vmax_reward_Nov30_2021.RDS") %>%
  mutate(id=as.numeric(id))

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


trial_set <- trial_df %>% 
  group_by(id, run) %>%
  mutate(v_max_wi_lag = lag(v_max_wi, order_by=run_trial)) %>%
  ungroup() %>%
  dplyr::select(id, run, run_trial, trial_neg_inv, rt_csv, rt_lag, v_entropy_wi, v_max_wi_lag, 
                rt_vmax_lag, last_outcome)

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
  left_join(labels_df, by="mask_value") %>%
  setDT()


#write.csv(to_plot %>% select(-rhs, -effect), file="fmri_brainbehavior_parcel_betas.csv", row.names=FALSE)
#write.csv(to_plot %>% select(-rhs, -effect), file="fmri_brainbehavior_parcel_betas_ptfce.csv", row.names=FALSE)
#write.csv(to_plot %>% select(-rhs, -effect), file="fmri_brainbehavior_parcel_betas_400.csv", row.names=FALSE)
write.csv(to_plot %>% select(-rhs, -effect), file="fmri_brainbehavior_parcel_entropy_betas_400.csv", row.names=FALSE)


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

                   split_on = splits, scale_predictors = c("trial_neg_inv", "rt_lag", "rt_vmax_lag", "v_max_wi_lag", "fmri_beta", "avg"),
                   tidy_args = list(effects = c("fixed", "ran_vals", "ran_pars", "ran_coefs"), conf.int = TRUE), 
                   calculate = c("parameter_estimates_reml"), ncores = 8, refit_on_nonconvergence = 5, padjust_by = "term",
                   emtrends_spec = list(
                     rt_lag.meg = list(outcome = "rt_csv", model_name = "main", var = "rt_lag", specs = c("last_outcome", "avg"), 
                                       at=list(avg = c(-2, 0, 2))),
                     rt_lag.fmri = list(outcome = "rt_csv", model_name = "main", var = "rt_lag", specs = c("last_outcome", "fmri_beta"), 
                                        at=list(fmri_beta = c(-2, 0, 2)))
                   )
)


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


plot_medusa(plot_df, x="plot_label", y="reg_region", color="statistic", term_filter="(avg|fmri_beta)",
            out_dir="parcel_betas_meg", p.value="p_FDR", plot_type="heat", width=15, height=7, )


emt <- meg_df$emtrends_list$rt_lag.meg
emt <- emt[,-3:-4]
emt <- emt %>% left_join(labels, by="mask_value") %>%
  mutate(avg=factor(avg, levels=c(-2, 0, 2), labels=c("low", "mean", "high")))
pdf("meg_emtrends_rt_lag.reward_with_fmri_betas.pdf", height = 15, width = 15)
ggplot(emt, aes(x=plot_label, y=rt_lag.trend, color=avg, ymin=rt_lag.trend - std.error, ymax = rt_lag.trend + std.error)) +
  geom_pointrange() + coord_flip() + facet_grid(reg_region~last_outcome)
dev.off()

emt1 <- meg_df$emtrends_list$rt_lag.fmri
emt1 <- emt1[,-3:-4]
emt1 <- emt1 %>% left_join(labels, by="mask_value") %>%
  mutate(fmri_beta=factor(fmri_beta, levels=c(-2, 0, 2), labels=c("low", "mean", "high")))

pdf("fmri_emtrends_rt_lag.reward_with_meg_betas.pdf", height = 15, width = 15)
ggplot(emt1, aes(x=plot_label, y=rt_lag.trend, color=fmri_beta, ymin=rt_lag.trend - std.error, ymax = rt_lag.trend + std.error)) +
  geom_pointrange() + coord_flip() + facet_grid(reg_region~last_outcome)
dev.off()


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
