# final parcel-wise analyses of DAN brain-to-behavior
library(data.table)
library(tidyverse)
library(afex)
library(lattice)
library(emmeans)
library(fmri.pipeline) # has mixed_By
setwd("~/Data_Analysis/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/betas_final")

source("../get_trial_data.R")
#source("~/Data_Analysis/r_packages/fmri.pipeline/R/mixed_by.R")

# location of whole brain betas for analysis
beta_dir <- "/Volumes/GoogleDrive/My Drive/SCEPTIC_fMRI/wholebrain_betas"

# MNH/AD internal labeling scheme
labels <- readxl::read_excel("~/Data_Analysis/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/MNH Dan Labels.xlsx") %>%
  dplyr::rename(mask_value=roinum) %>% select(mask_value, plot_label)

# BALSA parcel labels from whereami
labels_200 <- read.csv("/Volumes/GoogleDrive/My Drive/SCEPTIC_fMRI/wholebrain_betas/schaefer_200_whereami.csv") %>%
  dplyr::rename(mask_value = roi_num, plot_label = MNI_Glasser_HCP_v1.0) %>% dplyr::select(mask_value, plot_label) %>%
  mutate(plot_label = sub("Focus point:\\s+", "", plot_label, perl=TRUE))

trial_df <- get_trial_data(repo_directory = "~/Data_Analysis/clock_analysis") %>%
  group_by(id, run) %>%
  mutate(v_max_wi_lag = lag(v_max_wi, order_by=run_trial)) %>%
  ungroup() %>%
  dplyr::select(id, run, run_trial, trial_neg_inv, rt_csv, rt_lag, v_entropy_wi, v_max_wi_lag, 
                rt_vmax_lag, last_outcome)


#' run mixed_by model for each split within a csv file containing betas extracted from a GLM, esp. by extract_glm_betas_in_mask
#'   from fmri.pipeline
#' @param beta_csv The filename of the CSV containing betas/coefficients to be analyzed parcel-by-parcel
#' @param label_df A data.frame containing labels for each unique mask value.
#' @param trial_df A trial data.frame containing the data for the trial-level multilevel model
#' @param mask_file The name of the mask_file whose mask values correspond to extracted betas in \code{beta_csv}. Used
#'   to populate NIfTI images and surface/CIFTI files with the results of the mixed_by analyses
#' @param label_join_col The column name in \code{beta_csv} and \code{label_df} used to match-merge them
#' @param trial_join_col The column name in \code{beta_csv} and \code{trial_df} used to match-merge them
#' @param split_on The factors in \code{beta_csv} used to split the data into separate datasets for multilevel models.
#'   Default is 'mask_value', but if the beta_csv contains many contrasts, it may be more like, c('mask_value', 'l1_cope_name', 'l2_cope_name')
#' @param rhs_form The right-hand side formula used in the lmer() calls by mixed_by
#' @param ncores The number of cores to use for parallel computation of splits (passed through to mixed_by)
#' @param out_dir The output directory for statistic and image files
#'
#' @details Note that by default, this function also computes the robust mean of betas/coefficients in each mask_value and provides
#'   an overall map of means by parcel/mask value. This is useful to see in which parcels whole-brain activation is robust and significant
#'   and it provides a validation of the parcelwise beta extraction against the corresponding voxelwise maps from which the betas are drawn
#' @return The data.frame containing coefficients by splits (incl. mask value)
mixed_by_betas <- function(beta_csv, label_df, trial_df, mask_file = NULL, label_join_col = "mask_value", trial_join_col = "id", 
                           split_on = "mask_value", rhs_form=NULL, ncores = 16, out_dir = getwd()) {
  
  checkmate::assert_data_frame(label_df)
  checkmate::assert_file_exists(beta_csv)
  
  # the big betas file has all of the L2 contrasts, which are largely uninteresting and make the dataset massive
  # for now, subset down to overall contrast at L2.
  
  cope_df <- fread(beta_csv) %>%
    filter(l2_cope_name == "overall") %>%
    #filter(l1_cope_name == "EV_entropy_change_feedback" & l2_cope_name == "overall") %>%
    #filter(l1_cope_name == "EV_entropy_wiz_clock" & l2_cope_name == "overall") %>%
    dplyr::select(-feat_dir, -img, -mask_name, -session, -l1_cope_number, -l2_cope_number, -l2_model) %>%
    rename(fmri_beta=value) %>%
    merge(label_df, by=label_join_col, all.x = TRUE) 
  
  combo <- cope_df %>%
    merge(trial_df, by=trial_join_col, all.x = TRUE, allow.cartesian=TRUE) # trial_df and betas get crossed
    #merge(trial_df, by=.EACHI, all.x = TRUE)
  
  # %>%
  #   left_join(label_df, by=label_join_col) %>%
  #   left_join(trial_df, by=trial_join_col)
  
  # divide mixed_by by contrast
  # if ("l2_cope_name" %in% names(cope_df)) {
  #   cope_split <- split(cope_df, by=c("l1_cope_name", "l2_cope_name"))
  # } else {
  #   cope_split <- split(cope_df, by="l1_cope_name")
  # }
  # 
  # browser()
  # # run mixed_by for every contrast combination
  # for (df in cope_split) {
  
  # get overall maps for each parcel using a robust lm (one-sample t-test)
  # abc <- cope_df %>% filter(mask_value==1 & l1_cope_name == "EV_clock" & l2_cope_name == "overall")
  # 
  # zlm <- MASS::rlm(fmri_beta ~ 1, abc)
  # f.robftest(zlm, var=1)
  # 
  # cope_split <- split(cope_df, by=c("l1_cope_name", "l2_cope_name", "mask_value"))
  
  # use nest approach
  cope_test <- cope_df %>% group_by(l1_model, l1_cope_name, l2_cope_name, mask_value) %>% nest()
  
  # run RLM for every combination of contrasts and mask value
  # keep a data.frame that has the t and p values to build an afni dataset
  res <- cope_test %>% 
    mutate(model = map(data, function(df) {
      mobj <- MASS::rlm(fmri_beta ~ 1, df)
      ftest <- sfsmisc::f.robftest(mobj, 1)
      return(data.frame(t=summary(mobj)$coefficients[1,"t value"], p=ftest$p.value, logp=-1*log10(ftest$p.value)))
    })) %>% dplyr::select(-data) %>% unnest(model) %>% ungroup() %>% setDT()
  
  fill_mask_with_stats(mask_file, mask_col = "mask_value", stat_dt = res, subbrik_cols = c("t", "p", "logp"), 
                       split_on=c("l1_model", "l1_cope_name", "l2_cope_name"), afni_dir="~/afni", out_dir = file.path(out_dir, "parcel_mean"))
  
  ddf <- mixed_by(combo, outcomes = "rt_csv", rhs_model_formulae = list(main=rhs_form),
                  split_on = c("mask_value", "l1_cope_name", "l2_cope_name"), 
                  scale_predictors = c("trial_neg_inv", "rt_lag", "rt_vmax_lag", "v_max_wi_lag", "fmri_beta"),
                  tidy_args = list(effects = c("fixed"), conf.int = TRUE), 
                  calculate = c("parameter_estimates_reml"), ncores = ncores, refit_on_nonconvergence = 5, padjust_by = "term",
                  emtrends_spec = list(
                    rt_lag = list(outcome = "rt_csv", model_name = "main", var = "rt_lag", specs = c("last_outcome", "fmri_beta"), 
                                  at=list(fmri_beta = c(-2, 0, 2))) # z scores
                  ))
  
  saveRDS(ddf, file=file.path(out_dir, "test.rds"))
  
  to_plot <- ddf$coef_df_reml %>%
    filter(effect=="fixed") %>%
    group_by(term) %>%
    mutate(p_FDR=p.adjust(p.value, method="fdr")) %>%
    ungroup() %>% 
    #left_join(labels, by="mask_value") %>%
    left_join(label_df, by=!!label_join_col) %>%
    select(-rhs, -effect) %>%
    setDT()
  
  #fwrite(to_plot, file="fmri_brainbehavior_parcel_entropy_betas_200.csv", row.names=FALSE)
  return(to_plot)
}

fill_mask_with_stats <- function(mask_file, mask_col = "mask_value", stat_dt, subbrik_cols = c("t", "p"), subbrik_labels = NULL, split_on=NULL, 
                                 img_prefix = "maskfill", afni_dir = "~/abin", out_dir = getwd()) {
  require(RNifti)
  require(glue)
  
  checkmate::assert_file_exists(mask_file)
  vol <- RNifti::readNifti(mask_file)
  
  checkmate::assert_directory_exists(afni_dir)
  afni_dir <- normalizePath(afni_dir)
  
  checkmate::assert_data_table(stat_dt)
  checkmate::assert_subset(subbrik_cols, names(stat_dt))
  if (is.null(subbrik_labels)) {
    subbrik_labels <- subbrik_cols
  }
  
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  
  if (!is.null(split_on)) {
    stat_dt <- split(stat_dt, by=split_on)
    img_names <- paste0(names(stat_dt), ".nii.gz")
  } else {
    stat_dt <- list(stat_dt) # single element
    img_names <- img_prefix
  }
  
  # loop over splits (e.g., contrasts)
  for (dd in seq_along(stat_dt)) {
    data <- stat_dt[[dd]]
    to_combine <- sapply(seq_along(subbrik_cols), function(x) { tempfile(fileext=".nii.gz") })
    
    # loop over subbriks to create (e.g., t, p, and logp)
    for (ff in seq_along(subbrik_cols)) {
      new_stat <- vol # copy
      new_stat <- new_stat*0 # start with all zeros
      
      # loop over and populate each mask value
      for (val in unique(data$mask_value)) {
        new_stat[vol == val] <- data[[subbrik_cols[ff]]][data$mask_value == val] # the RHS should be one value... probably validate this if I extend the function
      }
      
      writeNifti(new_stat, file=to_combine[ff])
    }
    
    system(glue("{afni_dir}/3dTcat -overwrite -prefix {out_dir}/{img_names[dd]} {paste(to_combine, sep=' ')}"))
    system(glue("{afni_dir}/3drefit -relabel_all_str '{paste(subbrik_labels, sep=' ')}' {out_dir}/{img_names[dd]}"))
    unlink(to_combine)
  }
}

echange_l2_copes <- file.path(beta_dir, "L1m-echange/Schaefer2018_200Parcels_7Networks_order_fonov_2.3mm_ants_cope_l2.csv.gz")
model_base <- formula(~ (trial_neg_inv + rt_lag + v_max_wi_lag + v_entropy_wi + fmri_beta + last_outcome)^2 +
                        rt_lag:last_outcome:fmri_beta +
                        rt_vmax_lag*trial_neg_inv*fmri_beta +
                        (1 | id/run)
)


to_plot <- mixed_by_betas(echange_l2_copes, labels_200, trial_df, mask_file = "Schaefer2018_200Parcels_7Networks_order_fonov_1mm_ants.nii.gz",
                          rhs_form = model_base)




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
