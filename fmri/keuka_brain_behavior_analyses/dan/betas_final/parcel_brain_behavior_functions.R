
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
#'   an overall map of means by parcel/mask value. This is usezful to see in which parcels whole-brain activation is robust and significant
#'   and it provides a validation of the parcelwise beta extraction against the corresponding voxelwise maps from which the betas are drawn
#' @return The data.frame containing coefficients by splits (incl. mask value)
mixed_by_betas <- function(beta_csv, label_df, trial_df, mask_file = NULL, label_join_col = "mask_value", trial_join_col = "id", 
                           split_on = c("mask_value", "l1_cope_name"), rhs_form=NULL, ncores = 16, out_dir = NULL, out_prefix = NULL) {
  
  checkmate::assert_data_frame(label_df)
  checkmate::assert_file_exists(beta_csv)
  if (checkmate::test_formula(rhs_form)) {
    rhs_form <- list(main = rhs_form) # always work from list
  }
  
  if (is.null(out_dir)) {
    out_dir <- file.path(normalizePath(dirname(beta_csv)), "parcel_maps")
  }
  
  if (!dir.exists(out_dir)) { dir.create(out_dir) }
  
  if (is.null(out_prefix)) {
    out_prefix <- fmri.pipeline:::file_sans_ext(basename(beta_csv))
  }
  
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
  
  
  # run one-sample test statistics for each parcel to corroborate beta extraction against whole-brain voxelwise analysis
  onesamp_betas(cope_df, mask_file=mask_file, roi_column = "mask_value", 
                nest_by = c("l1_model", "l1_cope_name", "l2_cope_name"),
                out_dir = out_dir, img_prefix = "onesamp")
  
  out_file <- file.path(out_dir, paste0(out_prefix, "_mixed_by.rds"))
  
  if (file.exists(out_file)) {
    ddf <- readRDS(out_file)
  } else {
    # run mixed by across splits
    ddf <- mixed_by(combo, outcomes = "rt_csv", rhs_model_formulae = rhs_form,
                    split_on = split_on, 
                    scale_predictors = c("trial_neg_inv", "rt_lag", "rt_vmax_lag", "v_max_wi_lag", "fmri_beta"),
                    tidy_args = list(effects = c("fixed"), conf.int = TRUE), 
                    calculate = c("parameter_estimates_reml"), ncores = ncores, refit_on_nonconvergence = 5, padjust_by = "term",
                    emtrends_spec = list(
                      rt_lag = list(outcome = "rt_csv", model_name = "main", var = "rt_lag", specs = c("last_outcome", "fmri_beta"), 
                                    at=list(fmri_beta = c(-2, 0, 2))) # z scores
                    ))  
    
    rm(combo)
    saveRDS(ddf, file=out_file)
  }
  
  to_plot <- ddf$coef_df_reml %>%
    filter(effect=="fixed" & grepl("fmri_beta", term)) %>%
    dplyr::rename(p = p.value, t = statistic) %>%
    mutate(logp = -1*log10(p)) %>%
    group_by(term, model_name) %>%
    mutate(p_FDR=p.adjust(p, method="fdr")) %>%
    ungroup() %>% 
    left_join(label_df, by=label_join_col) %>%
    select(-rhs, -effect) %>%
    setDT()
  
  fill_mask_with_stats(mask_file, mask_col = "mask_value", stat_dt = to_plot, subbrik_cols = c("t", "p", "logp", "p_FDR"),
                       split_on=c("l1_cope_name", "l2_cope_name", "term", "model_name"), afni_dir="~/abin", out_dir = out_dir, img_prefix = NULL)
  
  #fwrite(to_plot, file="fmri_brainbehavior_parcel_entropy_betas_200.csv", row.names=FALSE)
  return(to_plot)
}

fill_mask_with_stats <- function(mask_file, mask_col = "mask_value", stat_dt, subbrik_cols = c("t", "p"), subbrik_labels = NULL, split_on=NULL, 
                                 img_prefix = "maskfill", afni_dir = "~/abin", out_dir = getwd(), overwrite = FALSE) {
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
  
  if (is.null(out_dir)) { out_dir <- getwd() }
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  
  if (!is.null(split_on)) {
    stat_dt <- split(stat_dt, by=split_on)
    if (is.null(img_prefix)) {
      img_names <- paste0(make.names(names(stat_dt)), ".nii.gz")  
    } else {
      img_names <- paste0(img_prefix, "_", make.names(names(stat_dt)), ".nii.gz")  
    }
  } else {
    stat_dt <- list(stat_dt) # single element
    img_names <- paste0(img_prefix, ".nii.gz")
  }
  
  # loop over splits (e.g., contrasts)
  for (dd in seq_along(stat_dt)) {
    this_out_file <- glue("{out_dir}/{img_names[dd]}")
    if (checkmate::test_file_exists(this_out_file) && isFALSE(overwrite)) {
      message(glue("Skipping existing file: {this_out_file}"))
      next
    }
    
    data <- stat_dt[[dd]]
    to_combine <- sapply(seq_along(subbrik_cols), function(x) { tempfile(fileext=".nii.gz") })
    
    # loop over subbriks to create (e.g., t, p, and logp)
    for (ff in seq_along(subbrik_cols)) {
      if (file.exists(to_combine[ff]) && isFALSE(overwrite)) {
        message(glue("Skipping existing file {to_combine[ff]}"))
        next
      }
      
      new_stat <- vol # copy
      new_stat <- new_stat*0 # start with all zeros
      
      # loop over and populate each mask value
      for (val in unique(data$mask_value)) {
        stat_val <- data[[subbrik_cols[ff]]][data$mask_value == val]
        if (length(stat_val) > 1L) { stop("more than one stat value found for a parcel. Problem with splits?")}
        new_stat[vol == val] <- stat_val # the RHS should be one value... probably validate this if I extend the function
      }
      
      writeNifti(new_stat, file=to_combine[ff])
    }
    
    # AFNI does not tolerate spaces in output file, so need to use a tmpfile and then move it
    tmpout <- tempfile(fileext = ".nii.gz")
    
    system(glue("{afni_dir}/3dTcat -overwrite -prefix '{tmpout}' {paste(to_combine, collapse=' ')}"))
    system(glue("{afni_dir}/3drefit -relabel_all_str '{paste(subbrik_labels, collapse=' ')}' '{tmpout}'"))
    
    file.copy(tmpout, this_out_file, overwrite = TRUE)
    unlink(c(tmpout, to_combine))
  }
}

onesamp_betas <- function(cope_df, mask_file, dv = "fmri_beta", roi_column = "mask_value", nest_by = NULL, 
                          out_dir = NULL, img_prefix = "onesamp") {
  require(tidyr)
  require(data.table)
  require(dplyr)
  
  # use nest approach
  cope_nest <- cope_df %>% group_by(across(all_of(unique(c(roi_column, nest_by))))) %>% nest()
  
  # run RLM for every combination of contrasts and mask value
  # keep a data.frame that has the t and p values to build an afni dataset
  res <- cope_nest %>% 
    mutate(model = map(data, function(df) {
      ff <- as.formula(paste(dv, "~ 1"))
      mobj <- MASS::rlm(ff, df)
      ftest <- sfsmisc::f.robftest(mobj, 1)
      return(data.frame(t=summary(mobj)$coefficients[1,"t value"], negp=1-ftest$p.value, logp=-1*log10(ftest$p.value)))
    })) %>% dplyr::select(-data) %>% unnest(model) %>% ungroup() %>% setDT()
  
  # write mean results to out_dir
  fill_mask_with_stats(mask_file, mask_col = "mask_value", stat_dt = res, subbrik_cols = c("t", "negp", "logp"), 
                       subbrik_labels = c("t", "1-p", "neglogp"),
                       split_on=nest_by, afni_dir="~/abin", out_dir = out_dir, img_prefix = img_prefix)
  
  return(res)
}



# LEFTOVERS

# get overall maps for each parcel using a robust lm (one-sample t-test)
# abc <- cope_df %>% filter(mask_value==1 & l1_cope_name == "EV_clock" & l2_cope_name == "overall")
# 
# zlm <- MASS::rlm(fmri_beta ~ 1, abc)
# f.robftest(zlm, var=1)
# 
# cope_split <- split(cope_df, by=c("l1_cope_name", "l2_cope_name", "mask_value"))

# use nest approach
# cope_nest <- cope_df %>% group_by(l1_model, l1_cope_name, l2_cope_name, mask_value) %>% nest()
# 
# # run RLM for every combination of contrasts and mask value
# # keep a data.frame that has the t and p values to build an afni dataset
# res <- cope_nest %>% 
#   mutate(model = map(data, function(df) {
#     mobj <- MASS::rlm(fmri_beta ~ 1, df)
#     ftest <- sfsmisc::f.robftest(mobj, 1)
#     return(data.frame(t=summary(mobj)$coefficients[1,"t value"], p=ftest$p.value, logp=-1*log10(ftest$p.value)))
#   })) %>% dplyr::select(-data) %>% unnest(model) %>% ungroup() %>% setDT()
# 
# rm(cope_nest)
# 
# # write mean results to out_dir
# fill_mask_with_stats(mask_file, mask_col = "mask_value", stat_dt = res, subbrik_cols = c("t", "p", "logp"), 
#                      split_on=c("l1_model", "l1_cope_name", "l2_cope_name"), afni_dir="~/afni", out_dir = out_dir, img_prefix = "onesamp")
