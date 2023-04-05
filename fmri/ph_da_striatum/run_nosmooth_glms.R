library(dplyr)
library(emmeans)
library(fmri.pipeline)
library(readr)

source("/proj/mnhallqlab/projects/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/get_trial_data.R")
trial_df <- get_trial_data(repo_directory = "/proj/mnhallqlab/projects/clock_analysis", dataset = "mmclock_fmri", groupfixed = TRUE) %>%
  dplyr::rename(run_number = "run")

subject_df <- readRDS("/proj/mnhallqlab/users/michael/fmri.pipeline/inst/example_files/mmclock_subject_data.rds") %>%
  mutate(mr_dir=paste0(mr_dir, "/mni_nosmooth_aroma")) #make sure we're looking in the right folder

# run_df <- readRDS("/proj/mnhallqlab/users/michael/fmri.pipeline/inst/example_files/mmclock_run_data.rds")

gpa <- setup_glm_pipeline(analysis_name="mmclock_nosmooth", scheduler="slurm",
  output_directory = "/proj/mnhallqlab/no_backup/mmclock_fmri",
  subject_data=subject_df, trial_data=trial_df, # run_data=run_df,
  tr=1.0, drop_volumes = 2,
  n_expected_runs=8,
  l1_models=NULL, l2_models=NULL, l3_models=NULL,
  fmri_file_regex="nfawuktm_clock[1-8]_vtavshipp_blur\\.nii\\.gz",
  fmri_path_regex="clock[0-9]",
  run_number_regex=".*clock([0-9]+)_vta.*",
  confound_settings=list(
    motion_params_file = "motion.par",
    confound_input_file="nuisance_regressors.txt",
    confound_input_colnames = c("csf", "dcsf", "wm", "dwm"),
    l1_confound_regressors = c("csf", "dcsf", "wm", "dwm"),
    exclude_run = "max(FD) > 5 | sum(FD > .9)/length(FD) > .10", #this must evaluate to a scalar per run
    exclude_subject = "n_good_runs < 4",
    truncate_run = "(FD > 0.9 & time > last_offset) | (time > last_offset + last_isi)",
    spike_volumes = NULL
  ),
  parallel = list(
    fsl = list(
      l1_feat_alljobs_time = "144:00:00"
    )
  )
)

gpa <- build_l1_models(gpa, from_spec_file = "/proj/mnhallqlab/projects/clock_analysis/fmri/ph_da_striatum/nosmooth_l1_models.yaml")
#gpa$parallel$fsl$l1_feat_alljobs_time <- "144:00:00"

gpa <- build_l2_models(gpa)
gpa <- build_l3_models(gpa)

saveRDS(gpa, file = "/proj/mnhallqlab/projects/clock_analysis/fmri/ph_da_striatum/nosmooth_glm_cache_3Apr2023.rds")

run_glm_pipeline(gpa)

# after GLMs complete
atlases <- c(
  "/proj/mnhallqlab/projects/clock_analysis/fmri/ph_da_striatum/masks/long_axis_l_12bins_2.3mm.nii.gz",
  "/proj/mnhallqlab/projects/clock_analysis/fmri/ph_da_striatum/masks/long_axis_r_12bins_2.3mm.nii.gz",
  "/proj/mnhallqlab/projects/clock_analysis/fmri/ph_da_striatum/masks/vta_vs_2.3mm.nii.gz"
)

betas <- extract_glm_betas_in_mask(
  gpa, atlases,
  extract_l3 = "none", extract_l1 = "prompt", extract_l2 = "prompt", what = c("cope", "zstat"),
  scheduler = "slurm", out_dir = file.path("/proj/mnhallqlab/no_backup/mmclock_fmri/mmclock_nosmooth/extracted_betas")
)


