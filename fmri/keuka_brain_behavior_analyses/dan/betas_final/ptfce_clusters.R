# hack betas of interest from group maps in entropy change
library(glue)
library(fmri.pipeline)
library(dplyr)
library(parallel)

atlas <- "/proj/mnhallqlab/projects/clock_analysis/fmri/pfc_entropy/original_masks/Schaefer2018_200Parcels_7Networks_order_fonov_2.3mm_ants.nii.gz"
basedir <- "/proj/mnhallqlab/users/michael/mmclock_entropy/mmclock_nov2021/feat_l3/L1m-echange/L2m-l2_l2c-overall"
echange_gfeats <- system(glue("find {basedir} -iname '*entropy_change_feedback.gfeat' -type d"), intern=TRUE)

y <- ptfce_spec$new(
  gfeat_dir = echange_gfeats,
  fwe_p = c(.05, .01)
)

#y$submit()

y$is_complete()
res <- y$get_clusters()

# res$cluster_obj[[1]]$get_call()
# res$cluster_obj[[1]]$is_complete()

echange_clusters <- res %>% filter(contrast_name == "overall" & grepl("age_sex", z))
obj <- echange_clusters$cluster_obj[[1]]
obj$get_clust_df()

obj$subset_atlas_against_clusters(atlas, minimum_overlap = 0.50, output_atlas = "default", mask_by_overlap=TRUE)
obj$get_outputs()
echange_rois <- read_csv(obj$get_outputs()$atlas_files[[1]]["summary"]) %>% rename(roi_num = roi_val)

schaefer_200_wami <- read.csv("/proj/mnhallqlab/projects/clock_analysis/fmri/pfc_entropy/original_masks/schaefer_200_whereami.csv")
echange_parcels <- schaefer_200_wami %>% inner_join(echange_rois) %>% filter(retained == TRUE)

write_csv(echange_parcels, file="schaefer_echange_parcel_overlap.csv")


# all schaefers that overlap overall echange FWE map at > .75
#gpa <- readRDS("/proj/mnhallqlab/users/michael/mmclock_entropy/mmclock_nov2021/mmclock_nov2021.rds")
gpa <- readRDS("/proj/mnhallqlab/users/michael/mmclock_pe/mmclock_nov2021/mmclock_nov2021.rds")

# WB entropy change betas
# mask_file <- "/proj/mnhallqlab/users/michael/mmclock_entropy/mmclock_nov2021/feat_l3/L1m-echange/L2m-l2_l2c-overall/L3m-age_sex/FEAT_l1c-EV_entropy_change_feedback.gfeat/cope1.feat/stats/zstat6_ptfce_Schaefer2018_200Parcels_7Networks_order_fonov_2.3mm_ants_overlap.nii.gz"
# betas <- extract_glm_betas_in_mask(gpa, mask_file, extract_l3="none", extract_l1="prompt", extract_l2="prompt", what=c("cope", "zstat"), ncores=8)

# all 200 parcels
betas <- extract_glm_betas_in_mask(
  gpa, atlas,
  extract_l3 = "none", extract_l1 = "prompt", extract_l2 = "prompt", what = c("cope", "zstat"),
  scheduler="slurm", out_dir = file.path("/proj/mnhallqlab/users/michael/mmclock_pe/mmclock_nov2021/extracted_betas_new")
)


# ENTROPY CHANGE MAPS
basedir <- "/proj/mnhallqlab/users/michael/mmclock_entropy/mmclock_nov2021/feat_l3/L1m-echange/L2m-l2_l2c-overall"
echange_gfeats <- system(glue("find {basedir} -iname '*entropy_change_feedback.gfeat' -type d"), intern=TRUE)

y <- ptfce_spec$new(
  gfeat_dir = entropy_gfeats,
  fwe_p = c(.05, .01)
)

#y$submit()

y$is_complete()
res <- y$get_clusters()

echange_overall <- res %>% filter(contrast_name == "overall" & grepl("age_sex", z))
obj <- echange_overall$cluster_obj[[1]]
abc <- obj$get_clust_df()

# generates 7 subclusters
obj$generate_subclusters(break_nvox = 1000, min_subclust_nvox = 30, max_subclust_nvox = 3200, print_progress = TRUE)
echange_df <- obj$get_clust_df()

echange_df %>% select(roi_num, subroi_num, Volume, CM.LR, CM.PA, CM.IS, Mean, Max.Int, MNI_Glasser_HCP_v1.0, Brainnetome_1.0, CA_ML_18_MNI)

first_roi <- echange_df %>%
  slice(1) %>%
  unnest(overlap) %>%
  select(roi_num, Volume, atlas, label, overlap_pct)





## PTFCE GIVES BIG MAPS for ECHANGE
# here are the acf results for 3dClustsim


# try 3dclustsim acf for entropy change
# these ACF results should apply to all group analyses since the L1 runs are identical
# feat_info <- read_gfeat_dir("/proj/mnhallqlab/users/michael/mmclock_pe/mmclock_nov2021/feat_l3/L1m-pe/L2m-l2_l2c-overall/L3m-age_sex/FEAT_l1c-EV_pe.gfeat")
gpa <- readRDS("/proj/mnhallqlab/users/michael/mmclock_entropy/mmclock_nov2021/mmclock_nov2021.rds")

# get feat dirs for echange
echange_info <- gpa$l1_model_setup$fsl %>% filter(l1_model == "echange")
echange_res4d <- file.path(echange_info$feat_dir, "stats", "res4d.nii.gz")
echange_mask <- file.path(echange_info$feat_dir, "mask.nii.gz")
table(file.exists(echange_res4d))
table(file.exists(echange_mask))

clustsim_acf <- afni_3dclustsim$new(
  fwhmx_input_files = echange_res4d, fwhmx_mask_files = echange_mask,
  scheduler = "slurm", out_dir = "/proj/mnhallqlab/users/michael/mmclock_entropy/mmclock_nov2021/feat_l3/L1m-echange",
  clustsim_mask = "/proj/mnhallqlab/lab_resources/standard/mni_icbm152_nlin_asym_09c/mni_icbm152_t1_tal_nlin_asym_09c_mask_2.3mm.nii", ncpus = 8, prefix = "echange"
)
clustsim_acf$is_complete()
#clustsim_acf$submit()
clustsim_acf$get_clustsim_df() %>% filter(nn == 1 & sided == "bi" & pthr == .001 & athr == .05)


# echange cluster with ACF option

# cobj$run() # execute 3dClusterize (use force = TRUE to force re-estimation)

cobj <- clustsim_acf$apply_clustsim(
  statistic_nifti = file.path(
    "/proj/mnhallqlab/users/michael/mmclock_entropy/mmclock_nov2021/feat_l3",
    "L1m-echange/L2m-l2_l2c-overall/L3m-age_sex/FEAT_l1c-EV_entropy_change_feedback.gfeat",
    "cope1.feat/stats/zstat6.nii.gz"
  ), NN = 1, sided = "bi", athr = .05, pthr = .001, output_thresholded_image = TRUE
)

cobj$generate_subclusters(break_nvox = 1000, min_subclust_nvox = 30, max_subclust_nvox = 1800, print_progress = TRUE)
cobj$get_clust_df()

# overlaps with mask
odat <- cobj$get_clust_df() %>%
  select(roi_num, subroi_num, Volume, CM.LR, CM.PA, CM.IS, overlap) %>%
  unnest(overlap)



# lookup of every parcel in Schaefer 200
#atlas <- "/proj/mnhallqlab/projects/clock_analysis/fmri/pfc_entropy/masks/Schaefer_DorsAttn_2.3mm.nii.gz"


wami_obj <- afni_whereami$new(
  coord_file = "/proj/mnhallqlab/projects/clock_analysis/fmri/pfc_entropy/original_masks/schaefer200_com_coords.txt",
  coord_orientation = "RAI", coord_space = "MNI", coord_file_columns = 0:2
)

wami_obj$run()
wami_obj$get_whereami_df()
