library(tidyverse)
library(lme4)
library(data.table)
library(readxl)

out_dir <- "/Volumes/GoogleDrive/My Drive/SCEPTIC_fMRI/bsocial_medusa"
#repo_directory <- "~/code/clock_analysis"
repo_directory <- "~/Data_Analysis/clock_analysis"
setwd("~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/dan")
source("get_trial_data.R")
source("medusa_final/plot_medusa.R")
source("~/code/fmri.pipeline/R/mixed_by.R")
out_dir <- "/Volumes/GoogleDrive/My Drive/SCEPTIC_fMRI/bsocial_medusa"

emm = T # extract EMMEANS estimates, e.g. for hi/lo abs(PE)
reprocess = F
alignment <- "rt"

# helper function to compile list of formulae
named_list <- function(...) {
  vnames <- as.character(match.call())[-1]
  return(setNames(list(...), vnames))
}

if (!reprocess) {
  d <- readRDS("/Volumes/GoogleDrive/My Drive/SCEPTIC_fMRI/explore_medusa/rt_visuomotor_long_dan_only_200_trial_df.rds")
} else {
  
  rt_visuomotor_long <- fread("~/OneDrive - University of Pittsburgh/Documents/SCEPTIC_fMRI/bsocial_medusa/transformed_schaefer_dan_3.125mm_rt_long_decon_aligned.csv.gz") %>%
    mutate(atlas_value = as.character(atlas_value),
           id = as.character(id))
  dan_labels <- setDT(read_excel("~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/MNH Dan Labels.xlsx")) %>%
    mutate(atlas_value = as.character(roinum),
           parcel = word(MNHLabel, 2, sep = "_"),
           hemi = substr(MNHLabel, 1,1)
    ) %>% select(atlas_value, MNHLabel, parcel, Stream, Visuomotor_Gradient, lobe, parcel, hemi)
  
  # add Schaefer 17 labels
  setwd("~/code/schaefer_wb_parcellation")
  schaefer_7 <- read.csv("labels/Schaefer2018_200Parcels_7Networks_order.csv") %>%
    mutate(network=factor(network), net_num = as.numeric(network)) %>%
    rename(network7=network, net_num7=net_num)
  
  # this has the spatial coordinate, spatial_roi_num
  schaefer_7_lookup <- read.csv("labels/Schaefer_200_7networks_labels.csv")
  
  schaefer_7 <- schaefer_7 %>% inner_join(schaefer_7_lookup, by="roi_num") %>%
    rename(roi_num7=roi_num, subregion7=subregion)
  
  schaefer_17 <- read.csv("labels/Schaefer2018_200Parcels_17Networks_order.csv") %>%
    mutate(network=factor(network), net_num = as.numeric(network)) %>%
    rename(network17=network, net_num17=net_num) %>%
    select(-hemi) # mirrored in 7
  
  # this has the spatial coordinate, spatial_roi_num
  schaefer_17_lookup <- read.csv("labels/Schaefer_200_17networks_labels.csv") %>%
    select(roi_num, spatial_roi_num) # x,y,z and labels already duplicated in 7-network lookup
  
  schaefer_17 <- schaefer_17 %>% inner_join(schaefer_17_lookup, by="roi_num") %>%
    rename(roi_num17=roi_num, subregion17=subregion)
  
  both <- inner_join(schaefer_7, schaefer_17, by="spatial_roi_num") %>%
    select(spatial_roi_num, roi_num7, roi_num17, network7, network17, net_num7, net_num17, subregion7, subregion17, everything())
  setDT(both)
  labels <- both %>% filter(net_num7==3 & (network17=="DorsAttnA" | network17=="DorsAttnB")) %>% 
    mutate(roi_num7 = as.factor(roi_num7)) %>% 
    # label lobes
    mutate(lobe = case_when(
      str_detect(subregion17, "Temp") ~ "temporal",
      str_detect(subregion17, "Par") | str_detect(subregion17, "SPL") | str_detect(subregion17, "PostC") |
        str_detect(subregion17, "IPS") | str_detect(subregion17, "IPL") | str_detect(subregion17, "pCun") ~ "parietal",
      str_detect(subregion17, "PFC") | str_detect(subregion17, "FEF") | str_detect(subregion17, "PrCv") ~ "frontal"),
      vm_gradient17 = case_when(
        lobe == "temporal" ~ "MT+",
        lobe == "parietal" & network17 == "DorsAttnA" ~ "PPCcaudal",
        lobe == "parietal" & network17 == "DorsAttnB" ~ "PPCrostral",
        lobe == "frontal" ~ "premotor",
        TRUE ~ as.character(network17)),
      plot_label = sub("Focus point:\\s+", "", MNI_Glasser_HCP_v1.0, perl=TRUE),
      mask_value = as.integer(as.character(roi_num7)),
      atlas_value = as.character(mask_value))
  
  dan_labels <-  dan_labels %>% select(atlas_value, MNHLabel) %>% inner_join(labels, by = "atlas_value")
  
  rt_visuomotor_long <- rt_visuomotor_long %>% inner_join(labels, by = "atlas_value") #%>% filter(!is.na(Stream))
  
  trial_df <- setDT(get_trial_data(dataset = "bsocial", repo_directory = out_dir))  
  
  trial_df <- trial_df %>% dplyr::select(id, run, run_trial, trial, trial_neg_inv_sc, rt_csv_sc, rt_lag_sc, pe_max, rew_om_c, abs_pe_c, abspexrew,
                                         v_entropy_wi, v_entropy_wi_change, kld3, v_max_wi, abs_pe, outcome) %>%
    mutate(log_kld3 = log(kld3 + .00001),
           id = as.character(id))
  # rt_visuomotor_long <- rt_visuomotor_long %>% inner_join(trial_df)
  # rt_visuomotor_long <- rt_visuomotor_long %>% filter(evt_time > -3 & evt_time < 6)
  
  #alignment <- "clock"
  #alignment <- "clock_online"
  
  # setDT(trial_df)
  
  message("Merging")
  if (alignment=="clock") {
    # subset to columns of interest
    trial_df <- trial_df %>%
      dplyr::select(id, run, run_trial, trial_neg_inv_sc, rt_csv_sc, rt_lag_sc, pe_max_lag,
                    v_entropy_wi, v_entropy_wi_change_lag, kld3_lag, v_max_wi, abs_pe_lag, outcome_lag) %>%
      mutate(log_kld3_lag = log(kld3_lag + .00001))
    
    d <- merge(trial_df, clock_visuomotor_long, by = c("id", "run", "trial"))
    d <- d %>% tidyr::separate(visuomotor_side, into=c("vm_gradient17", "side"), sep="_")  
  } else if (alignment == "clock_online") {
    # subset to columns of interest
    trial_df <- trial_df %>%
      dplyr::select(id, run, run_trial, trial_neg_inv_sc, rt_csv_sc, rt_lag_sc, pe_max_lag,
                    v_entropy_wi, v_entropy_wi_change_lag, kld3_lag, v_max_wi, abs_pe_lag, outcome_lag) %>%
      mutate(log_kld3_lag = log(kld3_lag + .00001))
    
    d <- merge(trial_df, clock_visuomotor_long_online, by = c("id", "run", "trial"))
    d <- d %>% tidyr::separate(visuomotor_side, into=c("vm_gradient17", "side"), sep="_") 
  } else if (alignment == "rt") {
    # subset to columns of interest
    trial_df <- trial_df %>%
      dplyr::select(id, run, run_trial, trial, trial_neg_inv_sc, rt_csv_sc, rt_lag_sc, pe_max, rew_om_c, abs_pe_c, abspexrew,
                    v_entropy_wi, v_entropy_wi_change, kld3, v_max_wi, abs_pe, outcome) %>%
      mutate(log_kld3 = log(kld3 + .00001))
    
    d <- merge(trial_df, rt_visuomotor_long, by = c("id", "run", "trial"))
    d <- d %>% rename(side = "hemi")
    # d <- d %>% tidyr::separate(visuomotor_side, into=c("vm_gradient", "side"), sep="_")
    saveRDS(d, "/Volumes/GoogleDrive/My Drive/SCEPTIC_fMRI/explore_medusa/rt_visuomotor_long_dan_only_200_trial_df.rds")
    }
}
# add subject characteristics
setwd(out_dir)
sub_df <- setDT(read_csv("b_social_demographics_May25_2022.csv", col_types = cols(ID = col_character())) %>% rename(id = ID))

d <- inner_join(d, sub_df, by = "id")

if (reprocess) {rm(rt_visuomotor_long)}
# rm(clock_visuomotor_long)
# rm(clock_visuomotor_long_online)
gc()

# cor(trial_df$v_entropy_wi, trial_df$v_entropy_wi_change, use="pairwise")
# vv <- trial_df %>% select(v_entropy_wi, v_entropy_wi_change) %>%
#   mutate(el1=lag(v_entropy_wi, 1), el2=lag(v_entropy_wi, 2)) %>% na.omit()
# ccf(vv$v_entropy_wi, vv$v_entropy_wi_change)
#
# cor(vv$v_entropy_wi_change, vv$el1)
# cor(vv$v_entropy_wi_change, vv$el2)
# acf(vv$v_entropy_wi)


if (alignment == "clock" || alignment == "clock_online") {
  # basal analysis
  # vmax_wi: targeted analysis to demonstrate MT+ (vmax-positive) versus DAN (vmax-negative)
  enc_clock_base <- formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + v_max_wi +
                              v_entropy_wi + v_entropy_wi_change_lag + abs_pe_lag + outcome_lag +
                              (1 | id) )
  
  enc_clock_rslope <- formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + v_max_wi +
                                v_entropy_wi + v_entropy_wi_change_lag + abs_pe_lag + outcome_lag +
                                (abs_pe_lag + v_entropy_wi + v_entropy_wi_change_lag + v_max_wi | id) )
  
  
  # clock-aligned using kld3_lag
  enc_clock_kld <- formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + v_max_wi + 
                             v_entropy_wi + v_entropy_wi_change_lag + abs_pe_lag + outcome_lag + kld3_lag +
                             (1 | id) )
  
  # log version of kld given the strong skew
  enc_clock_logkld <- formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + v_max_wi + 
                                v_entropy_wi + v_entropy_wi_change_lag + abs_pe_lag + outcome_lag + log_kld3_lag +
                                (1 | id) )
  
  # signed PE
  enc_clock_pe <- formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + v_max_wi + 
                            v_entropy_wi + v_entropy_wi_change_lag + pe_max_lag +
                            (1 | id) )
  
  # abs_pe x outcome interaction
  enc_clock_int <- formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + v_max_wi +
                             v_entropy_wi + v_entropy_wi_change_lag + abs_pe_lag * outcome_lag +
                             (1 | id) )
  
  # centered reward and abs pe to compute interaction
  # enc_clock_int_cent <- formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + v_max_wi +
  #                            v_entropy_wi + v_entropy_wi_change_lag + abs_pe_c + rew_om_c + abspexrew +
  #                            (1 | id) )
  
  # trial interactions
  enc_clock_trial <- formula(~ v_max_wi*run_trial + v_entropy_wi*run_trial + v_entropy_wi_change_lag*run_trial + abs_pe_lag*run_trial +
                               rt_csv_sc*run_trial + rt_lag_sc + outcome_lag +
                               (1 | id) )
  
  flist <- named_list(enc_clock_base, enc_clock_rslope, enc_clock_kld, enc_clock_logkld, enc_clock_pe, enc_clock_int, enc_clock_trial)
} else if (alignment == "rt") {
  # basal analysis
  enc_rt_base <- formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + v_max_wi +
                           v_entropy_wi + v_entropy_wi_change + abs_pe + outcome +
                           (1 | id) )
  enc_rt_base_grp <- formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + v_max_wi*groupLeth +
                               v_entropy_wi + v_entropy_wi_change*groupLeth + abs_pe*groupLeth + outcome +
                               (1 | id) )
  enc_rt_base_grp_age_edu <- formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + v_entropy_wi + outcome +
                                       v_max_wi*groupLeth + v_entropy_wi_change*groupLeth + abs_pe*groupLeth + 
                                       v_max_wi*age + v_entropy_wi_change*age + abs_pe*age + 
                                       v_max_wi*registration_edu + v_entropy_wi_change*registration_edu + abs_pe*registration_edu + 
                                       (1 | id) )
  
  # vmax_wi: targeted analysis to demonstrate MT+ (vmax-positive) versus DAN (vmax-negative)
  enc_rt_rslope <- formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + v_max_wi +
                             v_entropy_wi + v_entropy_wi_change + abs_pe + outcome +
                             (abs_pe + v_entropy_wi + v_entropy_wi_change + v_max_wi | id) )
  
  enc_rt_kld <- formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + v_max_wi + 
                          v_entropy_wi + v_entropy_wi_change + abs_pe + outcome + kld3 +
                          (1 | id) )
  
  enc_rt_logkld <- formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + v_max_wi + 
                             v_entropy_wi + v_entropy_wi_change + abs_pe + outcome + log_kld3 +
                             (1 | id) )
  
  # signed PE
  enc_rt_pe <- formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + v_max_wi + 
                         v_entropy_wi + v_entropy_wi_change + pe_max +
                         (1 | id) )
  # signed PE w/o covariates
  enc_rt_pe_only <- formula(~ pe_max +
                              (1 | id) )
  
  # abs_pe x outcome interaction
  enc_rt_int <- formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + v_max_wi +
                          v_entropy_wi + v_entropy_wi_change + abs_pe * outcome +
                          (1 | id) )
  
  enc_rt_int_cent <- formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + v_max_wi +
                               v_entropy_wi + v_entropy_wi_change + abs_pe_c + rew_om_c + abspexrew +
                               (1 | id) )
  
  # centered outcome approach to get absPE at average of reward/omissions
  enc_rt_int <- formula(~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + v_max_wi +
                          v_entropy_wi + v_entropy_wi_change + abs_pe * outcome +
                          (1 | id) )
  
  
  # trial interactions
  enc_rt_trial <- formula(~ rt_csv_sc*run_trial + rt_lag_sc + v_max_wi*run_trial +
                            v_entropy_wi*run_trial + v_entropy_wi_change*run_trial + abs_pe*run_trial + outcome +
                            (1 | id) )
  flist <- named_list(enc_rt_base, enc_rt_pe, enc_rt_base_grp, enc_rt_base_grp_age_edu)
  # flist <- named_list(enc_rt_rslope, enc_rt_kld, enc_rt_pe)
  # flist <- named_list(enc_rt_base, enc_rt_rslope, enc_rt_kld, enc_rt_logkld, enc_rt_pe, enc_rt_int, enc_rt_int_cent, enc_rt_trial)
}


splits <- c("vm_gradient17", "side", "evt_time")

message("Running mixed_by")
ddf <- mixed_by(d, outcomes = "decon_mean", rhs_model_formulae = flist,
                split_on = splits, scale_predictors = c("abs_pe", "abs_pe_lag", "pe_max", "pe_max_lag", "run_trial", "age", "registration_edu"),
                tidy_args = list(effects = c("fixed", "ran_vals", "ran_pars", "ran_coefs"), conf.int = TRUE), 
                calculate = c("parameter_estimates_reml"), ncores = 10, refit_on_nonconvergence = 5, padjust_by = "term"#,
                # emtrends_spec = list(
                #   abspe = list(outcome = "decon_mean", model_name = "enc_rt_int", var = "abs_pe", specs = c("outcome"))
                #)
)

saveRDS(ddf, file=file.path(out_dir, paste0(alignment, "_", splits[1], "by_200_roi_GRP_demo", "_encode_medusa_fmri.rds")))

dd <- ddf$emtrends_list$abspe
dd <- dd[,c(-4, -5)]
ggplot(dd, aes(x=evt_time, y=abs_pe.trend, ymin=abs_pe.trend-std.error, ymax=abs_pe.trend+std.error, color=vm_gradient)) + 
  geom_pointrange() + geom_line() + geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
  facet_grid(outcome...1 ~ side)
