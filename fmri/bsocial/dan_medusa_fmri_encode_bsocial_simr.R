library(tidyverse)
library(lme4)
library(data.table)
library(readxl)
library(fmri.pipeline)
library(simr)

out_dir <- "/Volumes/GoogleDrive/My Drive/SCEPTIC_fMRI/bsocial_medusa/"

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
  # d <- readRDS("G:/.shortcut-targets-by-id/1bycmTNiqFyrMhY3WKbteNgv36PjFKFIc/explore_medusa/rt_visuomotor_long_dan_only_200_trial_df.rds")
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
gc()

test_df <- d %>% filter(vm_gradient17 == "PPCcaudal" & side == "R" & evt_time == "2")
test_df <- na.omit(test_df) %>%#remove NAs for powerSim to work
  mutate(v_entropy_wi_change_scaled = scale(v_entropy_wi_change))
#full model
# enc_rt_base_grp <- lmer(decon_mean ~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + v_max_wi*groupLeth +
#                             v_entropy_wi + v_entropy_wi_change*groupLeth + abs_pe*groupLeth + outcome +
#                             (1 | id), data = test_df, na.action = na.omit)

#simple model
enc_rt_base_grp <- lmer(decon_mean ~ v_entropy_wi_change_scaled*groupLeth +
                          (v_entropy_wi_change_scaled | id), data = test_df, na.action = na.omit)
# (1 | id), data = test_df, na.action = na.omit)
enc_rt_base_grp <- lmer(decon_mean ~ trial_neg_inv_sc + rt_csv_sc + rt_lag_sc + v_entropy_wi + outcome +
                          v_max_wi*groupLeth + v_entropy_wi_change_scaled*groupLeth + abs_pe*groupLeth + 
                          v_max_wi*age + v_entropy_wi_change_scaled*age + abs_pe*age + 
                          v_max_wi*registration_edu + v_entropy_wi_change_scaled*registration_edu + abs_pe*registration_edu + 
                          (1 | id), data = test_df, na.action = na.omit)
summary(enc_rt_base_grp)



sample_size <- c(500)

# sample_size <- c(152, 275)
sequence <- seq(from = .1, to = .25, by = .01) # adjust by to increase precision


#toy example to check effects of alpha
# 
effect_size <- .2
# alpha_uncorrected <- .00001
model <- extend(enc_rt_base_grp, along="id", n=sample_size)
fixef(model)["v_entropy_wi_change_scaled:groupLethBPD_LL"] <- effect_size
alpha_uncorrected <- .05/(2*10*4) # 2 hemispheres, 10 timepoints, 4 streams
# alpha_uncorrected <- 10^-6
s <-  as.data.frame(summary(powerSim(model, nsim = 20, alpha = alpha_uncorrected,#increase to 10000
                                     test = fixed("v_entropy_wi_change_scaled:groupLethBPD_LL", "t"))))
s

f <- Sys.getenv('PBS_NODEFILE')
library(parallel)
library(foreach)
library(doParallel)
ncores <- detectCores()
nodelist <- if (nzchar(f)) readLines(f) else rep('localhost', ncores)

cat("Node list allocated to this job\n")
print(nodelist)

cl <- makePSOCKcluster(nodelist, outfile='')
print(cl) ##; print(unclass(cl))
registerDoParallel(cl)

out_list <- list()
for (num in sample_size) {
  print(paste0("Sample Size:", num))
  model <- extend(enc_rt_base_grp, along="id", n=num)
  simlist <- list() #list of power sim results, no alteration to sample size
  # parallelize simulations
  df <- foreach(i = 1:length(sequence), .packages=c("lme4", "tidyverse", "simr"), 
                .combine='rbind') %dopar% {
                  effect_size <- sequence[[i]]
                  # for (i in sequence) {
                  # j <- match(i, sequence) #can't index a list on decimals apparently
                  # print("Effect size of v_entropy_wi_change:groupLethBPD_LL")
                  print(effect_size)
                  fixef(model)["v_entropy_wi_change_scaled:groupLethBPD_LL"] <- effect_size
                  # alpha_uncorrected <- .05/(2*10*4) # 2 hemispheres, 10 timepoints, 4 streams
                  alpha_uncorrected <- 10^-5 # 2 hemispheres, 10 timepoints, 4 streams
                  s <-  as.data.frame(summary(powerSim(model, nsim = 100, alpha = alpha_uncorrected,#increase to 10000
                                                       test = fixed("v_entropy_wi_change_scaled:groupLethBPD_LL", "t"))))
                  s$effect_size <- effect_size
                  s$sample_size <- num
                  s$alpha = alpha_uncorrected
                  s}
  df <- as_tibble(df)
  
  k <- match(num, sample_size)
  out_list[[k]] <- df
}
stopCluster(cl)
sim_df_test <- rbindlist(out_list)
fdf_test <- sim_df_test %>% filter( sample_size < 300) %>% rename(Power = "mean") %>% mutate(sample_size_f = as.factor(sample_size))
ggplot(fdf_test, aes(effect_size, Power, lty = sample_size_f)) + geom_line() + geom_point()
ggplot(fdf1, aes(effect_size, Power, lty = sample_size_f)) + geom_line() + geom_point()
fdf1 <- fdf1 %>% mutate(interaction_effect_size = effect_size/.01)
ggplot(fdf1, aes(interaction_effect_size, Power, lty = sample_size_f)) + geom_line() + geom_point()
ggplot(fdf, aes(effect_size, Power, lty = sample_size_f)) + geom_line() + geom_point()

fdf <- fdf %>% mutate(interaction_effect_size = effect_size/0.011235)
setwd('~/OneDrive/collected_letters/grants/bsocial_renewal_2022/research_plan/figs/')

pdf("fmri_power_curve_aim3.pdf", width = 4, height = 3)
ggplot(fdf %>% filter (interaction_effect_size > .35), aes(interaction_effect_size, Power, lty = sample_size_f)) + geom_line(size = 1) + geom_point(size = 3) + 
  ggtitle("Figure 8. Power to detect group differences \nin dorsal stream responses") +
  geom_text(data = fdf %>% filter (interaction_effect_size > .35 & interaction_effect_size < .65 & sample_size_f == 275), aes(label = round(Power, digits = 2)), nudge_x = -.041, nudge_y = .02) + 
  geom_text(data = fdf %>% filter (interaction_effect_size > .35 & interaction_effect_size < .65 & sample_size_f == 152), aes(label = round(Power, digits = 2)), nudge_x = -.04, nudge_y = .02, color = "gray45") + 
  xlab("Interaction effect,\nproportion of the main effect of behavioral variable") + ylab("Statistical power") + 
  # theme(legend.title = "Sample size", legend.text = element_text (labels = c("Current, 152", "Proposed, 275"))) +
  scale_linetype_manual(name="Sample size",
                        labels=c("Current,\n152","Proposed,\n275"), values = c(1, 2)) + 
  geom_hline(yintercept=0.8, color="gray60") +
  theme_minimal() +
  theme(legend.key.size =  unit(0.5, "in"), legend.position = c(.8, .35), legend.text = element_text(size = 12), legend.title = element_text(size = 12))
dev.off()
# beepr::beep(sound = 2)
#element 1-3 of out_list is power analyses at n = 152, 275, 300; within each sample size, analyses are run over effect sizes from -02 to .01, iteration by .01

#e.g., to print the power analysis for an estimate of .01 with 152 people, you would do: 
# out_list[[1]][[4]]
