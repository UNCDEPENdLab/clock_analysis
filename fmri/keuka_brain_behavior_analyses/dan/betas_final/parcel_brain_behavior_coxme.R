# final parcel-wise analyses of DAN brain-to-behavior
library(data.table)
library(tidyverse)
library(afex)
library(lattice)
library(emmeans)
library(fmri.pipeline) # has mixed_By
library(readr)
library(purrr)
library(glue)
library(coxme)
library(viridis)
library(ggnewscale)
library(RColorBrewer)
library(ggrepel)
library(foreach)
library(doParallel)
source("~/code/Rhelpers/theme_black.R")

from_cache = T
beta_dir <- "/Volumes/GoogleDrive/My Drive/SCEPTIC_fMRI/wholebrain_betas"
source("~/code/fmri.pipeline/R/mixed_by.R")

for (study in c("meg", "fmri")) {
  for (censor_ends in c(F, T)) {
    for (decompose_within_between_trial in c(F, T)) {
      if (!from_cache) {
      if (Sys.getenv("USER")=="alexdombrovski") {
        load("~/code/clock_analysis/coxme/fMRI_MEG_coxme_objects_no_MEDUSA_Nov23_2020")
      } else {
        setwd("/Users/hallquist/Data_Analysis/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/betas_final")
        
        source("../get_trial_data.R")
        source("/Users/hallquist/Data_Analysis/r_packages/fmri.pipeline/R/mixed_by.R")
        source("../medusa_final/plot_medusa.R")
        labels <- readxl::read_excel("/Users/hallquist/Data_Analysis/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/MNH Dan Labels.xlsx") %>%
          dplyr::rename(mask_value=roinum) %>% select(mask_value, plot_label)
        trial_df <- get_trial_data(repo_directory = "/Users/hallquist/Data_Analysis/clock_analysis")
        
      }
      #analysis_dir <- "~/Data_Analysis/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/betas_final"
      analysis_dir <- "~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/betas_final"
      setwd(analysis_dir)
      
      source("../get_trial_data.R")
      source("parcel_brain_behavior_functions.R")
      # BALSA parcel labels from whereami
      setwd("~/code/schaefer_wb_parcellation")
      schaefer_7 <- read.csv("labels/Schaefer2018_400Parcels_7Networks_order.csv") %>%
        mutate(network=factor(network), net_num = as.numeric(network)) %>%
        rename(network7=network, net_num7=net_num)
      
      # this has the spatial coordinate, spatial_roi_num
      schaefer_7_lookup <- read.csv("labels/Schaefer_400_7networks_labels.csv")
      
      schaefer_7 <- schaefer_7 %>% inner_join(schaefer_7_lookup, by="roi_num") %>%
        rename(roi_num7=roi_num, subregion7=subregion)
      
      schaefer_17 <- read.csv("labels/Schaefer2018_400Parcels_17Networks_order.csv") %>%
        mutate(network=factor(network), net_num = as.numeric(network)) %>%
        rename(network17=network, net_num17=net_num) %>%
        select(-hemi) # mirrored in 7
      
      # this has the spatial coordinate, spatial_roi_num
      schaefer_17_lookup <- read.csv("labels/Schaefer_400_17networks_labels.csv") %>%
        select(roi_num, spatial_roi_num) # x,y,z and labels already duplicated in 7-network lookup
      
      schaefer_17 <- schaefer_17 %>% inner_join(schaefer_17_lookup, by="roi_num") %>%
        rename(roi_num17=roi_num, subregion17=subregion)
      
      both <- inner_join(schaefer_7, schaefer_17, by="spatial_roi_num") %>%
        select(spatial_roi_num, roi_num7, roi_num17, network7, network17, net_num7, net_num17, subregion7, subregion17, everything())
      setDT(both)
      labels_df <- both %>% 
        filter(net_num7==3 | (network17=="DorsAttnA" | network17=="DorsAttnB")) %>%
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
          mask_value = as.integer(as.character(roi_num7))
        ) %>% select(roi_num7, mask_value, vm_gradient17, hemi)
      
      
      #trial_df <- get_trial_data(repo_directory = "~/Data_Analysis/clock_analysis") %>%
      # trial_df <- get_trial_data(repo_directory = "~/code/clock_analysis") %>%
      if (study == "fmri" & !censor_ends) {
        trial_df <- bb
        out_dir <- file.path(beta_dir, "coxme")
      } else if (study == "fmri" & censor_ends) {
        trial_df <- fbb
        out_dir <- file.path(beta_dir, "coxme/censored")
      } else if (study == "meg" & !censor_ends){
        trial_df <- mbb
        out_dir <- file.path(beta_dir, "coxme")
      } else if (study == "meg" & censor_ends) {
        trial_df <- mfbb
        out_dir <- file.path(beta_dir, "coxme/censored")
      }
      
      trial_df <- trial_df %>% select(ID, run, trial, rewFunc, t2, t1, run_trial, omission_lag, bin, response, 
                                      value_wi, uncertainty_wi, value_wi_t, uncertainty_wi_t, value_b_t, uncertainty_b_t, value_b, uncertainty_b) 
      
      abspe_betas <- fread("/Volumes/GoogleDrive/My Drive/SCEPTIC_fMRI/wholebrain_betas/L1m-abspe_plus_rew/Schaefer_444_final_2009c_2.3mm_cope_l2.csv.gz") %>%
        filter(l2_cope_name == "overall" & !l1_cope_name  %in% c("EV_clock", "EV_feedback", "EV_rew_om")) %>% # only parametric modulators
        dplyr::select(-feat_dir, -img, -mask_name, -session, -l1_cope_number, -l2_cope_number, -l2_model) %>%
        rename(fmri_beta = value) %>%
        # merge(label_df, by = label_join_col, all.x = TRUE)
        merge(labels_df, by = "mask_value", all = FALSE)
      echange_betas <- fread("/Volumes/GoogleDrive/My Drive/SCEPTIC_fMRI/wholebrain_betas/L1m-echange/Schaefer_444_final_2009c_2.3mm_cope_l2.csv.gz") %>%
        filter(l2_cope_name == "overall" & !l1_cope_name  %in% c("EV_clock", "EV_feedback")) %>% # only parametric modulators
        dplyr::select(-feat_dir, -img, -mask_name, -session, -l1_cope_number, -l2_cope_number, -l2_model) %>%
        rename(fmri_beta = value) %>%
        # merge(label_df, by = label_join_col, all.x = TRUE)
        merge(labels_df, by = "mask_value", all = FALSE)
      
      betas <- rbind(abspe_betas, echange_betas) %>% select(mask_value, id, l1_cope_name, fmri_beta) %>% 
        mutate(l1_cope_name = str_remove(l1_cope_name, "EV_")) %>% mutate(l1_cope_name = str_remove(l1_cope_name, "_feedback")) %>%
        pivot_wider(names_from = c(l1_cope_name, mask_value), values_from = fmri_beta) %>% rename(ID = "id")
      labels <- names(betas %>% select(!ID))
      
      surv_df <- inner_join(trial_df, betas, by = "ID")
      # make cluster ----
      f <- Sys.getenv('PBS_NODEFILE')
      library(parallel)
      ncores <- detectCores()
      nodelist <- if (nzchar(f)) readLines(f) else rep('localhost', ncores)
      
      cat("Node list allocated to this job\n")
      print(nodelist)
      
      cl <- makePSOCKcluster(nodelist, outfile='')
      print(cl) ##; print(unclass(cl))
      registerDoParallel(cl)
      
      save.image(file="parcel_input_snapshot_coxme.RData")
      
      df <- foreach(i = 1:length(labels), .packages=c("lme4", "tidyverse", "broom", "coxme", "car"), 
                    .combine='rbind') %dopar% {
                      label <- labels[[i]]
                      # for (label in labels) {print(paste("Processing parcel", label))
                      # for (side in c("l", "r")) {
                      # for (t in -1:10) {
                      surv_df$h <- surv_df[[label]]
                      # form <- as.formula(paste0("Surv(t1,t2,response) ~ wvs1b1a1*", label, " + wvs2b1a1*", label, " + wvs3b1a1*", label, 
                      # " +  value_wi*", label," + uncertainty_wi*", label, " + (1|ID)"))
                      # alternative with within- vs. between-trials decomposition
                      if (decompose_within_between_trial) {
                        # this "full" model uses between- and within-trial predictors simultaneously, but 
                        # note that the interaction of h (decon) with between-trial uncertainty or value
                        # only means that they will respond faster when those are high
                        m  <- coxme(Surv(t1,t2,response) ~ value_wi_t*h + uncertainty_wi_t*h +
                                      (1|ID), surv_df)} else {
                                        # simplified model:
                                        m  <- coxme(Surv(t1,t2,response) ~ value_wi*h + uncertainty_wi*h + 
                                                      (1|ID), surv_df)
                                      }
                      stats <- as_tibble(insight::get_statistic(m))
                      # stats$p <- 2*(1-pnorm(stats$Statistic))
                      stats[3] <- insight::get_parameters(m)[2]
                      stats$p <- signif(1 - pchisq((stats$Statistic)^2, 1), 2)
                      stats$label <- label
                      # stats$side <- substr(as.character(str_match(label, "_R_|_L_")),2,2)
                      # stats$t <- as.numeric(gsub(".*_", "\\1", label))
                      # suffix <- paste0("_", stats$side[1], "_", stats$t[1])
                      # stats$region <- str_remove(label, suffix)
                      # newlist[[label]]<-stats
                      stats}
      df <- as_tibble(df)               
      stopCluster(cl)
      beepr::beep(sound = 2)
      
      ddf <- df
      
      ddf <- ddf %>% mutate(mask_value = parse_number(label),
                            fmri_beta = gsub("[::0-9::,_]","", label)) %>% rename(term = "Parameter") %>% filter(str_detect(term,"h")) %>%
        inner_join(labels_df, by = "mask_value")
      str(ddf)
      # gsub("[^[:alnum:] ]", "", str
      
      ddf <- ddf %>% group_by(term) %>% mutate(roi_num7 = as.factor(roi_num7),
                                               padj_BY_term = p.adjust(p, method = 'BH'),
                                               p_level_fdr = as.factor(case_when(
                                                 # p_fdr > .1 ~ '0',
                                                 # p_fdr < .1 & p_fdr > .05 ~ '1',
                                                 padj_BY_term > .05 ~ '1',
                                                 padj_BY_term < .05 & padj_BY_term > .01 ~ '2',
                                                 padj_BY_term < .01 & padj_BY_term > .001 ~ '3',
                                                 padj_BY_term <.001 & padj_BY_term > .0001 ~ '4',
                                                 padj_BY_term <.0001 & padj_BY_term > .00001 ~ '5',
                                                 padj_BY_term <.00001 ~ '6'
                                               ))
      ) %>% ungroup() %>% mutate(
        statistic_untrimmed = Statistic,
        Statistic = psych::winsor(Statistic, trim = .01))
      # filter(net_num7 ==3 & (net_num17==7 | net_num17==8))
      # filter(net_num7 ==3 ) 
      ddf$p_level_fdr <- factor(ddf$p_level_fdr, levels = c('1', '2', '3', '4', '5', '6'), labels = c("NS","p < .05", "p < .01", "p < .001", "p < .0001", "p < .00001")) 
      
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
          mask_value = as.integer(as.character(roi_num7))
        )
      
      ddf <- merge(ddf, labels)
      
      ############ Plot
      setwd(out_dir)
      
      if (decompose_within_between_trial) {method = "decomposed"
      } else {method = "non_decomposed"}
      saveRDS(ddf, paste0("beta_coxme_", method, "_", study, ".rds"))
      
      } else if (from_cache) {
        if (!censor_ends) {
          out_dir <- file.path(beta_dir, "coxme")
        } else if (censor_ends) {
          out_dir <- file.path(beta_dir, "coxme/censored")}
          
        if (decompose_within_between_trial) {method = "decomposed"
        } else {method = "non_decomposed"}
        setwd(out_dir)
        ddf <- readRDS(file.path(paste0("beta_coxme_", method, "_", study, ".rds")))
      }
      for (beta in c("abspe", "entropychange")) {
        pdf(paste0("value_by_network_axial_", beta, "_", method, "_", study, ".pdf"), height = 12, width = 16)
        print(ggplot(ddf %>% filter(fmri_beta == beta & str_detect(term, "value")),  aes(y, x, shape = network17)) + 
                geom_point(size = 24, color = "white", aes(alpha = p_level_fdr, fill = Statistic)) + scale_shape_manual(values = 21:22) + geom_point(size = 24, aes(color = network17)) +
                scale_fill_viridis(option = "inferno") + theme_black() + geom_text_repel(aes(label=subregion17), color = "white", force = 5, point.padding = 40, size = 8))
        dev.off()
        
        
        pdf(paste0("uncertainty_by_network_axial_", beta, "_", method, "_", study, ".pdf"), height = 12, width = 16)
        print(ggplot(ddf %>% filter(fmri_beta == beta & str_detect(term, "uncertainty")),  aes(y, x, shape = network17)) + 
                geom_point(size = 24, color = "white", aes(alpha = p_level_fdr, fill = Statistic)) + scale_shape_manual(values = 21:22) + geom_point(size = 24, aes(color = network17)) +
                scale_fill_viridis(option = "turbo") + theme_black() + geom_text_repel(aes(label=subregion17), color = "white", force = 5, point.padding = 40, size = 8))
        dev.off()
        
        pdf(paste0("value_by_network_saggital_",beta, "_",  method, "_", study, ".pdf"), height = 12, width = 16)
        print(ggplot(ddf %>% filter(fmri_beta == beta & str_detect(term, "value")),  aes(y, z, shape = network17)) + 
                geom_point(size = 32, color = "white", aes(alpha = p_level_fdr, 
                                                           fill = Statistic)) + 
                scale_shape_manual(values = 21:22) + geom_point(size = 32, color = "white") +
                scale_fill_viridis(option = "inferno") + theme_black() + geom_text_repel(aes(label=subregion17),  color="blue"))
        dev.off()
        
        pdf(paste0("uncertainty_by_network_saggital_",beta, "_",  method, "_", study, ".pdf"), height = 12, width = 16)
        print(ggplot(ddf %>% filter(fmri_beta == beta & str_detect(term, "uncertainty")),  aes(y, z, shape = network17)) + 
                geom_point(size = 32, color = "white", aes(alpha = p_level_fdr, 
                                                           fill = Statistic)) + 
                scale_shape_manual(values = 21:22) + geom_point(size = 32, color = "white") +
                scale_fill_viridis(option = "turbo") + theme_black() + geom_text_repel(aes(label=subregion17),  color="blue") )
        dev.off()
        
        pdf(paste0("value_by_network_", beta, "_", method, "_", study, ".pdf"), height = 8, width = 5)
        print(ggplot(ddf %>% filter(fmri_beta == beta & str_detect(term, "value"))) + 
                geom_jitter(size = 12, width = .1, height = 0,  aes(network17, Statistic, color = Statistic, alpha = p_level_fdr)) + 
                geom_violin(aes(network17, Statistic), alpha = .2) + scale_shape_manual(values = 21:22) +
                scale_color_viridis(option = "inferno") + theme_black() + 
                geom_text_repel(aes(network17, Statistic, alpha = p_level_fdr, label=subregion17),  color="blue"))
        dev.off()
        
        pdf(paste0("value_by_vm_gradient_", beta, "_", method, "_", study, ".pdf"), height = 8, width = 10)
        print(ggplot(ddf %>% filter(fmri_beta == beta & str_detect(term, "value"))) + 
                geom_jitter(size = 12, width = .1, height = 0,  aes(vm_gradient17, Statistic, color = Statistic, alpha = p_level_fdr)) + 
                geom_violin(aes(vm_gradient17, Statistic), alpha = .2) + scale_shape_manual(values = 21:22) +
                scale_color_viridis(option = "inferno") + theme_black() + 
                geom_text_repel(aes(vm_gradient17, Statistic, alpha = p_level_fdr, label=subregion17),  color="blue"))
        dev.off()
        
                
        pdf(paste0("uncertainty_by_network_", beta, "_", method, "_", study, ".pdf"), height = 8, width = 5)
        print(ggplot(ddf %>% filter(fmri_beta == beta & str_detect(term, "uncertainty"))) + 
                geom_jitter(size = 12, width = .1, height = 0,  aes(network17, Statistic, color = Statistic, alpha = p_level_fdr)) + 
                geom_violin(aes(network17, Statistic), alpha = .2) + scale_shape_manual(values = 21:22) +
                scale_color_viridis(option = "turbo") + theme_black() + 
                geom_text_repel(aes(network17, Statistic, alpha = p_level_fdr, label=subregion17),  color="blue"))
        dev.off()
        
        pdf(paste0("uncertainty_by_vm_gradient_", beta, "_", method, "_", study, ".pdf"), height = 8, width = 10)
        print(ggplot(ddf %>% filter(fmri_beta == beta & str_detect(term, "uncertainty"))) + 
                geom_jitter(size = 12, width = .1, height = 0,  aes(vm_gradient17, Statistic, color = Statistic, alpha = p_level_fdr)) + 
                geom_violin(aes(vm_gradient17, Statistic), alpha = .2) + scale_shape_manual(values = 21:22) +
                scale_color_viridis(option = "turbo") + theme_black() + 
                geom_text_repel(aes(vm_gradient17, Statistic, alpha = p_level_fdr, label=subregion17),  color="blue"))
        dev.off()
      }
    }
  }
}
