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
library(readxl)
source("~/code/Rhelpers/theme_black.R")

from_cache = F # uses previously saved coxme output ddfs
plots = F
verbose = F # prints coxme model output for each parcel
inspect = F # inspect model statistics for sanity checks
beta_dir <- "~/OneDrive - University of Pittsburgh/Documents/SCEPTIC_fMRI/wholebrain_betas"
source("~/code/fmri.pipeline/R/mixed_by.R")

studies = c("meg", "fmri")
# studies = c("meg", "fmri") # whether to include meg session replication

# sensitivity analyses for fMRI sample, recommend running first 2, one at a time:
censor = c(T)            # if "T", sensitivity analysis: censor first 500 ms and last second
decompose = c(F)            # if "F" sensitivity analyss: don't decompose U and V into within- vs between-trial components
omit_value = F              # do not change, needed to check once; not necessary ex-post
split_by_reward = F         # if "T", sensitivity analysis including interaction with last reward/omission; not necessary ex-post
signals = c("entropychange")
# signals = c("abspe", "echange")
ranefs = "rint"  # "rslope", "rint", testing random slopes for non-decomposed for now

# loop over model versions as defined above

if (omit_value) {rhs = "_no_value"} else {rhs = ""}
for (censor_ends in censor) {
  for (study in studies) {
    for (decompose_within_between_trial in decompose) {
      if (!from_cache) {
        if (Sys.getenv("USER")=="alexdombrovski") {
          load("~/code/clock_analysis/coxme/fMRI_MEG_coxme_objects_no_MEDUSA_Nov23_2020")
          setwd("~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/betas_final/")
          # source("../get_trial_data.R")
          source("parcel_brain_behavior_functions_alex.R")
        } else if (Sys.getenv("USER")=="Alex") {
          load("~/code/clock_analysis/coxme/fMRI_MEG_coxme_objects_Nov22_2020")
          setwd("~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/betas_final/")
          # source("../get_trial_data.R")
          source("parcel_brain_behavior_functions_alex.R")
        } else {
          setwd("/Users/hallquist/Data_Analysis/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/betas_final")
          source("../get_trial_data.R")
          source("/Users/hallquist/Data_Analysis/r_packages/fmri.pipeline/R/mixed_by.R")
          source("../medusa_final/plot_medusa.R")
          analysis_dir <- "~/Data_Analysis/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/betas_final"
        }
        analysis_dir <- "~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/betas_final"
        setwd(analysis_dir)
        
        
        
        labels_df <- setDT(read_excel("~/code/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/MNH DAN Labels 400 Good Only 47 parcels.xlsx")) %>%
          mutate(roi_num7 = as.factor(roi7_400), 
                 mask_value = as.integer(roi7_400),
                 plot_label = mnh_label_400, 
                 vm_gradient17 = parcel_group) %>% select(roi_num7, mask_value, plot_label, vm_gradient17, network17_400_DAN, hemi, x, y, z)
        
        
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
        rm(list = c("bb", "fbb", "mbb", "mfbb"))
        trial_df <- trial_df %>% select(ID, run, trial, rewFunc, t2, t1, run_trial, omission_lag, bin, response, trial_neg_inv_sc,
                                        value_wi, uncertainty_wi, value_wi_t, uncertainty_wi_t, value_b_t, uncertainty_b_t, value_b, uncertainty_b) 
        
        abspe_betas <- fread(file.path(paste0(beta_dir, "/L1m-abspe_plus_rew/Schaefer_444_final_2009c_2.3mm_cope_l2.csv.gz"))) %>%
          filter(l2_cope_name == "overall" & !l1_cope_name  %in% c("EV_clock", "EV_feedback", "EV_rew_om")) %>% # only parametric modulators
          dplyr::select(-feat_dir, -img, -mask_name, -session, -l1_cope_number, -l2_cope_number, -l2_model) %>%
          rename(fmri_beta = value) %>%
          # merge(label_df, by = label_join_col, all.x = TRUE)
          merge(labels_df, by = "mask_value", all = FALSE)
        echange_betas <- fread(file.path(paste0(beta_dir, "/L1m-echange/Schaefer_444_final_2009c_2.3mm_cope_l2.csv.gz"))) %>%
          filter(l2_cope_name == "overall" & !l1_cope_name  %in% c("EV_clock", "EV_feedback")) %>% # only parametric modulators
          dplyr::select(-feat_dir, -img, -mask_name, -session, -l1_cope_number, -l2_cope_number, -l2_model) %>%
          rename(fmri_beta = value) %>%
          # merge(label_df, by = label_join_col, all.x = TRUE)
          merge(labels_df, by = "mask_value", all = FALSE)
        if ("abspe" %in% signals & "entropychange" %in% signals) {
          betas <- rbind(abspe_betas, echange_betas)} else if ("abspe" %in% signals & !"entropychange" %in% signals) {betas <- abspe_betas
          } else if ((!"abspe" %in% signals) & "entropychange" %in% signals) {betas <- echange_betas}
        
        betas <-  betas %>% select(mask_value, id, l1_cope_name, fmri_beta) %>% 
          mutate(l1_cope_name = str_remove(l1_cope_name, "EV_")) %>% mutate(l1_cope_name = str_remove(l1_cope_name, "_feedback")) %>%
          pivot_wider(names_from = c(l1_cope_name, mask_value), values_from = fmri_beta) %>% rename(ID = "id")
        labels <- names(betas %>% select(!ID))
        
        surv_df <- inner_join(trial_df, betas, by = "ID")
        
        # parallel loop over parcels and signals ("labels")
        
        
        # make cluster ----
        f <- Sys.getenv('PBS_NODEFILE')
        library(parallel)
        ncores <- detectCores() - 1
        nodelist <- if (nzchar(f)) readLines(f) else rep('localhost', ncores)
        
        cat("Node list allocated to this job\n")
        print(nodelist)
        
        cl <- makePSOCKcluster(nodelist, outfile='')
        print(cl) ##; print(unclass(cl))
        registerDoParallel(cl)
        
        # save.image(file="parcel_input_snapshot_coxme.RData")
        print(paste("Running ", study, " decompose = ", decompose_within_between_trial, " censor = ", censor_ends))
        
        df <- foreach(i = 1:length(labels), .packages=c("lme4", "tidyverse", "broom", "coxme", "car"),
                      .combine='rbind') %dopar% {
                        label <- labels[[i]]
                        # print(paste(i," ", "Processing parcel", label))
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
                          if (omit_value) {
                            m  <- coxme(Surv(t1,t2,response) ~  uncertainty_wi_t * trial_neg_inv_sc * scale(h) +
                                          (1|ID), surv_df)} else {
                                            m  <- coxme(Surv(t1,t2,response) ~ value_wi_t * scale(h) + uncertainty_wi_t * trial_neg_inv_sc * scale(h) +
                                                          (1|ID), surv_df)} 
                        } else {
                          if (omit_value) {
                            m  <- coxme(Surv(t1,t2,response) ~  uncertainty_wi * trial_neg_inv_sc * scale(h) + 
                                          (1|ID), surv_df)} else { if ("rslope" %in% ranefs){
                                            m  <- coxme(Surv(t1,t2,response) ~ value_wi + uncertainty_wi + trial_neg_inv_sc + value_wi:scale(h) + uncertainty_wi:trial_neg_inv_sc + 
                                                          uncertainty_wi : scale(h) + uncertainty_wi : trial_neg_inv_sc : scale(h) + 
                                                          (1 + value_wi|ID), surv_df)
                                          } else {
                                            # simplified model:
                                            m  <- coxme(Surv(t1,t2,response) ~ value_wi * scale(h) + uncertainty_wi * trial_neg_inv_sc * scale(h) + 
                                                          (1|ID), surv_df)}}
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
                        if (verbose) {print(str(stats))}
                        print(paste(i," ", "Processed parcel", label))
                        stats}
        df <- as_tibble(df)               
        stopCluster(cl)
        beepr::beep(sound = 2)
        
        ddf <- df
        
        ddf <- ddf %>% mutate(mask_value = parse_number(label),
                              # save all terms, not only interactions, for validation
                              fmri_beta = gsub("[::0-9::,_]","", label)) %>% rename(term = "Parameter") %>% #filter(str_detect(term,"h")) %>%
          inner_join(labels_df, by = "mask_value")
        # str(ddf)
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
        ddf$study <- study
        ddf$method <- decompose_within_between_trial
        ddf$censor <- censor_ends
        ddf <- merge(ddf, labels_df)
        print(str(ddf))
        ############ Save
        # setwd(out_dir)
        
        if (decompose_within_between_trial) {method = "decomposed"
        value_term = "value_wi_t:scale(h)"
        uncertainty_term = "scale(h):uncertainty_wi_t"
        } else {method = "non_decomposed"
        value_term = "value_wi:scale(h)"
        uncertainty_term = "scale(h):uncertainty_wi"
        }
        if (!censor_ends) {
          out_dir <- file.path(beta_dir, "coxme")
        } else if (censor_ends) {
          out_dir <- file.path(beta_dir, "coxme/censored")}
        setwd(out_dir)
        saveRDS(ddf, paste0("beta_coxme_", method, "_", study, rhs, ranefs, "_trial.rds"))
        
      } else if (from_cache) {
        if (!censor_ends) {
          out_dir <- file.path(beta_dir, "coxme")
        } else if (censor_ends) {
          out_dir <- file.path(beta_dir, "coxme/censored")}
        
        if (decompose_within_between_trial) {method = "decomposed"
        value_term = "value_wi_t:scale(h)"
        uncertainty_term = "scale(h):uncertainty_wi_t"
        } else {method = "non_decomposed"
        value_term = "value_wi:scale(h)"
        uncertainty_term = "scale(h):uncertainty_wi"
        }
        setwd(out_dir)
        ddf <- readRDS(file.path(paste0("beta_coxme_", method, "_", study, rhs, ranefs, "_trial.rds")))
      }
      ############ Plot
      if (plots) {
        for (beta in signals) {
          
          pdf(paste0("value_axial_", beta, "_", method, "_", study, "_", ranefs, ".pdf"), height = 12, width = 16)
          print(ggplot(ddf %>% filter(fmri_beta == beta & term == value_term),  aes(y, x, shape = network17_400_DAN)) + 
                  geom_point(size = 24, color = "white", aes(alpha = p_level_fdr, fill = Statistic)) + scale_shape_manual(values = 21:22) + geom_point(size = 24, aes(color = network17_400_DAN)) +
                  scale_fill_viridis(option = "inferno") + theme_black() + geom_text_repel(aes(label=plot_label), color = "white", force = 5, point.padding = 40, size = 8))
          dev.off()
          
          pdf(paste0("uncertainty_axial_", beta, "_", method, "_", study, "_", ranefs, ".pdf"), height = 12, width = 16)
          print(ggplot(ddf %>% filter(fmri_beta == beta & term == uncertainty_term),  aes(y, x, shape = network17_400_DAN)) + 
                  geom_point(size = 24, color = "white", aes(alpha = p_level_fdr, fill = Statistic)) + scale_shape_manual(values = 21:22) + geom_point(size = 24, aes(color = network17_400_DAN)) +
                  scale_fill_viridis(option = "turbo") + theme_black() + geom_text_repel(aes(label=plot_label), color = "white", force = 5, point.padding = 40, size = 8))
          dev.off()
          
          pdf(paste0("value_saggital_",beta, "_",  method, "_", study, "_", ranefs, ".pdf"), height = 12, width = 16)
          print(ggplot(ddf %>% filter(fmri_beta == beta & term == value_term ),  aes(y, z, shape = network17_400_DAN)) + 
                  geom_point(size = 32, color = "white", aes(alpha = p_level_fdr, 
                                                             fill = Statistic)) + 
                  scale_shape_manual(values = 21:22) + geom_point(size = 32, color = "white") +
                  scale_fill_viridis(option = "inferno") + theme_black() + geom_text_repel(aes(label=plot_label),  color="blue"))
          dev.off()
          
          pdf(paste0("uncertainty_saggital_",beta, "_",  method, "_", study, "_", ranefs, ".pdf"), height = 12, width = 16)
          print(ggplot(ddf %>% filter(fmri_beta == beta & term == uncertainty_term),  aes(y, z, shape = network17_400_DAN)) + 
                  geom_point(size = 32, color = "white", aes(alpha = p_level_fdr, 
                                                             fill = Statistic)) + 
                  scale_shape_manual(values = 21:22) + geom_point(size = 32, color = "white") +
                  scale_fill_viridis(option = "turbo") + theme_black() + geom_text_repel(aes(label=plot_label),  color="blue") )
          dev.off()
          
          pdf(paste0("value_vioin_", beta, "_", method, "_", study, "_", ranefs, ".pdf"), height = 8, width = 5)
          print(ggplot(ddf %>% filter(fmri_beta == beta & term == value_term )) + 
                  geom_jitter(size = 12, width = .1, height = 0,  aes(network17_400_DAN, Statistic, color = Statistic, alpha = p_level_fdr)) + 
                  geom_violin(aes(network17_400_DAN, Statistic), alpha = .2) + scale_shape_manual(values = 21:22) +
                  scale_color_viridis(option = "inferno") + theme_black() + 
                  geom_text_repel(aes(network17_400_DAN, Statistic, alpha = p_level_fdr, label=plot_label),  color="blue"))
          dev.off()
          
          pdf(paste0("value_violin_by_vm_gradient_", beta, "_", method, "_", study, "_", ranefs, ".pdf"), height = 8, width = 10)
          print(ggplot(ddf %>% filter(fmri_beta == beta & term == value_term )) + 
                  geom_jitter(size = 12, width = .1, height = 0,  aes(vm_gradient17, Statistic, color = Statistic, alpha = p_level_fdr)) + 
                  geom_violin(aes(vm_gradient17, Statistic), alpha = .2) + scale_shape_manual(values = 21:22) +
                  scale_color_viridis(option = "inferno") + theme_black() + 
                  geom_text_repel(aes(vm_gradient17, Statistic, alpha = p_level_fdr, label=plot_label),  color="blue"))
          dev.off()
          
          
          pdf(paste0("uncertainty_violin_", beta, "_", method, "_", study, "_", ranefs, ".pdf"), height = 8, width = 5)
          print(ggplot(ddf %>% filter(fmri_beta == beta & term == uncertainty_term)) + 
                  geom_jitter(size = 12, width = .1, height = 0,  aes(network17_400_DAN, Statistic, color = Statistic, alpha = p_level_fdr)) + 
                  geom_violin(aes(network17_400_DAN, Statistic), alpha = .2) + scale_shape_manual(values = 21:22) +
                  scale_color_viridis(option = "turbo") + theme_black() + 
                  geom_text_repel(aes(network17_400_DAN, Statistic, alpha = p_level_fdr, label=plot_label),  color="blue"))
          dev.off()
          
          pdf(paste0("uncertainty_violin_by_vm_gradient_", beta, "_", method, "_", study, "_", ranefs, ".pdf"), height = 8, width = 10)
          print(ggplot(ddf %>% filter(fmri_beta == beta & term == uncertainty_term)) + 
                  geom_jitter(size = 12, width = .1, height = 0,  aes(vm_gradient17, Statistic, color = Statistic, alpha = p_level_fdr)) + 
                  geom_violin(aes(vm_gradient17, Statistic), alpha = .2) + scale_shape_manual(values = 21:22) +
                  scale_color_viridis(option = "turbo") + theme_black() + 
                  geom_text_repel(aes(vm_gradient17, Statistic, alpha = p_level_fdr, label=plot_label),  color="blue"))
          dev.off()
          
          pdf(paste0("uncertaintyBYtrial_violin_", beta, "_", method, "_", study, "_", ranefs, ".pdf"), height = 8, width = 5)
          print(ggplot(ddf %>% filter(fmri_beta == beta & term == "scale(h):uncertainty_wi_t:trial_neg_inv_sc")) + 
                  geom_jitter(size = 12, width = .1, height = 0,  aes(network17_400_DAN, Statistic, color = Statistic, alpha = p_level_fdr)) + 
                  geom_violin(aes(network17_400_DAN, Statistic), alpha = .2) + scale_shape_manual(values = 21:22) +
                  scale_color_viridis(option = "turbo") + theme_black() + 
                  geom_text_repel(aes(network17_400_DAN, Statistic, alpha = p_level_fdr, label=plot_label),  color="blue"))
          dev.off()
          
          pdf(paste0("uncertaintyBYtrial_violin_by_vm_gradient_", beta, "_", method, "_", study, "_", ranefs, ".pdf"), height = 8, width = 10)
          print(ggplot(ddf %>% filter(fmri_beta == beta & term == "scale(h):uncertainty_wi_t:trial_neg_inv_sc")) + 
                  geom_jitter(size = 12, width = .1, height = 0,  aes(vm_gradient17, Statistic, color = Statistic, alpha = p_level_fdr)) + 
                  geom_violin(aes(vm_gradient17, Statistic), alpha = .2) + scale_shape_manual(values = 21:22) +
                  scale_color_viridis(option = "turbo") + theme_black() + 
                  geom_text_repel(aes(vm_gradient17, Statistic, alpha = p_level_fdr, label=plot_label),  color="blue"))
          dev.off()
        }
      }
    }
  }
}


if (inspect) {
  
  summary(lm(uncertainty_wi_t ~ value_wi_t * trial_neg_inv_sc * omission_lag, trial_df))
  summary(lm(uncertainty_wi_t ~ value_wi_t * trial_neg_inv_sc * omission_lag, trial_df))
  summary(lm(uncertainty_wi ~ value_wi + trial_neg_inv_sc, trial_df))
  summary(lm(uncertainty_wi ~ value_wi, trial_df))
  
  m  <- coxme(Surv(t1,t2,response) ~ value_wi_t + uncertainty_wi_t * trial_neg_inv_sc + uncertainty_wi_t * omission_lag +
                (1|ID), surv_df)
  summary(m)
  
  m2  <- coxme(Surv(t1,t2,response) ~ value_wi_t * value_b_t + uncertainty_wi_t * uncertainty_b_t + trial_neg_inv_sc + uncertainty_wi_t * omission_lag +
                 (1|ID), surv_df)
  summary(m2)
  
  m0  <- coxme(Surv(t1,t2,response) ~  uncertainty_wi_t * trial_neg_inv_sc + uncertainty_wi_t * omission_lag +
                 (1|ID), surv_df)
  summary(m0)
  # inspect models to understand uncertainty sensitivity
  # look at the key "decomposed" models for fMRI and MEG
  ddf_fmri <- read_rds("/Volumes/GoogleDrive/My Drive/SCEPTIC_fMRI/wholebrain_betas/coxme/beta_coxme_decomposed_fmritrial.rds") %>% mutate(
    sample = "fmri",
    model_flavor = "decomposed"
  )
  ddf_meg <- read_rds("/Volumes/GoogleDrive/My Drive/SCEPTIC_fMRI/wholebrain_betas/coxme/beta_coxme_decomposed_megtrial.rds") %>% mutate(
    sample = "meg",
    model_flavor = "decomposed"
  )
  
  library(VIM)
  aggr_plot <- aggr(ddf %>% filter(term == "scale(h):uncertainty_wi_t:trial_neg_inv_sc"), col=c('navyblue','red'), numbers=TRUE, sortVars=TRUE, labels=names(data), cex.axis=.7, gap=3, ylab=c("Histogram of missing data","Pattern"))  
  
  # aggr_plot <- aggr(ddf_fmri, col=c('navyblue','red'), numbers=TRUE, sortVars=TRUE, labels=names(data), cex.axis=.7, gap=3, ylab=c("Histogram of missing data","Pattern"))
  # aggr_plot <- aggr(ddf_meg , col=c('navyblue','red'), numbers=TRUE, sortVars=TRUE, labels=names(data), cex.axis=.7, gap=3, ylab=c("Histogram of missing data","Pattern"))
  # 
  # 
  # ddf <- rbind(ddf_fmri, ddf_meg)
  # 
  toplot <- 
    ddf %>% filter(fmri_beta == "entropychange" & term == "uncertainty_wi_t")
  histogram(toplot$Statistic)
  ggplot() + 
    geom_jitter(size = 12, width = .1, height = 0,  aes(vm_gradient17, Statistic, color = Statistic, alpha = p_level_fdr)) + 
    geom_violin(aes(vm_gradient17, Statistic), alpha = .2) + scale_shape_manual(values = 21:22) +
    scale_color_viridis(option = "turbo") + theme_black() + 
    geom_text_repel(aes(vm_gradient17, Statistic, alpha = p_level_fdr, label=plot_label),  color="blue") #+ facet_wrap(~sample)
  
  ddf_meg %>% filter(fmri_beta == "abspe" & term == "h:uncertainty_wi_t") %>% View()
  
  m <- coxme(Surv(t1,t2,response) ~ value_wi_t  + uncertainty_wi_t + (1|ID), bb)
  summary(m)
  
  mt <- coxme(Surv(t1,t2,response) ~ value_wi_t  + uncertainty_wi_t * trial_neg_inv_sc + (1|ID), bb)
  summary(mt)
  
  
  mm <- coxme(Surv(t1,t2,response) ~ value_wi_t + uncertainty_wi_t + (1|ID), mbb)
  summary(mm)
  
  mmt <- coxme(Surv(t1,t2,response) ~ value_wi_t  + uncertainty_wi_t * trial_neg_inv_sc + (1|ID), mbb)
  summary(mmt)
  
}