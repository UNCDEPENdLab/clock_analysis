# runs "within-trial-invariant" mixed-effects Cox models 
# where response hazard is predicted by lagged deconvolved signal

# paths and packages ----
library(tidyverse)
library(lme4)
library(survival)
library(coxme)
library(survminer)
library(ggpubr)
library(broom)
library(broom.mixed)
# devtools::install_github('junkka/ehahelper') # requires gfortran, $ brew cask install gfortran
library(car)
library(foreach)
library(doParallel)
library(viridis)

reprocess = F # if running for the first time or need to reprocess MEDUSA data
if (reprocess) {censor_clock_post_rt = T # include ITI and clock epochs only, exclude post-RT epoch
}
# alignment = c("clock", "rt") # alignment of within-trial deconvolved signal
# alignment = c("rt")
alignment = c("clock")
uncensored = F # include first 1s and last 0.5s in survival models
decompose_within_between_trial = T # whether value and uncertainty are decomposed into within- and between-trial components


# basedir <- "~/Data_Analysis"
basedir <- "~/code"
coxme_dir <- file.path(basedir, "clock_analysis/coxme")
medusa_dir = "/Volumes/GoogleDrive/My Drive/SCEPTIC_fMRI/dan_medusa"
cache_dir = "/Volumes/GoogleDrive/My Drive/SCEPTIC_fMRI/dan_medusa/cache"
repo_directory <- file.path(basedir,"clock_analysis")

use_lagged_decons = c(F)
for (lagged_decon in use_lagged_decons) {
  # lagged_decon = T means use previous trial's signal to predict current RT
  if (lagged_decon) {
    rt_plot_dir <- "~/OneDrive/collected_letters/papers/sceptic_fmri/dan/plots/coxme/next_rt_rt_aligned"
    clock_plot_dir <- "~/OneDrive/collected_letters/papers/sceptic_fmri/dan/plots/coxme/next_rt_clock_aligned"
  } else {rt_plot_dir <- "~/OneDrive/collected_letters/papers/sceptic_fmri/dan/plots/coxme/this_rt_rt_aligned"
  clock_plot_dir <- "~/OneDrive/collected_letters/papers/sceptic_fmri/dan/plots/coxme/this_rt_clock_aligned"}
  cwd <- getwd()
  setwd(coxme_dir)
  
  # load survival and wide MEDUSA ----
  if (reprocess) {
    source(file.path(coxme_dir, "clock_coxme_prep_medusa_wtrial_inv.R"))
  } else {
    if (lagged_decon) {
      load(file.path(cache_dir, "fMRI_coxme_objects_with_wtrial_inv_medusa_lagged_Dec29_2020.RData"))  
    } else{
      load(file.path(cache_dir, "fMRI_coxme_objects_with_wtrial_inv_medusa_UNlagged_Dec29_2020.RData"))  
    }
  }
  
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
  
  
  # trial-invariant MEDUSA coxme ----
  # drop "_lag" suffix (kept in .RData to avoid confusion)
  if (lagged_decon) {
    if (uncensored) {
      clock_surv <- clock_lag_bb %>%
        rename_at(vars(contains("_R_") |  contains("_L_")), 
                  ~ str_remove(., "_lag"))
      rt_surv <- rt_lag_bb %>%
        rename_at(vars(contains("_R_") |  contains("_L_")), 
                  ~ str_remove(., "_lag"))
    } else {
      clock_surv <- clock_lag_fbb %>%
        rename_at(vars(contains("_R_") |  contains("_L_")), 
                  ~ str_remove(., "_lag"))
      rt_surv <- rt_lag_fbb %>%
        rename_at(vars(contains("_R_") |  contains("_L_")), 
                  ~ str_remove(., "_lag"))
    }
  } else {
    if (uncensored) {
      clock_surv <- clock_bb 
      rt_surv <- rt_bb 
    } else {
      clock_surv <- clock_fbb 
      rt_surv <- rt_fbb 
    }
    
  }
  # if (rt_tplus1) {
  #   rt_surv %>% mutate(evt_time = evt_time-1)
  # }
  rt_labels <- names(rt_surv[grepl("_L_|_R_", names(rt_surv))])
  clock_labels <- names(clock_surv[grepl("_L_|_R_", names(clock_surv))])
  
  # #test
  # labels <- labels[1:2]
  
  for (event in alignment) {
    if (event == 'rt') {surv_df <- rt_surv
    labels <- rt_labels} else {surv_df <- clock_surv
    labels <- clock_labels}
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
                      m  <- coxme(Surv(t1,t2,response) ~ omission_lag*wvs1b1a1*h + omission_lag2*wvs2b1a1*h +
                                    value_wi_t*h + uncertainty_wi_t*h +
                                    value_b_t*h + uncertainty_b_t*h +
                                    trial_neg_inv_sc*h +
                                    (1|ID), surv_df)} else {
                                      # simplified model:
                                      m  <- coxme(Surv(t1,t2,response) ~ omission_lag*wvs1b1a1*h + omission_lag2*wvs2b1a1*h + 
                                                    value_wi_t*h + uncertainty_wi_t*h + 
                                                    (1|ID), surv_df)
                                    }
                    stats <- as_tibble(insight::get_statistic(m))
                    stats$p <- 2*(1-pnorm(stats$Statistic))
                    stats$label <- label
                    stats$side <- substr(as.character(str_match(label, "_R_|_L_")),2,2)
                    stats$t <- as.numeric(gsub(".*_", "\\1", label))
                    suffix <- paste0("_", stats$side[1], "_", stats$t[1])
                    stats$region <- str_remove(label, suffix)
                    # newlist[[label]]<-stats
                    stats}
    df <- as_tibble(df)               
    # ddf$label <- as.factor(stringr::str_remove(pattern = paste0("_", gsub(".*_", "\\1", ddf$label), ""), ddf$label))
    terms <- unique(df$Parameter) 
    neural <- str_detect(terms,"h")
    terms <- terms[neural]
    # FDR ----
    df <- df  %>% group_by(Parameter) %>% mutate(p_fdr = p.adjust(p, method = 'fdr'),
                                                 p_level_fdr = as.factor(case_when(
                                                   # p_fdr > .1 ~ '0',
                                                   # p_fdr < .1 & p_fdr > .05 ~ '1',
                                                   p_fdr > .05 ~ '1',
                                                   p_fdr < .05 & p_fdr > .01 ~ '2',
                                                   p_fdr < .01 & p_fdr > .001 ~ '3',
                                                   p_fdr <.001 ~ '4')),
    ) %>% ungroup() 
    df$p_level_fdr <- factor(df$p_level_fdr, levels = c('1', '2', '3', '4'), labels = c("NS","p < .05", "p < .01", "p < .001"))
    df$`p, FDR-corrected` = df$p_level_fdr
    # plots ----
    if (event == "clock") {
      setwd(clock_plot_dir)  
      epoch_label = "Time relative to clock onset, seconds"
    } else if (event == "rt") {
      setwd(rt_plot_dir)
      epoch_label = "Time relative to outcome, seconds"
    }
    for (fe in terms) {
      edf <- df %>% filter(Parameter == paste(fe)) 
      termstr <- str_replace_all(fe, "[^[:alnum:]]", "_")
      if (decompose_within_between_trial) {
        pdf(paste(termstr, "_decomposed.pdf", sep = ""), width = 11, height = 6)
      } else {pdf(paste(termstr, ".pdf", sep = ""), width = 11, height = 6)}
      print(ggplot(edf, aes(t, region)) + geom_tile(aes(fill = Statistic, alpha = `p, FDR-corrected`), size = 1) +  
              geom_vline(xintercept = 0, lty = "dashed", color = "#FF0000", size = 2) + facet_wrap(~side) +
              scale_fill_viridis(option = "plasma") + scale_color_grey() + xlab(epoch_label) + ylab("Parcel") + 
              labs(alpha = expression(italic(p)[FDR])) + ggtitle(paste(termstr)))
      dev.off()
      ## save ----
      # save output for inspection
      save(file = "medusa_coxme_wtrial_inv_output.Rdata", df)} 
  }
}

