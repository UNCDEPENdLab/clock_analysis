#this script is intended to check the correlation of convolved PE and entropy regressors to respond to a reviewer request to check the robustness of
#fmri analyses to simultaneous inclusion of predictors
library(dplyr)
library(lme4)
library(emmeans)
library(brms)

setwd("/Users/michael/ics/clock_analysis/fmri")
useCache=TRUE
if (useCache) {
  pe_e_corr_df <- readRDS("pe_e_correlation_dataframe.rds")
} else {
  load("/gpfs/group/mnh5174/default/clock_analysis/fmri/fsl_pipeline/configuration_files/MMClock_aroma_preconvolve_fse_groupfixed.RData")

  #simult model is currently position 3
  odir <- "sceptic-clock-feedback-v_entropy-pe_max-preconvolve_fse_groupfixed"
  dmats <- file.path(fsl_model_arguments$subject_covariates$mr_dir,
    fsl_model_arguments$expectdir, odir, "designmatrix.RData")

  pe_e_corr <- lapply(dmats, function(subj_d) {
    load(subj_d) #will load a number of objects, with d being the most important (the design matrix)
    pec <- sapply(d$collin_convolve, function(run) {
      run$r["v_entropy", "pe_max"]
    })

    ret <- data.frame(id=subj_data$id[1], pe_entropy_corr=pec, run=names(pec), stringsAsFactors=FALSE)
    rm(d, subj_data) #just to be very safe across iterations since we're using load...
    return(ret)
  })

  pe_e_corr_df <- bind_rows(pe_e_corr) #assemble as a single data.frame
  saveRDS(pe_e_corr_df, file="pe_e_correlation_dataframe.rds") #save to disk
}

hist(pe_e_corr_df$pe_entropy_corr)

m1 <- lmer(pe_entropy_corr ~ run + (1 | id), pe_e_corr_df)
emmeans(m1, ~run)
emmeans(m1, ~1)

#allow heterogeneity in l1 variance by run
m2 <- brm(bf(pe_entropy_corr ~ run + (1 | id),
  sigma ~ run), data = pe_e_corr_df, chains=4, cores=4, iter=3000)

summary(m2)
emmeans(m2, ~1)
emg <- emmeans(m2, ~run)

#default ROPEs on coefficients
bayestestR::rope(m2)

#on the correlation estimates themselves, by run (run2 is a little funky)
bayestestR::rope(emg)

emg_overall <- emmeans(m2, ~1)
bayestestR::rope(emg_overall, range=c(-.1, .1)) #test on |r| < .1

#homogeneity (no run-level variation in L1 resid variance)
m3 <- brm(pe_entropy_corr ~ run + (1 | id), pe_e_corr_df, chains=4, cores=4, iter=3000)
emmeans(m3, ~1)
emmeans(m3, ~run)

