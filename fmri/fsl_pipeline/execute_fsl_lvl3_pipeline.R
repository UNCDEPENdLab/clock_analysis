# This script sets up the .fsf files to run a group analysis
## FSL Feat Level 2 analysis -- that is, fixed effects combinations of runs.

#load the master configuration file
to_run <- Sys.getenv("fsl_pipeline_file")

run_model_index <- as.numeric(Sys.getenv("run_model_index")) #which variant to execute
if (nchar(to_run) == 0L) { stop("Cannot locate environment variable fsl_pipeline_file") }
if (!file.exists(to_run)) { stop("Cannot locate configuration file", to_run) }
if (is.na(run_model_index)) { stop("Couldn't identify usable run_model_index variable.") }

load(to_run)

library(tidyverse)
library(dependlab)

#verify that mr_dir is present as expected
subinfo <- fsl_model_arguments$subject_covariates
feat_run_outdir <- fsl_model_arguments$outdir[run_model_index] #the name of the subfolder for the current run-level model
feat_lvl3_outdir <- file.path(fsl_model_arguments$group_output_dir, feat_run_outdir) #output directory for this run-level model
ncopes <- fsl_model_arguments$n_l1_copes[run_model_index] #number of l1 copes determines number of FEAT LVL3 analyses to run (1 per LVL1 cope)

subinfo$dir_found <- file.exists(subinfo$mr_dir)

cat("The following subjects were in the covariate file, but not the processed MRI data\n")
print(subset(subinfo, dir_found==FALSE))

dir.create(feat_lvl3_outdir, showWarnings=FALSE, recursive=TRUE)
setwd(feat_lvl3_outdir)

#cope structure for preconvolve models
#1 = clock_onset
#2 = feedback_onset
#3 = regressor of interest (in single-param models)

feat_lvl2_dirname <- "FEAT_LVL2_runtrend.gfeat" #should populate this to the structure at some point
models <- fsl_model_arguments$group_model_variants #different covariate models for the current run-level model (run_model_index)

##rework using subinfo structure as the authoritative guide (rather than repeated searches)
copedf <- c()
for (s in 1:nrow(subinfo)) {
  for (cope in 1:ncopes) {
    expectdir <- file.path(subinfo[s,"mr_dir"], fsl_model_arguments$expectdir, feat_run_outdir, feat_lvl2_dirname, paste0("cope", cope, ".feat"))
    if (dir.exists(expectdir)) {
      copedf <- rbind(copedf, data.frame(ID=subinfo[s,"ID"], model=feat_run_outdir, cope=cope, fsldir=expectdir))
    } else {
      message("could not find expected directory: ", expectdir)
    }
  }
}

mdf <- merge(subinfo, copedf, by="ID", all.y=TRUE)
badids <- c(11335, #low IQ, ADHD Hx, loss of consciousness
            11332, #should be excluded, but scan was terminated early due to repeated movement
            11282, #RTs at the floor for essentially all runs. Not appropriate
            11246, #huge movement and RTs at floor
            #10637, #large and many movements in later runs (need to revisit to confirm) ### OCT2018: 6 of 8 runs pass our algorithmic thresholds for motion
            10662  #I think there are reconstruction problems here -- need to revisit
            ) 

mdf <- mdf %>% filter(!ID %in% badids)
mdf <- arrange(mdf, ID, model, cope)

##fsl constructs models by cope
bycope <- lapply(split(mdf, mdf$cope), droplevels)

#loop over group-level models, setup the design matrix and spawn a FSL Level 3 job
l3template <- readLines(file.path(getMainDir(), "clock_analysis", "fmri", "fsf_templates", "feat_lvl3_sceptic_template.fsf"))

#loop over copes and group models, setting up .fsf files for each combination
for (cope in 1:length(bycope)) {
  if (is.null(bycope[[cope]]$Intercept)) { bycope[[cope]]$Intercept <- 1 } #add the column of ones

  #cope-level subfolder
  model_output_dir <- file.path(feat_lvl3_outdir, paste0("cope", cope))
  dir.create(model_output_dir, showWarnings=FALSE)
  
  for (this_model in models) {

    model_df <- bycope[[cope]]
    fsf_syntax <- l3template #copy shared ingredients
    fsf_syntax <- gsub(".OUTPUTDIR.", file.path(model_output_dir, paste(this_model, collapse="-")), fsf_syntax, fixed=TRUE)
    if (!"Intercept" %in% this_model) { this_model <- c("Intercept", this_model) } #at present, force an intercept column
    
    if (fsl_model_arguments$center_l3_predictors) {
      for (p in this_model) {
        if (p != "Intercept" && is.numeric(model_df[[p]])) {
          model_df[[p]] <- model_df[[p]] - mean(model_df[[p]], na.rm=TRUE)
        }
      }
    }

    model_df$dummy_ <- rnorm(nrow(model_df))
    mform <- as.formula(paste("dummy_ ~ -1 + ", paste(this_model, collapse=" + ")))
    fit_lm <- lm(mform, model_df)
    dmat <- model.matrix(fit_lm) #eventually allow interactions and so on??
    
    #add design matrix
    fsf_syntax <- c(fsf_syntax, generate_fsf_ev_syntax(inputs=model_df$fsldir, dmat=dmat))

    #generate diagonal contrast matrix, one per EV
    cmat <- diag(length(this_model))
    rownames(cmat) <- this_model #just name the contrasts after the EVs themselves

    fsf_syntax <- c(fsf_syntax, generate_fsf_contrast_syntax(cmat))
    
    #write the FSF to file
    out_fsf <- file.path(model_output_dir, paste0(paste(this_model, collapse="-"), ".fsf"))
    writeLines(fsf_syntax, con=out_fsf)

    #run the L3 analysis in parallel using qsub
    qsub_file(script=file.path(fsl_model_arguments$pipeline_home, "qsub_feat_lvl3.bash"), env_variables=c(torun=out_fsf))
    
  }

}
