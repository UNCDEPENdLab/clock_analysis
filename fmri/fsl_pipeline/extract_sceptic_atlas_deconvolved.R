#this script is not dependent on run_model_index because it reads raw data at l1
Sys.setenv(fsl_pipeline_file="/gpfs/group/mnh5174/default/clock_analysis/fmri/fsl_pipeline/configuration_files/MMClock_aroma_preconvolve_fse_groupfixed.RData")
Sys.setenv(run_model_index=1)

to_run <- Sys.getenv("fsl_pipeline_file")

run_model_index <- as.numeric(Sys.getenv("run_model_index")) #which variant to execute
if (nchar(to_run) == 0L) { stop("Cannot locate environment variable fsl_pipeline_file") }
if (!file.exists(to_run)) { stop("Cannot locate configuration file", to_run) }
if (is.na(run_model_index)) { stop("Couldn't identify usable run_model_index variable.") }

load(to_run)

source(file.path(fsl_model_arguments$pipeline_home, "functions", "glm_helper_functions.R"))
source(file.path(fsl_model_arguments$pipeline_home, "functions", "deconvolve_funcs.R"))
source(file.path(fsl_model_arguments$pipeline_home, "functions", "spm_funcs.R"))

library(tidyverse)
library(abind)
library(oro.nifti)
library(reshape2)
library(dependlab)
library(oro.nifti)
library(parallel)
library(foreach)
library(doParallel)

#verify that mr_dir is present as expected
subinfo <- fsl_model_arguments$subject_covariates
feat_run_outdir <- fsl_model_arguments$outdir[run_model_index] #the name of the subfolder for the current run-level model
feat_lvl3_outdir <- file.path(fsl_model_arguments$group_output_dir, feat_run_outdir) #output directory for this run-level model

#used for reading l1 data
load(file.path(fsl_model_arguments$pipeline_home, "configuration_files", paste0(paste(fsl_model_arguments$analysis_name, feat_run_outdir, "lvl2_inputs", sep="_"), ".RData")))

#registerDoSEQ()
cl <- makeCluster(18) #hard code for now
registerDoParallel(cl)
clusterExport(cl, c("sigmoid", "spm_hrf", "generate_feature", "dsigmoid", "deconvolve_nlreg")) #make sure functions are available

subinfo$dir_found <- file.exists(subinfo$mr_dir)

atlas_files <- c("/gpfs/group/mnh5174/default/clock_analysis/fmri/hippo_voxelwise/long_axis_l_2.3mm.nii.gz",
  "/gpfs/group/mnh5174/default/clock_analysis/fmri/hippo_voxelwise/long_axis_r_2.3mm.nii.gz")

atlas_imgs <- lapply(atlas_files, readNIfTI, reorient=FALSE)

out_dir <- "/gpfs/group/mnh5174/default/clock_analysis/fmri/hippo_voxelwise/deconvolved_timeseries"

#use this for 1st-level results
l1_inputs <- feat_l2_inputs_df$feat_dir

TR <- 1.0 #seconds

#determine nifti files for atlas
l1_niftis <- sapply(l1_inputs, function(x) {
  fsf <- readLines(file.path(x, "design.fsf"))
  nifti <- grep("^set feat_files\\(1\\)", fsf, perl=TRUE, value=TRUE)
  stopifnot(length(nifti)==1L)
  nifti <- paste0(sub("set feat_files\\(1\\) \"([^\"]+)\"", "\\1", nifti, perl=TRUE), ".nii.gz")
  return(nifti)
})

#Bush deconvolution settings
nev_lr <- .01 #neural events learning rate (default in algorithm)
epsilon <- .005 #convergence criterion (default)
kernel <- spm_hrf(TR)$hrf #canonical SPM difference of gammas

#for testing
#l1_niftis <- l1_niftis[1:5]

#loop over atlas files
for (ai in 1:length(atlas_files)) {
  cat("Working on atlas: ", atlas_files[ai], "\n")
  aimg <- atlas_imgs[[ai]]
  a_indices <- which(aimg != 0, arr.ind=TRUE)
  a_coordinates <- cbind(a_indices, t(apply(a_indices, 1, function(r) { translateCoordinate(i=r, nim=aimg, verbose=FALSE) })))
  a_coordinates <- as.data.frame(a_coordinates) %>% setNames(c("i", "j", "k", "x", "y", "z")) %>%
    mutate(vnum=1:n(), atlas_value=aimg[a_indices], atlas_name=basename(atlas_files[ai])) %>%
    mutate_at(vars(x,y,z), round, 2) %>% select(vnum, atlas_value, everything())

  #loop over niftis
  for (si in 1:length(l1_niftis)) {
    this_subj <- feat_l2_inputs_df %>% select(subid, run_num, contingency, emotion, drop_volumes) %>% dplyr::slice(si)
    out_name <- file.path(out_dir, with(this_subj[1,,drop=F], paste0(subid, "_run", run_num, "_", basename(atlas_files[ai]), "_deconvolved.csv.gz")))

    if (file.exists(out_name)) {
      message("Deconvolved file already exists: ", out_name)
      next
    }
    
    cat("  Deconvolving subject: ", l1_niftis[si], "\n")
    dump_out <- tempfile()
    runAFNICommand(paste0("3dmaskdump -mask ", atlas_files[ai], " -o ", dump_out, " ", l1_niftis[si]))
    ts_out <- data.table::fread(dump_out) #read time series
    #ijk <- ts_out[, 1:3] + 1 #keep in 1-based indexing in outputs (AFNI uses zero-based)
    #names(ijk) <- c("i", "j", "k")
    #ijk$vnum <- 1:nrow(ijk)

    #to_deconvolve is a voxels x time matrix
    to_deconvolve <- as.matrix(ts_out[, -1:-3]) #remove ijk
    to_deconvolve <- t(apply(to_deconvolve, 1, scale)) #need to unit normalize for algorithm not to choke on near-constant 100-normed data

    deconv_mat <- foreach(vox_ts=iter(to_deconvolve, by="row"), .combine="rbind") %dopar% {
      reg <- tryCatch(deconvolve_nlreg(as.vector(vox_ts), kernel=kernel, nev_lr=nev_lr, epsilon=epsilon),
        error=function(e) { cat("Problem deconvolving: ", l1_niftis[si], "\n", file="deconvolve_errors", append=TRUE); return(rep(NA, length(vox_ts))) })
      return(reg)
    }
    
    #deconv_mat <- t(parApply(cl=cl, X=to_deconvolve, MARGIN=1, FUN=function(vox_ts) { #transpose to maintain voxels x time
    #  deconvolve_nlreg(vox_ts, kernel=kernel, nev_lr=nev_lr, epsilon=epsilon)
    #  #res <- tryCatch(deconvolve_nlreg(vox_ts, kernel=kernel, nev_lr=nev_lr, epsilon=epsilon),
    #  #  error=function(e) { browser(); return(rep(NA, length(ts))) }) #return NA in event of failure
    #  #return(res)
    #}))

    ## deconv_mat <- t(apply(X=to_deconvolve, 1, FUN=function(ts) { #transpose to maintain voxels x time
    ##   browser()
    ##   tryCatch(deconvolve_nlreg(ts, kernel=kernel, nev_lr=nev_lr, epsilon=epsilon),
    ##     error=function(e) { print(e); return(rep(NA, length(ts))) }) #return NA in event of failure
    ## }))

    #melt this for combination
    deconv_melt <- reshape2::melt(deconv_mat, value.name="decon", varnames=c("vnum", "time"))
    to_deconvolve_melt <- reshape2::melt(to_deconvolve, value.name="BOLD_z", varnames=c("vnum", "time"))
    
    deconv_df <- deconv_melt %>% mutate(vnum=as.numeric(vnum)) %>% left_join(a_coordinates, by="vnum") %>%
      mutate(nifti=sub("/gpfs/group/mnh5174/default/MMClock/MR_Proc/", "", l1_niftis[si], fixed=TRUE)) %>%
      cbind(feat_l2_inputs_df %>% select(subid, run_num, contingency, emotion, drop_volumes) %>% dplyr::slice(si)) %>% #add metadata
      mutate(time=time+drop_volumes) %>% select(-drop_volumes, -nifti) %>% #dropping nifti for now to save on file size
      select(subid, run_num, contingency, emotion, time, atlas_name, atlas_value, vnum, x, y, z, decon) #omitting i,j,k for now

    orig_df <- to_deconvolve_melt %>% mutate(vnum=as.numeric(vnum)) %>% left_join(a_coordinates, by="vnum") %>%
      mutate(nifti=sub("/gpfs/group/mnh5174/default/MMClock/MR_Proc/", "", l1_niftis[si], fixed=TRUE)) %>%
      cbind(feat_l2_inputs_df %>% select(subid, run_num, contingency, emotion, drop_volumes) %>% dplyr::slice(si)) %>% #add metadata
      mutate(time=time+drop_volumes) %>% select(-drop_volumes, -nifti) %>% #dropping nifti for now to save on file size
      select(subid, run_num, contingency, emotion, time, atlas_name, atlas_value, vnum, x, y, z, BOLD_z)

    out_name <- file.path(out_dir, with(deconv_df[1,,drop=F], paste0(subid, "_run", run_num, "_", basename(atlas_name), "_deconvolved.csv.gz")))
    write_csv(deconv_df, path=out_name)

    out_name <- file.path(out_dir, with(orig_df[1,,drop=F], paste0(subid, "_run", run_num, "_", basename(atlas_name), "_original_bold.csv.gz")))
    write_csv(orig_df, path=out_name)
    
  }
}
