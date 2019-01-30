# This script saves the significant clusters for each map in a SCEPTIC group analysis

#used for testing: group fixed entropy
#Sys.setenv(fsl_pipeline_file="/gpfs/group/mnh5174/default/clock_analysis/fmri/fsl_pipeline/configuration_files/MMClock_aroma_preconvolve_fse_groupfixed.RData")
#Sys.setenv(run_model_index=2)

#load the master configuration file

to_run <- Sys.getenv("fsl_pipeline_file")

run_model_index <- as.numeric(Sys.getenv("run_model_index")) #which variant to execute
if (nchar(to_run) == 0L) { stop("Cannot locate environment variable fsl_pipeline_file") }
if (!file.exists(to_run)) { stop("Cannot locate configuration file", to_run) }
if (is.na(run_model_index)) { stop("Couldn't identify usable run_model_index variable.") }

load(to_run)

#load("/gpfs/group/mnh5174/default/clock_analysis/fmri/fsl_pipeline/configuration_files/MMClock_aroma_preconvolve_fse_groupfixed.RData") #current arguments
#run_model_index <- 1

source(file.path(fsl_model_arguments$pipeline_home, "functions", "glm_helper_functions.R"))

#1) load spatial maps RData object
#2) rebuild into 4d cube (where fourth dimension is run/subject)
#3) clusterize each effect of interest in stats outputs using 3dclust and generating mask
#library(ggplot2)
library(tidyverse)
library(abind)
library(oro.nifti)
library(reshape2)
library(robust)
library(car)
library(dependlab)
library(oro.nifti)
library(parallel)
library(foreach)
library(doParallel)

#verify that mr_dir is present as expected
subinfo <- fsl_model_arguments$subject_covariates
feat_run_outdir <- fsl_model_arguments$outdir[run_model_index] #the name of the subfolder for the current run-level model
feat_lvl3_outdir <- file.path(fsl_model_arguments$group_output_dir, feat_run_outdir) #output directory for this run-level model
n_l1_copes <- fsl_model_arguments$n_l1_copes[run_model_index] #number of l1 copes determines number of FEAT LVL3 analyses to run (1 per LVL1 cope)
l1_cope_names <- fsl_model_arguments$l1_cope_names[[run_model_index]] #names of l1 copes (used for folder naming)
zthresh <- fsl_model_arguments$zthresh #3.09
clustsize <- fsl_model_arguments$clustsize #34

#registerDoSEQ()
cl <- makeCluster(fsl_model_arguments$n_cluster_beta_cpus)
registerDoParallel(cl)

subinfo$dir_found <- file.exists(subinfo$mr_dir)

feat_lvl2_dirname <- "FEAT_LVL2_runtrend.gfeat" #should populate this to the structure at some point
models <- sapply(fsl_model_arguments$group_model_variants, function(x) { paste(x, collapse="-") }) #different covariate models for the current run-level model (run_model_index)

for (l1 in 1:n_l1_copes) {
  l1_contrast_name <- l1_cope_names[l1]
  model_output_dir <- file.path(feat_lvl3_outdir, l1_contrast_name)

  #generate separate files for each l1 contrast (reset here)
  all_metadata <- list()
  all_rois <- list()
  all_beta_series <- list()
  
  for (this_model in models) {
    expect_gfeat <- file.path(model_output_dir, paste0(l1_contrast_name, "-", this_model, ".gfeat"))

    if (!file.exists(expect_gfeat)) {
      message("Could not locate expected .gfeat directory for group analysis: ", expect_gfeat)
      next
    }

    design_fsf <- file.path(expect_gfeat, "design.fsf") #expected design file
    stopifnot(file.exists(design_fsf))
    design_txt <- readLines(design_fsf)
    subject_inputs <- grep("^\\s*set feat_files\\(\\d+\\).*", design_txt, value=TRUE, perl=TRUE)
    subject_inputs <- sub("\\s*set feat_files\\(\\d+\\)\\s*\"?([^\"]+)\"?", "\\1", subject_inputs, perl=TRUE) #just keep the directory itself
    
    l3_ev_txt <- grep("\\s*set fmri\\(evtitle\\d+\\).*", design_txt, value=TRUE, perl=TRUE)
    l3_ev_names <- sub("\\s*set fmri\\(evtitle\\d+\\)\\s*\"?([^\"]+)\"?", "\\1", l3_ev_txt, perl=TRUE)
    ev_title_nums <- as.numeric(sub("\\s*set fmri\\(evtitle(\\d+)\\).*", "\\1", l3_ev_txt, perl=TRUE))
    n_l3_copes <- max(ev_title_nums) #NB. Should come back here and use the copes, not evs!
    l3_ev_names <- l3_ev_names[order(ev_title_nums)] #order l3 copes in ascending order to match l3 loop
    
    evs <- grep("\\s*fmri\\(evg[0-9.]+\\)", design_txt, value=TRUE, perl=TRUE)
    evnums <- as.numeric(sub("\\s*set fmri\\(evg[0-9]+\\.(\\d+)\\).*", "\\1", evs, perl=TRUE))
    subnums <- as.numeric(sub("\\s*set fmri\\(evg([0-9]+)\\.\\d+\\).*", "\\1", evs, perl=TRUE))
    ev_values <- as.numeric(sub("\\s*set fmri\\(evg[0-9]+\\.\\d+\\)\\s+([-\\d+.]+)", "\\1", evs, perl=TRUE))
    
    dmat <- matrix(NA_real_, nrow=max(subnums), ncol=max(evnums))
    dmat[cbind(subnums,evnums)] <- ev_values
    colnames(dmat) <- make.names(gsub("\"", "", l3_ev_names[order(ev_title_nums)], fixed=TRUE))

    design_df <- data.frame(dmat) %>% mutate(feat_input_id=1:n())

    #figure out what the l2 contrasts are and read relevant statistics for each
    l2fsf <- file.path(subject_inputs[1], "design.fsf")
    stopifnot(file.exists(l2fsf))
    l2_syntax <- readLines(l2fsf)

    l2_contrast_info <- grep("\\s*set fmri\\(conname_real\\.\\d+\\).*", l2_syntax, value=TRUE, perl=TRUE)
    l2_contrast_nums <- as.numeric(sub("\\s*set fmri\\(conname_real\\.(\\d+)\\).*", "\\1", l2_contrast_info, perl=TRUE))
    l2_contrast_names <- sub("\\s*set fmri\\(conname_real\\.\\d+\\)\\s*\"?([^\"]+)\"?.*", "\\1", l2_contrast_info, perl=TRUE)

    n_l2_contrasts <- max(l2_contrast_nums)
    l2_contrast_names <- l2_contrast_names[order(l2_contrast_nums)] #order l2 contrast names in ascending order to match l2 loop below

    #to get L1 betas (per run), we need to extract the inputs to the L2 analysis. These are embedded in the L2 FSF file
    #for () {
    #}

    #hard coding location of beta series analysis for now
    beta_series_inputs <- sub(paste0("^(", fsl_model_arguments$fmri_dir, "/", fsl_model_arguments$idregex, "/", fsl_model_arguments$expectdir,
      ")/.*"), "\\1/sceptic-clock_bs-feedback-preconvolve_fse_groupfixed", subject_inputs)
    
    #loop over l2 contrasts
    #l2_loop_outputs <- foreach(l2=iter(1:n_l2_contrasts), .packages=c("oro.nifti", "dplyr")) %do% {
    l2_loop_outputs <- list()
    for (l2 in 1:n_l2_contrasts) {
      l2_loop_cluster_metadata <- list()
      l2_loop_rois <- list()
      l2_loop_bs <- list()
      
      l2_contrast_name <- l2_contrast_names[l2] #current l2 contrast
      copefiles <- file.path(subject_inputs, "stats", paste0("cope", l2, ".nii.gz"))
      imgdims <- dim(oro.nifti::readNIfTI(copefiles[1], read_data=FALSE))

      #generate concatenated cope file image of l2 images (one per subject)
      copeconcat <- array(0, dim=c(imgdims, length(copefiles)))
      for (i in 1:length(copefiles)) { copeconcat[,,,i] <- readNIfTI(copefiles[i], reorient=FALSE)@.Data }

      #clusterize the current l1 cope for a given l3 covariate and a given l2 contrast
      
      for (l3 in 1:n_l3_copes) {
        l3_contrast_name <- l3_ev_names[l3] #current covariate
        groupmap <- file.path(expect_gfeat, paste0("cope", l2, ".feat"), "stats", paste0("zstat", l3, ".nii.gz")) #this is for numeric naming of .gfeat dirs
        
        #gdat <- readNIfTI(groupmap, reorient=FALSE)
        #generate cluster mask
        clust_1d <- paste0(tempfile(), "_tmpclust.1D")
        clust_brik <- paste0(tempfile(), "_tmpclust")
        runAFNICommand(paste0("3dclust -overwrite -1Dformat -nosum -1dindex 0 -1tindex 0",
          " -1thresh ", zthresh, " -dxyz=1 -savemask ", clust_brik, " 1.01 ", clustsize, " ", groupmap), 
          stdout=clust_1d)

        #get coordinates and names of regions
        lookup <- runAFNICommand(paste0("whereami -coord_file ", clust_1d, "'[1,2,3]' -space MNI -lpi -atlas CA_ML_18_MNIA"),
          stderr="/dev/null", intern=TRUE)
        
        exitstatus <- attr(lookup, "status")  
        if (!is.null(exitstatus) && exitstatus != 0) next #whereami failed, which occurs when there are no clusters. Skip to next tbrik

        #get voxel sizes of clusters
        vsizes <- read.table(clust_1d)$V1

        section_map <- grep("+++++++ nearby Atlas structures +++++++", lookup, fixed=TRUE)
        section_split <- rep.int(seq_along(section_map), times=diff(c(section_map, length(lookup) + 1)))
        lookup_split <- split(lookup, section_split)
        bestguess <- sapply(lookup_split, function(sec) {
          atlaslines <- grep("Atlas CA_ML_18_MNIA: Macro Labels (N27)", sec, fixed=TRUE)
          nomatch <- grep("***** Not near any region stored in databases *****", sec, fixed=TRUE)
          if (length(nomatch) > 0L) {
            return("Unable to identify label")
          } else {
            return(sub("(^\\s*|\\s*$)", "", sec[atlaslines+1], perl=TRUE)) #first match after atlas for each cluster          
          }
        })
        
        coordlines <- grep("Focus point (LPI)=", lookup, fixed=TRUE)
        coords <- lookup[coordlines+2] #first line after header is TLRC, second is MNI
        #coords <- sub("<a href=.*$", "", coords, perl=TRUE)
        coords <- sub("^\\s*(-?\\d+\\s*mm.*\\{MNI\\})\\s*<a href=.*$", "\\1", coords, perl=TRUE)

        coords_l <- as.numeric(sub("^\\s*(-*\\d+) mm.*", "\\1", coords, perl=TRUE))
        coords_p <- as.numeric(sub("^\\s*(-*\\d+) mm \\[(?:L|R)\\],\\s+(-*\\d+) mm.*", "\\2", coords, perl=TRUE))
        coords_i <- as.numeric(sub("^\\s*(-*\\d+) mm \\[(?:L|R)\\],\\s+(-*\\d+) mm \\[(?:A|P)\\],\\s+(-*\\d+) mm.*", "\\3", coords, perl=TRUE))
        
        cluster_metadata <- data.frame(model=this_model, l1_contrast=l1_contrast_name, l2_contrast=l2_contrast_name, l3_contrast=l3_contrast_name,
          cluster_number=1:length(coords_l), cluster_size=vsizes, cluster_threshold=clustsize, z_threshold=zthresh,
          x=coords_l, y=coords_p, z=coords_i, label=bestguess, stringsAsFactors=FALSE)
        
        roimask <- readAFNI(paste0(clust_brik, "+tlrc.HEAD"), vol=1)
        #afni masks tend to read in as 4D matrix with singleton 4th dimension. Fix this
        if (length(dim(roimask)) == 4L) { roimask@.Data <- roimask[,,,,drop=T] }
        
        maskvals <- sort(unique(as.vector(roimask)))
        maskvals <- maskvals[!maskvals == 0]
        
        #generate a matrix of roi averages across subjects
        #this should be subjects x clusters in size
        roimats <- get_cluster_means(roimask, copeconcat)
        
        roi_df <- reshape2::melt(roimats, value.name="cope_value", varnames=c("feat_input_id", "cluster_number")) %>%
          mutate(model=this_model, l1_contrast=l1_contrast_name, l2_contrast=l2_contrast_name, l3_contrast=l3_contrast_name) %>%
          #full_join(design_df %>% select(subject, !!l3_contrast_name), by="subject") #merge with relevant covariate
          full_join(design_df, by="feat_input_id") #merge with all covariates

        #handle beta series extraction (NB. beta_series_inputs should be in same order as subject_inputs based on use of sub above)
        beta_series_df <- get_beta_series(beta_series_inputs, roimask, n_bs=50)

        #for identification, add cluster information to beta series from ROI data.frame
        beta_series_df <- beta_series_df %>% left_join(roi_df, by=c("feat_input_id", "cluster_number")) %>%
          select(feat_input_id, run, trial, cluster_number, everything())

        #beta_series_df %>% group_by(feat_input_id) %>% summarize(mean(bs_value), mean(cope_value))
        
        coords <- lapply(maskvals, function(v) {
          mi <- which(roimask==v, arr.ind=TRUE)
          return(mi)
        })

        names(coords) <- make.names(bestguess, unique=TRUE)
        
        ##get correlation matrix for each ROI (note that the dimension of the latent subspace is constrained by the number of subjects -- eigenvalues become 0 after N)
        #corrmats <- lapply(roimats, function(r) { return(cor(r)) })
        
        #thiscope <- list(cluster_metadata=cluster_metadata, roivals=roimats, coords=coords, corrmats=corrmats)
        #curoutdir <- file.path(getwd(), l1copes[l1], l2copes[l2])
        #dir.create(curoutdir, showWarnings=FALSE, recursive=TRUE)
        #save(cluster_metadata, roimats, coords, corrmats, file=file.path(curoutdir, paste0(paste(l1copes[l1], l2copes[l2], l3copes[l3], sep="_"), "_betas.RData")))

        #all_metadata[[paste(l1, l2, l3, sep=".")]] <- cluster_metadata
        #all_rois[[paste(l1, l2, l3, sep=".")]] <- roi_df
        
        l2_loop_cluster_metadata[[paste(l1, l2, l3, sep=".")]] <- cluster_metadata
        l2_loop_rois[[paste(l1, l2, l3, sep=".")]] <- roi_df
        l2_loop_bs[[paste(l1, l2, l3, sep=".")]] <- beta_series_df
      }

      #return(list(cluster_metadata=l2_loop_cluster_metadata, rois=l2_loop_rois))
      l2_loop_outputs[[l2]] <- list(cluster_metadata=l2_loop_cluster_metadata, rois=l2_loop_rois, beta_series=l2_loop_bs)
    }

    #tack on metadata and roi betas from this l2 contrast to the broader set
    all_metadata <- bind_rows(all_metadata, rlang::flatten(lapply(l2_loop_outputs, "[[", "cluster_metadata")))
    all_rois <- bind_rows(all_rois, rlang::flatten(lapply(l2_loop_outputs, "[[", "rois")))
    all_beta_series <- bind_rows(all_beta_series, rlang::flatten(lapply(l2_loop_outputs, "[[", "beta_series")))
  }

  #organize models intelligently
  all_metadata <- all_metadata %>% arrange(model, l1_contrast, l2_contrast, l3_contrast)
  all_rois <- all_rois %>% arrange(model, l1_contrast, l2_contrast, l3_contrast, cluster_number, feat_input_id)
  all_beta_series <- all_beta_series %>% arrange(model, l1_contrast, l2_contrast, l3_contrast, cluster_number, feat_input_id, run, trial)

  readr::write_csv(x=all_metadata, file.path(model_output_dir, paste0(l1_contrast_name, "_cluster_metadata.csv")))
  readr::write_csv(x=all_rois, file.path(model_output_dir, paste0(l1_contrast_name, "_roi_betas.csv")))
  readr::write_csv(x=all_beta_series, file.path(model_output_dir, paste0(l1_contrast_name, "_roi_beta_series.csv.gz")))

  #not uniquely useful at present (CSVs have it all)
  #save(all_metadata, all_rois, dmat, file=file.path(model_output_dir, "sceptic_clusters.RData"))

}

try(stopCluster(cl)) #cleanup pool upon exit of this script
