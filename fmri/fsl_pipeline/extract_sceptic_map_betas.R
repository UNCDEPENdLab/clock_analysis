# This script saves the significant clusters for each map in a SCEPTIC group analysis

#load the master configuration file
to_run <- Sys.getenv("fsl_pipeline_file")

run_model_index <- as.numeric(Sys.getenv("run_model_index")) #which variant to execute
if (nchar(to_run) == 0L) { stop("Cannot locate environment variable fsl_pipeline_file") }
if (!file.exists(to_run)) { stop("Cannot locate configuration file", to_run) }
if (is.na(run_model_index)) { stop("Couldn't identify usable run_model_index variable.") }

load(to_run)

source(file.path(fsl_model_arguments$pipeline_home, "functions", "glm_helper_functions.R"))

library(tidyverse)
library(dependlab)
library(oro.nifti)
library(parallel)

#1) load spatial maps RData object
#2) rebuild into 4d cube (where fourth dimension is run/subject)
#3) clusterize each effect of interest in stats outputs using 3dclust and generating mask
#library(ggplot2)
library(abind)
library(oro.nifti)
library(reshape2)
library(robust)
library(car)
library(dplyr)
library(doParallel)

#verify that mr_dir is present as expected
subinfo <- fsl_model_arguments$subject_covariates
feat_run_outdir <- fsl_model_arguments$outdir[run_model_index] #the name of the subfolder for the current run-level model
feat_lvl3_outdir <- file.path(fsl_model_arguments$group_output_dir, feat_run_outdir) #output directory for this run-level model
n_l1_copes <- fsl_model_arguments$n_l1_copes[run_model_index] #number of l1 copes determines number of FEAT LVL3 analyses to run (1 per LVL1 cope)
l1_cope_names <- fsl_model_arguments$l1_cope_names[[run_model_index]] #names of l1 copes (used for folder naming)
zthresh <- fsl_model_arguments$zthresh #3.09
clustsize <- fsl_model_arguments$clustsize #34

cl <- makeCluster(fsl_model_arguments$n_cluster_beta_cpus)
registerDoParallel(cl)

subinfo$dir_found <- file.exists(subinfo$mr_dir)

feat_lvl2_dirname <- "FEAT_LVL2_runtrend.gfeat" #should populate this to the structure at some point
models <- fsl_model_arguments$group_model_variants #different covariate models for the current run-level model (run_model_index)

all_metadata <- list()
all_rois <- list()

for (l1 in 1:n_l1_copes) {
  l1_contrast_name <- l1_cope_names[l1]
  model_output_dir <- file.path(feat_lvl3_outdir, l1_contrast_name)
  
  for (this_model in models) {
    expect_gfeat <- file.path(model_output_dir, paste0(l1_contrast_name, "-", paste(this_model, collapse="-"), ".gfeat"))

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

    design_df <- data.frame(dmat) %>% mutate(subject=1:n())

    #figure out what the l2 contrasts are and read relevant statistics for each
    l2fsf <- file.path(subject_inputs[1], "design.fsf")
    stopifnot(file.exists(l2fsf))
    l2_syntax <- readLines(l2fsf)

    l2_contrast_info <- grep("\\s*set fmri\\(conname_real\\.\\d+\\).*", l2_syntax, value=TRUE, perl=TRUE)
    l2_contrast_nums <- as.numeric(sub("\\s*set fmri\\(conname_real\\.(\\d+)\\).*", "\\1", l2_contrast_info, perl=TRUE))
    l2_contrast_names <- sub("\\s*set fmri\\(conname_real\\.\\d+\\)\\s*\"?([^\"]+)\"?.*", "\\1", l2_contrast_info, perl=TRUE)

    n_l2_contrasts <- max(l2_contrast_nums)
    l2_contrast_names <- l2_contrast_names[order(l2_contrast_nums)] #order l2 contrast names in ascending order to match l2 loop below
    
    #loop over l2 contrasts
    l2_loop_outputs <- foreach(l2=iter(1:n_l2_contrasts), .packages=c("oro.nifti", "dplyr")) %dopar% {
      #for (l2 in 1:n_l2_contrasts) {
      l2_loop_cluster_metadata <- list()
      l2_loop_rois <- list()
      
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
        
        cluster_metadata <- data.frame(l1_contrast=l1_contrast_name, l2_contrast=l2_contrast_name, l3_contrast=l3_contrast_name,
          cluster_number=1:length(coords_l), cluster_size=vsizes, cluster_threshold=clustsize, z_threshold=zthresh,
          x=coords_l, y=coords_p, z=coords_i, labels=bestguess, stringsAsFactors=FALSE)
        
        roimask <- readAFNI(paste0(clust_brik, "+tlrc.HEAD"), vol=1)
        #afni masks tend to read in as 4D matrix with singleton 4th dimension. Fix this
        if (length(dim(roimask)) == 4L) { roimask@.Data <- roimask[,,,,drop=T] }
        
        maskvals <- sort(unique(as.vector(roimask)))
        maskvals <- maskvals[!maskvals == 0]
        
        #generate a matrix of roi averages across subjects
        #this should be subjects x clusters in size
        roimats <- sapply(maskvals, function(v) {
          mi <- which(roimask==v, arr.ind=TRUE)
          nsubj <- length(copefiles)
          nvox <- nrow(mi)
          mi4d <- cbind(pracma::repmat(mi, nsubj, 1), rep(1:nsubj, each=nvox))
          
          mat <- matrix(copeconcat[mi4d], nrow=nvox, ncol=nsubj) #need to manually reshape into matrix from vector
          #for each subject, compute huber m-estimator of location/center Winsorizing at 2SD across voxels (similar to voxel mean)
          #clusavg <- apply(mat, 2, function(x) { MASS::huber(x, k=2)$mu })
          clusavg <- apply(mat, 2, mean)
          return(clusavg)
          
          #return(t(mat)) #transpose to get subjects x voxels matrix            
        })
        
        roi_df <- reshape2::melt(roimats, value.name="cope_value", varnames=c("subject", "cluster_number")) %>%
          mutate(l1_contrast=l1_contrast_name, l2_contrast=l2_contrast_name, l3_contrast=l3_contrast_name) %>%
          #full_join(design_df %>% select(subject, !!l3_contrast_name), by="subject") #merge with relevant covariate
          full_join(design_df, by="subject") #merge with all covariates

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
      }

      return(list(cluster_metadata=l2_loop_cluster_metadata, rois=l2_loop_rois))
    }

    all_metadata <- bind_rows(rlang::flatten(lapply(l2_loop_outputs, "[[", "cluster_metadata")))
    all_rois <- bind_rows(rlang::flatten(lapply(l2_loop_outputs, "[[", "rois")))
  }  
}

readr::write_csv(x=all_metadata, file.path(model_output_dir, "cluster_metadata.csv"))
readr::write_csv(x=all_rois, file.path(model_output_dir, "roi_betas.csv"))

save(all_metadata, all_rois, dmat,  file=file.path(model_output_dir,"sceptic_clusters.RData"))

try(stopCluster(cl)) #cleanup pool upon exit of this script
