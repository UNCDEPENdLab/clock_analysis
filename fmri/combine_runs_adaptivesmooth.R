setwd("/Volumes/Serena/MMClock/MR_Proc/10711_20140826/native_nosmooth/rglm_tc_carry")
mrfiles <- list.files(path="/Volumes/Serena/MMClock/MR_Proc/10711_20140826/native_nosmooth/rglm_tc_carry", pattern="^run.*me.*.nii.gz")

source(file.path(getMainDir(), "clock_analysis", "fmri", "glm_helper_functions.R"))

conditions <- factor(sub("run[0-9]+_me\\.(\\w+)\\.nii\\.gz", "\\1", mrfiles, perl=TRUE))
runs <- factor(sub("run([0-9])+_me\\.\\w+\\.nii\\.gz", "\\1", mrfiles, perl=TRUE))

mrfiles_bycondition <- split(mrfiles, conditions)

library(plyr)
library(oro.nifti)
library(foreach)
library(doSNOW)
library(abind)

setDefaultClusterOptions(master="localhost")
clusterobj <- makeSOCKcluster(8)
registerDoSNOW(clusterobj)

for (cond in 1:length(mrfiles_bycondition)) {
    ##transform into standard space
    for (f in 1:length(mrfiles_bycondition[[cond]])) {
        run <- sub("run([0-9])+_me\\.\\w+\\.nii\\.gz", "\\1", mrfiles_bycondition[[cond]][f], perl=TRUE)
        runFSLCommand(paste0("applywarp --in=", mrfiles_bycondition[[cond]][f], " --out=warp_", mrfiles_bycondition[[cond]][f],
                             " --interp=trilinear --warp=../clock", run, "/func_to_MNI_2.3mm_warpfield ",
                             "--ref=/Users/michael/standard/mni_icbm152_nlin_asym_09c/mni_icbm152_t1_tal_nlin_asym_09c_2.3mm.nii ",
                             "--mask=/Users/michael/standard/mni_icbm152_nlin_asym_09c/mni_icbm152_t1_tal_nlin_asym_09c_mask_2.3mm"), fsldir="/opt/ni_tools/fsl")
    }

    allruns <- lapply(mrfiles_bycondition[[cond]], function(f) {
        dat <- readNIfTI(paste0("warp_", f), reorient = FALSE)
    })

    betas <- abind(lapply(allruns, function(r) { r@.Data[,,,1] }), along=4)
    variances <- abind(lapply(allruns, function(r) { r@.Data[,,,2] }), along=4)
    mrout <- allruns[[1]] #copy first run as template for output
    rm(allruns)

    any_zero <- apply(betas, 4, function(timepoint) { as.numeric(timepoint==0.0) })
    any_zero <- array(apply(any_zero, 1, sum), dim=dim(betas)[1:3])
    analyze_indices <- which(any_zero==0.0, arr.ind=TRUE) #only keep voxels that are not missing across runs

    #voxels x runs matrices for betas and variances
    beta_analyze <- apply(betas, 4, "[", analyze_indices)
    var_analyze <- apply(variances, 4, "[", analyze_indices)
      
    ##build a weighted least squares model from first level analyses where weighted estimates of the betas are combined across runs
    combruns <- foreach(row=1:nrow(var_analyze), .combine=rbind, .inorder=TRUE) %dopar% {
        b_vec <- beta_analyze[row,]
        w_vec <- 1/var_analyze[row,]
        b_wls <- sum(b_vec*w_vec)/sum(w_vec)
        v_wls <- sqrt(1/sum(w_vec))
        t_wls <- b_wls / v_wls
        c(b=b_wls, v=v_wls, t=t_wls)
    }

    wls_data <- array(NA_real_, dim=c(dim(betas)[1:3], ncol(combruns)))

    ind4d <- cbind(pracma::repmat(analyze_indices, ncol(combruns), 1), rep(1:ncol(combruns), each=nrow(combruns)))
    wls_data[ind4d] <- combruns
    
    mrout@.Data <- wls_data
    mrout@dim_[5] <- ncol(combruns)
    mrout@cal_min <- min(mrout, na.rm=TRUE) #need to have min and max vals in header that match data for oro.nifti to succeed
    mrout@cal_max <- max(mrout, na.rm=TRUE)
    writeNIfTI(mrout, filename=paste0("combined_", names(mrfiles_bycondition)[cond]))
    
}

stopCluster(clusterobj)
