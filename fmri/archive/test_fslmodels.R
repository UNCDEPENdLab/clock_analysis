# TODO: Add comment
# 
# Author: michael
###############################################################################
setwd(file.path(getMainDir(), "clock_analysis/fmri/10873_testbed"))

#quick test on 10873 as model subject who gives reasonable EV maps in early runs
#trying to sort out fixed versus estimated learning rates for rewards and whether to carryover value from one run to the next

mrfiles <- c(
    "/Volumes/Serena/MMClock/MR_Proc/10873_20140918/mni_5mm_wavelet/clock1/nfswudktm_clock1_5_drop6_trunc276.nii.gz",
    "/Volumes/Serena/MMClock/MR_Proc/10873_20140918/mni_5mm_wavelet/clock2/nfswudktm_clock2_5_drop6_trunc264.nii.gz",
    "/Volumes/Serena/MMClock/MR_Proc/10873_20140918/mni_5mm_wavelet/clock3/nfswudktm_clock3_5_drop6_trunc276.nii.gz",
    "/Volumes/Serena/MMClock/MR_Proc/10873_20140918/mni_5mm_wavelet/clock4/nfswudktm_clock4_5_drop6_trunc287.nii.gz",
    "/Volumes/Serena/MMClock/MR_Proc/10873_20140918/mni_5mm_wavelet/clock5/nfswudktm_clock5_5_drop6_trunc255.nii.gz",
    "/Volumes/Serena/MMClock/MR_Proc/10873_20140918/mni_5mm_wavelet/clock6/nfswudktm_clock6_5_drop6_trunc298.nii.gz",
    "/Volumes/Serena/MMClock/MR_Proc/10873_20140918/mni_5mm_wavelet/clock7/nfswudktm_clock7_5_drop6_trunc296.nii.gz",
    "/Volumes/Serena/MMClock/MR_Proc/10873_20140918/mni_5mm_wavelet/clock8/nfswudktm_clock8_5_drop6_trunc294.nii.gz"
)

mrrunnums <- 1:8

mrrunlengths <- c(276, 264, 276, 287, 255, 298, 296, 294)

#models

library(fitclock)

source("../fslValueModel.R")
source("../glm_helper_functions.R")
##Sys.setenv(FSLDIR="/usr/local/ni_tools/fsl")
Sys.setenv(FSLDIR="/opt/ni_tools/fsl")
if (!file.exists("10873_fits.RData")) {
  s <- clockdata_subject(subject_ID='10873', dataset="fMRIEmoClock_10873_tc_tcExport.csv")
  vm_2rates <- deltavalue_model(clock_data=s, alphaV=0.2, betaV=0.1, carryover_value=TRUE)

  f_fixed <- vm_2rates$predict(returnFit=TRUE)
  fslValueModel(f_fixed, mrfiles, mrrunlengths, mrrunnums, run=TRUE, force=FALSE, dropVolumes=6, outdir="fsl_fixed_carryover")
    
  f_estim <- vm_2rates$fit() #estimate learning rate as a free parameter
  fslValueModel(f_estim, mrfiles, mrrunlengths, mrrunnums, run=TRUE, force=FALSE, dropVolumes=6, outdir="fsl_estim_carryover")

  #now model without carryover of ev from one run to the next
  vm_2rates$carryover_value <- FALSE

  vm_2rates$theta["alphaV", "cur"] <- 0.2
  vm_2rates$theta["betaV", "cur"] <- 0.1
  f_fixed_nocarryover <- vm_2rates$predict(returnFit=TRUE) #estimate learning rate as a free parameter
  fslValueModel(f_fixed_nocarryover, mrfiles, mrrunlengths, mrrunnums, run=TRUE, force=FALSE, dropVolumes=6, outdir="fsl_fixed_nocarryover")  
  
  f_estim_nocarryover <- vm_2rates$fit() #estimate learning rate as a free parameter
  fslValueModel(f_estim_nocarryover, mrfiles, mrrunlengths, mrrunnums, run=TRUE, force=FALSE, dropVolumes=6, outdir="fsl_estim_nocarryover")  
  
}

