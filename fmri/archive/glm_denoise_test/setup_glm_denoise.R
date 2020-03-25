#try to upsample fMRI data to 0.1s by linear interpolation, and also lay down a design matrix on a 0.1s grid
#then run glmdenoise

setwd(file.path(getMainDir(), "clock_analysis", "fmri"))
source("glm_helper_functions.R")
library(oro.nifti)
library(R.matlab)
library(rhdf5)
library(plyr)

#go down to 3x3x3 voxels for memory constraints
for (run in 1:8) {
  runAFNICommand(paste0("3dresample -overwrite -rmode Li -dxyz 3.0 3.0 3.0 ",
          "-inset /Volumes/Serena/MMClock/MR_Proc/10711_20140826/mni_5mm_wavelet/clock", run, "/fswudktm_clock", run, "_5.nii.gz ",
          "-prefix /Volumes/Serena/MMClock/MR_Proc/10711_20140826/mni_5mm_wavelet/clock", run, "/fswudktm_clock", run, "_5_3mm.nii.gz"), afnidir="/opt/ni_tools/afni")
}


#glmdenoise mentions that it's important not to mean-normalize the data, so we need to go back to fswudktm files
#only retain data matrix, not header info
upfactor=5
tr=1.0

starts <- rep(7, 8) #always start with 7th volume
ends <- c(281, 283, 278, 281, 280, 280, 288, 275)

for (run in 1:8) {
  r_orig <- readNIfTI(paste0("/Volumes/Serena/MMClock/MR_Proc/10711_20140826/mni_5mm_wavelet/clock", run, "/fswudktm_clock", run, "_5_3mm.nii.gz"), reorient=FALSE)@.Data[,,,starts[run]:ends[run]]
  acqtimes <- 0:(dim(r_orig)[4]*tr - tr)
  
  #linear interpolation for each voxel
  acqtimes_interp <- seq(0,max(acqtimes), by=tr/upfactor)

  system.time(r_up <- apply(r_orig, c(1,2,3), function(v) {
        approx(acqtimes, v, xout=acqtimes_interp)$y
      }))
  
  #returns time as first dimension
  r_up <- aperm(r_up, c(2,3,4,1))
  
  #assign(paste0("r", run), r_up)
  #writeMat(con=paste0("10711_r", run, ".mat"), rdat=r_up)
  
  outdir <- "/Volumes/Serena/MMClock/MR_Proc/10711_20140826/mni_5mm_wavelet/glmdenoise_test"
  h5out <- file.path(outdir, paste0("r", run, ".h5"))
  dir.create(outdir, showWarnings=FALSE)
  
  #use hdf5 format to read into matlab (writeMat above was failing due to large size)
  h5createFile(h5out)
  h5createDataset(file=h5out, dataset=paste0("run", run), dims=dim(r_up), level=4) #level of compression runs 0-9, use 4 for speed
  h5write(r_up, file=h5out, name=paste0("run", run))
  H5close() #close any handles
  
}



tdir <- "/Volumes/Serena/MMClock/MR_Proc/10711_20140826/mni_5mm_wavelet/fsl_tc_nocarry/run_timing_tc/"
#"run1_clock_FSL3col.txt"
#r1clock=read.table('/Volumes/Serena/MMClock/MR_Proc/10711_20140826/mni_5mm_wavelet/fsl_tc_nocarry/run_timing_tc/run1_clock_FSL3col.txt', header=FALSE, col.names=c("onset", "duration", "value"))
designs <- list()
designs_zerodur <- list()
conditions <- c("clock", "feedback")
for (run in 1:8) {
  events <- list()
  events_zerodur <- list()
  for(cond in conditions) {
    timing <- read.table(paste0(tdir, "run", run, "_", cond, "_FSL3col.txt"), header=FALSE, col.names=c("onset", "duration", "value"))
    rlength <- ends[run] - starts[run] + 1 #number of volumes
    acqtimes <- 0:(rlength-tr)*tr
    grid <- seq(0,max(acqtimes), by=tr/upfactor)
    onsets <- sapply(timing$onset, round_any, tr/upfactor)
    durations <- sapply(timing$duration, round_any, tr/upfactor)

    #resample timing onto 0.2s grid
    events_grid <- do.call(c, apply(cbind(onsets, durations), 1, function(row) {
              #events are all 0.2 seconds long, so always subtract off the event length when making sequence so that
              #final event spans the last 0.2s interval
              round_any(seq(row[1], row[1] + row[2] - tr/upfactor, by=tr/upfactor), tr/upfactor) #need to use round_any again to avoid tiny differences in grid vs. events_grid comparison below 
            }))

   
    stick_function <- rep(0, length(grid))
    stick_function[grid %in% events_grid] <- 1.0

    #onsets only (no duration -- for FIR)
    stick_function_zerodur <- rep(0, length(grid))
    stick_function_zerodur[grid %in% onsets] <- 1.0
    
    events[[cond]] <- stick_function
    events_zerodur[[cond]] <- stick_function_zerodur
  }
  
  designs[[paste0("run", run)]] <- do.call(cbind, events)
  designs_zerodur[[paste0("run", run)]] <- do.call(cbind, events_zerodur)
}

#all designs as a structure
writeMat(con="10711_designmats_0p2s.mat", designs=designs, designs_zerodur=designs_zerodur)


tdir <- "/Volumes/Serena/MMClock/MR_Proc/10711_20140826/mni_5mm_wavelet/afni_tc"
designconv <- local({load(file.path(tdir, "afni_tc_design.RData")); environment()})$d
mats <- lapply(designconv$design.convolve, function(run) {
      as.matrix(run[,!grepl("base[0-9]+", names(run))])
    })

writeMat(con="10711_designmats_convolved.mat", designs=mats, regnames=dimnames(mats[[1]])[[2]])