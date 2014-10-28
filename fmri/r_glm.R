library(fmri)
library(fitclock)
#load a subject
setwd(file.path(getMainDir(), "clock_analysis", "fmri", "fmri_fits"))
load(file.path(getMainDir(), "clock_analysis", "fmri", "fmri_fits", "10811_fitinfo.RData"))

# This is a demonstration of how mean centering a parametric regressor prior to convolution separates the neural activation into the eventness and
# parametric pieces. Summed together, these recover the direct convolution of the parametric regressor with the HRF.
onsets1 <- f_value$clock_onset[1,]
durations1 <- f_value$RTraw[1,]
evRun1 <- f_value$ev[1,]

pdf("demonstration of mean-centered parametric before convolution.pdf", width=10, height=8)
par(mfrow=c(4,1))
evtOnly <- fmri.stimulus(scans=283, times=onsets1, durations=durations1/1000, values=1.0, mean=FALSE, rt=1.0)
plot(1:length(evtOnly), evtOnly, type="l", main="event convolved", xlab="")

meanPreConvolve <- fmri.stimulus(scans=283, times=onsets1, durations=durations1/1000, values=(evRun1-mean(evRun1)), mean=FALSE, rt=1.0)
plot(1:length(noMean), meanPreConvolve, type="l", col="blue", main="mean-centered EV pre-convolution", xlab="")

noMean <- fmri.stimulus(scans=283, times=onsets1, durations=durations1/1000, values=evRun1, mean=FALSE, rt=1.0)
plot(1:length(noMean), noMean, type="l", main="ev convolved", col="red", xlab="")

m <- lm(noMean ~ meanPreConvolve + evtOnly)
summary(m)

f <- fitted(m)

plot(1:length(f), f, type="l", main="fitted values of event (1) + mean-cent EV (2) = EV convolved regressor (3)", col="orange", xlab="Time")
dev.off()
fmri.stimulus(scans=last_fmri_volume[i], values=reg[,"value"], times=reg[,"onset"], durations=reg[,"duration"], rt=1.0, mean=FALSE) #hard-coded 1.0s TR for now

build_design_matrix(fitobj=f_value, regressors=c("clock", "feedback", "ev", "rpe_neg", "rpe_pos"), 
    event_onsets=c("clock_onset", "feedback_onset", "clock_onset", "feedback_onset", "feedback_onset"), 
    durations=c("clock_duration", "feedback_duration", "clock_duration", "feedback_duration", "feedback_duration"), 
    normalizations=c("none", "none", "none", "none", "none"),
    baselineCoefOrder=2, writeTimingFiles=c("AFNI"), output_directory="~/10811_timingtest_meanonly")


build_design_matrix(fitobj=f_value, regressors=c("clock", "feedback", "ev", "rpe_neg", "rpe_pos"), 
    event_onsets=c("clock_onset", "feedback_onset", "clock_onset", "feedback_onset", "feedback_onset"), 
    durations=c("clock_duration", "feedback_duration", "clock_duration", "feedback_duration", "feedback_duration"), 
    normalizations=c("durmax_1.0", "durmax_1.0", "evtmax_1.0", "evtmax_1.0", "evtmax_1.0"),
    baselineCoefOrder=2, writeTimingFiles=c("AFNI"), output_directory="~/10811_timingtest")

build_design_matrix(fitobj=f_value, regressors=c("clock", "feedback", "ev", "rpe_neg", "rpe_pos"), 
    event_onsets=c("clock_onset", "feedback_onset", "clock_onset", "feedback_onset", "feedback_onset"), 
    durations=c("clock_duration", "feedback_duration", "clock_duration", "feedback_duration", "feedback_duration"), 
    baselineCoefOrder=2, writeTimingFiles=c("AFNI"), output_directory="~/10811_timingtest")

#runVolumes=runlengths, runsToOutput=mrrunnums, output_directory=timingdir)


#read in convolved regressors from above
clock_onset <- read.table("/Volumes/Serena/MMClock/MR_Raw/11178_20140310/MBclock_recon//afni_tc/run_timing_tc/run8_clock.1D")$V1
rel_uncertainty <- read.table("/Volumes/Serena/MMClock/MR_Raw/11178_20140310/MBclock_recon/clock8_nosmooth/afni_tc/run_timing_tc/run8_rel_uncertainty.1D")$V1
mean_uncertainty <- read.table("/Volumes/Serena/MMClock/MR_Raw/11178_20140310/MBclock_recon/clock8_nosmooth/afni_tc/run_timing_tc/run8_mean_uncertainty.1D")$V1
pos_rpe <- read.table("/Volumes/Serena/MMClock/MR_Raw/11178_20140310/MBclock_recon/clock8_nosmooth/afni_tc/run_timing_tc/run8_rpe_pos.1D")$V1
neg_rpe <- read.table("/Volumes/Serena/MMClock/MR_Raw/11178_20140310/MBclock_recon/clock8_nosmooth/afni_tc/run_timing_tc/run8_rpe_neg.1D")$V1

dmat <- cbind(clock_onset, rel_uncertainty, mean_uncertainty, pos_rpe, neg_rpe)
dmat_fmri <- fmri.design(dmat, order=3)

mrdat <- read.NIFTI("/Volumes/Serena/MMClock/MR_Raw/11178_20140310/MBclock_recon/clock8_nosmooth/nfwudktm_clock8_trunc303.nii", level = 0.75,setmask=TRUE)
anat <- read.NIFTI("/Volumes/Serena/MMClock/MR_Raw/11178_20140310/MBclock_recon/clock8_nosmooth/template_brain.nii", level = 0.75,setmask=FALSE)
#result <- fmri.lm(mrdat, dmat_fmri, actype="smooth", contrast=diag(9), vvector=rep(1,9))

#just contrast for first beta
result <- fmri.lm(mrdat, dmat_fmri, actype="smooth", contrast=1, vvector=1)

sm <- fmri.smooth(result, adaptation="fullaws")
summary.fmridata(sm)

pv <- fmri.pvalue(sm, mode="local")

plot(pv, anatomic=anat, maxpvalue=.05)

#the $segm field has values 0 (null), 1 (positive activation), and -1 (negative activation) for each voxel basd on segmentation.
sm_seg <- fmri.smooth(result, adaptation="segment", alpha=.05)

save(result, sm, file="fmritest.RData")

library(abind)
b_se_t_p_seg <- abind(sm$cbeta, sm$var, sm$cbeta/sqrt(sm$var), 1 - pv$pvalue, sm_seg$segm, along=4)
niftiHead <- sm$header
niftiHead$dimension[5] <- 5 #b, se, t, p, seg
write.NIFTI(b_se_t_p_seg, niftiHead, filename="beta_test")

d_demean <- build_design_matrix(fitobj=f_value, regressors=c("clock", "feedback", "ev", "rpe_neg", "rpe_pos"), 
    event_onsets=c("clock_onset", "feedback_onset", "clock_onset", "feedback_onset", "feedback_onset"), 
    durations=c("clock_duration", "feedback_duration", "clock_duration", "feedback_duration", "feedback_duration"),
    normalizations=c("durmax_1.0", "durmax_1.0", "evtmax_1.0", "evtmax_1.0", "evtmax_1.0"),
    baselineCoefOrder=2)

visualizeDesignMatrix <- function(df) {
  require(ggplot2)
  require(reshape2)
  print(cor(df))
  df$volume <- 1:nrow(df)
  df.m <- melt(df, id.vars="volume")
  ggplot(df.m, aes(x=volume, y=value)) + geom_line(size=1.5) + theme_bw(base_size=15) + facet_grid(variable ~ ., scales="free_y") 
}

#\
r_valueModel <- function(f_value, mrfiles, runlengths, mrrunnums) {
  library(orthopolynom)
  unnormalized.p.list <- legendre.polynomials( 3, normalized=FALSE )
  polynomial.values(polynomials=unnormalized.p.list, x=seq(-1,1, 0.1))
  
  
  d_value <- build_design_matrix(fitobj=f_value, regressors=c("clock", "ev", "feedback", "rpe_neg", "rpe_pos"), 
      event_onsets=c("clock_onset", "clock_onset", "feedback_onset", "feedback_onset", "feedback_onset"), 
      durations=c("clock_duration", "clock_duration", "feedback_duration", "feedback_duration", "feedback_duration"),
      normalizations=c("durmax_1.0", "evtmax_1.0", "durmax_1.0", "evtmax_1.0", "evtmax_1.0"),
      baselineCoefOrder=2, writeTimingFiles=c("AFNI"),
      runVolumes=runlengths, runsToOutput=mrrunnums, output_directory=timingdir)
  
  for (f in mrfiles) {
    
  }
  #basic main effect matrix
  regressornames <- c("clock", "feedback", "ev", "rpe_pos", "rpe_neg")
  
  contrastmat <- c()
  
  for (r in 1:length(regressornames)) {
    thisCon <- rep(0.0, length(regressornames))
    thisCon[r] <- 1.0
    contrastmat <- rbind(contrastmat,thisCon)
    rownames(contrastmat)[r] <- paste("me", regressornames[r], sep=".") 
  }
  
#need to get an emotion-modulated matrix at some point, although this is between runs
  
}

r_glm <- 
    
    




contrastMat <- cbind(
    c()
    )



#library(oro.nifti)
#o <- fmri2oro(mrdat)
