source(file.path(getMainDir(), "fitclock", "R", "utility_functions.R"))
setwd(file.path(getMainDir(), "clock_analysis", "fmri"))
onsets <- c(10,30,50,70,90,110,130,150,170,190,210,230)
durations <- rep(c(0.5, 10), length(onsets)/2)
#durations <- c(1:10,13,15)
#durations <- rep(10, length(onsets))
height1 <- rep(c(5, 10, 15), length(onsets)/3)
height2 <- rep(c(20,-20), length(onsets)/2)


h1.c <- height1 - mean(height1)
h2.c <- height2 - mean(height2)

#spm friston 1998 double gamma
#evt_boxcar <- fmri.stimulus(scans=300, values=1, times=10, durations=150, rt=1.0, a1=6, a2=16, b1=1, b2=1, c=1/6, mean=FALSE) #hard-coded 1.0s TR for now

#glover 1999 double gamma
evt_boxcar <- fmri.stimulus(scans=300, values=1, times=10, durations=150, rt=1.0, a1=6, a2=12, b1=0.9, b2=0.9, c=0.35, mean=FALSE) #hard-coded 1.0s TR for now

evt_scaleup <- fmri.stimulus(scans=300, values=1, times=onsets, durations=durations, rt=1.0, a1=6, a2=12, b1=0.9, b2=0.9, c=0.0, mean=FALSE) #hard-coded 1.0s TR for now



evt_rescale <- evt_scaleup/max(evt_scaleup)
#evt_boxcar <- fmri.stimulus(scans=300, values=1, times=10, durations=150, rt=1.0, a1=6, a2=16, b1=1, b2=1, c=0, mean=FALSE) #hard-coded 1.0s TR for now
plot(1:length(evt_boxcar), evt_rescale, type="l")

evt <- fmri.stimulus(scans=300, values=rep(1, length(onsets)), times=onsets, durations=durations, rt=1.0, cc=0) #hard-coded 1.0s TR for now
h1.convolve <- fmri.stimulus(scans=300, values=h1.c, times=onsets, durations=durations, rt=1.0, cc=0) #hard-coded 1.0s TR for now
h2.convolve <- fmri.stimulus(scans=300, values=h2.c, times=onsets, durations=durations, rt=1.0, cc=0) #hard-coded 1.0s TR for now

#AFNI dmUBLOCK(0): longer durations lead to larger modulations with the max being the height of the parametric regressor (after 10 seconds)
#AFNI dmUBLOCK(1): each event is convolved with gamma maxing at height=1.0 (but with proper duration)

dmblock <- read.table("nodata.xmat.1D")[,c(-1,-2)]
dmublock0 <- read.table("dmblocktest_x1d_dmUBLOCK0.xmat.1D")[,c(-1,-2)]
dmublock1 <- read.table("dmblocktest_x1d_dmUBLOCK1.xmat.1D")[,c(-1,-2)]

names(dmublock0) <- names(dmublock1) <- c("afni_evt", "afni_h1", "afni_h2")

round(cor(cbind(dmublock0, evt, h1.convolve, h2.convolve)), 4)

dev.new()
par(mfrow=c(3,1))
evt.rescale <- plotrix::rescale(evt, c(0,1))
plot(1:length(evt), evt, type="l", col="blue")
plot(1:length(h1.convolve), h1.convolve, type="l", col="red")
plot(1:length(h2.convolve), h2.convolve, type="l", col="black")

dev.new()
par(mfrow=c(3,1))
plot(1:nrow(dmublock0), dmublock0[,1], type="l", col="blue")
plot(1:nrow(dmublock0), dmublock0[,2], type="l", col="red")
plot(1:nrow(dmublock0), dmublock0[,3], type="l", col="black")

dev.new()
par(mfrow=c(3,1))
plot(1:nrow(dmublock1), dmublock1[,1], type="l", col="blue")
plot(1:nrow(dmublock1), dmublock1[,2], type="l", col="red")
plot(1:nrow(dmublock1), dmublock1[,3], type="l", col="black")


#so, the fmri.stimulus with mean centering recovers the dmUBLOCK(0) coding of activation increasing in proportion to both duration and regressor values
#the scaling difference comes from the fact that the amplitude is not constrained between 0 and 1 in the canonical convolution

#should we try to rescale these to have bounds defined by the parameter?

#rescale for 0-1 range (per AFNI standard)
afniblock <- function(t, q=4, rescale=5.05){#19285) {
    ret <- t^q*exp(-t)/(4^4*exp(-4))/rescale
    ret[which(ret < 0.0 | ret > 1.0)] <- 0.0
    ret
}


boxcar <- rep(c(0,1), each=50)
#boxcar <- c(rep(0,length(boxcar)), boxcar, rep(0,length(boxcar)))
test <- convolve(x=boxcar, y=afniblock(1:length(boxcar)), type="open")
plot(1:length(test), test, type="l")

time <- 300
boxcar <- rep(0, time)

for (i in 1:length(onsets)) {
    boxcar[(onsets[i]):(onsets[i]+durations[i] - 1)] <- 1.0
}

#zeropad
boxcar <- c(rep(0,length(boxcar)), boxcar, rep(0,length(boxcar)))

#can't remember why, but need to reverse kernel convolution in R to get it right...
test <- convolve(x=boxcar, y=rev(afniblock(1:length(boxcar))))
#test <- convolve(x=boxcar, y=afniblock(1:length(boxcar)), conj=FALSE)
#test <- convolve(x=boxcar, y=afniblock(1:length(boxcar)), conj=FALSE)
test <- test[(length(boxcar)/3+1):(length(boxcar)/3*2)]
#plot(1:length(test), test, type="l")
#lines(1:length(test), dmublock0[,1], type="l", col="blue")
#cor(test, dmublock0[,1])
#round(cbind(test, dmublock0[,1]), 4)

#dev.new()


#compare to block
#plot(1:length(test), test, type="l")
#lines(1:length(test), dmblock, type="l", col="blue")

padconv <- c(0,0,test[1:(length(test)-2)])
plot(1:length(test), padconv, type="l")
lines(1:length(test), dmblock, type="l", col="blue")


cor(padconv, dmblock)
cbind(padconv, dmblock)
round(cbind(padconv, dmblock), 3)

#function convolves regressor with normalized HRF. Normalization options are:
# 1) pre-convolution HRF max=1.0 normalization of each stimulus regardless of duration: identical to dmUBLOCK(1)
# 2) pre-convolution HRF max=1.0 normalization for long events (15+ sec) -- height of HRF is modulated by duration of event: identical to dmUBLOCK(0) 
hrf_convolve_normalize <- function(totalscans, onsets, durations, values, rt=1.0, normeach=FALSE, demean=FALSE, ...) {
  
  #this is my hacky way to figure out the peak value 
  #obtain an estimate of the peak HRF height of a long event convolved with the HRF at these settings of a1, b1, etc.
  hrf_boxcar <- fmri.stimulus(scans=300, values=1.0, times=100, durations=100, rt=1.0, mean=FALSE, ...) #don't mean center for computing max height
  hrf_max <- max(hrf_boxcar)
  
  #for each event, convolve it with hrf, normalize, then sum convolved events to get full timecourse
  normedEvents <- sapply(1:length(onsets), function(i) {
        #obtain unit-convolved duration-modulated regressor to define HRF prior to modulation by parametric regressor
        stim_conv <- fmri.stimulus(scans=totalscans, values=1.0, times=onsets[i], durations=durations[i], rt=rt, mean=FALSE, ...)
        if (normeach) {
          stim_conv <- stim_conv/max(stim_conv) #rescale HRF to a max of 1.0 for each event, regardless of duration -- EQUIVALENT TO dmUBLOCK(1)
        } else {
          stim_conv <- stim_conv/hrf_max #rescale HRF to a max of 1.0 for long event -- EQUIVALENT TO dmUBLOCK(0)
        }
        
        stim_conv <- stim_conv*values[i] #for each event, multiply by parametric regressor value
      })
  
  tc_conv <- apply(normedEvents, 1, sum)
  
  #grand mean center convolved regressor
  if (demean) { tc_conv <- tc_conv - mean(tc_conv) }
  
  tc_conv
}

#afniNormConvolve <- function(totalscans, onsets, durations, values, rt=1.0, normeach=FALSE) {
#  require(plotrix)
#  
#  if (normeach) {
#    #THIS IS EQUIVALENT TO dmUBLOCK(1)
#    #obtain unit-convolved duration-modulated regressor to define HRF prior to modulation by parametric regressor
#    normedEvents <- sapply(1:length(onsets), function(i) {
#          stim <- fmri.stimulus(scans=totalscans, values=1.0, times=onsets[i], durations=durations[i], rt=rt, cc=0.0) #no undershoot
#          stim <- stim/max(stim) #rescale HRF to a max of 1.0 for each event
#          stim <- stim*values[i] #for each event, multiply by parametric regressor value
#        })
#
#    conv <- apply(normedEvents, 1, sum)
#  } else {
#    #THIS IS EQUIVALENT TO dmUBLOCK(0)
#    #obtain a duration-modulated HRF for the whole run (longer events lead to larger durations)
#    
#    #a hack for now, but find the max of a long convolved boxcar to scale the hrf to max = 1.0 for long events (~15 seconds or more)
#    hrf_boxcar <- fmri.stimulus(scans=300, values=1.0, times=100, durations=100, rt=1.0, cc=0.0, mean=FALSE) #don't mean center
#    hrf_max <- max(hrf_boxcar)
#    
#    hrf <- fmri.stimulus(scans=totalscans, values=values, times=onsets, durations=durations, rt=rt, cc=0.0) #no undershoot
#    if (all(values==1.0)) {
#      conv <- rescale(hrf, c(0,1)) #unit regressor -- rescale to 0-1
#    } else { 
#      conv <- rescale(hrf, c(min(values), max(values))) #rescale to min and max of parametric regressor
#    }
#  }
#  
#  conv
#}

#VALIDATE DMUBLOCK(0)
normed_evt <- hrf_convolve_normalize(300, onsets, durations, values=rep(1, length(onsets)), normeach=FALSE, cc=0)
normed_h1 <- hrf_convolve_normalize(300, onsets, durations, values=h1.c, cc=0)
normed_h2 <- hrf_convolve_normalize(300, onsets, durations, values=h2.c, cc=0)

#compare to AFNI 
#event regressor
par(mfrow=c(3,1))
plot(1:length(normed_evt), normed_evt, type="l")
lines(1:length(normed_evt), dmublock0[,1], type="l", col="blue")
cor(dmublock0[,1], normed_evt)
#round(cbind(dmublock0[,1], normed_evt), 4)

#+/- 5
plot(1:length(normed_h1), normed_h1, type="l")
lines(1:length(normed_h1), dmublock0[,2], type="l", col="blue")
cor(dmublock0[,2], normed_h1)

#+/- 20
plot(1:length(normed_h2), normed_h2, type="l")
lines(1:length(normed_h2), dmublock0[,3], type="l", col="blue")
cor(dmublock0[,3], normed_h2)


#VALIDATE DMUBLOCK(1)
normed_evt <- hrf_convolve_normalize(300, onsets, durations, values=rep(1, length(onsets)), normeach=TRUE)
normed_h1 <- hrf_convolve_normalize(300, onsets, durations, values=h1.c, normeach=TRUE)
normed_h2 <- hrf_convolve_normalize(300, onsets, durations, values=h2.c, normeach=TRUE)

#compare to AFNI 
#event regressor
par(mfrow=c(3,1))
plot(1:length(normed_evt), normed_evt, type="l")
lines(1:length(normed_evt), dmublock1[,1], type="l", col="blue")
cor(dmublock1[,1], normed_evt)
round(cbind(dmublock0[,1], normed_evt), 4)

#+/- 5
plot(1:length(normed_h1), normed_h1, type="l")
lines(1:length(normed_h1), dmublock1[,2], type="l", col="blue")
cor(dmublock1[,2], normed_h1)

#+/- 20
plot(1:length(normed_h2), normed_h2, type="l")
lines(1:length(normed_h2), dmublock1[,3], type="l", col="blue")
cor(dmublock1[,3], normed_h2)

# Notes from discussion with Alex: the "Eventness" regressor where the height of the hrf varies by duration makes good sense because we will be getting a time on task effect where regions
# that are involved in stimulus processing should show more activity for longer trials.

# For the parametric regressor, allowing duration to modulate the height of the HRF leads to a confounding of the parametric values with duration. For example, a large value of the regressor
# with a short duration could look identical to a small value with a long duration.

# Thus, for cognitive regressors, makes sense for the parametric value either to be 1) stick-convolved (zero duration) or 2) RT convolved, but normalized to 1.0 for each stimulus, which is equivalent to dmUBLOCK(1).
# The nice part about the per-stimulus normalization is that the bounds of the regressor follow the max and min of the parametric values independent of duration, whereas for the dmUBLOCK(0) approach where only the entire
# HRF (all trials) is normalized to 0-1, then only for long events with large parametric values does the HRF approach the max/min of the parametric regressor.


# provisional solution:
# 1) separate eventness regressors for clock and feedback. regressors will be duration convolved with no per-stimulus normalization, but a grand normalization of 0-1.0 (although with undershoot, it should be negative, right?).
# 2) parametric regressors with per-stimulus 0-1 normalization 

# things to compare across models
# 1) model with whole-trial eventness regressor versus separate events
# 2) model with eventness regressor where duration modulates height versus per stimulus normalization



load(file.path(getMainDir(), "clock_analysis", "fmri", "fmri_fits", "10811_fitinfo.RData"))

# This is a demonstration of how mean centering a parametric regressor prior to convolution separates the neural activation into the eventness and
# parametric pieces. Summed together, these recover the direct convolution of the parametric regressor with the HRF.
onsets1 <- f_value$clock_onset[1,]
durations1 <- f_value$RTraw[1,]
evRun1 <- f_value$ev[1,]
rpes <- f_value$rpe[1,]
par(mfrow=c(2,1))
plot(1:length(rpes), rpes, type="l")

plot(1:length(evRun1), evRun1, type="l")

par(mfrow=c(3,1))
evtOnly <- fmri.stimulus(scans=283, times=onsets1, durations=durations1/1000, values=1.0, mean=FALSE, rt=1.0)
plot(1:length(evtOnly), evtOnly, type="l", main="event convolved", xlab="")

rpeCent <- fmri.stimulus(scans=283, times=onsets1, durations=durations1/1000, values=(rpes-mean(rpes)), mean=FALSE, rt=1.0)
plot(1:length(rpeCent), rpeCent, type="l", main="rpe centered convolved", xlab="")

evCent <- fmri.stimulus(scans=283, times=onsets1, durations=durations1/1000, values=(evRun1-mean(evRun1)), mean=FALSE, rt=1.0)
plot(1:length(evCent), evCent, type="l", main="ev centered convolved", xlab="")

cor(cbind(evtOnly, rpeCent, evCent))


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

d <- build_design_matrix(fitobj=f_value, regressors=c("clock", "feedback", "ev", "rpe_neg", "rpe_pos"), 
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
clock_onset <- read.table("/Volumes/Serena/MMClock/MR_Proc/11178_20140310/mni_5mm_wavelet/afni_value/run_timing_deltavalue/run8_clock.1D")$V1
#rel_uncertainty <- read.table("/Volumes/Serena/MMClock/MR_Raw/11178_20140310/MBclock_recon/clock8_nosmooth/afni_tc/run_timing_tc/run8_rel_uncertainty.1D")$V1
#mean_uncertainty <- read.table("/Volumes/Serena/MMClock/MR_Raw/11178_20140310/MBclock_recon/clock8_nosmooth/afni_tc/run_timing_tc/run8_mean_uncertainty.1D")$V1
pos_rpe <- read.table("/Volumes/Serena/MMClock/MR_Proc/11178_20140310/mni_5mm_wavelet/afni_value/run_timing_deltavalue/run8_rpe_pos.1D")$V1
neg_rpe <- read.table("/Volumes/Serena/MMClock/MR_Proc/11178_20140310/mni_5mm_wavelet/afni_value/run_timing_deltavalue/run8_rpe_neg.1D")$V1

#dmat <- cbind(clock_onset, rel_uncertainty, mean_uncertainty, pos_rpe, neg_rpe)
dmat <- cbind(clock_onset, pos_rpe, neg_rpe)
dmat_fmri <- fmri.design(dmat, order=3)

library(oro.nifti)
library(fmri)
mrdat <- oro2fmri(readNIfTI("/Volumes/Serena/MMClock/MR_Proc/11178_20140310/native_nosmooth/clock8/nfudktm_clock8_trunc303.nii.gz", reorient=FALSE), level = 0.75,setmask=TRUE)
#anat <- read.NIFTI("/Volumes/Serena/MMClock/MR_Proc/11178_20140310/native_nosmooth/clock8/template_brain.nii", level = 0.75,setmask=FALSE)
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
