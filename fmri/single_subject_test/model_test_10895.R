#setwd(file.path(getMainDir(), "clock_analysis", "fmri"))
library(fitclock)
#load(file.path(getMainDir(), "clock_analysis", "fmri", "fmri_fits", "10895_fitinfo.RData"))
load(file.path(getMainDir(), "clock_analysis", "fmri", "fmri_fits", "11262_fitinfo.RData"))
source(file.path(getMainDir(), "rs-fcMRI_Motion", "rs-fcMRI_Motion_Functions.R"))
source(file.path(getMainDir(), "clock_analysis", "fmri", "glm_model_selection.R"))
setwd("/Volumes/bek/clock_10895_singlesubject_modeling")

#these are the number of good volumes per run for this subject (10895)
subj_runlengths <- c(300, 286, 276, 307, 278, 274, 298, 299)

#d <- build_design_matrix(fitobj=f_value, regressors=c("clock", "feedback"), 
#    event_onsets=c("clock_onset", "feedback_onset"), 
#    durations=c("clock_duration", 0), 
#    normalizations=c("none", "none"),
#    baselineCoefOrder=2, writeTimingFiles=c("FSL", "AFNI"), runVolumes=subj_runlengths, output_directory="~/10895_model_tests/10895_clockRT_feedback0")
#
#d <- build_design_matrix(fitobj=f_value, regressors=c("clock", "feedback"), 
#    event_onsets=c("clock_onset", "feedback_onset"), 
#    durations=c(0, 0), 
#    normalizations=c("none", "none"),
#    baselineCoefOrder=2, writeTimingFiles=c("AFNI"), runVolumes=subj_runlengths, output_directory="~/10895_model_tests/10895_clock0_feedback0")
#
#d <- build_design_matrix(fitobj=f_value, regressors=c("clock", "feedback", "rpe_pos", "rpe_neg"), 
#    event_onsets=c("clock_onset", "feedback_onset", "feedback_onset", "feedback_onset"), 
#    durations=c("clock_duration", 0, 0, 0), 
#    normalizations=c("none", "none", "none", "none"),
#    baselineCoefOrder=2, writeTimingFiles=c("AFNI"), runVolumes=subj_runlengths, output_directory="~/10895_model_tests/10895_clockRT_feedback0_rpepos0_rpeneg0_nonorm")
#
#d <- build_design_matrix(fitobj=f_value, regressors=c("clock", "feedback", "rpe_pos", "rpe_neg"), 
#    event_onsets=c("clock_onset", "feedback_onset", "feedback_onset", "feedback_onset"), 
#    durations=c("clock_duration", 0, 0, 0), 
#    normalizations=c("durmax1.0", "durmax1.0", "evtmax1.0", "evtmax1.0"),
#    baselineCoefOrder=2, writeTimingFiles=c("AFNI"), runVolumes=subj_runlengths, output_directory="~/10895_model_tests/10895_clockRT_feedback0_rpepos0_rpeneg0_evtmaxRPE")
#
#
#d <- build_design_matrix(fitobj=f_value, regressors=c("clock", "feedback", "rpe_pos", "rpe_neg"), 
#    event_onsets=c("clock_onset", "feedback_onset", "feedback_onset", "feedback_onset"), 
#    durations=c("clock_duration", "feedback_duration", 0, 0), 
#    normalizations=c("none", "none", "none", "none"),
#    baselineCoefOrder=2, writeTimingFiles=c("AFNI"), runVolumes=subj_runlengths, output_directory="~/10895_model_tests/10895_clockRT_feedbackDur_rpepos0_rpeneg0_nonorm")


#as one would predict, using the feedback duration in the convolution is almost identical to a 0-duration feedback event
#thorndike:10895_model_tests michael$ 1dCorrelate 10895_clockRT_feedbackDur_rpepos0_rpeneg0_nonorm/feedback_concat.1D 10895_clockRT_feedback0_rpepos0_rpeneg0_nonorm/feedback_concat.1D
# Pearson correlation [n=2318 #col=2]
# Name                   Name                    Value   BiasCorr   2.50%   97.50%  N: 2.50% N:97.50%
# ---------------------  ---------------------  -------- -------- -------- -------- -------- --------
# feedback_concat.1D[0]  feedback_concat.1D[0]  +0.98289 +0.98287 +0.98143 +0.98416 +0.98145 +0.98422

#this makes sense given that feedback is a fixed and relatively short duration event: .9s

#the clock regressor, however, depends much more on the duration convolution
#thorndike:10895_model_tests michael$ 1dCorrelate 10895_clock0_feedback0/clock_concat.1D 10895_clockRT_feedback0/clock_concat.1D
# Pearson correlation [n=2318 #col=2]
# Name                Name                 Value   BiasCorr   2.50%   97.50%  N: 2.50% N:97.50%
# ------------------  ------------------  -------- -------- -------- -------- -------- --------
# clock_concat.1D[0]  clock_concat.1D[0]  +0.48130 +0.48101 +0.45004 +0.51138 +0.44939 +0.51198

#this is also unsurprising because clock durations vary a lot


#d <- build_design_matrix(fitobj=f_value, regressors=c("clock", "feedback", "rpe_pos", "rpe_neg", "ev"), 
#    event_onsets=c("clock_onset", "feedback_onset", "feedback_onset", "feedback_onset", "clock_onset"), 
#    durations=c("clock_duration", "feedback_duration", 0, 0, "clock_duration"), 
#    normalizations=c("none", "none", "none", "none"),
#    baselineCoefOrder=2, writeTimingFiles=c("AFNI"), runVolumes=subj_runlengths, output_directory="~/10895_model_tests/10895_clockRT_feedbackDur_rpepos0_rpeneg0_nonorm")
#
#
##look at high-pass filtering design
#spec.pgram(d$design.convolve[[1]], plot=TRUE, detrend=TRUE, main="fMRI periodogram", log="no")

#spectrum(d$design.convolve[[1]], plot=TRUE, detrend=TRUE, taper=0, main="full.corrupt MR ~ motion.bp, then bandpass BP", log="yes")

##try out ar(1) modeling and AIC checking

library(oro.nifti)
runfiles <- sort(list.files(pattern="nfswudktm.*drop6_3mm.nii.gz"))
runmask <- readNIfTI("runmask_3mm.nii.gz", reorient=FALSE)@.Data
##allruns <- foreach(f=1:length(runfiles), .packages="oro.nifti") %dopar% {
voxAnalyze <- which(runmask > 0, arr.ind=TRUE) #vector of voxel indices corresponding to good voxels

if (!file.exists("allRunsRaw_Sparse_3mm.RData")) {
    
    library(abind)

    library(Matrix)
    library(slam)
    library(plyr)
    library(iterators)
    library(doSNOW)

    allruns <- list()

    ##generate a sparse storage approach that can be used to analyze voxels in parallel
    ##giving up on Matrix because it does not allow row-wise iterator and I don't want to export huge dataset to each worked
    ##giving up on slam because it takes 5 seconds just to access one row!

    for (f in 1:length(runfiles)) {
        mr <- readNIfTI(runfiles[f], reorient=FALSE)@.Data
        mr <- apply(mr, 4, function(vol) { vol[voxAnalyze] })
        
        ##mr <- Matrix(apply(mr, 4, function(vol) { vol * runmask }), sparse=TRUE) #mask by runmask and convert to voxel x time sparse matrix
        allruns[[f]] <- mr

        ##  mr3 <- aaply(mr, 4, function(vol) { vol * runmask }) #mask by runmask
        ##  mr3 <- aperm(mr3, c(2,3,4,1)) #aaply always puts the split dimension first
        ##  mr4 <- Matrix(mr3, sparse=TRUE)

        rm(mr)
    }


    ##prewhiten each run separately (whitening concatenated data makes no sense given transitions between runs)
    ##skipping for the moment given that it's slow and we really want to whiten the residuals from the lm calls, not the time series only
    ##white <- foreach(r=iter(allruns), .noexport=c("concatMR", "allruns"), .inorder=TRUE) %dopar% {
    ##  w <- apply(r, 1, function(vox) {
    ##        #based on playing around with a few time series, looks at ARMA(2,2) is a good compromise for most time series
    ##        #and doesn't take too long to compute. use faster CSS-ML method unless an error occurs, in which case the safer, slower ML is used.
    ##        m <- tryCatch(arima(vox, c(2, 0, 2), method="CSS-ML"), error=function(e) { arima(vox, c(2, 0, 2), method="ML") } )
    ##        as.vector(residuals(m))
    ##      })
    ##}

    ##bind into a voxels x time concatenated matrix
    concatMR <- do.call(cbind, allruns)
    save(concatMR, file="allRunsRaw_Sparse_3mm.RData")

    ##concatWhite <- do.call(cbind, white)
    ##save(concatWhite, file="allRunsWhite_Sparse.RData")

}

#load voxel x time matrix of all runs
load("allRunsRaw_Sparse_3mm.RData")

#ev + rpe with parametric values centered pre-convolution
#d <- build_design_matrix(fitobj=f_value, regressors=c("clock", "feedback", "ev", "rpe_pos", "rpe_neg"), 
#    event_onsets=c("clock_onset", "feedback_onset", "clock_onset", "feedback_onset", "feedback_onset"), 
#    durations=c("clock_duration", 0, "clock_duration", 0, 0), 
#    normalizations=c("none", "none", "none", "none", "none"), center_values=TRUE,
#    baselineCoefOrder=2, runVolumes=subj_runlengths, convolve_wi_run=FALSE, dropVolumes = 6)
#
#d_allruns_center <- concatDesign(d, hpass=.01)
#g <- visualizeDesignMatrix(d_allruns_center, outfile="dmat_center.pdf", includeBaseline = FALSE)
#
#d <- build_design_matrix(fitobj=f_value, regressors=c("clock", "feedback", "ev", "rpe_pos", "rpe_neg"), 
#    event_onsets=c("clock_onset", "feedback_onset", "clock_onset", "feedback_onset", "feedback_onset"), 
#    durations=c("clock_duration", 0, "clock_duration", 0, 0), 
#    normalizations=c("none", "none", "none", "none", "none"), center_values=FALSE,
#    baselineCoefOrder=2, runVolumes=subj_runlengths, convolve_wi_run=FALSE, dropVolumes = 6)
#
#d_allruns_nocenter <- concatDesign(d, hpass=.01)
#g <- visualizeDesignMatrix(d_allruns_nocenter, outfile="dmat_uncenter.pdf", includeBaseline = FALSE)
#
#d <- build_design_matrix(fitobj=f_value, regressors=c("clock", "feedback", "ev", "rpe_pos", "rpe_neg"), 
#    event_onsets=c("clock_onset", "feedback_onset", "clock_onset", "feedback_onset", "feedback_onset"), 
#    durations=c("clock_duration", 0, "clock_duration", 0, 0), 
#    normalizations=c("durmax_1.0", "evtmax_1.0", "durmax_1.0", "evtmax_1.0", "evtmax_1.0"), center_values=TRUE,
#    baselineCoefOrder=2, runVolumes=subj_runlengths, convolve_wi_run=FALSE, dropVolumes = 6)
#
#d_allruns_evtmax <- concatDesign(d, hpass=.01)
#g <- visualizeDesignMatrix(d_allruns_evtmax, outfile="dmat_center_evtmaxPar.pdf", includeBaseline = FALSE)

##conclusion: mean centering parametric regressors has major effect on collinearity
##normalization largely affects scaling. in the case of evtmax_1.0 (each hrf impulse is rescaled to a height of 1.0 before convolution),
##the convolved regressor follows the units of regressor. in the case of durmax_1.0, the convolved regressor would only have the same units
##as the regressor for long events, but captures the modulation of the hrf by duration.
##cor(cbind(d_allruns_center[,1:5], d_allruns_evtmax[,1:5]))

##add high pass to match fmri
#d <- build_design_matrix(fitobj=f_value, regressors=c("clock", "feedback", "ev", "rpe_pos", "rpe_neg"), 
#    event_onsets=c("clock_onset", "feedback_onset", "clock_onset", "feedback_onset", "feedback_onset"), 
#    durations=c("clock_duration", 0, "clock_duration", 0, 0), 
#    normalizations=c("durmax_1.0", "evtmax_1.0", "durmax_1.0", "evtmax_1.0", "evtmax_1.0"), center_values=TRUE, convolve_wi_run=TRUE,
#    baselineCoefOrder=2, runVolumes=subj_runlengths, dropVolumes = 6, high_pass=.01)


#MODEL 1: 
##clock (RT durmax) + feedback (0 evtmax) sticks
##ev@clock: RT durmax
##rpe_pos@feedback: 0 evtmax
##rpe_neg@feedback: 0 evtmax
d <- build_design_matrix(fitobj=f_value, regressors=c("clock", "feedback", "ev", "rpe_pos", "rpe_neg"), 
    event_onsets=c("clock_onset", "feedback_onset", "clock_onset", "feedback_onset", "feedback_onset"), 
    durations=c("clock_duration", 0, "clock_duration", 0, 0), 
    normalizations=c("durmax_1.0", "evtmax_1.0", "durmax_1.0", "evtmax_1.0", "evtmax_1.0"), center_values=TRUE, convolve_wi_run=TRUE,
    baselineCoefOrder=2, runVolumes=subj_runlengths, dropVolumes = 6)

res <- incr_fit(concatDesign(d, hpass=.01), concatMR, mask3d=voxAnalyze, modelname="01_ev_rpepos_rpeneg_evtmax_wirunscaling", outputPvals=FALSE, njobs=20)

##MODEL 2:
##Update model 1: do not mean center parametric regressor prior to convolution
d <- build_design_matrix(fitobj=f_value, regressors=c("clock", "feedback", "ev", "rpe_pos", "rpe_neg"), 
    event_onsets=c("clock_onset", "feedback_onset", "clock_onset", "feedback_onset", "feedback_onset"), 
    durations=c("clock_duration", 0, "clock_duration", 0, 0), 
    normalizations=c("durmax_1.0", "evtmax_1.0", "durmax_1.0", "evtmax_1.0", "evtmax_1.0"), center_values=FALSE, convolve_wi_run=TRUE,
    baselineCoefOrder=2, runVolumes=subj_runlengths, dropVolumes = 6)

res <- incr_fit(concatDesign(d, hpass=.01), concatMR, mask3d=voxAnalyze, modelname="02_ev_rpepos_rpeneg_evtmax_nocenter_wirunscaling", outputPvals=FALSE, njobs=20)

##MODEL 3: (EV aligned at feedback)
##clock (RT durmax) + feedback (0 evtmax) sticks
##ev@feedback: 0 evtmax
##rpe_pos@feedback: 0 evtmax
##rpe_neg@feedback: 0 evtmax
d <- build_design_matrix(fitobj=f_value, regressors=c("clock", "feedback", "ev", "rpe_pos", "rpe_neg"), 
    event_onsets=c("clock_onset", "feedback_onset", "feedback_onset", "feedback_onset", "feedback_onset"), 
    durations=c("clock_duration", 0, 0, 0, 0), 
    normalizations=c("durmax_1.0", "evtmax_1.0", "evtmax_1.0", "evtmax_1.0", "evtmax_1.0"), center_values=TRUE, convolve_wi_run=TRUE,
    baselineCoefOrder=2, runVolumes=subj_runlengths, dropVolumes = 6)

res <- incr_fit(concatDesign(d, hpass=.01), concatMR, mask3d=voxAnalyze, modelname="03_evfeedback_rpepos_rpeneg_evtmax_wirunscaling", outputPvals=FALSE, njobs=20)

##MODEL 4:
##look at sqrt RPE
d <- build_design_matrix(fitobj=f_value, regressors=c("clock", "feedback", "ev", "rpe_pos_sqrt", "rpe_neg_sqrt"), 
    event_onsets=c("clock_onset", "feedback_onset", "clock_onset", "feedback_onset", "feedback_onset"), 
    durations=c("clock_duration", 0, "clock_duration", 0, 0), 
    normalizations=c("durmax_1.0", "evtmax_1.0", "durmax_1.0", "evtmax_1.0", "evtmax_1.0"), center_values=TRUE, convolve_wi_run=TRUE,
    baselineCoefOrder=2, runVolumes=subj_runlengths, dropVolumes = 6)

res <- incr_fit(concatDesign(d, hpass=.01), concatMR, mask3d=voxAnalyze, modelname="04_ev_rpepossqrt_rpenegsqrt_evtmax_wirunscaling", outputPvals=FALSE, njobs=20)

##MODEL 5:
##switch EV to evtmax_1.0 scaling (duration does not increase height relative to other events)
d <- build_design_matrix(fitobj=f_value, regressors=c("clock", "feedback", "ev", "rpe_pos", "rpe_neg"), 
    event_onsets=c("clock_onset", "feedback_onset", "clock_onset", "feedback_onset", "feedback_onset"), 
    durations=c("clock_duration", 0, "clock_duration", 0, 0), 
    normalizations=c("durmax_1.0", "evtmax_1.0", "evtmax_1.0", "evtmax_1.0", "evtmax_1.0"), center_values=TRUE, convolve_wi_run=TRUE,
    baselineCoefOrder=2, runVolumes=subj_runlengths, dropVolumes = 6)

res <- incr_fit(concatDesign(d, hpass=.01), concatMR, mask3d=voxAnalyze, modelname="05_evevtmax_rpepos_rpeneg_wirunscaling", outputPvals=FALSE, njobs=20)

##MODEL 6:
##all 0 duration events
d <- build_design_matrix(fitobj=f_value, regressors=c("clock", "feedback", "ev", "rpe_pos", "rpe_neg"), 
    event_onsets=c("clock_onset", "feedback_onset", "clock_onset", "feedback_onset", "feedback_onset"), 
    durations=c(0, 0, 0, 0, 0), 
    normalizations=c("evtmax_1.0", "evtmax_1.0", "evtmax_1.0", "evtmax_1.0", "evtmax_1.0"), center_values=TRUE, convolve_wi_run=TRUE,
    baselineCoefOrder=2, runVolumes=subj_runlengths, dropVolumes = 6)

res <- incr_fit(concatDesign(d, hpass=.01), concatMR, mask3d=voxAnalyze, modelname="06_ev_rpepos_rpeneg_wirunscaling_all0dur", outputPvals=FALSE, njobs=20)

##MODEL 7:
#additional variance of EV at feedback after adding EV at clock
#clock + feedback nuisance
#ev@clock: evtmax 1.0
#rpe_pos@feedback: 0dur
#rpe_neg@feedback: 0dur
#ev@feedback: 0dur
d <- build_design_matrix(fitobj=f_value, regressors=c("clock", "feedback", "ev", "rpe_pos", "rpe_neg", "ev"), 
    event_onsets=c("clock_onset", "feedback_onset", "clock_onset", "feedback_onset", "feedback_onset", "feedback_onset"), 
    durations=c("clock_duration", 0, "clock_duration", 0, 0, 0), 
    normalizations=c("durmax_1.0", "none", "evtmax_1.0", "none", "none", "none"), center_values=TRUE, convolve_wi_run=TRUE,
    baselineCoefOrder=2, runVolumes=subj_runlengths, dropVolumes = 6)

res <- incr_fit(concatDesign(d, hpass=.01), concatMR, mask3d=voxAnalyze, modelname="07_ev_rpepos_rpeneg_evfeedback_wirunscaling", outputPvals=FALSE, njobs=20)

##MODEL 8:
##badre model
##clock: RT durmax 1.0
##mean_uncertainty@clock: RT evtmax 1.0
##rel_uncertainty@clock: RT evtmax 1.0
##ev@clock: RT evtmax 1.0
##feedback: RT durmax 1.0
##rpe_pos@feedback: RT evtmax 1.0
##rpe_neg@feedback: RT evtmax 1.0

d <- build_design_matrix(fitobj=f, regressors=c("clock", "mean_uncertainty", "rel_uncertainty", "ev", "feedback", "rpe_pos", "rpe_neg"), 
    event_onsets=c("clock_onset", "clock_onset", "clock_onset", "clock_onset", "feedback_onset", "feedback_onset", "feedback_onset"), 
    durations=c("clock_duration", "clock_duration", "clock_duration", "clock_duration", 0, 0, 0), 
    normalizations=c("durmax_1.0", "evtmax_1.0", "evtmax_1.0", "evtmax_1.0", "durmax_1.0", "evtmax_1.0", "evtmax_1.0"), center_values=TRUE, convolve_wi_run=TRUE,
    baselineCoefOrder=2, runVolumes=subj_runlengths, dropVolumes = 6)

res <- incr_fit(concatDesign(d, hpass=.01), concatMR, mask3d=voxAnalyze, modelname="08_badre_evtmax", outputPvals=FALSE, njobs=20)

##MODEL 9:
##badre model with durmax 1.0 clock parametric regressors (i.e., longer events have bigger scaling)
##clock: RT durmax 1.0
##mean_uncertainty@clock: RT durmax 1.0
##rel_uncertainty@clock: RT durmax 1.0
##ev@clock: RT durmax 1.0
##feedback: RT evtmax 1.0
##rpe_pos@feedback: RT evtmax 1.0
##rpe_neg@feedback: RT evtmax 1.0

d <- build_design_matrix(fitobj=f, regressors=c("clock", "mean_uncertainty", "rel_uncertainty", "ev", "feedback", "rpe_pos", "rpe_neg"), 
    event_onsets=c("clock_onset", "clock_onset", "clock_onset", "clock_onset", "feedback_onset", "feedback_onset", "feedback_onset"), 
    durations=c("clock_duration", "clock_duration", "clock_duration", "clock_duration", 0, 0, 0), 
    normalizations=c("durmax_1.0", "durmax_1.0", "durmax_1.0", "durmax_1.0", "evtmax_1.0", "evtmax_1.0", "evtmax_1.0"), center_values=TRUE, convolve_wi_run=TRUE,
    baselineCoefOrder=2, runVolumes=subj_runlengths, dropVolumes = 6)

res <- incr_fit(concatDesign(d, hpass=.01), concatMR, mask3d=voxAnalyze, modelname="09_badre_durmax_clockregressors_NEW", outputPvals=FALSE, outputDf=FALSE, njobs=22)

#13Nov2014: visual comparison of models 8 and 9 indicates that the clock-aligned model regressors are typically much stronger for the evtmax normalization than durmax.
#This is consistent with the idea that the 


##MODEL 10:
##As with model 9 (badre durmax), but add temporal derivatives
d <- build_design_matrix(fitobj=f, regressors=c("clock", "mean_uncertainty", "rel_uncertainty", "ev", "feedback", "rpe_pos", "rpe_neg"), 
    event_onsets=c("clock_onset", "clock_onset", "clock_onset", "clock_onset", "feedback_onset", "feedback_onset", "feedback_onset"), 
    durations=c("clock_duration", "clock_duration", "clock_duration", "clock_duration", 0, 0, 0), 
    normalizations=c("durmax_1.0", "durmax_1.0", "durmax_1.0", "durmax_1.0", "evtmax_1.0", "evtmax_1.0", "evtmax_1.0"), center_values=TRUE, convolve_wi_run=TRUE,
    baselineCoefOrder=2, runVolumes=subj_runlengths, dropVolumes = 6)

res <- incr_fit(concatDesign(d, hpass=.01), concatMR, mask3d=voxAnalyze, modelname="10_badre_durmax_clockregressors_tempderiv", outputPvals=FALSE, outputDf=FALSE, njobs=20, add_derivs=TRUE)

##MODEL 11:
##As with model 8 (badre evtmax), but add temporal derivatives
d <- build_design_matrix(fitobj=f, regressors=c("clock", "mean_uncertainty", "rel_uncertainty", "ev", "feedback", "rpe_pos", "rpe_neg"), 
    event_onsets=c("clock_onset", "clock_onset", "clock_onset", "clock_onset", "feedback_onset", "feedback_onset", "feedback_onset"), 
    durations=c("clock_duration", "clock_duration", "clock_duration", "clock_duration", 0, 0, 0), 
    normalizations=c("durmax_1.0", "evtmax_1.0", "evtmax_1.0", "evtmax_1.0", "durmax_1.0", "evtmax_1.0", "evtmax_1.0"), center_values=TRUE, convolve_wi_run=TRUE,
    baselineCoefOrder=2, runVolumes=subj_runlengths, dropVolumes = 6)

res <- incr_fit(concatDesign(d, hpass=.01), concatMR, mask3d=voxAnalyze, modelname="11_badre_evtmax_tempderiv", outputPvals=FALSE, njobs=22, add_derivs=TRUE)

#looking at 8 vs. 11, there is only slight improvement in rsq across a number of regions, and activation for the primary beta is weaker in a number of areas
#given the complexity of doing statistics on the combination of the regressor and derivative betas, going with model 8 for now for group analysis.
#but interested in pursuing the "derivative boost" idea from Lindquist, which necessitates nonparametric stats at the group level












####TEST USE OF DERIVATIVES WITHIN BUILD DESIGN
#d <- build_design_matrix(fitobj=f, regressors=c("clock", "mean_uncertainty", "rel_uncertainty", "ev", "feedback", "rpe_pos", "rpe_neg"), 
#    event_onsets=c("clock_onset", "clock_onset", "clock_onset", "clock_onset", "feedback_onset", "feedback_onset", "feedback_onset"), 
#    durations=c("clock_duration", "clock_duration", "clock_duration", "clock_duration", 0, 0, 0), 
#    normalizations=c("durmax_1.0", "durmax_1.0", "durmax_1.0", "durmax_1.0", "evtmax_1.0", "evtmax_1.0", "evtmax_1.0"), center_values=TRUE, convolve_wi_run=TRUE,
#    baselineCoefOrder=2, runVolumes=subj_runlengths, dropVolumes = 6, add_derivs=TRUE)
#
#g <- visualizeDesignMatrix(concatDesign(d, hpass=.01), outfile="dmat_deriv_test.pdf", includeBaseline = FALSE)

###
#try fitting model without value carryover
s <- clockdata_subject(subject_ID='10895', dataset='fMRIEmoClock_10895_tc_tcExport.csv')
vm1 <- deltavalue_model(clock_data=s, alphaV=0.3, betaV=0.3, carryover_value=FALSE) #N.B. This matches V matrix from full time-clock algorithm fit.
f_value1 <- vm1$fit() #estimate learning rate as a free parameter

vm2 <- deltavalue_model(clock_data=s, alphaV=0.3, betaV=NULL, carryover_value=FALSE)
f_value2 <- vm2$fit() #estimate learning rate as a free parameter

d <- build_design_matrix(fitobj=f_value2, regressors=c("clock", "feedback", "ev", "rpe_pos", "rpe_neg"), 
    event_onsets=c("clock_onset", "feedback_onset", "clock_onset", "feedback_onset", "feedback_onset"), 
    durations=c("clock_duration", 0, "clock_duration", 0, 0), 
    normalizations=c("durmax_1.0", "evtmax_1.0", "durmax_1.0", "evtmax_1.0", "evtmax_1.0"), center_values=TRUE, convolve_wi_run=TRUE,
    baselineCoefOrder=2, runVolumes=subj_runlengths, dropVolumes = 6)

res <- incr_fit(concatDesign(d, hpass=.01), concatMR, mask3d=voxAnalyze, modelname="10_ev_rpepos_rpeneg_evtmax_wirunscaling_nocarryover", outputPvals=FALSE, njobs=20)


##look at model with temporal derivatives added in each step



####LOOK AT ICA
#try badre model
#center parametric regressors, but place ev at feedback onset as impulse event
#d <- build_design_matrix(fitobj=f, regressors=c("clock", "mean_uncertainty", "rel_uncertainty", "ev", "feedback", "rpe_pos", "rpe_neg"), 
#    event_onsets=c("clock_onset", "clock_onset", "clock_onset", "clock_onset", "feedback_onset", "feedback_onset", "feedback_onset"), 
#    durations=c("clock_duration", "clock_duration", "clock_duration", "clock_duration", 0, 0, 0), 
#    normalizations=c("none", "none", "none", "none", "none", "none", "none"), center_values=FALSE, convolve_wi_run=TRUE,
#    baselineCoefOrder=2, runVolumes=subj_runlengths)

#each event gets a 1.0 scaling prior to convolution (such that longer events are not thought to have a bigger effect on neural activation
d <- build_design_matrix(fitobj=f, regressors=c("clock", "mean_uncertainty", "rel_uncertainty", "ev", "feedback", "rpe_pos", "rpe_neg"), 
    event_onsets=c("clock_onset", "clock_onset", "clock_onset", "clock_onset", "feedback_onset", "feedback_onset", "feedback_onset"), 
    durations=c("clock_duration", "clock_duration", "clock_duration", "clock_duration", 0, 0, 0), 
    normalizations=c("durmax_1.0", "evtmax_1.0", "evtmax_1.0", "evtmax_1.0", "durmax_1.0", "none", "none"), center_values=FALSE, convolve_wi_run=TRUE,
    baselineCoefOrder=2, runVolumes=subj_runlengths)

d_allruns <- concatDesign(d, hpass=.01)

corrs_evtmax <- corwithtarget(cbind(d_allruns[,!grepl("run[0-9]base.*", dimnames(d_allruns)[[2]]) ], ic_time[,goodics]), omit=NULL, 
    target=c("clock", "feedback", "ev", "mean_uncertainty", "rel_uncertainty", "rpe_pos", "rpe_neg"), absrmin=.1)

#same but with all-run (between) scaling
d <- build_design_matrix(fitobj=f, regressors=c("clock", "mean_uncertainty", "rel_uncertainty", "ev", "feedback", "rpe_pos", "rpe_neg"), 
    event_onsets=c("clock_onset", "clock_onset", "clock_onset", "clock_onset", "feedback_onset", "feedback_onset", "feedback_onset"), 
    durations=c("clock_duration", "clock_duration", "clock_duration", "clock_duration", 0, 0, 0), 
    normalizations=c("durmax_1.0", "evtmax_1.0", "evtmax_1.0", "evtmax_1.0", "durmax_1.0", "none", "none"), center_values=FALSE, convolve_wi_run=FALSE,
    baselineCoefOrder=2, runVolumes=subj_runlengths)

d_allruns <- concatDesign(d, hpass=.01)

corrs_evtmax_bw <- corwithtarget(cbind(d_allruns[,!grepl("run[0-9]base.*", dimnames(d_allruns)[[2]]) ], ic_time[,goodics]), omit=NULL, 
    target=c("clock", "feedback", "ev", "mean_uncertainty", "rel_uncertainty", "rpe_pos", "rpe_neg"), absrmin=.1)

#unified (signed) RPE regressor
d <- build_design_matrix(fitobj=f, regressors=c("clock", "mean_uncertainty", "rel_uncertainty", "ev", "feedback", "rpe"), 
    event_onsets=c("clock_onset", "clock_onset", "clock_onset", "clock_onset", "feedback_onset", "feedback_onset"), 
    durations=c("clock_duration", "clock_duration", "clock_duration", "clock_duration", 0, 0), 
    normalizations=c("durmax_1.0", "evtmax_1.0", "evtmax_1.0", "evtmax_1.0", "durmax_1.0", "none"), center_values=FALSE, convolve_wi_run=TRUE,
    baselineCoefOrder=2, runVolumes=subj_runlengths)

d_allruns <- concatDesign(d, hpass=.01)
print(visualizeDesignMatrix(d_allruns, includeBaseline = FALSE))

corrs_evtmax_rpeunified <- corwithtarget(cbind(d_allruns[,!grepl("run[0-9]base.*", dimnames(d_allruns)[[2]]) ], ic_time[,goodics]), omit=NULL, 
    target=c("clock", "feedback", "ev", "mean_uncertainty", "rel_uncertainty", "rpe"), absrmin=.01)

#look at 1/0 prediction error
d <- build_design_matrix(fitobj=f, regressors=c("clock", "mean_uncertainty", "rel_uncertainty", "ev", "feedback", "rpe_pos_bin", "rpe_neg_bin"), 
    event_onsets=c("clock_onset", "clock_onset", "clock_onset", "clock_onset", "feedback_onset", "feedback_onset", "feedback_onset"), 
    durations=c("clock_duration", "clock_duration", "clock_duration", "clock_duration", 0, 0, 0), 
    normalizations=c("durmax_1.0", "evtmax_1.0", "evtmax_1.0", "evtmax_1.0", "none", "none", "none"), center_values=FALSE, convolve_wi_run=TRUE,
    baselineCoefOrder=2, runVolumes=subj_runlengths)

d_allruns <- concatDesign(d, hpass=.01)

corrs_evtmax_rpebin <- corwithtarget(cbind(d_allruns[,!grepl("run[0-9]base.*", dimnames(d_allruns)[[2]]) ], ic_time[,goodics]), omit=NULL, 
    target=c("clock", "feedback", "ev", "mean_uncertainty", "rel_uncertainty", "rpe_pos_bin", "rpe_neg_bin"), absrmin=.2)

whitened<- corwithtarget(cbind(d_allruns[,!grepl("run[0-9]base.*", dimnames(d_allruns)[[2]]) ], ic_time), omit=NULL, 
    target=c("clock", "feedback", "ev", "mean_uncertainty", "rel_uncertainty", "rpe_pos_bin", "rpe_neg_bin"), absrmin=.3, prewhiten=TRUE)


ccf(d_allruns[,"rpe_pos_bin"], ic_time[,20])
ccf(d_allruns[,"rpe_pos_bin"], ic_time[,50])
ccf(d_allruns[,"rpe_pos_bin"], ic_time[,41])



#put ev at feedback
d <- build_design_matrix(fitobj=f, regressors=c("clock", "mean_uncertainty", "rel_uncertainty", "feedback", "rpe_pos_sqrt", "rpe_neg_sqrt"), 
    event_onsets=c("clock_onset", "clock_onset", "clock_onset", "feedback_onset", "feedback_onset", "feedback_onset"), 
    durations=c("clock_duration", "clock_duration", "clock_duration", 0, 0, 0), 
    normalizations=c("durmax_1.0", "evtmax_1.0", "evtmax_1.0", "none", "none", "none"), center_values=TRUE, convolve_wi_run=TRUE,
    baselineCoefOrder=2, runVolumes=subj_runlengths)

d_allruns <- concatDesign(d, hpass=.01)

corrs_evfeedback <- corwithtarget(cbind(d_allruns[,!grepl("run[0-9]base.*", dimnames(d_allruns)[[2]]) ], ic_time), omit=NULL, 
    target=c("rpe_pos_sqrt", "rpe_neg_sqrt"), partial="feedback", absrmin=.1)

corrs_nonorm$ev
corrs_evtmax$ev
corrs_evtmax_bw$ev

corrs_nonorm$rpe_pos
corrs_evtmax$rpe_pos
corrs_evtmax_bw$rpe_pos




d_dropvols <- build_design_matrix(fitobj=f, regressors=c("clock", "mean_uncertainty", "rel_uncertainty", "ev", "feedback", "rpe_pos_sqrt", "rpe_neg_sqrt"), 
    event_onsets=c("clock_onset", "clock_onset", "clock_onset", "clock_onset", "feedback_onset", "feedback_onset", "feedback_onset"), 
    durations=c("clock_duration", "clock_duration", "clock_duration", "clock_duration", 0, 0, 0), 
    normalizations=c("none", "none", "none", "none", "none", "none", "none"), center_values=FALSE, convolve_wi_run=TRUE,
    baselineCoefOrder=2, runVolumes=subj_runlengths, dropVolumes=6)

d_allruns_drop <- concatDesign(d, hpass=.01)


#ran 50 component ICA Melodic on concatenated runs
#here are the time courses -- correlate with variables of interest
ic_time <- read.table("melodic_concat/melodic_Tmodes")
names(ic_time) <- paste0("IC", 1:ncol(ic_time))
goodics <- c(1,2,6,7,8,10,11,13,16,18,20,24,25,28,31,35,37,38,40,41,42,43,45,50)

v <- ccf(d_allruns[,"clock"], ic_time[,38])
v2 <- ccf(d_allruns[,"rel_uncertainty"], ic_time[,38])

corwithtarget(cbind(d_allruns[,!grepl("run[0-9]base.*", dimnames(d_allruns)[[2]]) ], ic_time), omit=NULL, 
    target=c("clock", "feedback", "ev", "mean_uncertainty", "rel_uncertainty", "rpe_pos_sqrt", "rpe_neg_sqrt"), partial=c("clock", "feedback"), absrmin=.1)

corwithtarget(cbind(d_allruns[,!grepl("run[0-9]base.*", dimnames(d_allruns)[[2]]) ], ic_time), omit=NULL, 
    target=c("ev", "mean_uncertainty", "rel_uncertainty"), partial=c("clock"), absrmin=.1)

corwithtarget(cbind(d_allruns[,!grepl("run[0-9]base.*", dimnames(d_allruns)[[2]]) ], ic_time), omit=NULL, 
    target=c("rpe_pos", "rpe_neg"), partial=c("feedback"), absrmin=.1)


#ev at feedback


#g <- visualizeDesignMatrix(cbind(d_allruns, ic_time[,goodics]), outfile="mimic_badre_icscorr.pdf", runboundaries=subj_runlengths, events=d$concat_onsets[c("clock", "feedback")], includeBaseline = FALSE)
limit_vols <- 300 #nrow(ic_time) #1500
tr <- 1.0

runboundaries <- cumsum(subj_runlengths)*tr #timing of run onsets in seconds
#events <- d$concat_onsets[c("clock", "feedback")]
events <- d$concat_onsets[c("rpe_pos", "rpe_neg")]
#events <- d$concat_onsets["clock"]

runboundaries <- runboundaries[which(runboundaries < tr*limit_vols)]
events <- lapply(events, function(v) { v[which(v < tr*limit_vols)] })

#g <- visualizeDesignMatrix(cbind(d_allruns, ic_time[,c(1,2,6,30)])[1:limit_vols,], outfile="mimic_badre_icscorr.pdf", runboundaries=runboundaries, events=events, includeBaseline = FALSE)
g <- visualizeDesignMatrix(ic_time[1:limit_vols,goodics], outfile="ica_eventmap.pdf", runboundaries=runboundaries, events=events, includeBaseline = FALSE)

ic_time[as.vector(sapply(runboundaries, function(x) { (x-10):(x+10) })),30]


#junk
#fitARMA is slower than arima (despite advertising its speed)

#choose ARMA(3,0,3)
#  arm <- FitARMA(as.vector(vox), order = c(1, 0, 0), demean = TRUE, MeanMLEQ = FALSE, pApprox = 30, MaxLag = 30)
#  arm <- FitARMA(as.vector(vox), order = c(1, 0, 1), demean = TRUE, MeanMLEQ = FALSE, pApprox = 30, MaxLag = 30) 
#  arm <- FitARMA(as.vector(vox), order = c(2, 0, 1), demean = TRUE, MeanMLEQ = FALSE, pApprox = 30, MaxLag = 30)
#  arm <- FitARMA(as.vector(vox), order = c(3, 0, 1), demean = TRUE, MeanMLEQ = FALSE, pApprox = 30, MaxLag = 30)
#  arm <- FitARMA(as.vector(vox), order = c(6, 0, 1), demean = TRUE, MeanMLEQ = FALSE, pApprox = 30, MaxLag = 30)
#  arm <- FitARMA(as.vector(vox), order = c(6, 0, 2), demean = TRUE, MeanMLEQ = FALSE, pApprox = 30, MaxLag = 30)
#  arm <- FitARMA(as.vector(vox), order = c(10, 1, 0), demean = TRUE, MeanMLEQ = FALSE, pApprox = 30, MaxLag = 30)
#  arm <- FitARMA(as.vector(vox), order = c(13, 1, 0), demean = TRUE, MeanMLEQ = FALSE, pApprox = 30, MaxLag = 50)
#  arm <- FitARMA(as.vector(vox), order = c(20, 1, 0), demean = TRUE, MeanMLEQ = FALSE, pApprox = 30, MaxLag = 50)
#  arm <- FitARMA(as.vector(vox), order = c(15, 1, 0), demean = TRUE, MeanMLEQ = FALSE, pApprox = 30, MaxLag = 50)
#  arm <- FitARMA(as.vector(vox), order = c(16, 1, 0), demean = TRUE, MeanMLEQ = FALSE, pApprox = 30, MaxLag = 50)
#  white <- residuals(arm)

#  system.time(arm <- FitARMA(as.vector(vox), order = c(5, 0, 4), demean = TRUE, MeanMLEQ = FALSE, pApprox = 30, MaxLag = 30))
#  vdemean <- as.vector(vox) - mean(vox)
#  system.time(arm <- GetFitARMA(vdemean, p=2, q=2))
#  system.time(arima(vox[1,], c(2, 0, 2)))
#  
#  test <- ar.mle(as.vector(vox), aic=TRUE)
#  system.time(auto <- auto.arima(as.vector(vox)))
#  system.time(autod <- auto.arima(vdemean))
#  arima(vox[1,], c(5, 0, 4))


## x11()
## visualizeDesignMatrix(d$design.convolve[[1]])

## dmatFilter <- as.data.frame(lapply(d$design.convolve[[1]], function(r) lmBandpass(r, dt=1.0, .01, 10) ))

## x11()
## visualizeDesignMatrix(dmatFilter)


##
load("11275_fitinfo.RData")

plotRTs=function(f) {
  rtDf <- data.frame(
      trial=rep(1:ncol(f$RTobs), nrow(f$RTobs)),
      run=rep(1:nrow(f$RTobs), each=ncol(f$RTobs)),
      run_condition=rep(f$run_condition, each=ncol(f$RTobs)), 
      rew_function=rep(f$rew_function, each=ncol(f$RTobs)),
      reward=gdata::unmatrix(f$Reward, byrow=TRUE),
      rt=c(gdata::unmatrix(f$RTobs, byrow=TRUE), gdata::unmatrix(f$RTpred, byrow=TRUE)),
      rt_type=rep(c("observed", "predicted"), each=length(f$RTobs))
  )
  
  rtPlot <- ggplot(rtDf, aes(x=trial, y=rt, color=rt_type))
  
  rtPlot <- rtPlot + geom_line(size=1.1) + ggtitle(paste(rtDf$rew_function, rtDf$run_condition, sep=", ")) +
      theme_bw(base_size=14) + facet_grid(run_condition ~ rew_function)
  
  print(rtPlot)
  return(invisible(rtPlot))
  
}

plotRTs(f)

#look at problems with collinearity in fit
d <- build_design_matrix(fitobj=f, regressors=c("clock", "mean_uncertainty", "rel_uncertainty", "ev", "feedback", "rpe_pos", "rpe_neg"), 
    event_onsets=c("clock_onset", "clock_onset", "clock_onset", "clock_onset", "feedback_onset", "feedback_onset", "feedback_onset"), 
    durations=c("clock_duration", "clock_duration", "clock_duration", "clock_duration", 0, 0, 0), 
    normalizations=c("durmax_1.0", "evtmax_1.0", "evtmax_1.0", "evtmax_1.0", "durmax_1.0", "evtmax_1.0", "evtmax_1.0"),
    center_values=TRUE, convolve_wi_run=TRUE, baselineCoefOrder=2,  runsToOutput=1:8,
    dropVolumes = 6, high_pass=.01)

#look at problems with collinearity in fit
d <- build_design_matrix(fitobj=f, regressors=c("clock", "mean_uncertainty", "rel_uncertainty", "ev", "feedback", "rpe_pos", "rpe_neg"), 
    event_onsets=c("clock_onset", "clock_onset", "clock_onset", "clock_onset", "feedback_onset", "feedback_onset", "feedback_onset"), 
    durations=c("clock_duration", "clock_duration", "clock_duration", "clock_duration", 0, 0, 0), 
    normalizations=rep("none", 7),
    center_values=TRUE, convolve_wi_run=TRUE, baselineCoefOrder=2,  runsToOutput=1:8,
    dropVolumes = 6, high_pass=.01)

d$collin.convolve$run8
d$collin.raw$run8

#cbind mean and relative uncertainty parametric values before convolution
rel_mean <- cbind(d$design[["run8", "mean_uncertainty"]][,"value"], d$design[["run8", "rel_uncertainty"]][,"value"])
par(mfrow=c(2,1))
plot(1:50, rel_mean[,1], type="l", main="relative uncertainty")
plot(1:50, rel_mean[,2], type="l", main="mean uncertainty")



#regressors without mean centering prior to convolution (horrid!)
d <- build_design_matrix(fitobj=f, regressors=c("clock", "mean_uncertainty", "rel_uncertainty", "ev", "feedback", "rpe_pos", "rpe_neg"), 
    event_onsets=c("clock_onset", "clock_onset", "clock_onset", "clock_onset", "feedback_onset", "feedback_onset", "feedback_onset"), 
    durations=c("clock_duration", "clock_duration", "clock_duration", "clock_duration", 0, 0, 0), 
    normalizations=rep("none", 7),
    center_values=FALSE, convolve_wi_run=TRUE, baselineCoefOrder=2,  runsToOutput=1:8,
    dropVolumes = 6, high_pass=.01)

d$collin.convolve$run8

#overall centering of parametric values prior to convolution
d <- build_design_matrix(fitobj=f, regressors=c("clock", "mean_uncertainty", "rel_uncertainty", "ev", "feedback", "rpe_pos", "rpe_neg"), 
    event_onsets=c("clock_onset", "clock_onset", "clock_onset", "clock_onset", "feedback_onset", "feedback_onset", "feedback_onset"), 
    durations=c("clock_duration", "clock_duration", "clock_duration", "clock_duration", 0, 0, 0), 
    normalizations=rep("none", 7),
    center_values=TRUE, convolve_wi_run=FALSE, baselineCoefOrder=2,  runsToOutput=1:8,
    dropVolumes = 6, high_pass=.01)

d$collin.convolve$run8



#duration max 1.0 normalization
d <- build_design_matrix(fitobj=f, regressors=c("clock", "mean_uncertainty", "rel_uncertainty", "ev", "feedback", "rpe_pos", "rpe_neg"), 
    event_onsets=c("clock_onset", "clock_onset", "clock_onset", "clock_onset", "feedback_onset", "feedback_onset", "feedback_onset"), 
    durations=c("clock_duration", "clock_duration", "clock_duration", "clock_duration", 0, 0, 0), 
    normalizations=rep("durmax_1.0", 7),
    center_values=TRUE, convolve_wi_run=TRUE, baselineCoefOrder=2,  runsToOutput=1:8,
    dropVolumes = 6, high_pass=.01)

d$collin.convolve$run8
