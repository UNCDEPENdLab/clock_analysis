#script to fit behavior using fitclock package.
#also compute more basic statistics on the response times outside of a formal model
library(abind)
library(fitclock)

#help with Tomasz fitting
s <- clockdata_subject(subject_ID='test', dataset='/Users/michael/Downloads/o1_formatted.csv')
f <- clock_fit()
f$populate_fit(s)

options(warn=2)
poseps <- clock_model()
poseps$add_params(
    meanRT(max_value=4000),
    autocorrPrevRT(),
    goForGold(),
    go(),
    noGo(),
    meanSlowFast(),
    exploreBeta()
)

poseps$set_data(s)
poseps$carryover_value <- TRUE

#incr_fit_poseps <- poseps$incremental_fit(njobs=7, plot=FALSE)

f_poseps <- poseps$fit(random_starts=4)

negeps <- clock_model()
negeps$add_params(
    meanRT(max_value=4000),
    stickyChoice(),
    go(),
    noGo(),
    meanSlowFast(),
    exploreBeta(min_value=c(epsilon=-100000)) #allow negative epsilon
)

negeps$set_data(s)
f_negeps <- negeps$fit(random_starts=4)

expDiff_model <- clock_model(fit_RT_diffs=TRUE)
expDiff_model$add_params(
    meanRT(max_value=4000),
    autocorrPrevRT(),
    goForGold(),
    go(),
    noGo(),
    meanSlowFast(),
    exploreBeta()
)

behavDir=file.path(getMainDir(), "clock_analysis", "fmri", "behavior_files")
setwd(file.path(getMainDir(), "clock_analysis", "fmri", "fmri_fits"))

behavFiles <- list.files(path=behavDir, pattern="*tcExport.csv", full.names=TRUE, recursive=TRUE)

#prototype of single subject clock_fit without parameter fitting
library(fitclock)
library(ggplot2)
s <- clockdata_subject(subject_ID='test', dataset=behavFiles[1])
f <- clock_fit()
f$populate_fit(s)

#create design matrix object d that contains convolved and unconvolved variants of the design (also write out FSL and AFNI-style stimulus timing) 
d <- build_design_matrix(fitobj=f, regressors=c("clock", "feedback"), 
    event_onsets=c("clock_onset", "feedback_onset"), 
    durations=c("clock_duration", "feedback_duration"),
    normalizations=c("durmax_1.0", "durmax_1.0"),
    baselineCoefOrder=2, center_values=TRUE, convolve_wi_run = TRUE,
    writeTimingFiles=c("AFNI", "FSL"), output_directory="timing_files", high_pass=.01, tr=1.0)

#create an AFNI-style concatenated design matrix with unique baselines per run (for visualization)
d_concat <- concatDesignRuns(d)
#g <- visualizeDesignMatrix(d_concat, runboundaries=cumsum(d$runVolumes), events=d$concat_onsets[c("clock", "feedback")], includeBaseline = FALSE) #include clock and feedback onsets
g <- visualizeDesignMatrix(d_concat, runboundaries=cumsum(d$runVolumes), events=d$concat_onsets["clock"], includeBaseline = FALSE) #just clock onsets
ggsave(filename="design_visualization.pdf", plot=g, width=50, height=9, limitsize=FALSE) #save a (large) pdf of the design matrix

## 1) Fit subject behavior using Frank TC model with positive and negative epsilon variants
##				- Also fit variants where expected value does not carry over across run
##				- Also fit each run individually

poseps <- clock_model()
poseps$add_params(
    meanRT(max_value=4000),
    autocorrPrevRT(),
    goForGold(),
    go(),
    noGo(),
    meanSlowFast(),
    exploreBeta()
)

negeps <- clock_model()
negeps$add_params(
    meanRT(max_value=4000),
    stickyChoice(),
    go(),
    noGo(),
    meanSlowFast(),
    exploreBeta(min_value=c(epsilon=-100000)) #allow negative epsilon
)


for (b in behavFiles) {
  #example location of file on bea_res, which contains scan date
  #/Volumes/bea_res/Data/Tasks/EmoClockfMRI/Basic/11229/20140521/Raw/fMRIEmoClock_11229_tc_tcExport.csv
  subid <- sub("^.*fMRIEmoClock_(\\d+)_tc_tcExport.csv$", "\\1", b, perl=TRUE)
  
  ##setup clock data subject object for fitting 
  s <- clockdata_subject(subject_ID=subid, dataset=b)
  
  if (file.exists(paste0(subid, "_fitinfo.RData"))) { 
    cat("Fit data already present for: ", subid, "\n")
    load(paste0(subid, "_fitinfo.RData"))
  } else {           
    cat("Fitting behavioral data for subject: ", subid, "\n")
    
    ##set data for model fit
    poseps$set_data(s)
    poseps$carryover_value <- TRUE
    negeps$set_data(s)
    negeps$carryover_value <- TRUE
    
    incr_fit_poseps <- poseps$incremental_fit(njobs=7, plot=FALSE)
    incr_fit_negeps <- negeps$incremental_fit(njobs=7, plot=FALSE)
    
    png(file.path(paste0(subid, "_incrfit.png")), width=9, height=6, units="in", res=300)
    print(incr_fit_poseps$AICplot)
    dev.off()
    
    png(file.path(paste0(subid, "_incrfit_negeps.png")), width=9, height=6, units="in", res=300)
    print(incr_fit_negeps$AICplot)
    dev.off()
    
    f_poseps <- poseps$fit(random_starts=40)
    f_negeps <- negeps$fit(random_starts=40)
    
    #run-wise fits
    f_runs <- lapply(s$runs, function(r) {
          poseps$fit(toFit=r, random_starts=5)
        })
    
    #refit without value carryover
    poseps$carryover_value <- FALSE
    negeps$carryover_value <- FALSE
    f_nocarryover <- poseps$fit(random_starts=40)
    f_nocarryover_negeps <- negeps$fit(random_starts=40)
    
    ##delta rule value model (simple)
    vm <- deltavalue_model(clock_data=s, alphaV=0.3, betaV=0.3) #N.B. This matches V matrix from full time-clock algorithm fit.
    f_value <- vm$fit() #estimate learning rate as a free parameter
    
    vm$carryover_value <- FALSE
    f_value_nocarryover <- vm$fit()
    
    save(s, f_value, f_value_nocarryover, f_poseps, f_nocarryover, f_negeps, f_nocarryover_negeps, f_runs, s, incr_fit_poseps, incr_fit_negeps, file=paste0(subid, "_fitinfo.RData"))
  }
  
}

## new model with an increment parameter, tau, that is positive for happy blocks, 0 for scrambled, and negative for fear.
epstau <- clock_model()
epstau$add_params(
    meanRT(max_value=4000),
    stickyChoice(),
    go(),
    noGo(),
    meanSlowFast(),
    exploreBeta(min_value=c(epsilon=0, tau=0))
)

for (b in behavFiles) {
  #example location of file on bea_res, which contains scan date
  #/Volumes/bea_res/Data/Tasks/EmoClockfMRI/Basic/11229/20140521/Raw/fMRIEmoClock_11229_tc_tcExport.csv
  subid <- sub("^.*fMRIEmoClock_(\\d+)_tc_tcExport.csv$", "\\1", b, perl=TRUE)
  
  ##setup clock data subject object for fitting 
  s <- clockdata_subject(subject_ID=subid, dataset=b)
  
  if (file.exists(paste0(subid, "_taufitinfo.RData"))) { 
    cat("Fit data already present for: ", subid, "\n")
    load(paste0(subid, "_taufitinfo.RData"))
  } else {           
    cat("Fitting behavioral data for subject: ", subid, "\n")
    
    ##set data for model fit
    epstau$set_data(s)
    epstau$carryover_value <- TRUE
    
    #incr_fit <- poseps$incremental_fit(njobs=7, plot=FALSE)  
    f <- epstau$fit(random_starts=40)
    
    save(s, f, file=paste0(subid, "_taufitinfo.RData"))
  }
}




## new model with an increment parameter, tau, that is positive for happy blocks, 0 for scrambled, and negative for fear.
sticky <- clock_model()
sticky$add_params(
    meanRT(max_value=4000),
    stickyChoice(),
    go(),
    noGo(),
    meanSlowFast(),
    exploreBeta(min_value=c(epsilon=0))
)

for (b in behavFiles) {
  #example location of file on bea_res, which contains scan date
  #/Volumes/bea_res/Data/Tasks/EmoClockfMRI/Basic/11229/20140521/Raw/fMRIEmoClock_11229_tc_tcExport.csv
  subid <- sub("^.*fMRIEmoClock_(\\d+)_tc_tcExport.csv$", "\\1", b, perl=TRUE)
  
  ##setup clock data subject object for fitting 
  s <- clockdata_subject(subject_ID=subid, dataset=b)
  
  if (file.exists(paste0(subid, "_stickyfitinfo.RData"))) { 
    cat("Fit data already present for: ", subid, "\n")
    load(paste0(subid, "_stickyfitinfo.RData"))
  } else {           
    cat("Fitting behavioral data for subject: ", subid, "\n")
    
    ##set data for model fit
    sticky$set_data(s)
    sticky$carryover_value <- TRUE
    
    #incr_fit <- poseps$incremental_fit(njobs=7, plot=FALSE)  
    f <- sticky$fit(random_starts=40)
    
    save(s, f, file=paste0(subid, "_stickyfitinfo.RData"))
  }
}





#get AIC params for tau model for comparison to standard model
tauobjs <- list.files(pattern=".*_taufitinfo.RData")

options(error=recover)
taudf <- c()
for (o in tauobjs) {
  load(o)
  id <- sub("(\\d+)_taufitinfo\\.RData", "\\1", o, perl=TRUE)
  thisguy <- data.frame(id=factor(id), tauAIC=f$AIC, tauSSE=f$SSE, as.data.frame(t(f$theta[,"cur_value"])))
  taudf <- rbind(taudf, thisguy)
}

#regular fits (with sticky choice)
regobjs <- list.files(pattern=".*_stickyfitinfo.RData")

options(error=recover)
regdf <- c()
for (o in regobjs) {
  load(o)
  id <- sub("(\\d+)_stickyfitinfo\\.RData", "\\1", o, perl=TRUE)
  thisguy <- data.frame(id=factor(id), regAIC=f$AIC, regSSE=f$SSE, as.data.frame(t(f$theta[,"cur_value"])))
  regdf <- rbind(regdf, thisguy)
}

both <- merge(taudf, regdf, by="id")

t.test(both$tauAIC, both$regAIC)
t.test(both$tauSSE, both$regSSE)

lattice::histogram(both$tau)

cor.test(taudf$tau, taudf$epsilonBeta)
cbind(both$tau, both$epsilonBeta.x)

psych::phi(table(taudf$epsilonBeta > 0, taudf$tau > 0))
prop.table(table(taudf$epsilonBeta > 0, taudf$tau > 0), margin=1)
table(taudf$epsilonBeta > 0)

chisq.test(table(taudf$epsilonBeta > 0, taudf$tau > 0))

both$aicDiff <- both$regAIC - both$tauAIC
cor.test(~tau + aicDiff, subset(both, tau>0))


cor.test(~tau + aicDiff, subset(both, tau>0 & tau < 8000)) #8000 is a bit of an outlier

xyplot(tau ~ aicDiff, subset(both, tau>0 & tau < 8000))

cor.test(~tau + tauAIC, subset(taudf, tau > 0))


#SPM BMS approach: need to generate a subjects x models matrix
library(R.matlab)
AICmat <- -1*both[,c("regAIC", "tauAIC")]

#writeMat(con="AICmatrix_n36.mat", AICmat=AICmat, mnames=levels(allM$model))
writeMat(con="AICmatrix_tauvsregular.mat", AICmat=as.matrix(AICmat), mnames=c("posEpsSticky", "tauSticky"))

##run SPM BMS in MATLAB
#system("matlab -nodisplay < computeBMSprobs.m")

options(matlab="/Applications/MATLAB_R2015a.app/bin/matlab")

Matlab$startServer()
matlab <- Matlab()
isOpen <- open(matlab)
if (!isOpen) { throw("MATLAB server is not running: waited 30 seconds.") }
print(matlab)

setVariable(matlab, AICmat=as.matrix(AICmat))
evaluate(matlab, "addpath('/Users/michael/Data_Analysis/clock_analysis/spm_BMS'); \
[alpha, expost_r, xp] = spm_BMS(AICmat, 1e6, 1);")

alpha <- getVariable(matlab, "alpha")[[1L]]
expost_r <- getVariable(matlab, "expost_r")[[1L]]
xp <- getVariable(matlab, "xp")[[1L]]

close(matlab)






## 2) Allow for emotion-varying exploration, fitting controls only

fitobjs <- list.files(pattern="\\d+_*fitinfo.RData")
#cw behavior looks pretty ugly... increasing AIC for more parameters, no go for gold or beyond
#fitobjs <- fitobjs[!which(fitobjs=="15_fitinfo.RData")]
subids <- as.integer(sub("^(\\d+)_fitinfo.RData", "\\1", fitobjs, perl=TRUE))
bpdsub <- as.integer(subids < 10000)
control_fitobjs <- fitobjs[subids > 10000] #only keep controls from MMY3

if (!file.exists("fits_epsilon_emo_controls_5May2015.RData")) {
  posEps <- clock_model()
  posEps$add_params(
      meanRT(max_value=4000),
      autocorrPrevRT(),
      goForGold(),
      go(),
      noGo(),
      meanSlowFast(),
      exploreBeta(by="run_condition")   #allow for emotion-varying exploration
  )
  
  epsemo_controls <- list()
  for (fpos in 1:length(control_fitobjs)) {
    load(control_fitobjs[fpos])
    posEps$set_data(s)
    
    fitobj <- posEps$fit(random_starts=5)
    epsemo_controls[[fpos]] <- fitobj 
  }
  
  save(epsemo_controls, file="fits_epsilon_emo_controls_5May2015.RData")
}

#NB. Getting non-convergence for several subjects
if (!file.exists("fits_negepsilon_emo_controls_1Jun2015.RData")) {
  
  #allow for exploration by emotion and negative values
  negEps <- clock_model()
  negEps$add_params(
      meanRT(max_value=4000),
      stickyChoice(),
      go(),
      noGo(),
      meanSlowFast(),
      exploreBeta(min_value=-100000, by="run_condition")
  )
  
  negepsemo_controls <- list()
  for (fpos in 1:length(control_fitobjs)) {
    load(control_fitobjs[fpos])
    negEps$set_data(s)
    
    fitobj <- negEps$fit(random_starts=5)
    negepsemo_controls[[fpos]] <- fitobj 
  }
  
}



allfits <- get_fit_array(fitobjs)

## 3) Examine effects of a) max value, b) prior RT change, and c) prior omission on RTs.
##
##		- General idea is to extract important features of response time behavior without a heavy duty
##		  computational model. 
##		- From basis analyses, it's clear that people tend to speed up if their prior RT
##      was slower than average (ala Moustafa) and slow down if prior RT was faster than average.
##		- Likewise, prior reward omissions (here treated as 0 points, as opposed to an NPE based on model)
##		  tend to drive faster RTs, whereas positive RPEs lead to slower subsequent RTs.


###seems like we should look at behavior on trial t + 1 after RPE + versus RPE -
#how do these reward prediction errors influence subsequent behavior?

alldf <- c()
ids <- as.integer(sub("(\\d+)_fitinfo.RData", "\\1", fitobjs, perl=TRUE))
#evaluate Gaussian function (mean=mu, sd=sigma) for a vector x 
gausseval <- function(x, mu, sigma) {
  exp(-(x-mu)^2/(2*(sigma*length(x))^2))    
}

for (i in 1:length(control_fitobjs)) {
  #cat("Loading object: ", fitobjs[i], "\n")
  loc <- local({load(control_fitobjs[i]); environment()})$f #time-clock fit object
  
  rewards <- gdata::unmatrix(loc$Reward, byrow=TRUE)
  
  #compute a basic value function from the trialwise RTs
  #this is done by tracking EV for each trial and timestep (so, 50 x 4000 for a run)
  #we take the obtained reward and spread its effect on EV using a Gaussian generalization function (ala Skeptic)
  valuefunc <- lapply(1:nrow(loc$RTraw), function(run) {
        v_it <- array(0, dim=c(dim(loc$RTraw)[2], 4000), dimnames=list(trials=NULL, timesteps=NULL)) #runs=NULL,         
        for (t in 1:length(loc$RTraw[run,])) {
          elig <- gausseval(1:dim(v_it)[2], mu=loc$RTraw[run, t], sigma=.05)
          
          thisRew <- elig*loc$Reward[run, t]
          if (t > 1) v_it[t,] <- 0.9*v_it[(t-1),] + thisRew #0.9 is the decay (10%) rate for learned value
          else v_it[t,] <- thisRew
          #cat("t: ", t, "\n")
          #plot(1:4000, v_vec, type="l")          
        }
        
        v_it
      })
  
  valuefunc <- do.call(abind, list(along=0, valuefunc))
  maxv <- apply(valuefunc, c(1,2), which.max)
  maxv[which(maxv == 1)] <- NA 
  
  #same idea here as above, but use a decaying uncertainty function ala the old logistic operator
  #in practice, this is not working well at all. need to debug/improve
  uncfunc <- lapply(1:nrow(loc$RTraw), function(run) {
        u_it <- array(0, dim=c(dim(loc$RTraw)[2], 4000), dimnames=list(trials=NULL, timesteps=NULL)) #runs=NULL,         
        for (t in 1:length(loc$RTraw[run,])) {
          elig <- gausseval(1:dim(u_it)[2], mu=loc$RTraw[run, t], sigma=.01)
          #thisRew <- smth(thisRew, window=0.25, method="gaussian", tails=TRUE) #SLOW!
          if (t > 1) u_it[t,] <- 0.9*u_it[(t-1),] - elig #0.9 is the decay (10%) rate for learned value
          else u_it[t,] <- elig
          #cat("t: ", t, "\n")
          u_it[t,1:200] <- u_it[t,201]
          u_it[t,3800:4000] <- u_it[t,3799]
          #plot(1:4000, u_it[t,], type="l")
          #browser()
          #Sys.sleep(0.1)
          
        }
        #stop("die")
        u_it
      })
  
  uncfunc <- do.call(abind, list(along=0, uncfunc))
  maxu <- apply(uncfunc[,,200:3800], c(1,2), which.max) + 200
  #maxu[which(maxu == 1)] <- NA 
  
  #what if we used ecdf of uncertainty to quantify how uncertain the choice was? (rather than hard max)
  #ecdf function will ignore constant tails in computation, which is nice
  ulookup <- matrix(0, nrow(loc$RTraw), ncol(loc$RTraw))
  for (j in 1:nrow(loc$RTraw)) {
    for (k in 1:ncol(loc$RTraw)) {
      e <- ecdf(uncfunc[j,k,])
      ulookup[j,j] <- e(uncfunc[j,k,ifelse(loc$RTraw[j,k] > 4000, 4000, loc$RTraw[j,k]) ])      
    }
  }
  
  #should first RT be considered maximally uncertain? (1.0)
  #ef <- apply(uncfunc, c(1,2), ecdf)
  #rtunc <- apply(loc$RTraw, c(1,2), ef)
  
  ev <- gdata::unmatrix(loc$ev, byrow=TRUE)
  rpes <- gdata::unmatrix(loc$rpe, byrow=TRUE)
  rpes_lag <- gdata::unmatrix(t(apply(loc$rpe, 1, Hmisc::Lag, shift=1)), byrow=TRUE)
  rts <- gdata::unmatrix(loc$RTraw, byrow=TRUE)
  
  #compute deviation of current RT from max expected value RT 
  vdev <- t(sapply(1:nrow(maxv), function(r) {
            loc$RTraw[r,] - Hmisc::Lag(maxv[r,]) #maxv[r,] #lag the maxv here because we want to use max value that does not include this trial's outcome.
          }))
#  vdev2 <- t(sapply(1:nrow(maxv), function(r) {
#            loc$RTraw[r,] - maxv[r,] #old version where maxv(i) includes the reward(i)
#          }))
  
  #lag this variable to give us ability to look at how deviation from max value on prior trial
  #predicts RT adaptation on current trial
  vdevlag <- gdata::unmatrix(t(apply(vdev, 1, Hmisc::Lag, shift=1)), byrow=TRUE)
  
  udev <- t(sapply(1:nrow(maxu), function(r) {
            loc$RTraw[r,] - maxu[r,]
          }))
  
  udevlag <- gdata::unmatrix(t(apply(udev, 1, Hmisc::Lag, shift=1)), byrow=TRUE)
  maxu <- gdata::unmatrix(maxu, byrow=TRUE) #make into vector for including with data.frame
  maxv <- gdata::unmatrix(maxv, byrow=TRUE) 
  
  #replace low RTs with something more reasonable
  rts[which(rts < 100)] <- 100
  rt_lag <- gdata::unmatrix(t(apply(loc$RTraw, 1, Hmisc::Lag, shift=1)), byrow=TRUE)
  rtchange <- gdata::unmatrix(loc$RTraw - t(apply(loc$RTraw, 1, Hmisc::Lag, shift=1)), byrow=TRUE)
  overallrt <- mean(rts, na.rm=TRUE)
  runningrt <- t(apply(loc$RTraw, 1, function(run) {
            sapply(1:length(run), function(x) { mean(run[1:(x-1)], na.rm=TRUE) })
          }))
  runningrt[,1] <- NA #running RT undefined on first trial
  priorrtfast <- factor(as.integer(gdata::unmatrix(t(apply(loc$RTraw, 1, Hmisc::Lag, shift=1)) > overallrt, byrow=TRUE)), levels=c(0,1), labels=c("Prior RT < Avg", "Prior RT > Avg"))
  priorrtfast_running <- factor(as.integer(gdata::unmatrix(t(apply(loc$RTraw, 1, Hmisc::Lag, shift=1)) > runningrt, byrow=TRUE)), levels=c(0,1), labels=c("Prior RT < Avg", "Prior RT > Avg"))
  runningrtdev <- t(sapply(1:nrow(runningrt), function(r) {
            loc$RTraw[r,] - runningrt[r,] 
          }))
  
  runningrtdev_lag <- gdata::unmatrix(t(apply(runningrtdev, 1, Hmisc::Lag, shift=1)), byrow=TRUE)
  runningrtdev <- gdata::unmatrix(runningrtdev, byrow=TRUE) #unmatrix for data.frame
  
  trial <- 1:length(rts)
  cond <- rep(loc$run_condition, each=50)
  contingency <- rep(loc$rew_function, each=50)
  omission <- factor(as.integer(rewards==0), levels=c(0,1), labels=c("Rew delivered", "Rew omitted")) 
  omission_lag <- factor(as.integer(gdata::unmatrix(t(apply(loc$Reward==0, 1, Hmisc::Lag, shift=1)), byrow=TRUE)), levels=c(0,1), labels=c("PriorRew", "PriorOmission"))
  df <- data.frame(lunaid=ids[i], trial, contingency, cond, rts, rt_lag, rtchange, rpes, rpes_lag, rewards, ev, omission, omission_lag, priorrtfast, priorrtfast_running, vdevlag, udevlag, maxu, maxv,
      runningrtdev, runningrtdev_lag)
  alldf <- rbind(alldf, df)
}

alldf$rpes_lag_pos <- factor(as.numeric(alldf$rpes_lag > 0), levels=c(0,1), labels=c("Prior RPE-", "Prior RPE+"))
alldf$rpes_lag_abs <- abs(alldf$rpes_lag)


save(alldf, file="allRTs_withLaggedVars_1Jun2015.RData")
