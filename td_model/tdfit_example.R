#play around with fitting TD using example subject
setwd(file.path(getMainDir(), "clock_analysis", "td_model")) #alex: you'll need to run setwd(file.path("<directory where you did the git clone>", "clock_analysis", "td_model"))
source("tdrt_model.R") #bring td_fit function into workspace

#testing temporal difference RT model

if (!require(fitclock)) {
  if (!require(devtools)) { install.packages("devtools"); require(devtools) }
  install_github("fitclock", "LabNeuroCogDevel", args="--byte-compile")
  library(fitclock)
}

#example clock data object (part of fitclock)
exSubj = clockdata_subject(subject_ID="008", dataset=clocksubject_fMRI_008jh)

#print names of contingencies across runs
sapply(exSubj$runs, "[[", "rew_function")

#print names of emotions across runs
sapply(exSubj$runs, "[[", "run_condition")

fit_all <- list()
#loop over and fit runs
for (r in exSubj$runs) {
  rts = round(r$RTraw/10) #round to 10ms bins
  rewards = r$Reward
  condition = paste(r$rew_function, r$run_condition, sep=".") #name of condition
  #fit_all[[condition]]= td_fit(rewards, rts, cs_times=c(1), ntimesteps=400) #one CS at onset
  fit_all[[condition]]= td_fit(rewards, rts, cs_times=c(1, 100, 200, 300), ntimesteps=400) #add second CS at 1000ms
  
}

#names of fitted objects
names(fit_all)

fit_all$DEV.scram$ms_basisplot #plot basis functions


replay_fit_plot(fit_all$DEV.scram, fps=2.5) #replay DEV scrambled fit at 2.5 frames per second
replay_fit_plot(fit_all$IEV.scram, fps=3.5) #replay IEV scrambled fit at 3.5 frames per second

replay_fit_plot(fit_all$IEV.scram, fps=3.5, display="value") #replay value for IEV scrambled fit at 3.5 frames per second
replay_fit_plot(fit_all$IEV.scram, fps=2.5, display="delta") #replay td errors for IEV scrambled fit at 3.5 frames per second
replay_fit_plot(fit_all$IEV.scram, fps=2.5, display="action") #replay action for IEV scrambled fit at 3.5 frames per second

replay_fit_plot(fit_all$IEV.scram, fps=2.5, display="uncertainty") #replay action for IEV scrambled fit at 3.5 frames per second





### simulation

rts <- runif(200, 150, 3950) #100 samples between 150 and 3950 ms
rewards <- sapply(rts, getScore, scrfunc="IEV")
rts <- round(rts/10) #need to round to 10ms bins for current td fit

fitobj= td_fit(rewards, rts, cs_times=c(1, 100, 200, 300), ntimesteps=400) #add second CS at 1000ms
replay_fit_plot(fitobj, fps=5, display="action_value") #replay value for IEV scrambled fit at 3.5 frames per second
replay_fit_plot(fitobj, fps=10, display="value") #replay value for IEV scrambled fit at 3.5 frames per second

rts <- runif(300, 150, 3900) #100 samples between 150 and 3950 ms
rewards <- sapply(rts, getScore, scrfunc="DEV")
rts <- round(rts/10) #need to round to 10ms bins for current td fit

fitobj= td_fit(rewards, rts, cs_times=c(1, 100, 200, 300), ntimesteps=400) #add second CS at 1000ms
replay_fit_plot(fitobj, fps=5, display="action_value") #replay value for IEV scrambled fit at 3.5 frames per second
replay_fit_plot(fitobj, fps=10, display="value") #replay value for IEV scrambled fit at 3.5 frames per second



rts <- rep(3950, 200) #100 samples between 150 and 3950 ms
rewards <- sapply(rts, getScore, scrfunc="IEV")
rts <- round(rts/10) #need to round to 10ms bins for current td fit

fitobj= td_fit(rewards, rts, cs_times=c(1, 100, 200, 300), ntimesteps=400) #add second CS at 1000ms
replay_fit_plot(fitobj, fps=3.5, display="value") #replay value for IEV scrambled fit at 3.5 frames per second
replay_fit_plot(fitobj, fps=3.5, display="action_value") #replay value for IEV scrambled fit at 3.5 frames per second
