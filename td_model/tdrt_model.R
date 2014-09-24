replay_fit_plot <- function(outlist, fps=3, display="value") { #}, anim_loop=FALSE) {
  require(animation)
  
  #plot sequential fit across trials
  ani.options(interval=1/fps) #, loop=anim_loop) #nmax=anim_loops)
  if (display=="value") {
    ani.replay(outlist$vplot)  
  } else if (display=="action") {
    ani.replay(outlist$aplot)
  } else if (display=="delta") {
    ani.replay(outlist$dplot)
  } else if (display=="uncertainty") {
    ani.replay(outlist$uplot)
  } else if (display=="action_value") {
    ani.replay(outlist$avplot)
  }
  
  
}

td_fit <- function( 
    rewards=NULL,       #vector of reward magnitudes (US) received for each trial
    us_times=NULL,      #vector of US stimulus times (in time step units)
    cs_times=NULL,      #either vector or matrix of CS times (in time step units). If vector, then assume each elemtn is a consistent CS at a given time.
    ntimesteps=NULL, #number of time steps for estimating TD
    nbasis=12,          #number of microstimuli basis functions per stimulus
    alpha = 0.05,       #step-size    
    decay = 0.995,      #memory trace decay rate
    gamma = 0.99,       #discount factor
    lambda = 0.95,      #eligibility trace decay rate
    sigma = 0.04,       #microstimulus width
    theta = 0.25,       #action threshold
    upsilon = 0.9,      #action decay
    presence = 0.2      #presence microstimulus level
) {
  
  #adapted ludvig code
  #Simple Code for CSC/PR/MS TD Model with Actions.
  
  require(reshape2)
  require(ggplot2)
  require(animation)
  
  tracklist = list() #stuff to return
  
  ntrials = length(rewards)
  stopifnot(ntrials == length(us_times)) #rewards and US times must be the same length
  
  if (inherits(cs_times, "matrix")) {
    if (nrow(cs_times) == 1L) { #col vector (not suggested): replicate for i trials
      cs_times = matrix(rep(cs_times, ntrials), nrow=ntrials, byrow=TRUE)
    }
  } else {
    cs_times = matrix(rep(cs_times, ntrials), nrow=ntrials, byrow=TRUE) #replicate CS timing vector for number of trials
  }
  
  stim_times = cbind(cs_times, us_times) #us time is last column
  colnames(stim_times) = c(paste0("CS", 1:ncol(cs_times)), "US")
  
  nstimuli=ncol(stim_times) #how many stimuli are being tracked
  
  ## Initialize tracking matrices:
  x = matrix(data=0, nrow=nstimuli, ncol=nbasis)        #Stimulus representation matrix (nstimuli x nbasis)
  w = rep(0, nstimuli*nbasis)                           #Weight vector (length nstimuli x nbasis)
  w_uncertainty = rep(0, nstimuli*nbasis)               #Weight vector for uncertainty (length nstimuli x nbasis)
  e = rep(0, nstimuli*nbasis)                           #Eligibility Traces (length nstimuli x nbasis)
  
  #older vector-based representation (from Ludvig)
  #note the nbasis+1 is a carryover from the MATLAB code, where the last n elements are used for the presence representation
  #for ms basis, ncol=nbasis*nstimuli is correct.
  #x = matrix(data=0, nrow=1, ncol=(nbasis+1)*nstimuli)
  #w = matrix(data=0, nrow=1, ncol=(nbasis+1)*nstimuli)
  #e = matrix(data=0, nrow=1, ncol=(nbasis+1)*nstimuli)
  
  delta = matrix(data=0, nrow=ntrials, ncol=ntimesteps)    #TD Errors
  value = matrix(data=0, nrow=ntrials, ncol=ntimesteps)    #Value functions
  action = matrix(data=0, nrow=ntrials, ncol=ntimesteps)   #Response Levels
  uncertainty = matrix(data=0, nrow=ntrials, ncol=ntimesteps) #new matrix that is updated when subject samples
  ms = matrix(data=0, nrow=ntimesteps, ncol=nbasis)        #Microstimulus levels
  #maxaction = rep(0, ntrials)
  
  ## MS is the master matrix of basis functions for any given stimulus
  ## Microstimulus Generation (NIPS2009 Version)
  trace = 1
  for (timestep in 1:ntimesteps) {
    for (micro in 1:nbasis) {
      ms[timestep, micro] = trace * ((1/sqrt(2*pi))*exp(-((trace-(micro/nbasis))^2)/(2*(sigma^2))));    
    }
    trace = trace * decay;
  }
  
  #plot microstimuli
  ms_melt = melt(ms, varnames=c("timestep", "msnumber"), value.name="msheight")
  ms_melt$msnumber = factor(ms_melt$msnumber)
  tracklist$ms_basisplot = ggplot(ms_melt, aes(x=timestep, y=msheight, color=msnumber)) + geom_line()
  
  ##fit data
  
  rewdeliv = 0 #number of rewards delivered in fit (debugging)
  
  #matrices for tracking weights and eligibility traces for each timestep of a given trial
  wtrial = matrix(0, nrow=ntimesteps, ncol=nstimuli*nbasis)
  etrial = matrix(0, nrow=ntimesteps, ncol=nstimuli*nbasis)
  
  x_w_perTrial <- array(0, dim=c(ntrials, ntimesteps, nstimuli))
  
  #i trials
  #t timesteps
  #s stimuli
  for (i in 1:ntrials) {
    oldvalue = 0               # Reset on every trial
    oldvalue_uncertainty = 0
    #oldreward = 0;
    oldaction = 0
    e = rep(0, nstimuli*nbasis) #reset eligibility trace
    rewardtime = rts[i] #time of reward (RT in clock)
    
    msfigure <- list()
    
    for (t in 1:ntimesteps) {
      if (t == rewardtime) {
        rewdeliv = rewdeliv + 1
        reward = rewards[i]
        cat("rew count: ", rewdeliv, "\n")
        timesamp_reward = 1 #sampling gives a unit-magnitude reward
      } else {    
        reward = 0
        timesamp_reward = 0
      }
      
      #par(mfrow=c(5,1))
      
      #microstimulus representation of stimulus
      #x is updated for each stimulus at each timestep
      #so, it always represents the MS heights of a given stimulus relative to its onset time: stim_times[i, s]
      for (s in 1:nstimuli) {
        if (t >= stim_times[i, s]) {
          #noise = rand*0.02 -.01; (unused)
          noise = 0;
          x [s, ] = ms[t - stim_times[i, s] + 1, ] + noise
          #x [((s-1)*nbasis +1):((s-1)*nbasis +nbasis)] = ms[t-stim_times[i, s]+1,] + noise; #old vector-based update     
        } else {
          x [s, ] = 0 #0 stimulus representation prior to its onset time (may need to modify for symmetric update in time)
          #x [((s-1)*nbasis +1):((s-1)*nbasis +nbasis)] = 0;
          #x [nstimuli * nbasis + s] = 0; #this is updating the numms + 1 initial specification of x (19:21 in the 3 stim x 6 basis setup...). Don't get it. -- this is used by presence representation! 
        }
        #plot(1:ncol(x), x[s,], xlab="basis function", ylab=colnames(stim_times)[s])
        #text(x=ncol(x), y=max(x[s,]), paste0("t=", t)) #add trial number at
      }
      
      xvec = as.vector(x) #flatten x for computation of value, w, and e. (x as mat for ease of indexing/representation above)
      
      #if (i == 1 && t == 20) browser()
      
      value[i, t] = crossprod(xvec,w) #inner product across all stimuli
      
      w_mat <- matrix(w, nrow=nstimuli)
      x_w <- w_mat*x
#      for (s in 1:nstimuli) {
#        plot(1:ncol(x_w), x_w[s,], xlab="basis function", ylab=paste0("w*", colnames(stim_times)[s]))
#        text(x=ncol(x_w), y=max(x_w[s,]), paste0("t=", t)) #add trial number at
#      }
      
      x_w_perTrial[i,t,] <- rowSums(x_w)
            
      uncertainty[i, t] = crossprod(xvec, w_uncertainty)
      
      #Action Selection:
      action[i, t] = upsilon * oldaction + max (0, oldvalue - theta)
      oldaction=action[i, t]
      
      #Learning Algorithm:
#      if (t > rewardtime) {
#        t=ntimesteps + 1 #exit t
#        next
#        #delta[i, t] = 0
#        #delta_uncertainty=0
#      } else {
        delta[i, t] = reward + (gamma * value[i, t]) - oldvalue #TD Learning
        delta_uncertainty = timesamp_reward + (gamma * uncertainty[i, t]) - oldvalue_uncertainty
#      }
      
      w = w + (alpha * delta[i, t] * e)       #weight update
      w_uncertainty = w_uncertainty + (alpha * delta_uncertainty * e)
      e = xvec + (gamma * lambda * e)         #update eligibility trace
      wtrial[t, ] = w                         #DEBUG: track weight vector at this timestep
      etrial[t, ] = e                         #DEBUG: track eligibility trace at this timestep
      
      oldvalue = crossprod(xvec,w) #Or oldvalue = value. Difference in which weights are used.
      oldvalue_uncertainty = crossprod(xvec,w_uncertainty) #Or oldvalue = value. Difference in which weights are used.
      
      #cat("range of value: ", range(value), "\n")
      #oldreward = reward
    
      #msfigure[[t]] <- recordPlot()
    } #end timestep loop    
    
#    for(s in 1:nstimuli) {
#      plot(1:length(ntimesteps), sum(x_w[s,]))  
#    }
    
    
    #browser()
    #rewardtime
    #ani.replay(msfigure)
    
  } #end trial loop
  
  
  #generate animation of value function across trials
  vplot <- list()
  for (v in 1:nrow(value)) {
    plot(1:ntimesteps, value[v,], type="l", xlab="timestep", ylab="expected value", ylim=c(floor(min(value)), ceiling(max(value))))
    text(x=ntimesteps, y=max(value), paste0("i=", v)) #add trial number at
    text(x=us_times[v], y=(max(value) - min(value))*0.5, us_times[v])
    abline(v=c(stim_times[1,-ncol(stim_times)])) #plot cs_times using first row of stim_times (minus us time, the last col)
    vplot[[v]] <- recordPlot()   #ani.record()
  }
  
  #action function
  aplot <- list()
  for (a in 1:nrow(action)) {
    plot(1:ntimesteps, action[a,], type="l", xlab="timestep", ylab="action", ylim=c(floor(min(action)), ceiling(max(action))))
    text(x=ntimesteps, y=max(action), paste0("i=", a)) #add trial number at
    text(x=us_times[a], y=(max(action) - min(action))*0.5, us_times[a])
    abline(v=c(stim_times[1,-ncol(stim_times)])) #plot cs_times using first row of stim_times (minus us time, the last col)
    aplot[[a]] <- recordPlot()   #ani.record()
  }
  
  #action + value
  both <- list()
  par(mfrow=c(2,1))
  for (a in 1:nrow(action)) {
    plot(1:ntimesteps, action[a,], type="l", xlab="timestep", ylab="action", ylim=c(floor(min(action)), ceiling(max(action))))
    text(x=ntimesteps, y=max(action), paste0("i=", a)) #add trial number at
    text(x=us_times[a], y=(max(action) - min(action))*0.5, us_times[a])
    abline(v=c(stim_times[1,-ncol(stim_times)])) #plot cs_times using first row of stim_times (minus us time, the last col)
    
    plot(1:ntimesteps, value[a,], type="l", xlab="timestep", ylab="expected value", ylim=c(floor(min(value)), ceiling(max(value))))
    text(x=ntimesteps, y=max(value), paste0("i=", a)) #add trial number at
    text(x=us_times[a], y=(max(value) - min(value))*0.5, us_times[a])
    abline(v=c(stim_times[1,-ncol(stim_times)])) #plot cs_times using first row of stim_times (minus us time, the last col)
    
    both[[a]] <- recordPlot()   #ani.record()
  }
  par(mfrow=c(1,1))
  
  #td error function
  dplot <- list()
  for (d in 1:nrow(delta)) {
    plot(1:ntimesteps, delta[d,], type="l", xlab="timestep", ylab="delta (td error)", ylim=c(floor(min(delta)), ceiling(max(delta))))
    text(x=ntimesteps, y=max(delta), paste0("i=", d)) #add trial number at
    text(x=us_times[d], y=(max(delta) - min(delta))*0.5, us_times[d])
    abline(v=c(stim_times[1,-ncol(stim_times)])) #plot cs_times using first row of stim_times (minus us time, the last col)
    dplot[[d]] <- recordPlot()   #ani.record()
  }
  
  #uncertainty value function
  uplot <- list()
  for (u in 1:nrow(uncertainty)) {
    plot(1:ntimesteps, uncertainty[u,], type="l", xlab="timestep", ylab="uncertainty (value of sampling)", ylim=c(floor(min(uncertainty)), ceiling(max(uncertainty))))
    text(x=ntimesteps, y=max(uncertainty), paste0("i=", u)) #add trial number at
    text(x=us_times[u], y=(max(uncertainty) - min(uncertainty))*0.5, us_times[u])
    abline(v=c(stim_times[1,-ncol(stim_times)])) #plot cs_times using first row of stim_times (minus us time, the last col)
    uplot[[u]] <- recordPlot()   #ani.record()
  }
  
  tracklist$value <- value
  tracklist$delta <- delta
  tracklist$action <- action
  tracklist$vplot <- vplot
  tracklist$aplot <- aplot
  tracklist$dplot <- dplot
  tracklist$uplot <- uplot
  tracklist$avplot <- both 
  
  return(tracklist)  
}