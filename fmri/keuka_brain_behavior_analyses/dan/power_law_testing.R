#power law goofing
ntimesteps=500

ev = 10*sin(2*pi*(1:ntimesteps)*1/ntimesteps) + 2.5*sin(2*pi*(1:ntimesteps)*2/ntimesteps) + 2.0*cos(2*pi*(1:ntimesteps)*4/ntimesteps)
ev = ev + abs(min(ev)) + 10;
prb = 25*cos(2*pi*(1:ntimesteps)*1/ntimesteps) + 10*cos(2*pi*(1:ntimesteps)*3/ntimesteps) + 6*sin(2*pi*(1:ntimesteps)*5/ntimesteps)
prb_max=0.7
prb_min=0.3
prb = (prb - min(prb))*(prb_max-prb_min)/(max(prb)-min(prb)) + prb_min

allshift = array(NA_real_, dim=c(ntimesteps, ntimesteps, 3))

for (i in 1:ntimesteps) {
  if (i > 1) {
    shift = c(i:ntimesteps, 1:(i-1))
  } else { shift <- 1:ntimesteps }
  evi = ev[shift]
  prbi = prb[shift]
  
  allshift[i,,1] = evi
  allshift[i,,2] = prbi
  allshift[i,,3] = evi/prbi
  
}

EV=allshift[,1,1]

sigmoid <- function(x, beta=1) {
  y <- 1/(1+exp(-beta*x))
  return(y)
}


softmax <- function(m, beta=1) {
  
  m <- (m - max(m))/beta #avoid floating point overflow
  
  #p_choice = (exp((v_func-max(v_func))/beta) / (sum(exp((v_func-max(v_func))/beta))); %Divide by temperature
  
  p <- exp(m)/sum(exp(m))
  return(p)
}


#EV1 <- sigmoid(EV, beta=0.5)
#plot(EV^-1, type="l")

plot(softmax(EV, beta=.001), type="l")



plot(EV^1, type="l")

plot(exp(EV/2), type="l")
plot(exp(EV/5), type="l")
plot(exp(EV/10), type="l")
plot(exp(EV/200), type="l")

exp_renorm <- function(x, beta=1) {
  trans <- exp((x-max(x))/beta) #avoid floating point overflow

  #https://www.ibm.com/support/pages/transforming-different-likert-scales-common-scale
  trans <- (max(x) - min(x))*(trans - min(trans))/(max(trans) - min(trans)) + min(x)
  
  # minshift <- min(trans) - min(x)
  # trans <- trans - minshift
  # rpeak <- (max(trans) - min(trans))/(max(x) - min(x))
  return(trans)
}

plot(exp_renorm(EV), type="l")
plot(EV, type="l")
exp_renorm(c(rnorm(100)/1e6, 0), beta=1)

exp_renorm(rep(10, 100), beta=1)

#plot(exp_renorm(EV, beta=.1), type="l", main="SCEPTIC value compression", xlab = "Time", ylab = "Value")
plot(EV, type="l", main="SCEPTIC value compression", xlab = "Time", ylab = "Value")
lines(exp_renorm(EV, beta=1), type="l", col="red")
lines(exp_renorm(EV, beta=3), type="l", col="green")
lines(exp_renorm(EV, beta=10), type="l", col="gray")
lines(exp_renorm(EV, beta=20), type="l", col="blue")
lines(exp_renorm(EV, beta=30), type="l", col="orange")
lines(exp_renorm(EV, beta=1000), type="l", col="orange")


range(EV)
range(exp_renorm(EV, beta=25))

plot(EV1, type="l")

vartrans <- function(x, target=0.5) {
  sdx <- sd(x, na.rm=T)
  sdt <- sdx * target
}



#u model
udf <- read.csv("/Users/hallquist/Data_Analysis/clock_analysis/fmri/keuka_brain_behavior_analyses/dan/vba_originals/mmclock_fmri_fixed_uv_ureset_mfx_trial_statistics.csv.gz")
