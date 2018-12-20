#reverse engineer rewfunc
#rather than solve for rt, just do a brute force computation... :-)
gen_rt_objects <- function(maxrt=4050) { #a bit above trial length since the flip time can pass the trial end slightly
  rtvec <- 1:maxrt
  
  Shift = 700
  rt_extended = 7000
  
  CEV_frq <- sapply(rtvec, function(rt) { 1-((rt+ Shift)/rt_extended) })
  CEV_mag <- sapply(rtvec, function(rt) { rt_extended*37/(rt_extended-(rt+Shift)) })
  
  IEV_frq <- sapply(rtvec, function(rt) { CEV_frq[rt] + CEV_frq[rt]*(0.25*sin((rt*pi)/5000)) })
  IEV_mag <- sapply(rtvec, function(rt) { 2*CEV_mag[rt] -(10*log(rt+ Shift)) })
  
  CEVR_frq <- sapply(rtvec, function(rt) { (rt_extended)*37/(200*(rt_extended-(rt+Shift))) })
  CEVR_mag <- sapply(rtvec, function(rt) { 200*(1-((rt+Shift)/ rt_extended)) })
  
  DEV_frq <- sapply(rtvec, function(rt) { 2*CEV_frq[rt] -  IEV_frq[rt] })
  DEV_mag <- sapply(rtvec, function(rt) { 10*log(rt+Shift) })
  
  return(list(RT=rtvec, CEV_frq=CEV_frq, CEV_mag=CEV_mag, IEV_frq=IEV_frq, IEV_mag=IEV_mag,
              CEVR_frq=CEVR_frq, CEVR_mag=CEVR_mag, DEV_frq=DEV_frq, DEV_mag=DEV_mag))
}

lookup <- gen_rt_objects()

RewFunction_Reversed <- function(frq, cond, L) { #mag, 
  if (cond=='CEV') {
    matchv <- abs(frq - L$CEV_frq)
    #matchv <- abs(mag - CEV_mag)
  } else if (cond=="CEVR") {
    matchv <- abs(frq - L$CEVR_frq)
  } else if (cond=="DEV") {
    matchv <- abs(frq - L$DEV_frq)
  } else if (cond=="IEV") {
    matchv <- abs(frq - L$IEV_frq)
  } else {
    stop("unknown cond:", cond)
  }

  #stopifnot(sum(matchv < .0001) == 1)
  rt <- which.min(matchv)
  trustworthy <- sum(matchv < .000001) == 1 #is the match singular?

  return(data.frame(rt_backcalc=rt, backcalc_trustworthy=trustworthy))  
}
