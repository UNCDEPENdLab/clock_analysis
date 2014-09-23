getScore <- function(RT, scrfunc, noise=5) {
	
	fm = getMagFreq(RT, scrfunc)
	
	if (noise > 0) {
		a = -noise
		b = noise
		r = a + (b-a) * runif(1, 0, 1)
		#noise is an integer from -5 to 5
		r = round(r)
	} else {
		r <- 0
	}
	
	actualMag = fm["Mag"] + r
	ev = actualMag*fm["Freq"]
	actualMag = round(actualMag);
	
	rd = runif(1, 0, 1)
	
	if (fm["Freq"] > rd && RT > 0 && RT <= 4000) {
		score = unname(actualMag)
	} else {
		score = 0
	}
	
	return(score)
}



getMagFreq <- function(RT,scrfunc) { 
	k = 37;
	Shift = 700;
	rt_extended = 7000;
	DEV_factor = 10;
	DEV_factor2= 1;
	sin_factor = 0.25;
	
	if (scrfunc == "CEV") {
		Mag = (k*rt_extended)/(rt_extended-(RT+Shift))
		Freq = 1-((RT+Shift)/rt_extended)
	} else if (scrfunc=="DEV") {
		Mag = DEV_factor*log(DEV_factor2*(RT+Shift))
		CEV_x = 1-((RT+Shift)/rt_extended)
		IEV_x = CEV_x + (CEV_x*(sin_factor*sin((RT*pi)/5000)))
		Freq = (2*CEV_x)-IEV_x	
	} else if (scrfunc=="IEV") {
		CEV_x = (k*rt_extended)/(rt_extended-(RT+Shift))
		DEV_x = DEV_factor*log(DEV_factor2*(RT+Shift))
		Mag = (2*CEV_x)-(DEV_x)
		CEV_x2 = 1-((RT+Shift)/rt_extended)
		Freq = CEV_x2 + (CEV_x2*(sin_factor*sin((RT*pi)/5000)))
	} else if (scrfunc=="CEVR") {
		Mag = 1-((RT+Shift)/rt_extended)
		Mag = Mag*200
		Freq = (k*rt_extended)/(rt_extended-(RT+Shift))
		Freq = Freq/200
	} else { stop("unknown scrfunc: ", scrfunc) }
	
	return(c(Mag=Mag, Freq=Freq))
	
}

