: this model is built-in to neuron with suffix syn_g

COMMENT
synaptic current with exponential rise and decay conductance defined by
        i = g * (v - e)      i(nanoamps), g(micromhos);
        where
         g = 0 for t < onset and
         g=amp*((1-exp(-(t-onset)/tau0))-(1-exp(-(t-onset)/tau1)))
          for t > onset
ENDCOMMENT
					       
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	POINT_PROCESS syn_g
	RANGE onset, tau0, tau1, gmax, e, iampa, gampa
	NONSPECIFIC_CURRENT iampa
}
UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(umho) = (micromho)
}

PARAMETER {
	onset=0  (ms)
	tau0=0.3 (ms)
	tau1=3.0 (ms)
	gmax=0	 (umho)
	e=0	 (mV)
	v	 (mV)
}

ASSIGNED { iampa (nA)  gampa (umho) }

LOCAL   a[2]
LOCAL   tpeak
LOCAL   adjust
LOCAL   amp

BREAKPOINT {
        gampa = cond(t)
	iampa = gampa*(v - e)
}

FUNCTION cond(x) {
	tpeak=tau0*tau1*log(tau0/tau1)/(tau0-tau1)
	adjust=1/((1-exp(-tpeak/tau0))-(1-exp(-tpeak/tau1)))
	amp=adjust*gmax
	if (x < onset) {
		cond = 0
	}else{
		a[0]=1-exp(-(x-onset)/tau0)
		a[1]=1-exp(-(x-onset)/tau1)
		cond = amp*(a[0]-a[1])
	}
}
