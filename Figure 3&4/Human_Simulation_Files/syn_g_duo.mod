: this model is built-in to neuron with suffix syn_g_duo
: implemented by Suraj Honnuraiah and Greg Stuart
: contact hs@ini.phys.ethz.ch

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
	POINT_PROCESS syn_g_duo
	RANGE onset, tau0, tau1, tau2, tau3, gmax, ntar, e, i, gampa, gnmda, iampa, inmda
	NONSPECIFIC_CURRENT i, iampa, inmda
}
UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(umho) = (micromho)
}

PARAMETER {
	onset=0  (ms)
	: ampa
	tau0=0.3 (ms)
	tau1=3.0 (ms)	:2 ms 80% dexp from Spruston et al 1995;  (1.7ms @ 35)
	: nmda from Larkum 2009
	tau2=5.0 (ms) 	tau3=200 (ms) 	:200ms 80% dexp Spruston et al 1995; Hausser et al 1997 JNS
	gmax=0	 (umho)
	ntar = 2
	: nmda to amp ratio from Larkum 2009
	e=0	 (mV)
	v	 (mV)
}

ASSIGNED { i (nA) iampa (nA)  inmda (nA) g (umho) gampa (umho) gnmda (umho) }

LOCAL   a[4]
LOCAL   tpeak
LOCAL   adjust
LOCAL   amp

BREAKPOINT {
        gampa = conda(t)
        gnmda = condn(t)
	g = gampa+gnmda
	iampa= (gampa*(v - e))
	inmda= (gnmda*(v - e) *mgblock(v))
	i = iampa + inmda
}

FUNCTION conda(x) {
	tpeak=tau0*tau1*log(tau0/tau1)/(tau0-tau1)
	adjust=1/((1-exp(-tpeak/tau0))-(1-exp(-tpeak/tau1)))
	amp=adjust*gmax
	if (x < onset) {
		conda = 0
	}else{
		a[0]=1-exp(-(x-onset)/tau0)
		a[1]=1-exp(-(x-onset)/tau1)
		conda = amp*(a[0]-a[1])
	}
}

FUNCTION condn(x) {
	tpeak=tau2*tau3*log(tau2/tau3)/(tau2-tau3)
	adjust=1/((1-exp(-tpeak/tau2))-(1-exp(-tpeak/tau3)))
	amp=adjust*gmax*ntar
	if (x < onset) {
		condn = 0
	}else{
		a[2]=1-exp(-(x-onset)/tau2)
		a[3]=1-exp(-(x-onset)/tau3)
		condn = amp*(a[2]-a[3])
	}
}

FUNCTION mgblock(v) {	: from Jahr & Stevens 1990	: Slope of the voltage dependence is changed from 0.062 to 0.08 from Rhodes 2006
	: Extracellular Magnesium concentration changed from 2 to 0.7 mM
	mgblock = 1 / (1 + exp(0.08 * -v) * (0.7 / 3.57))}
