TITLE hyperpolarization-activated current (H-current) 

COMMENT
Based on Williams and Stuart J. Neurophysiol 83:3177,2000
ENDCOMMENT

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX iHalmog
	USEION h READ eh WRITE ih VALENCE 1
	RANGE gbar, h_inf, tau, ih, t0,t1,off, slo, offt1, offt2, slot1, slot2, q10, temp, tadj, vmin,vmax
	
}

UNITS {
	(mA) 	= (milliamp)
	(mV) 	= (millivolt)
	(molar)	= (1/liter)
	(mM) 	= (millimolar)
	(pS) = (picosiemens)
	(um) = (micron)

}

PARAMETER {
	v		(mV)
	
	celsius		(degC)
	eh	   (mV)     
	gbar	= 0.0	(pS/um2)

	off = -91       (mV)   		: V1/2 of activation	
	slo=6		(mV)	 	: slope of activation
	
	t0 = 2542.5883549	(ms) 	: parameters for time constant of activation    
	t1 = 11.40250855	(ms)     
	offt1 = 0		(mV)
	offt2 = 0		(mV)
	slot1 = 40.1606426		(mV)     
	slot2 = 16.1290323		(mV)
			
	temp = 21	(degC)		: original temp 
	q10  = 2.3			: temperature sensitivity
	vmin = -120 (mV)
	vmax = 100 (mV)     
	
}


ASSIGNED {
	ih		(mA/cm2)
        h_inf
        tau        (ms)
	tadj
	
}

STATE { h }


INITIAL {
	rates(v)
      	h = h_inf
}

BREAKPOINT { 
	SOLVE states METHOD cnexp
      	ih = (1e-4) * gbar * h * (v-eh)
}


DERIVATIVE states  { 

	rates(v) 
	h' = (h_inf-h)/tau  
}


PROCEDURE rates( v (mV)) {

	tadj= q10^((celsius-22)/10)
	h_inf = 1/(1+exp((v-off)/slo))	
        tau = 1/(tadj*(exp(-(v-offt1)/slot1)/t0+exp((v-offt2)/slot2)/t1))
}	

