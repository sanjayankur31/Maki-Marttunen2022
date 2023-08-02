TITLE K-fast channel from Korngreen and Sakmann 2000
: M.Migliore June 2001

NEURON {
	SUFFIX iA
	USEION k READ ek WRITE ik
        RANGE gbar,ninf,linf,taul,taun,g
        GLOBAL tq,qq, q10,vmin, vmax,tadj,temp
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(um) = (micron)

}

PARAMETER {
	v (mV)
	celsius		(degC)
	Tscale = 10	(degC)
	gbar=0.0 (pS/um2)
        offn=-47   (mV)
        offl=-66   (mV)
        slon=29   (mV)
        slol=10   (mV)
	qq=5
	tq=-55
	ek      (mV)
	vmin = -120	(mV)
	vmax = 100	(mV)
        temp = 21       (degC)          : original temp
        q10  = 2.3

        offmt = -71 (mV)
        slomt = 59 (mV)
        taummin = 0.34 (ms)
        taumdiff = 0.92 (ms)
        offht = -73 (mV)
        sloht = 23 (mV)
        tauhmin = 8 (ms)
        tauhdiff = 49 (ms)
	
}



STATE {
	n
        l
}

ASSIGNED {
	ik (mA/cm2)
        g  (pS/um2)
        ninf
        linf      
        taul  (ms)
        taun   (ms)
	tadj
}

INITIAL {
	rates(v)
	n=ninf
	l=linf
}


BREAKPOINT {
	SOLVE states METHOD cnexp
	g = gbar*n^4*l
	ik = (1e-4)*g*(v-ek)
}


DERIVATIVE states {     : exact when v held constant; integrates over dt step
        rates(v)
        n' = (ninf - n)/taun
        l' =  (linf - l)/taul
}

PROCEDURE rates(v (mV)) { :callable from hoc

        tadj= q10^((celsius-temp)/Tscale)
        ninf = 1/(1 + exp(-(v-offn)/slon))
        linf = 1/(1 + exp((v-offl)/slol))
        taun = (taummin+taumdiff*exp(-((offmt-v)/slomt)^2))/tadj
        taul = (tauhmin+tauhdiff*exp(-((offht-v)/sloht)^2))/tadj
}














