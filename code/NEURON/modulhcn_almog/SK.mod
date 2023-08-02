:modified 1/7/2007 by Chris Deister for the GP neuron (to remove some of the background current that existed in Mercer 2007)

NEURON {
	SUFFIX sk
	USEION k READ ek WRITE ik
        USEION ca READ cai
        RANGE  gbar,gkahp,ik, inf,tau,g
        GLOBAL Cq10
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(molar) = (1/liter)
	(mM) = (millimolar)
	(pS) = (picosiemens)
	(um) = (micron)
}

PARAMETER {
	gbar = 0	(pS/um2)
        n = 4
        cai = 50.e-6	(mM)
        b0inv = 16.6666667	(ms)			:1/b0
	celsius = 37	(degC)
	offc = 0.04635023	(mM)                    :(b0/a0)^4
	sloc = 4.0
	Cq10 = 3
}

STATE {	w }

ASSIGNED {
	ik	(mA/cm2)
        g	(pS/um2)
        inf
        tau	(ms)
	a	(1/ms)
        v	(mV)
        ek	(mV)
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	g = gbar*w
	ik = (1e-4)* g*(v-ek)
}

INITIAL {
	rate(cai)
	w=inf
}

DERIVATIVE state {
	rate(cai)
	w' = (inf - w)/tau
}

PROCEDURE rate(cai (mM)) {
	LOCAL q10
	q10 = Cq10^((celsius - 22 (degC))/10 (degC) )
	tau = q10*b0inv/(1+(cai/offc)^sloc)
	inf = 1/(1+(offc/cai)^sloc)
}


