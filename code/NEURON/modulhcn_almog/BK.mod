TITLE BK-type Purkinje calcium-activated potassium current

COMMENT

NEURON implementation of a BK-channel in Purkinje cells
Kinetical Scheme: Hodgkin-Huxley (m^3*z^2*h)

Modified from Khaliq et al., J.Neurosci. 23(2003)4899
 
Laboratory for Neuronal Circuit Dynamics
RIKEN Brain Science Institute, Wako City, Japan
http://www.neurodynamics.brain.riken.jp

Reference: Akemann and Knoepfel, J.Neurosci. 26 (2006) 4602
Date of Implementation: May 2005
Contact: akemann@brain.riken.jp

Modified by Tuomo Maki-Marttunen: Moved CONSTANT block contents
to PARAMETER block to allow mutation-specific changes to ion
channels

ENDCOMMENT

NEURON {
       SUFFIX bk
       USEION k READ ek WRITE ik
       USEION ca READ cai
       RANGE gbar, gk,  ik, minf, taum, hinf, tauh, zinf, tauz
       GLOBAL zhalf
}

UNITS {
	(mV) = (millivolt)
	(mA) = (milliamp)
	(nA) = (nanoamp)
	(pA) = (picoamp)
	(S)  = (siemens)
	(nS) = (nanosiemens)
	(pS) = (picosiemens)
	(um) = (micron)
	(molar) = (1/liter)
	(mM) = (millimolar)		
}

PARAMETER {
	v (mV)
	celsius (degC)

	gbar = 40 (pS/um2)

	ek (mV)
	cai (mM)

	zhalf = 0.01 (mM)

	q10 = 3
	
	offm = -28.9 (mV)
	slom = 6.2 (mV)

	ctm = 0.000505 (s)
	ctmmax = 1.0 (s)
	offmt1 = -86.4 (mV)
	slomt1 = 10.1 (mV)
	offmt2 = 33.3 (mV)
	slomt2 = 10 (mV)

	ctauz = 1 (ms)

	ch = 0.085
	offh = -32 (mV)
	sloh = 5.8 (mV)
	cth = 0.0019 (s)
	cthmax = 1.0 (s)
	offht1 = -48.5 (mV)
	sloht1 = 5.2 (mV)
	offht2 = 54.2 (mV)
	sloht2 = 12.9 (mV)
}

ASSIGNED {
	ik (mA/cm2)
	qt
	gk (pS/um2)   
	minf
	taum (ms)
	hinf
	tauh (ms)
	zinf
	tauz (ms)
}

STATE {
	m   FROM 0 TO 1
	z   FROM 0 TO 1
	h   FROM 0 TO 1
}

INITIAL {
	qt = q10^((celsius-22 (degC))/10 (degC))
	rates(v)
	m = minf
	z = zinf
	h = hinf
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	gk = gbar * m^3 * z^2 * h      
	ik = (1e-4)* gk * (v - ek)
}

DERIVATIVE states {
	rates(v)
	m' = (minf-m)/taum
	z' = (zinf-z)/tauz
	h' = (hinf-h)/tauh
}

PROCEDURE rates( v (mV) ) {
	v = v + 5 (mV)
	minf = 1 / ( 1+exp((offm-v)/slom) )
	taum = (1e3) * ( ctm + ctmmax / ( exp(-(offmt1-v)/slomt1) + exp((offmt2-v)/slomt2) ) ) / qt
	
	zinf = 1 /(1 + zhalf/cai)
	tauz = ctauz/qt

	hinf = ch + (1-ch) / ( 1+exp(-(offh-v)/sloh) )
	tauh = (1e3) * ( cth + cthmax / ( exp(-(offht1-v)/sloht1) + exp((offht2-v)/sloht2) ) ) / qt
}

