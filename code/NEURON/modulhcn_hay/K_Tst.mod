:Comment : The transient component of the K current
:Reference : :		Voltage-gated K+ channels in layer 5 neocortical pyramidal neurones from young rats:subtypes and gradients,Korngreen and Sakmann, J. Physiology, 2000
:Comment : shifted -10 mv to correct for junction potential
:Comment: corrected rates using q10 = 2.3, target temperature 34, orginal 21

NEURON	{
	SUFFIX K_Tst
	USEION k READ ek WRITE ik
	RANGE gK_Tstbar, gK_Tst, ik, offm, slom, offh, sloh, offmt, slomt, taummin, taumdiff,offht, sloht, tauhmin, tauhdiff
}

UNITS	{
	(S) = (siemens)
	(mV) = (millivolt)
	(mA) = (milliamp)
}

PARAMETER	{
	gK_Tstbar = 0.00001 (S/cm2)
        offm = -10 (mV)
        slom = 19 (mV)
        offh = -76 (mV)
        sloh = 10 (mV)
        offmt = -81 (mV)
        slomt = 59 (mV)
	taummin = 0.34 (ms)
	taumdiff = 0.92 (ms)
        offht = -83 (mV)
        sloht = 23 (mV)
	tauhmin = 8 (ms)
	tauhdiff = 49 (ms)
}

ASSIGNED	{
	v	(mV)
	ek	(mV)
	ik	(mA/cm2)
	gK_Tst	(S/cm2)
	mInf
	mTau
	hInf
	hTau
}

STATE	{
	m
	h
}

BREAKPOINT	{
	SOLVE states METHOD cnexp
	gK_Tst = gK_Tstbar*(m^4)*h
	ik = gK_Tst*(v-ek)
}

DERIVATIVE states	{
	rates()
	m' = (mInf-m)/mTau
	h' = (hInf-h)/hTau
}

INITIAL{
	rates()
	m = mInf
	h = hInf
}

PROCEDURE rates(){
  LOCAL qt
  qt = 2.3^((34-21)/10)

	UNITSOFF
		mInf =  1/(1 + exp((offm-v)/slom))
		hInf =  1/(1 + exp(-(offh-v)/sloh))
		mTau =  (taummin+taumdiff*exp(-((offmt-v)/slomt)^2))/qt
		hTau =  (tauhmin+tauhdiff*exp(-((offht-v)/sloht)^2))/qt
	UNITSON
}
