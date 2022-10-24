:Comment :
:Reference : :		Kole,Hallermann,and Stuart, J. Neurosci. 2006

NEURON	{
	SUFFIX Ih
	NONSPECIFIC_CURRENT ihcn
	RANGE gIhbar, gIh, ihcn, ehcn, offma, sloma, tauma, offmb, slomb, taumb
}

UNITS	{
	(S) = (siemens)
	(mV) = (millivolt)
	(mA) = (milliamp)
}

PARAMETER	{
	gIhbar = 0.00001 (S/cm2) 
	ehcn =  -45.0 (mV)
	offma = -154.9 (mV)
	sloma = 11.9 (mV)
	tauma = 155.521 (ms)
	offmb = 0.0 (mV)
	slomb = 33.1 (mV)
	taumb = 5.18135 (ms)
}

ASSIGNED	{
	v	(mV)
	ihcn	(mA/cm2)
	gIh	(S/cm2)
	mInf
	mTau
	mAlpha
	mBeta
}

STATE	{ 
	m
}

BREAKPOINT	{
	SOLVE states METHOD cnexp
	gIh = gIhbar*m
	ihcn = gIh*(v-ehcn)
}

DERIVATIVE states	{
	rates()
	m' = (mInf-m)/mTau
}

INITIAL{
	rates()
	m = mInf
}

PROCEDURE rates(){
	UNITSOFF
        if(v == offma){
            v = v + 0.0001
        }
		mAlpha = -(offma-v)/tauma/(exp(-(offma-v)/sloma)-1)
		mBeta  = exp(-(offmb-v)/slomb)/taumb
		mInf = mAlpha/(mAlpha + mBeta)
		mTau = 1/(mAlpha + mBeta)
	UNITSON
}
