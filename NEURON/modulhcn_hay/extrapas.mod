:Comment :

NEURON	{
	SUFFIX extrapas
	NONSPECIFIC_CURRENT i
	RANGE g, i, e
}

UNITS	{
	(S) = (siemens)
	(mV) = (millivolt)
	(mA) = (milliamp)
}

PARAMETER	{
	g = 0.000000 (S/cm2) 
	e =  -45.0 (mV)
}

ASSIGNED	{
	v	(mV)
	i	(mA/cm2)
}

STATE	{ 
}

BREAKPOINT	{
	i = g*(v-e)
}

INITIAL{
}

