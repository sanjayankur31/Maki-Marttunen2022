TITLE GABA receptor


COMMENT
Conductance-based GABA synaptic current. Copied from AMPANMDA.mod, retained only the AMPA-form current
ENDCOMMENT


NEURON {

        POINT_PROCESS GABA
        RANGE gGABAmax
	RANGE E_Cl, tau_sGABA
        RANGE i, i_GABA, g_GABA, sGABA
        NONSPECIFIC_CURRENT i, i_GABA

}

PARAMETER {

	gGABAmax = 0.01  (uS)

	E_Cl = -80       (mV)
	tau_sGABA = 2   (ms)

}

ASSIGNED {

        v (mV)
        i (nA)
	i_GABA (nA)
        g_GABA (uS)
}

STATE {

        sGABA       : GABA state variable to construct the single-exponential profile - decays with conductance tau_sGABA
}

INITIAL{

	sGABA = 0
        
}

BREAKPOINT {

        SOLVE state METHOD cnexp
        g_GABA = gGABAmax*sGABA          :compute time varying conductance
        i_GABA = g_GABA*(v-E_Cl) :compute the GABA driving force based on the time varying conductance, membrane potential, and GABA reversal
	i = i_GABA
}

DERIVATIVE state{

        sGABA' = -sGABA/tau_sGABA
}


:NET_RECEIVE (weight, Pv, Pr, u, tsyn (ms)){
NET_RECEIVE (weight){
	
        sGABA = sGABA + 1
        if (sGABA > 1) { :Do not allow larger values than 1
          sGABA = 1
        }
}

