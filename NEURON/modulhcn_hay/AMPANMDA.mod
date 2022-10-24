TITLE AMPA and NMDA receptor with short-term plasticity 


COMMENT
Conductance-based AMPA-NMDA synaptic current with Wang-type of short-term synaptic depression
Implementation by Tuomo Maki-Marttunen, 2016
Tuomo 2021: removed depression
ENDCOMMENT


NEURON {

        POINT_PROCESS AMPANMDA  
        RANGE gAMPAmax, gNMDAmax, MgCon
	RANGE E_Glu, tau_sAMPA, tau_sNMDA, tau_xNMDA, alphas
        RANGE i, i_AMPA, i_NMDA, g_AMPA, g_NMDA, sAMPA, sNMDA, xNMDA
        NONSPECIFIC_CURRENT i, i_AMPA,i_NMDA

}

PARAMETER {

	gAMPAmax = 0.01  (uS)
	gNMDAmax = 0.007 (uS)
	MgCon = 0.69     
	mggate

	E_Glu = 0       (mV)
	tau_sAMPA = 2   (ms)
	tau_sNMDA = 100 (ms)
	tau_xNMDA = 2   (ms)
	alphas = 0.5    (kHz)

}

ASSIGNED {

        v (mV)
        i (nA)
	i_AMPA (nA)
	i_NMDA (nA)
        g_AMPA (uS)
	g_NMDA (uS)
}

STATE {

        sAMPA       : AMPA state variable to construct the single-exponential profile - decays with conductance tau_sAMPA
	sNMDA       : NMDA state variable to construct the dual-exponential profile - decays with conductance tau_sAMPA
        xNMDA       : NMDA state variable to construct the dual-exponential profile - decays with conductance tau_xAMPA
}

INITIAL{

	sAMPA = 0
	sNMDA = 0
	xNMDA = 0
        
}

BREAKPOINT {

        SOLVE state METHOD cnexp
        if (sNMDA > 1) { :Do not allow larger values than 1
          sNMDA = 1
        }
	mggate = 1 / (1 + exp(0.062 (/mV) * -(v)) * (MgCon / 3.57 (mM))) :mggate kinetics - Jahr & Stevens 1990
        g_AMPA = gAMPAmax*sAMPA          :compute time varying conductance
	g_NMDA = gNMDAmax*sNMDA * mggate :compute time varying conductance using mggate kinetics
        i_AMPA = g_AMPA*(v-E_Glu) :compute the AMPA driving force based on the time varying conductance, membrane potential, and AMPA reversal
	i_NMDA = g_NMDA*(v-E_Glu) :compute the NMDA driving force based on the time varying conductance, membrane potential, and NMDA reversal
	i = i_AMPA + i_NMDA
}

DERIVATIVE state{

        sAMPA' = -sAMPA/tau_sAMPA
	sNMDA' = -sNMDA/tau_sNMDA + alphas*xNMDA*(1-sNMDA)
        xNMDA' = -xNMDA/tau_xNMDA
}


:NET_RECEIVE (weight, Pv, Pr, u, tsyn (ms)){
NET_RECEIVE (weight){
	
        sAMPA = sAMPA + 1
	xNMDA = xNMDA + 1
        if (sAMPA > 1) { :Do not allow larger values than 1
          sAMPA = 1
        }
}

