COMMENT
Copied from IClamp, use amp and e instead of just amp
ENDCOMMENT

NEURON {
        POINT_PROCESS GClamp
        RANGE del, dur, amp, e, i
        ELECTRODE_CURRENT i
}
UNITS {
        (nA) = (nanoamp)
}

PARAMETER {
        del (ms)
        dur (ms)        <0,1e9>
        amp (uS)
	e (mV)
}
ASSIGNED {
         v (mV)
         i (nA)
}

INITIAL {
        i = 0
}

BREAKPOINT {
        at_time(del)
        at_time(del+dur)

        if (t < del + dur && t >= del) {
                i = amp*(e-v)
        }else{
                i = 0
        }
}
