COMMENT
Copied from IClamp, use amp and e instead of just amp
ENDCOMMENT

NEURON {
        POINT_PROCESS GClampCa
        USEION ca READ eca WRITE ica
        RANGE del, dur, amp, ica
}
UNITS {
        (nA) = (nanoamp)
}

PARAMETER {
        del (ms)
        dur (ms)        <0,1e9>
        amp (uS)
}
ASSIGNED {
         v (mV)
         eca     (mV)
         ica (nA)
}

INITIAL {
        ica = 0
}

BREAKPOINT {
        at_time(del)
        at_time(del+dur)

        if (t < del + dur && t >= del) {
                ica = amp*(eca-v)
        }else{
                ica = 0
        }
}
