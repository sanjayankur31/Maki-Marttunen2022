#!/usr/bin/env python3
"""
Simulations for figure 1 for the HL5PC cell with ScZ variants.

File: NeuroML/experiments/figure_01_runner_HL5PC_ScZ.py

Copyright 2022 Ankur Sinha
Author: Ankur Sinha <sanjay DOT ankur AT gmail DOT com>
"""

import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(".")))
sys.path.append(os.path.dirname(os.path.abspath("..")))

import numpy as np
from figure_01_experiment import runner
from multiprocessing import set_start_method


def run():
    """Run simulations"""
    return runner(
        cellname="HL5PC", celldir="SkinnerLabHL5PC", num_data_points=10,
        step_sim=True,
        if_curve=True,
        sim_current_na="0.5nA",
        ifcurve_custom_amps=list(np.arange(0, 150E-3, 10E-3)),
        scz=True)


if __name__ == "__main__":
    set_start_method("spawn")
    simlist = run()
    print(simlist)
