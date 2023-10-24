#!/usr/bin/env python3
"""
Simulations for figure 2 for the HL5PC cell.

File: NeuroML/experiments/figure_02/figure_02_runner_HL5PC.py

Copyright 2023 Ankur Sinha
Author: Ankur Sinha <sanjay DOT ankur AT gmail DOT com>
"""

import os
import sys
# add experiments folder to path to ensure things can be imported
# needs to be done in all runner scripts
sys.path.append(os.path.dirname(os.path.abspath(".")))
sys.path.append(os.path.dirname(os.path.abspath("..")))


from figure_02_experiment import runner
from multiprocessing import set_start_method


def run():
    """Run simulations """
    return runner(cellname="HL5PC", celldir="SkinnerLabHL5PC", num_processes=8)


if __name__ == "__main__":
    set_start_method("spawn")
    simlist = run()
    print(simlist)
