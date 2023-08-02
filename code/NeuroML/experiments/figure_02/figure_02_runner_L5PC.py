#!/usr/bin/env python3
"""
Simulations for figure 2 for the L5PC cell.

File: NeuroML/experiments/figure_02/figure_02_runner_L5PC.py

Copyright 2023 Ankur Sinha
Author: Ankur Sinha <sanjay DOT ankur AT gmail DOT com>
"""

import os
import sys
# add experiments folder to path to ensure things can be imported
# needs to be done in all runner scripts
sys.path.append(os.path.dirname(os.path.abspath(".")))

import numpy as np
import matplotlib

from figure_02_experiment import (create_model, simulate_model,get_segments_at_distances)
from common import (delete_neuron_special_dir, get_timestamp)


# increase plot size
matplotlib.rcParams['figure.figsize'] = [19.2, 10.8]

def run():
    """Run simulations """
    simlist = []
    cellname = "L5PC"
    celldir = "../../../cells/HayEtAlL5PC/"

    # run simulations in new directory so that they don't conflict with
    # existing files.
    currenttime = get_timestamp()
    simdir = f"{currenttime}_{cellname}"
    os.mkdir(simdir)
    os.chdir(simdir)

    dist_vs_seg_dict = get_segments_at_distances(celldir, cellname,
                                                 "apical_dendrites_group", 50., [500])

    # g: 0 is blocked, 1 is unchanged
    # log scale starting at 10, but include 30, 100 for figure 2a
    current_range = list(np.logspace(start=1, stop=2.5, num=20))
    current_range.extend([30, 100])
    current_range = sorted(current_range)
    lems_file = ""
    for g in [0.001, 1.0]:
        for current in current_range:
            model_file_name = create_model(
                cellname=cellname,
                celldir=celldir,
                dist_vs_seg_dict=dist_vs_seg_dict,
                current_nA=f"{current} nA",
                g_Ih_multiplier=f"{g}",
            )
            lems_file = simulate_model(model_file_name, cellname, duration_ms=50.,
                                       skip_run=False, plot=True)
            simlist.append(lems_file)

    # delete compiled mod files
    delete_neuron_special_dir()

    os.chdir("../")
    print(f"Simdir: {simdir}")
    return simlist


if __name__ == "__main__":
    simlist = run()
    print(simlist)
