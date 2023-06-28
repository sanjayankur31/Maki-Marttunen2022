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
import json

from pyneuroml.analysis import generate_current_vs_frequency_curve
from pyneuroml.plot import generate_plot
from pyneuroml.pynml import reload_saved_data

from figure_02_experiment import (create_model, simulate_model,
                                  delete_neuron_special_dir,
                                  get_timestamp)


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

    # { g: { i: {data} } }
    simdata = {}
    # 0 is blocked, 1 is unchanged
    lems_file = ""
    for g in [0.0, 1.0]:
        simdata[g] = {}
        for current in [10, 30, 100]:
            model_file_name = create_model(
                cellname=cellname,
                celldir=celldir,
                current_nA=f"{current} nA",
                g_Ih_multiplier=f"{g}",
                distance_from_soma_inc=50
            )
            lems_file = simulate_model(model_file_name, cellname, duration_ms=10.,
                                       skip_run=False, plot=True)
            simlist.append(lems_file)
            simdata[g][current] = reload_saved_data(lems_file,
                                                    show_plot_already=False)

    # TODO: move to different plotter file so that you don't have to repeat
    # simulations
    # generate figure 2A
    # get segment indices: same for all sims, so just use the last one
    with open(lems_file.replace("LEMS_", "").replace(".xml", ".segs.json"), 'r') as f:
        segs_data = json.load(f)

    i = 0
    for dist, seg in segs_data:
        for current in [10, 30, 100]:
            generate_plot(
                xvalues=[simdata[0.0][current]['t']] * 2,
                yvalues=[simdata[0.0][current][f'L5PC_pop[{i}]/v'], simdata[1.0][30][f'L5PC_pop[{i}]/v']],
                labels=["g*0", "g*1"],
                title=f"Membrane potential ({current} nA at {dist})",
                show_plot_already=False,
                save_figure_to=f"2a_{current}_nA_{dist}-v.png",
                xaxis="time (s)",
                yaxis="membrane potential (V)",
                ylim=[-0.085, 0.065]
            )
        i = i + 1

    # delete compiled mod files
    delete_neuron_special_dir()

    os.chdir("../")
    print(f"Simdir: {simdir}")
    return simlist


if __name__ == "__main__":
    simlist = run()
    print(simlist)
