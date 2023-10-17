#!/usr/bin/env python3
"""
Simulations for figure 1 for the HL5PC cell.

File: NeuroML/experiments/figure_01_runner_HL5PC.py

Copyright 2022 Ankur Sinha
Author: Ankur Sinha <sanjay DOT ankur AT gmail DOT com>
"""


import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(".")))
sys.path.append(os.path.dirname(os.path.abspath("..")))

import numpy as np
import matplotlib

from pyneuroml.analysis import generate_current_vs_frequency_curve

from figure_01_experiment import (create_model, simulate_model)
from common import (delete_neuron_special_dir, get_run_dir,
                    get_abs_celldir)


# increase plot size
matplotlib.rcParams['figure.figsize'] = [19.2, 10.8]


def run():
    """Run simulations"""
    simlist = []
    cellname = "HL5PC"
    celldir = get_abs_celldir("SkinnerLabHL5PC")
    simdir = get_run_dir(cellname, "figure_01")

    os.mkdir(simdir)
    os.chdir(simdir)

    for g in [0., 1.0, 2.0]:
        model_file_name = create_model(
            cellname=cellname,
            celldir=celldir,
            current_nA="0.050 nA",
            g_Ih_multiplier=f"{g}",
        )
        simlist.append(model_file_name)
        simulate_model(model_file_name, cellname)

        # delete compiled mod files
        delete_neuron_special_dir()

        generate_current_vs_frequency_curve(
            nml2_file=model_file_name,
            cell_id="HL5PC",
            custom_amps_nA=list(np.arange(0E-3, 150E-3, 5E-4)),
            analysis_duration=2000,
            analysis_delay=200,
            temperature="34 degC",
            simulator="jNeuroML_NEURON",
            plot_if=True,
            plot_iv=True,
            pre_zero_pulse=300,
            post_zero_pulse=300,
            plot_voltage_traces=True,
            save_iv_figure_to=f"{model_file_name}_iv.png",
            save_iv_data_to=f"{model_file_name}_iv.dat",
            save_if_figure_to=f"{model_file_name}_if.png",
            save_if_data_to=f"{model_file_name}_if.dat",
            save_voltage_traces_to=f"{model_file_name}_v.png",
            show_plot_already=False,
            num_processors=8,
        )

    os.chdir("../")


if __name__ == "__main__":
    simlist = run()
    print(simlist)