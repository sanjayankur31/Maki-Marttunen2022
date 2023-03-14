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

from pyneuroml.analysis import generate_current_vs_frequency_curve

from figure_02_experiment import (create_model, simulate_model,
                                  delete_neuron_special_dir,
                                  get_timestamp)


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

    # 0 is blocked, 1 is unchanged
    for g in [0.0, 1.0]:
        model_file_name = create_model(
            cellname=cellname,
            celldir=celldir,
            current_nA="0.5 nA",
            g_Ih_multiplier=f"{g}",
            distance_from_soma_inc=50
        )
        simlist.append(model_file_name)
        simulate_model(model_file_name, cellname, duration_ms=10.)

        """
        # delete compiled mod files
        delete_neuron_special_dir()

        generate_current_vs_frequency_curve(
            nml2_file=model_file_name,
            cell_id=cellname,
            custom_amps_nA=list(np.arange(0, 1.6, 0.1)),
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
        """

    os.chdir("../")
    print(f"Simdir: {simdir}")
    return simlist


if __name__ == "__main__":
    simlist = run()
    print(simlist)
