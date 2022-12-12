#!/usr/bin/env python3
"""
Simulations for figure 1 for the L5PC cell.

File: NeuroML/experiments/figure_01_runner_L5PC.py

Copyright 2022 Ankur Sinha
Author: Ankur Sinha <sanjay DOT ankur AT gmail DOT com>
"""

import os
import numpy as np

from pyneuroml.analysis import generate_current_vs_frequency_curve

from figure_01_experiment import (create_model, simulate_model,
                                  delete_neuron_special_dir, plot,
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

    # 2 is "normal"
    for g in [0.0, 2.0, 4.0]:
        model_file_name = create_model(
            cellname=cellname,
            celldir=celldir,
            current_nA="0.5 nA",
            gIh_S_per_m2=f"{g} S_per_m2",
        )
        simlist.append(model_file_name)
        simulate_model(model_file_name, cellname)

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

    os.chdir("../")


if __name__ == "__main__":
    simlist = run()
    print(simlist)
    plot(simlist)
