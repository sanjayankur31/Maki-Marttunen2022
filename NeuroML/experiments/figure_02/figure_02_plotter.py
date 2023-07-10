#!/usr/bin/env python3
"""
Plot graphs from simulation data

File: figure_02_plotter.py

Copyright 2023 Ankur Sinha
Author: Ankur Sinha <sanjay DOT ankur AT gmail DOT com>
"""


import os
import sys
import textwrap
from pathlib import Path
import re
import logging
import numpy
from pyneuroml.plot.Plot import generate_plot
from pyneuroml.pynml import reload_saved_data
import matplotlib
from matplotlib import pyplot as plt
import json
import numpy as np

logger = logging.Logger(__name__)
logger.setLevel(logging.INFO)


# increase plot size
matplotlib.rcParams['figure.figsize'] = [19.2, 10.8]
matplotlib.rcParams['font.weight'] = "bold"
matplotlib.rcParams['axes.labelweight'] = "bold"
matplotlib.rcParams['lines.linewidth'] = 3

def plot(simfolder: str) -> None:
    """Plot various figures


    :param simfolder: name of folder containing simulation and generated data
    :type simfolder: str
    :raises ValueError: if a dir named `simfolder` does not exist (or cannot be accessed for some reason)
    :returns: None

    """
    simdir = Path(simfolder)
    if not simdir.exists() or not simdir.is_dir():
        ValueError(f"A directory named {simfolder} does not exist or cannot be accessed")

    # generate figure 2A
    # get segment indices:
    cellname = simfolder.split("_")[1].replace("/", "")
    print(cellname)
    with open(f"{simfolder}/{cellname}.segs.json", 'r') as f:
        segs_data = json.load(f)

    lems_files = list(simdir.glob("LEMS_*xml"))

    # { g: { i: {data} } }
    simdata = {}
    current_range = list(np.logspace(start=1, stop=2.5, num=20))
    current_range.extend([30, 100])
    current_range = sorted(current_range)
    for g in [0.0, 1.0]:
        simdata[g] = {}
        for current in current_range:
            suffix = f"{g}_{current} nA".replace(".", "_").replace(" ", "_")
            for lems_file in lems_files:
                if suffix in lems_file.__str__():
                    simdata[g][current] = reload_saved_data(lems_file.__str__(),
                                                            show_plot_already=False)

    os.chdir(simdir)
    i = 0
    for dist, seg in segs_data.items():
        for current in current_range:
            generate_plot(
                xvalues=[simdata[0.0][current]['t']] * 2,
                yvalues=[simdata[0.0][current][f'{cellname}_pop[{i}]/v'], simdata[1.0][current][f'{cellname}_pop[{i}]/v']],
                linewidths=[5, 5],
                labels=["g*0", "g*1"],
                title=f"Membrane potential ({current} nA at {dist})",
                show_plot_already=False,
                save_figure_to=f"2a_{current}_nA_{dist}-v.png",
                xaxis="time (s)",
                yaxis="membrane potential (V)",
                ylim=[-0.085, 0.085],
                xlim=[0, 50e-3]
            )
        i = i + 1

    # plt.show()
    os.chdir("../")


def usage():
    """Print usage

    :returns: None
    """
    print(textwrap.dedent(
        f"""
        Usage: {__file__.split("/")[-1]} <simulation folder>
        """.strip()
    ))


if __name__ == "__main__":
    if len(sys.argv) != 2:
        logging.error("Wrong arguments received")
        usage()
        sys.exit(-1)

    plot(sys.argv[1])
