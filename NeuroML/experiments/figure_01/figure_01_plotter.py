#!/usr/bin/env python3
"""
Plot graphs from simulation data

File: figure_01_plotter.py

Copyright 2023 Ankur Sinha
Author: Ankur Sinha <sanjay DOT ankur AT gmail DOT com>
"""


import sys
import textwrap
from pathlib import Path
import re
import logging
import numpy
from pyneuroml.plot.Plot import generate_plot
import matplotlib
from matplotlib import pyplot as plt

logger = logging.Logger(__name__)
logger.setLevel(logging.INFO)

# increase plot size
matplotlib.rcParams['figure.figsize'] = [19.2, 10.8]
matplotlib.rcParams['font.weight'] = "bold"
matplotlib.rcParams['axes.labelweight'] = "bold"
matplotlib.rcParams['lines.linewidth'] = 3


def plot_if(simfolder: str) -> None:
    """Plot I-F curves from data in simfolder.

    The datafiles must end in `_if.dat`. This will read all such files and
    generate a line for each.

    :param simfolder: name of folder containing simulation and generated data
    :type simfolder: str
    :raises ValueError: if a dir named `simfolder` does not exist (or cannot be accessed for some reason)
    :returns: None

    """
    simdir = Path(simfolder)
    if not simdir.exists() or not simdir.is_dir():
        ValueError(f"A directory named {simfolder} does not exist or cannot be accessed")

    datafiles = list(simdir.glob("net*_if.dat"))

    xvalues = []
    yvalues = []
    labels = []
    thresholds = []
    gmuls = []
    control_val = 0.0
    camuls = []
    for afile in datafiles:
        print(afile)
        threshold_file = afile.parent.__str__() + f"/threshold_i_{afile.name}"
        logger.debug(f"Processing {afile}")
        # can be improved: should ideally store tags in model and get values
        # from there instead of parsing filenames
        if "ScZ" not in afile.__str__():
            label = re.sub(r"_(\d+)_(\d+)_(\d+)_(\d+)_", r" \1.\2 \3.\4 ", afile.name.split(".")[0]).split("_")[2:]
            cellname = label[0].split(" ")[0]
            gmul = float(label[0].split(" ")[1])  # type: float

            flabel = f"{cellname}, gmul = {gmul}"
        else:
            label = re.sub(r"_(\d+)_(\d+)_(\d+)_(\d+)_(\d+)_(\d+)_", r" \1.\2 \3.\4 \5.\6 ", afile.name.split(".")[0]).split("_")[2:]
            cellname = label[0].split(" ")[0]
            gmul = float(label[0].split(" ")[1])  # type: float
            camul = float(label[0].split(" ")[2])  # type: float
            camuls.append(camul)

            flabel = f"{cellname}, g_Ih * {gmul:.3f}, g_CaLVast * {camul:.3f}"

        logger.debug(f"Label: {flabel}")
        labels.append(flabel)
        data = (numpy.loadtxt(afile))
        threshold_data = float(numpy.loadtxt(threshold_file))
        thresholds.append(threshold_data)
        gmuls.append(gmul)
        if gmul == 1:
            control_val = threshold_data
        # convert nA to pA
        xvalues.append(data[:, 0] / 1000)
        yvalues.append(data[:, 1])

    logger.debug(thresholds)
    logger.debug(gmuls)

    generate_plot(xvalues,
                  yvalues,
                  linewidths=[5] * len(xvalues),
                  title="F-I curve for different Ih/CaLVAst conductances",
                  xaxis="I(nA)", yaxis="f(spikes/s)",
                  show_plot_already=False, labels=labels,
                  bottom_left_spines_only=True, close_plot=False)

    # add inset with threshold values for non ScZ sims
    # https://matplotlib.org/stable/gallery/subplots_axes_and_figures/axes_demo.html#sphx-glr-gallery-subplots-axes-and-figures-axes-demo-py
    # https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.gcf.html
    fig = plt.gcf()
    if len(camuls) == 0:
        inset = fig.add_axes([0.2, 0.4, 0.1, 0.4])
        # to get the same random colours that matplotlib uses
        for i in range(0, len(gmuls)):
            barlabel = ((thresholds[i] - control_val) / control_val) * 100
            labelstr = f"{barlabel:.2f}%"
            print(labelstr)
            inset.bar(gmuls[i], thresholds[i])
            if barlabel != 0:
                inset.annotate(text=labelstr, xy=(gmuls[i], thresholds[i]),
                               xytext=(gmuls[i] + 0.5, thresholds[i] - 0.01))

        inset.spines[['right', 'top']].set_visible(False)
        inset.set_xlabel("g mul")
        inset.set_ylabel("I (nA)")
        inset.set_yticks([0, 0.05])
        inset.set_xticks(gmuls)
    fig.savefig(f"{simdir}-F-I.png")
    plt.show()


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

    plot_if(sys.argv[1])
