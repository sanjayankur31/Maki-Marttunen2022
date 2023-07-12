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
matplotlib.rcParams['font.size'] = 35
matplotlib.rcParams['axes.labelweight'] = "bold"
matplotlib.rcParams['axes.titlesize'] = "xx-large"
matplotlib.rcParams['axes.titleweight'] = "bold"
matplotlib.rcParams['lines.linewidth'] = 3
matplotlib.rcParams['legend.fontsize'] = 25
matplotlib.rcParams['axes.labelsize'] = 35
matplotlib.rcParams['xtick.labelsize'] = 30
matplotlib.rcParams['ytick.labelsize'] = 30


def plot_if(simfolder: str) -> None:
    """Plot I-F curves from data in simfolder.

    The datafiles must end in `_if.dat`. This will read all such files and
    generate a line for each.

    :param simfolder: name of folder containing simulation and generated data
    :type simfolder: str
    :raises ValueError: if a dir named `simfolder` does not exist (or cannot be
        accessed for some reason)
    :returns: None

    """
    simdir = Path(simfolder)
    if not simdir.exists() or not simdir.is_dir():
        raise ValueError(f"A directory named {simfolder} does not exist or cannot be accessed")

    datafiles = sorted(list(simdir.glob("net*_if.dat")))

    # temporary hack to make the control for the healthy exps be plotted first
    if "ScZ" not in datafiles[0].parent.name:
        temp = datafiles[0]
        datafiles[0] = datafiles[1]
        datafiles[1] = temp

    xvalues = []
    yvalues = []
    labels = []
    thresholds = []
    gmuls = []
    control_val = 0.0
    camuls = []
    for afile in datafiles:
        print(afile)
        data = numpy.loadtxt(afile)
        threshold_file = str(afile.parent) + f"/threshold_i_{afile.name}"
        threshold_data = float(numpy.loadtxt(threshold_file))
        logger.debug("Processing %s", afile)
        # can be improved: should ideally store tags in model and get values
        # from there instead of parsing filenames
        if "ScZ" not in str(afile):
            label = re.sub(r"_(\d+)_(\d+)_(\d+)_(\d+)_", r" \1.\2 \3.\4 ", afile.name.split(".")[0]).split("_")[2:]
            cellname = label[0].split(" ")[0]
            gmul = float(label[0].split(" ")[1])  # type: float
            if gmul == 1:
                control_val = threshold_data

            flabel = f"g_Ih * {gmul:.3f}"
        else:
            label = re.sub(r"_(\d+)_(\d+)_(\d+)_(\d+)_(\d+)_(\d+)_", r" \1.\2 \3.\4 \5.\6 ", afile.name.split(".")[0]).split("_")[2:]
            cellname = label[0].split(" ")[0]
            gmul = float(label[0].split(" ")[1])  # type: float
            camul = float(label[0].split(" ")[2])  # type: float
            camuls.append(camul)

            if gmul == 1 and camul == 1:
                control_val = threshold_data
            flabel = f"g_Ih * {gmul:.3f}, g_CaLVast * {camul:.3f}"

        logger.debug("Label: %s", flabel)
        labels.append(flabel)
        thresholds.append(threshold_data)
        gmuls.append(gmul)
        # convert nA to pA
        xvalues.append(data[:, 0] / 1000)
        yvalues.append(data[:, 1])

    logger.debug(thresholds)
    logger.debug(gmuls)

    # tweak title
    title = "    "
    if cellname == "L5PC":
        title += "Rodent "
    else:
        title += "Human "

    if "ScZ" in datafiles[0].parent.name:
        title += "(SCZ)"
    else:
        title += "(Health)"

    generate_plot(xvalues,
                  yvalues,
                  linewidths=[5] * len(xvalues),
                  title=None,
                  xaxis="I(nA)", yaxis="f(spikes/s)",
                  show_plot_already=False, labels=labels,
                  bottom_left_spines_only=True, close_plot=False,
                  cols_in_legend_box=1, legend_position="lower right",
                  title_above_plot=True, show_xticklabels=False)

    # skip tick labels for clarity
    ticks = len(xvalues[0])
    xtick_labels = []
    for i in range(ticks):
        if i % int(ticks / 4) == 0:
            xtick_labels.append(f"{xvalues[0][i]}")
        else:
            xtick_labels.append("")
    plt.xticks(xvalues[0], labels=xtick_labels)
    plt.title(title, fontsize=60)

    # add inset with threshold values for non ScZ sims
    # https://matplotlib.org/stable/gallery/subplots_axes_and_figures/axes_demo.html#sphx-glr-gallery-subplots-axes-and-figures-axes-demo-py
    # https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.gcf.html
    fig = plt.gcf()
    inset = fig.add_axes([0.18, 0.56, 0.1, 0.4])
    # if only g_Ih was changed, different inset
    if len(camuls) == 0:
        # to get the same random colours that matplotlib uses
        for i in range(0, len(gmuls)):
            barlabel = ((thresholds[i] - control_val) / control_val) * 100
            labelstr = f"{barlabel:+.2f}%"
            print(labelstr)
            inset.bar(gmuls[i], thresholds[i])
            if barlabel != 0:
                inset.annotate(text=labelstr, xy=(gmuls[i], thresholds[i]),
                               xytext=(i * 0.9, 1.04 * thresholds[i]),
                               fontsize=25)

        inset.spines[['right', 'top']].set_visible(False)
        inset.set_xlabel("g_Ih*")
        inset.set_ylabel("I (nA)")
        inset.set_xticks(gmuls)
    else:
        xtics = []
        for i in range(0, len(gmuls)):
            barlabel = ((thresholds[i] - control_val) / control_val) * 100
            labelstr = f"{barlabel:+.2f}%"
            inset.bar(i, thresholds[i])
            if gmuls[i] == 1 and camuls[i] == 1:
                xtics.append("Control")
            else:
                xtics.append(f"g_Ih*{gmuls[i]:.2f}, g_Ca*{camuls[i]:.2f}")
            if barlabel != 0:
                inset.annotate(text=labelstr, xy=(i, thresholds[i]),
                               xytext=(i * 0.9, 1.04 * thresholds[i]),
                               fontsize=25)
        inset.set_xlabel("")
        inset.set_xticks(ticks=list(range(len(gmuls))), labels=xtics, rotation="vertical")

    inset.spines[['right', 'top']].set_visible(False)
    inset.set_ylabel("I (nA)")
    plt.tight_layout()
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
