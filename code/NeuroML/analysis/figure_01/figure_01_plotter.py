#!/usr/bin/env python3
"""
Plot graphs from simulation data

File: figure_01_plotter.py

Copyright 2023 Ankur Sinha
Author: Ankur Sinha <sanjay DOT ankur AT gmail DOT com>
"""


from pathlib import Path
from natsort import realsorted
import logging
import numpy
from pyneuroml.plot.Plot import generate_plot
import matplotlib
from matplotlib import pyplot as plt

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

# increase plot size
matplotlib.rcParams['figure.figsize'] = [19.2, 10.8]


def __get_gmul_camul(afile, scz=True):
    """Get values from file name"""
    if scz is False:
        label = afile.name.split(".")[0].split("_")
        logger.debug(f"Split {afile.name} into labels: {label}")
        cellname = label[2]
        gmul = float(f"{label[3]}.{label[4]}")  # type: float
        camul = None
    else:
        label = afile.name.split(".")[0].split("_")
        logger.debug(f"Split {afile.name} into labels: {label}")
        cellname = label[2]
        gmul = float(f"{label[3]}.{label[4]}")  # type: float
        camul = float(f"{label[5]}.{label[6]}")  # type: float
    return (cellname, gmul, camul)


def __comparison_key(afile):
    keys = __get_gmul_camul(afile)
    logger.debug(f"Sorting by: ({keys[1]}, {keys[2]})")
    return (keys[1], keys[2])


def plot_if() -> None:
    """Plot I-F curves from data.

    The datafiles must end in `_if.dat`. This will read all such files and
    generate a line for each.

    :param simfolder: name of folder containing simulation and generated data
    :type simfolder: str
    :raises ValueError: if a dir named `simfolder` does not exist (or cannot be
        accessed for some reason)
    :returns: None

    """
    simdir = Path(".")
    simdir_name = simdir.absolute().name

    # sort by gmul: can't get it to sort by both
    datafiles = realsorted(list(simdir.glob("**/net*_if.dat")), key=__comparison_key)

    logger.debug(f"Sorted order is: {datafiles}")

    xvalues = []
    yvalues = []
    labels = []
    thresholds = []
    gmuls = []
    control_val = 0.0
    camuls = []
    flabel = ""
    for afile in datafiles:
        print(afile.name)
        data = numpy.loadtxt(afile)
        threshold_file = str(afile.parent) + f"/threshold_i_{afile.name}"
        threshold_data = float(numpy.loadtxt(threshold_file))
        logger.debug("Processing %s", afile)
        (cellname, gmul, camul) = __get_gmul_camul(afile, "ScZ" in simdir_name)

        if camul is None:
            flabel = f"g_Ih * {gmul:.3f}"
            if gmul == 1:
                control_val = threshold_data

        else:
            # current = f"{label[5]}.{label[6]}"
            camuls.append(camul)
            flabel = f"g_Ih * {gmul:.3f}, g_CaLVast * {camul:.3f}"
            if gmul == 1 and camul == 1:
                control_val = threshold_data

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
    title = "F-I curve"
    if cellname == "L5PC":
        title += ": rodent "
    else:
        title += ":  human "

    if "ScZ" in str(simdir_name):
        title += "(SCZ)"
    else:
        title += "(Health)"

    generate_plot(xvalues,
                  yvalues,
                  title=title,
                  xaxis="I(nA)", yaxis="f(spikes/s)",
                  show_plot_already=False, labels=labels,
                  bottom_left_spines_only=True, close_plot=False,
                  cols_in_legend_box=1, legend_position="lower right",
                  title_above_plot=True)

    # add inset with threshold values for non ScZ sims
    # https://matplotlib.org/stable/gallery/subplots_axes_and_figures/axes_demo.html#sphx-glr-gallery-subplots-axes-and-figures-axes-demo-py
    # https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.gcf.html
    fig = plt.gcf()
    inset = fig.add_axes([0.15, 0.5, 0.1, 0.4])
    # track labels, only print them once to prevent cluttering of plot
    labels = []
    # if only g_Ih was changed, different inset
    if len(camuls) == 0:
        # to get the same random colours that matplotlib uses
        for i in range(0, len(gmuls)):
            barlabel = ((thresholds[i] - control_val) / control_val) * 100
            labelstr = f"{barlabel:+.2f}%"
            logger.debug(f"labelstr is: {labelstr}")
            inset.bar(gmuls[i], thresholds[i])
            if labelstr not in labels:
                labels.append(labelstr)
                inset.annotate(text=labelstr, xy=(gmuls[i], thresholds[i]),
                               xytext=(gmuls[i] + 0.5, thresholds[i] - 0.01))
        inset.set_xlabel("g_Ih *")
        inset.set_xticks(gmuls)
    else:
        xtics = []
        for i in range(0, len(gmuls)):
            barlabel = ((thresholds[i] - control_val) / control_val) * 100
            labelstr = f"{barlabel:+.2f}%"
            logger.debug(f"labelstr is: {labelstr}")
            inset.bar(i, thresholds[i])
            if gmuls[i] == 1 and camuls[i] == 1:
                xtics.append("Control")
            else:
                xtics.append(f"g_Ih * {gmuls[i]:.2f}, g_Ca_LVAst * {camuls[i]:.2f}")
            if labelstr not in labels:
                labels.append(labelstr)
                inset.annotate(text=labelstr, xy=(i, thresholds[i]),
                               xytext=(i * 1.2, 0.9 * thresholds[i]))
        inset.set_xlabel("")
        inset.set_xticks(ticks=list(range(len(gmuls))), labels=xtics, rotation="vertical")

    inset.spines[['right', 'top']].set_visible(False)
    inset.set_ylabel("I (nA)")
    plt.tight_layout()
    fig.savefig(f"{simdir_name}-F-I.png")
    plt.show()


if __name__ == "__main__":
    plot_if()
