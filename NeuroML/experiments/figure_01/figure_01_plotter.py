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
from matplotlib import pyplot as plt

logger = logging.Logger(__name__)
logger.setLevel(logging.INFO)


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

    datafiles = list(simdir.glob("*_if.dat"))

    xvalues = []
    yvalues = []
    labels = []
    for afile in datafiles:
        # can be improved
        label = "_".join(re.sub(r"_(\d+)_(\d+)_(\d+)_(\d+)_", r" \1.\2 \3.\4 ", (re.sub(r"_m2.*", "_m2", afile.name.split(".")[0]))).split("_")[2:])
        labels.append(label)
        data = (numpy.loadtxt(afile))
        # convert nA to pA
        xvalues.append(data[:, 0] / 1000)
        yvalues.append(data[:, 1])

    generate_plot(xvalues,
                  yvalues,
                  title="F-I curve for different Ih/CaLVAst conductances",
                  xaxis="I(nA)", yaxis="f(spikes/s)",
                  show_plot_already=False, labels=labels)
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
