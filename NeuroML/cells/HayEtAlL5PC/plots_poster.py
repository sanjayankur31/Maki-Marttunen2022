#!/usr/bin/env python3
"""
Enter one line description here.

File:

Copyright 2023 Ankur Sinha
Author: Ankur Sinha <sanjay DOT ankur AT gmail DOT com>
"""

import matplotlib
from matplotlib import pyplot as plt
from pyneuroml.plot.PlotMorphology import plot_2D

matplotlib.rcParams['font.weight'] = "bold"
matplotlib.rcParams['axes.labelweight'] = "bold"
matplotlib.rcParams['axes.titlesize'] = "large"
matplotlib.rcParams['axes.titleweight'] = "bold"


inputfile = "L5PC.morph.cell.nml"
title = "Rodent model\n(Hay)"
plane = "xy"
plot_2D(inputfile, min_width=6, nogui=True, plane2d=plane,
        close_plot=False, title=title, square=False)
fig = plt.gcf()
plt.savefig(inputfile.replace(".cell.nml", f"{plane}.poster.png"), dpi=200, bbox_inches="tight")
