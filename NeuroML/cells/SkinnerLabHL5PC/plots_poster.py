#!/usr/bin/env python3
"""
Flip cell to make it upright

File: flip.py

Copyright 2023 Ankur Sinha
Author: Ankur Sinha <sanjay DOT ankur AT gmail DOT com>
"""


import math
import neuroml
from neuroml.loaders import read_neuroml2_file
from neuroml.writers import NeuroMLWriter
try:
    from pyneuroml.utils import rotate_cell
except ImportError:
    pass
from pyneuroml.plot.PlotMorphology import plot_2D

import matplotlib
from matplotlib import pyplot as plt

matplotlib.rcParams['font.weight'] = "bold"
matplotlib.rcParams['axes.labelweight'] = "bold"
matplotlib.rcParams['axes.titlesize'] = "large"
matplotlib.rcParams['axes.titleweight'] = "bold"


inputfile = "HL5PC.cell.nml"
title = "Human\n(Rich)"
plane = "xy"
nml_doc = read_neuroml2_file("HL5PC.cell.nml")
thecell = nml_doc.cells[0]
rotated_cell = rotate_cell(thecell, x=0, y=0, z=math.pi, relative_to_soma=True)

plot_2D(rotated_cell, min_width=6, nogui=True, plane2d=plane,
        close_plot=False, title=title, square=False)
fig = plt.gcf()
plt.title(title, fontsize=18)
plt.savefig(inputfile.replace(".cell.nml", f".{plane}.poster.png"), dpi=200, bbox_inches="tight")
