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

nml_doc = read_neuroml2_file("HL5PC.cell.nml")
thecell = nml_doc.cells[0]
rotated_cell = rotate_cell(thecell, x=0, y=0, z=math.pi, relative_to_soma=True)

plot_2D(rotated_cell, save_to_file="HL5PC.rotated.png", min_width=6)
