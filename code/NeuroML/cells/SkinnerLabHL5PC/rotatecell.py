#!/usr/bin/env python3
"""
Enter one line description here.

File:

Copyright 2023 Ankur Sinha
Author: Ankur Sinha <sanjay DOT ankur AT gmail DOT com>
"""


from pyneuroml.pynml import read_neuroml2_file, write_neuroml2_file
from pyneuroml.utils import rotate_cell
from neuroml.utils import component_factory
import math


hl5doc = read_neuroml2_file("./HL5PC.cell.nml")
hl5cell = hl5doc.cells[0]
rotated_cell = rotate_cell(hl5cell, x=0, y=0, z=math.pi, relative_to_soma=True)
hl5doc.cells[0] = rotated_cell

write_neuroml2_file(hl5doc, "./HL5PC.cell.nml", validate=False)
