#!/usr/bin/env python3
"""
Translate cell such that the soma is at (0, 0, 0)

File: translatecell.py

Copyright 2023 Ankur Sinha
Author: Ankur Sinha <sanjay DOT ankur AT gmail DOT com>
"""

from pyneuroml.pynml import read_neuroml2_file, write_neuroml2_file
import neuroml


hl5doc = read_neuroml2_file("./HL5PC.cell.nml")
hl5cell = hl5doc.cells[0]  # type: neuroml.Cell

soma_seg = hl5cell.get_segment(0)
x, y, z = soma_seg.proximal.x, soma_seg.proximal.y, soma_seg.proximal.z
for seg in hl5cell.morphology.segments:
    try:
        seg.proximal.x -= x
        seg.proximal.y -= y
        seg.proximal.z -= z
    except AttributeError:
        pass

    try:
        seg.distal.x -= x
        seg.distal.y -= y
        seg.distal.z -= z
    except AttributeError:
        pass

write_neuroml2_file(hl5doc, "./HL5PC.translated.cell.nml", validate=False)
