#!/usr/bin/env python3
"""
Mark segments where inputs are being placed.

File: figure_02_mark_segments.py

Copyright 2024 Ankur Sinha
Author: Ankur Sinha <sanjay DOT ankur AT gmail DOT com>
"""


import sys
from pyneuroml.plot.PlotMorphology import plot_2D_cell_morphology
from pyneuroml.io import read_neuroml2_file
import json


def run(
    cellfile: str,
    segment_dict_file: str,
    marker_color: str = "blue",
    marker_size: int = 10,
):
    """Run plotter

    :param cellfile: path to cell file
    :type cellfile: str
    :param segment_dict: dict geneated by experiment
    :type segment_dict: dict
    """
    cell_doc = read_neuroml2_file(cellfile)
    cell = cell_doc.cells[0]

    segment_dict = None
    with open(segment_dict_file, "r") as f:
        segment_dict = json.load(f)

    segment_spec = {}
    for vals in segment_dict.values():
        segment_id = vals[0]
        segment_spec[segment_id] = {
            "marker_color": marker_color,
            "marker_size": marker_size,
        }

    plot_2D_cell_morphology(cell=cell, plot_type="constant", highlight_spec={f"{cell.id}": segment_spec})


if __name__ == "__main__":
    if len(sys.argv) > 2:
        print("Usage: script <path to cell file> <path to segments json>")
    run(sys.argv[0], sys.argv[1])
