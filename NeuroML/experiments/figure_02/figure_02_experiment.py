#!/usr/bin/env python3
"""
Methods used for experiments in Figure 2
These are used by the individual scripts, they are not run using this file.

File: NeuroML/experiments/figure_02/figure_02_experiment.py

Copyright 2023 Ankur Sinha
Author: Ankur Sinha <sanjay DOT ankur AT gmail DOT com>
"""

import numpy as np

import neuroml
from neuroml.utils import component_factory
from neuroml.loaders import read_neuroml2_file
from pyneuroml.pynml import write_neuroml2_file
from pyneuroml import pynml
from pyneuroml.lems import LEMSSimulation
from pyneuroml.plot import generate_plot

from figure_01.figure_01_experiment import (create_modified_cell,
                                            delete_neuron_special_dir,
                                            get_timestamp)


def create_model(
    cellname: str,
    celldir: str,
    current_nA: str = "0.5nA",
    g_Ih_multiplier: str = None,
    g_Ca_LVAst_multiplier: str = None,
    distance_from_soma_inc: int = 0
) -> str:
    """Create model with different cells receiving inputs at different parts,
    for a particular set of Ih and CA_LVast maximal conductances.

    :param cellname: name of cell
    :param celldir: directory of cell nml file
    :param current_nA: current to apply in nA
    :param g_Ih_multiplier: multiplier for Ih channels
    :param g_Ca_LVAst_multiplier: multiplier for Ca_LVAst
    :param distance_from_soma_inc: incremental distance from soma at which to
        place inputs (inputs will be placed at N, 2N, 3N, and so on, use 0 to
        disable)
    :returns: network model file

    """
    timestamp = get_timestamp()
    nml_doc_name = f"net_{timestamp}_{cellname}"
    if g_Ih_multiplier is not None:
        nml_doc_name += f"_{g_Ih_multiplier}"
    if g_Ca_LVAst_multiplier is not None:
        nml_doc_name += f"_{g_Ca_LVAst_multiplier}"
    nml_doc_name += f"_{current_nA}"

    nml_doc_name = nml_doc_name.replace(".", "_").replace(" ", "_")
    nml_doc = component_factory(neuroml.NeuroMLDocument, id=nml_doc_name)

    cell_doc = read_neuroml2_file(f"{celldir}/{cellname}.cell.nml")
    # we could even store the cell in a separate file and include that, but
    # since we're only working with one file, we can include its NeuroML in our
    # network file here
    cell = cell_doc.cells[0]  # type: neuroml.Cell
    nml_doc.add(cell)

    create_modified_cell(cell, g_Ih_multiplier, g_Ca_LVAst_multiplier)

    # include channel file definitions
    for inc in cell_doc.includes:
        updated_path = f"{celldir}/{inc.href}"
        nml_doc.add(neuroml.IncludeType, href=updated_path)

    # create a population of N cells, each with an input at different location
    network = nml_doc.add(neuroml.Network, id=f"{cellname}_net", validate=False)
    apical_segments = cell.get_all_segments_in_group("apical_dendrite_group")

    # an input at 500 from soma
    segments_at_500 = cell.get_segments_at_distance(500)
    apical_segments_at_500 = {}
    for seg, frac_along in segments_at_500.items():
        if seg in apical_segments:
            apical_segments_at_500[seg] = frac_along

    a_segment_at_500 = list(segments_at_500.keys())[0]
    fraction_along_segment_at_500 = apical_segments_at_500[a_segment_at_500]

    if distance_from_soma_inc > 0:
        extremeties = cell.get_extremeties()

        furthest_tip = max(list(extremeties.values()))
        # create 25 points
        distances_for_inputs = np.arange(
            0, furthest_tip, distance_from_soma_inc,
            dtype=float
        ).round(3)
        # drop first point, linspace is [start, stop]
        distances_for_inputs = distances_for_inputs[1:]

    # create new population
    popsize = (1 + len(distances_for_inputs))
    print(f"Population size is {popsize}")
    population = network.add(
        neuroml.Population,
        id=f"{cellname}_pop",
        component=f"{cellname}",
        size=popsize
    )

    # add at 500, first
    pg = nml_doc.add(
        neuroml.PulseGenerator,
        id="pulseGen_%i" % 0,
        delay="5ms",
        duration="0.2ms",
        amplitude=current_nA,
    )
    network.add(
        neuroml.ExplicitInput, target=f"{population.id}[0]/{a_segment_at_500}", input=pg.id
    )

    # also use the same pulse generator to provide input to other cells at
    # different distances from soma
    for p in range(1, popsize):
        s_at_d = cell.get_segments_at_distance(distances_for_inputs[p - 1])
        apical_segments_at_d = {}
        for seg, frac_along in s_at_d.items():
            if seg in apical_segments:
                apical_segments_at_d[seg] = frac_along

        a_segment_at_d = list(apical_segments_at_d.keys())[0]
        if a_segment_at_d is None:
            print(f"No segment found at distance {distances_for_inputs[p - 1]}")
            continue
        print(f"{distances_for_inputs[p - 1]}: {a_segment_at_d}")
        print(cell.get_segment_location_info(a_segment_at_d))

        network.add(
            neuroml.ExplicitInput, target=f"{population.id}[{p}]/{a_segment_at_d}", input=pg.id
        )

    nml_doc_name += ".net.nml"
    write_neuroml2_file(nml_doc, nml_doc_name)
    print("Written network file to: " + nml_doc_name)
    return nml_doc_name


def simulate_model(model_file_name: str, cellname: str, num_cells: int = 1,
                   plot: bool = True,
                   skip_run=True, duration_ms=10):
    """Simulate the model, generating current plots if required.

    :param model_file_name: name of model file
    :type model_file_name: str
    :param cellname: name of cell
    :type cellname: str
    :param plot: toggle plotting
    :type plot: bool
    :returns: None
    """
    delete_neuron_special_dir()

    network_id = f"{cellname}_net"
    simulation_id = model_file_name.split(".")[0]
    simulation = LEMSSimulation(
        sim_id=simulation_id, duration=duration_ms, dt=0.1, simulation_seed=123
    )
    simulation.assign_simulation_target(network_id)
    simulation.include_neuroml2_file(model_file_name)

    simulation.create_output_file("output0", "%s.v.dat" % simulation_id)
    for i in range(0, num_cells):
        simulation.add_column_to_output_file(
            "output0", f"{cellname}_pop", f"{cellname}_pop[{i}]/v"
        )

    lems_simulation_file = simulation.save_to_file()

    pynml.run_lems_with_jneuroml_neuron(
        lems_simulation_file,
        max_memory="2G",
        nogui=True,
        plot=False,
        verbose=True,
        compile_mods=True,
        skip_run=skip_run
    )
    if not skip_run and plot:
        data_array = np.loadtxt("%s.v.dat" % simulation_id)
        generate_plot(
            [data_array[:, 0]],
            [data_array[:, 1]],
            "Membrane potential",
            show_plot_already=False,
            save_figure_to="%s-v.png" % simulation_id,
            xaxis="time (s)",
            yaxis="membrane potential (V)",
        )
