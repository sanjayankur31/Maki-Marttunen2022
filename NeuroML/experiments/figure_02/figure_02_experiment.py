#!/usr/bin/env python3
"""
Methods used for experiments in Figure 2
These are used by the individual scripts, they are not run using this file.

File: NeuroML/experiments/figure_02/figure_02_experiment.py

Copyright 2023 Ankur Sinha
Author: Ankur Sinha <sanjay DOT ankur AT gmail DOT com>
"""

import numpy as np
import json

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
    first_apical_segment_at_500 = None
    fraction_along_segment_at_500 = None
    for seg, frac_along in segments_at_500.items():
        # just take the first one
        if seg in apical_segments:
            first_apical_segment_at_500 = seg
            fraction_along_segment_at_500 = frac_along
            break

    # get distances from soma to get segments at
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

    # remove 500 if it's also in the distances
    distances_for_inputs = distances_for_inputs[distances_for_inputs != 500]

    # create new population: one cell each with an input at a different
    # location, so we need to know the number of cells (number of points for
    # inputs, calculated above)
    popsize = (1 + len(distances_for_inputs))
    print(f"Population size is {popsize}")
    population = network.add(
        neuroml.Population,
        id=f"{cellname}_pop",
        component=f"{cellname}",
        size=popsize
    )

    # input component
    pg = nml_doc.add(
        neuroml.PulseGenerator,
        id="pulseGen_%i" % 0,
        delay="5ms",
        duration="0.2ms",
        amplitude=current_nA,
    )

    # input list
    inputlist = network.add(
        neuroml.InputList, id="i1", component=pg.id, populations=population.id,
        validate=False
    )

    # add input at 500 first
    input_id = 0
    inputlist.add(
        neuroml.Input, id=input_id, target=f"../{population.id}[0]/{cellname}",
        segment_id=first_apical_segment_at_500,
        fraction_along=fraction_along_segment_at_500,
        destination="synapses"
    )

    # also use the same pulse generator to provide input to other cells at
    # different distances from soma
    dist_vs_seg_dict = {}
    dist_vs_seg_dict["500.0"] = first_apical_segment_at_500

    for p in range(1, popsize):
        s_at_d = cell.get_segments_at_distance(distances_for_inputs[p - 1])
        first_apical_segment_at_d = None
        frac_along_at_d = None
        for seg, frac_along in s_at_d.items():
            if seg in apical_segments:
                first_apical_segment_at_d = seg
                fraction_along_at_d = frac_along
                break

        if first_apical_segment_at_d is None:
            print(f"No segment found at distance {distances_for_inputs[p - 1]}")
            continue

        print(f"{distances_for_inputs[p - 1]}: {first_apical_segment_at_d}")
        dist_vs_seg_dict[distances_for_inputs[p - 1]] = first_apical_segment_at_d

        input_id += 1
        inputlist.add(
            neuroml.Input, id=input_id, target=f"../{population.id}[{p}]/{cellname}",
            segment_id=first_apical_segment_at_d,
            fraction_along=fraction_along_at_d,
            destination="synapses"
        )

    # write segments at different distances to a file, useful later when
    # plotting and so on
    # should be identical for all, since we're using only the one cell
    with open(f"{nml_doc_name}.segs.json", 'w') as f:
        f.write(json.dumps(dist_vs_seg_dict))

    nml_doc_name += ".net.nml"
    write_neuroml2_file(nml_doc, nml_doc_name)
    print("Written network file to: " + nml_doc_name)
    return nml_doc_name


def simulate_model(model_file_name: str, cellname: str,
                   plot: bool = True,
                   skip_run=True, duration_ms=10) -> str:
    """Simulate the model, generating current plots if required.

    :param model_file_name: name of model file
    :type model_file_name: str
    :param cellname: name of cell
    :type cellname: str
    :param plot: toggle plotting
    :type plot: bool
    :returns: name of LEMS simulation file
    :rtype: str
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

    model = read_neuroml2_file(model_file_name)
    inputlist = model.networks[0].input_lists[0].input
    num_cells = 0
    for aninput in inputlist:
        simulation.add_column_to_output_file(
            "output0", f"{cellname}_pop_{aninput.segment_id}", f"{cellname}_pop[{num_cells}]/v"
        )
        num_cells = num_cells + 1

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
        segs_data = {}
        with open(f"{simulation_id}.segs.json", 'r') as f:
            segs_data = json.load(f)

        x_vals = [data_array[:, 0]] * num_cells
        y_vals = []

        for i in range(0, num_cells):
            y_vals.append(data_array[:, i])

        generate_plot(
            xvalues=x_vals,
            yvalues=y_vals,
            labels=list(segs_data.keys()),
            title="Membrane potential",
            show_plot_already=False,
            save_figure_to="%s-v.png" % simulation_id,
            xaxis="time (s)",
            yaxis="membrane potential (V)",
            ylim=[-0.085, 0.065]
        )

    return lems_simulation_file
