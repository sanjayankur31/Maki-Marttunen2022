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
import matplotlib

import neuroml
from neuroml.utils import component_factory
from neuroml.loaders import read_neuroml2_file
from pyneuroml.pynml import write_neuroml2_file
from pyneuroml import pynml
from pyneuroml.lems import LEMSSimulation
from pyneuroml.plot import generate_plot

from figure_01.figure_01_experiment import (create_modified_cell)
from common import (delete_neuron_special_dir, get_timestamp)


# increase plot size
matplotlib.rcParams['figure.figsize'] = [19.2, 10.8]

def get_segments_at_distances(celldir: str, cellname: str, segment_group: str,
                              distance_from_soma_inc: float, extra_points =
                              list[float]):
    """Get dict of thickest segment at various distances from the soma, and
    their fraction along values

    :param celldir: cell directory
    :type celldir: str
    :param cellname: cell name
    :type cellname: str
    :param segment_group: segment group to limit to
    :type segment_group: str
    :param distance_from_soma_inc: incremental value to find points at
    :type distance_from_soma_inc: float
    :param extra_points: extra points to include at starting of list
    :type extra_points: list of floats
    :returns: dict with key as distance, value as [segment id, fraction along]

    """
    cell_doc = read_neuroml2_file(f"{celldir}/{cellname}.cell.nml")
    # we could even store the cell in a separate file and include that, but
    # since we're only working with one file, we can include its NeuroML in our
    # network file here
    cell = cell_doc.cells[0]  # type: neuroml.Cell

    apical_segments = cell.get_all_segments_in_group("apical_dendrite_group")
    dist_vs_seg_dict = {}

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

    distances = []
    # complete list
    distances.extend(extra_points)
    distances.extend(distances_for_inputs)

    for d in distances:
        s_at_d = cell.get_segments_at_distance(d)
        thickest_segment_at_d = None
        fraction_along_at_d = None
        segment_width_at_d = 0
        for seg, frac_along in s_at_d.items():
            thisseg_width_at_d = 0
            thisseg = cell.get_segment(seg)
            # proximal diam + frac_along * delta(diam)
            thisseg_width_at_d = (cell.get_actual_proximal(seg).diameter + ((abs(thisseg.distal.diameter - cell.get_actual_proximal(seg).diameter)/cell.get_segment_length(seg)) * frac_along))
            if (seg in apical_segments) and (thisseg_width_at_d > segment_width_at_d):
                thickest_segment_at_d = seg
                fraction_along_at_d = frac_along
                segment_width_at_d = thisseg_width_at_d
                break

        if thickest_segment_at_d is None:
            print(f"No segment found at distance {d}")
            continue

        print(f"{d}: {thickest_segment_at_d}, {fraction_along_at_d}")
        dist_vs_seg_dict[d] = [thickest_segment_at_d, fraction_along_at_d]

    # write segments at different distances to a file, useful later when
    # plotting
    with open(f"{cellname}.segs.json", 'w') as f:
        f.write(json.dumps(dist_vs_seg_dict))

    return dist_vs_seg_dict


def create_model(
    cellname: str,
    celldir: str,
    dist_vs_seg_dict: dict,
    current_nA: str = "0.5nA",
    g_Ih_multiplier: str = None,
    g_Ca_LVAst_multiplier: str = None,
) -> str:
    """Create model with different cells receiving inputs at different parts,
    for a particular set of Ih and CA_LVast maximal conductances.

    :param cellname: name of cell
    :param celldir: directory of cell nml file
    :param current_nA: current to apply in nA
    :param g_Ih_multiplier: multiplier for Ih channels
    :param g_Ca_LVAst_multiplier: multiplier for Ca_LVAst
    :param dist_vs_seg_dict: dict with keys as distance from soma, and value as
        segment id
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
    # create new population: one cell each with an input at a different
    # location, so we need to know the number of cells (number of points for
    # inputs, calculated above)
    popsize = len(dist_vs_seg_dict.keys())
    print(f"Population size is {popsize}")
    population = network.add(
        neuroml.Population,
        id=f"{cellname}_pop",
        component=f"{cellname}",
        size=popsize
    )

    # input component
    print(f"Current input is {current_nA}")
    pg = nml_doc.add(
        neuroml.PulseGenerator,
        id="pulseGen_%i" % 0,
        delay="10ms",
        duration="0.2ms",  # from code
        amplitude=current_nA,
    )

    # input list
    inputlist = network.add(
        neuroml.InputList, id="i1", component=pg.id, populations=population.id,
        validate=False
    )

    # input id, and cell id
    ctr = 0
    for d, val in dist_vs_seg_dict.items():
        inputlist.add(
            neuroml.Input, id=ctr, target=f"../{population.id}[{ctr}]/{cellname}",
            segment_id=val[0],
            fraction_along=val[1],
            destination="synapses"
        )
        ctr += 1

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
        sim_id=simulation_id, duration=duration_ms, dt=0.01, simulation_seed=123
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
        max_memory="8G",
        nogui=True,
        plot=False,
        verbose=True,
        compile_mods=True,
        skip_run=skip_run
    )
    if not skip_run and plot:
        data_array = np.loadtxt("%s.v.dat" % simulation_id)
        segs_data = {}
        with open(f"{cellname}.segs.json", 'r') as f:
            segs_data = json.load(f)

        x_vals = [data_array[:, 0]] * num_cells
        y_vals = []

        for i in range(1, num_cells + 1):
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
