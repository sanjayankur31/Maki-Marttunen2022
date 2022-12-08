#!/usr/bin/env python3
"""
Simulations to generate figure 1

File: sim_fig1.py

Copyright 2022 Ankur Sinha
Author: Ankur Sinha <sanjay DOT ankur AT gmail DOT com>
"""


import datetime
import numpy as np

import neuroml
from neuroml.utils import component_factory
from neuroml.loaders import read_neuroml2_file
from pyneuroml.pynml import write_neuroml2_file
from pyneuroml import pynml
from pyneuroml.lems import LEMSSimulation
from pyneuroml.plot import generate_plot
from pyneuroml.analysis import generate_current_vs_frequency_curve


def create_model(
    cellname: str, celldir: str, current_nA: str = "0.5nA", gIh_S_per_m2: str = None
) -> str:
    """Create model with different values of conductance for the Ih channel.

    :param cellname: TODO
    :param gIh_S_per_m2: conductance for Ih channels in gIh_S_per_m2
    :returns: network model file

    """
    timestamp = datetime.datetime.now().strftime("%Y%m%d%H%M%S")
    nml_doc_name = f"net_{timestamp}_{cellname}_{gIh_S_per_m2}_{current_nA}".replace(
        ".", "_"
    ).replace(" ", "_")
    nml_doc = component_factory(neuroml.NeuroMLDocument, id=nml_doc_name)

    cell_doc = read_neuroml2_file(f"{celldir}/{cellname}.cell.nml")
    # we could even store the cell in a separate file and include that, but
    # since we're only working with one file, we can include its NeuroML in our
    # network file here
    cell = cell_doc.cells[0]  # type: neuroml.Cell
    nml_doc.add(cell)

    # modify Ih conductance
    if gIh_S_per_m2 is not None:
        ih_somatic = None  # type: neuroml.ChannelDensity
        ih_basal = None  # type: neuroml.ChannelDensity
        ih_apical = None  # type: neuroml.ChannelDensityNonUniform

        # get Ih channel densities
        for cd in (
            cell.biophysical_properties.membrane_properties.channel_densities
            + cell.biophysical_properties.membrane_properties.channel_density_non_uniforms
        ):
            if cd.id == "Ih_somatic":
                ih_somatic = cd
            elif cd.id == "Ih_basal":
                ih_basal = cd
            elif cd.id == "Ih_apical":
                ih_apical = cd

        ih_somatic.cond_density = gIh_S_per_m2
        ih_basal.cond_density = gIh_S_per_m2

        # ih_apical is a non homogenous param, so needs a little more tinkering
        # the expression is of the form: gbat * <exp etc dependent on distance>
        value_exp = ih_apical.variable_parameters[
            0
        ].inhomogeneous_value.value  # type: str
        [gbar, expression] = value_exp.split("*", maxsplit=1)
        new_value_exp = f"{gIh_S_per_m2.replace('S_per_m2', '')} * {expression}"
        ih_apical.variable_parameters[0].inhomogeneous_value.value = new_value_exp

    # include channel file definitions
    for inc in cell_doc.includes:
        updated_path = f"{celldir}/{inc.href}"
        nml_doc.add(neuroml.IncludeType, href=updated_path)

    network = nml_doc.add(neuroml.Network, id=f"{cellname}_net", validate=False)
    population = network.add(
        neuroml.Population,
        id=f"{cellname}_pop",
        component=f"{cellname}",
        size=1,
    )

    pg = nml_doc.add(
        neuroml.PulseGenerator,
        id="pulseGen_%i" % 0,
        delay="50ms",
        duration="300ms",
        amplitude=current_nA,
    )
    exp_input = network.add(
        neuroml.ExplicitInput, target=f"{population.id}[0]", input=pg.id
    )

    nml_doc_name += ".net.nml"
    write_neuroml2_file(nml_doc, nml_doc_name)
    print("Written network file to: " + nml_doc_name)
    return nml_doc_name


def simulate_model(model_file_name: str, cellname: str, plot: bool = True):
    """Simulate the model, generating current plots if required.

    :param model_file_name: name of model file
    :type model_file_name: str
    :param cellname: name of cell
    :type cellname: str
    :param plot: toggle plotting
    :type plot: bool
    :returns: None
    """
    network_id = f"{cellname}_net"
    simulation_id = model_file_name.split(".")[0]
    simulation = LEMSSimulation(
        sim_id=simulation_id, duration=500, dt=0.1, simulation_seed=123
    )
    simulation.assign_simulation_target(network_id)
    simulation.include_neuroml2_file(model_file_name)

    simulation.create_output_file("output0", "%s.v.dat" % simulation_id)
    simulation.add_column_to_output_file(
        "output0", f"{cellname}_pop", f"{cellname}_pop[0]/v"
    )

    lems_simulation_file = simulation.save_to_file()

    pynml.run_lems_with_jneuroml_neuron(
        lems_simulation_file,
        max_memory="2G",
        nogui=True,
        plot=False,
        verbose=True,
        compile_mods=True,
    )
    if plot:
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


if __name__ == "__main__":
    # 2 is "normal"
    for g in [0, 2, 4]:
        model_file_name = create_model(
            cellname="L5PC",
            celldir="../cells/HayEtAlL5PC/",
            current_nA="0.5 nA",
            gIh_S_per_m2=f"{g} S_per_m2",
        )
        simulate_model(model_file_name, "L5PC")

        generate_current_vs_frequency_curve(
            nml2_file=model_file_name,
            cell_id="L5PC",
            custom_amps_nA=list(np.arange(0, 1.6, 0.1)),
            analysis_duration=2000,
            analysis_delay=200,
            temperature="34 degC",
            simulator="jNeuroML_NEURON",
            plot_if=True,
            save_if_figure_to=f"{model_file_name}_iv.png",
            save_if_data_to=f"{model_file_name}_iv.dat",
            show_plot_already=False,
            num_processors=8,
        )
