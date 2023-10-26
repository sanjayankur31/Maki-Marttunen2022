#!/usr/bin/env python3
"""
Methods used for experiments in Figure 1
These are used by the individual scripts, they are not run using this file.

File: NeuroML/experiments/figure_01_experiment.py

Copyright 2022 Ankur Sinha
Author: Ankur Sinha <sanjay DOT ankur AT gmail DOT com>
"""

import sys
import typing
import os
sys.path.append(os.path.dirname(os.path.abspath(".")))
sys.path.append(os.path.dirname(os.path.abspath("..")))

import random
import numpy as np
import matplotlib
from multiprocessing import Pool, Process

import neuroml
from neuroml.utils import component_factory
from neuroml.loaders import read_neuroml2_file
from pyneuroml.pynml import write_neuroml2_file
from pyneuroml.analysis import generate_current_vs_frequency_curve
from pyneuroml import pynml
from pyneuroml.lems import LEMSSimulation
from pyneuroml.plot import generate_plot
from common import (get_timestamp, get_relative_dir, delete_neuron_special_dir, get_run_dir, get_abs_celldir, data_means_cz)

try:
    from common.CM import scz_data
except ImportError:
    print("Data set not found, using defaults")
    scz_data = []

# increase plot size
matplotlib.rcParams['figure.figsize'] = [19.2, 10.8]


def create_modified_cell(cell, g_Ih_multiplier, g_Ca_LVAst_multiplier):
    """Modifies a cell object to multiply the Ih and Ca_LVAst maximal
    conductances with given multipliers

    :param cell: cell object to modify
    :param g_Ih_multiplier: multiplier for Ih channels
    :param g_Ca_LVAst_multiplier: multiplier for Ca_LVAst
    :returns: None

    """
    # modify Ih conductance
    if g_Ih_multiplier is not None:
        print(f"Ih multiplier: {g_Ih_multiplier}")
        value_exp = ""
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

        cond_den, unit = pynml.split_nml2_quantity(ih_somatic.cond_density)
        new_cond_den = float(g_Ih_multiplier) * float(cond_den)
        print(f"Ih somatic: {ih_somatic.cond_density} -> {new_cond_den}")
        ih_somatic.cond_density = f"{new_cond_den} {unit}"
        cond_den, unit = pynml.split_nml2_quantity(ih_basal.cond_density)
        new_cond_den = float(g_Ih_multiplier) * float(cond_den)
        print(f"Ih basal: {ih_basal.cond_density} -> {new_cond_den}")
        ih_basal.cond_density = f"{new_cond_den} {unit}"

        # ih_apical is a non homogenous param, so needs a little more tinkering
        # the expression is of the form: gbat * <exp etc dependent on distance>
        value_exp = ih_apical.variable_parameters[
            0
        ].inhomogeneous_value.value  # type: str
        # here we can just include the new multiplier in the expression
        new_value_exp = f"{g_Ih_multiplier} * {value_exp}"
        print(f"Ih apical: {ih_apical.variable_parameters[0].inhomogeneous_value.value} -> {new_value_exp}")
        ih_apical.variable_parameters[0].inhomogeneous_value.value = new_value_exp

    # Not present in basal dendrites
    if g_Ca_LVAst_multiplier is not None:
        print(f"Ca_LVAst multiplier: {g_Ca_LVAst_multiplier}")
        value_exp = ""
        Ca_LVAst_somatic = None  # type: neuroml.ChannelDensity
        Ca_LVAst_apical = None  # type: neuroml.ChannelDensity

        # get Ca_LVAst channel densities
        for cd in (
            cell.biophysical_properties.membrane_properties.channel_density_nernsts
            + cell.biophysical_properties.membrane_properties.channel_density_non_uniform_nernsts
            +
            cell.biophysical_properties.membrane_properties.channel_density_non_uniforms
        ):
            print(f"Looking at {cd.id}")
            if cd.ion_channel == "Ca_LVAst":
                if "somatic" in cd.id:
                    Ca_LVAst_somatic = cd
                elif "apic" in cd.id:
                    Ca_LVAst_apical = cd

        assert (Ca_LVAst_somatic is not None)
        assert (Ca_LVAst_apical is not None)

        cond_den, unit = pynml.split_nml2_quantity(Ca_LVAst_somatic.cond_density)
        new_cond_den = float(g_Ca_LVAst_multiplier) * float(cond_den)
        print(f"Ca_LVAst somatic: {Ca_LVAst_somatic.cond_density} -> {new_cond_den}")
        Ca_LVAst_somatic.cond_density = f"{new_cond_den} {unit}"

        # Ca_LVAst_apical is a non homogenous param, so needs a little more tinkering
        # the expression is of the form: mul1 * gbar * mul2 * <exp etc dependent on distance>
        value_exp = Ca_LVAst_apical.variable_parameters[
            0
        ].inhomogeneous_value.value  # type: str
        # here we can directly include the multiplier in the expression
        new_value_exp = f"{g_Ca_LVAst_multiplier} * {value_exp}"
        print(f"Ca_LVAst apical: {Ca_LVAst_apical.variable_parameters[0].inhomogeneous_value.value} -> {new_value_exp}")
        Ca_LVAst_apical.variable_parameters[0].inhomogeneous_value.value = new_value_exp


def create_model(
    cellname: str, celldir: str, current_nA: str = "0.5nA", g_Ih_multiplier:
    typing.Optional[str] = None, g_Ca_LVAst_multiplier: typing.Optional[str] = None, cwd=None
) -> typing.Tuple[str, str]:
    """Create model with different values of conductance for the Ih channel.

    :param cellname: TODO
    :param g_Ih_multiplier: multiplier for Ih channels
    :param g_Ca_LVAst_multiplier: multiplier for Ca_LVAst
    :returns: network model file

    """
    if cwd is not None:
        os.mkdir(cwd)
        os.chdir(cwd)
    timestamp = get_timestamp()

    nml_doc_name = f"net_{timestamp}_{cellname}"
    cell_doc_name = f"{cellname}"
    if g_Ih_multiplier is not None:
        nml_doc_name += f"_{g_Ih_multiplier}"
        cell_doc_name += f"_{g_Ih_multiplier}"
    if g_Ca_LVAst_multiplier is not None:
        nml_doc_name += f"_{g_Ca_LVAst_multiplier}"
        cell_doc_name += f"_{g_Ca_LVAst_multiplier}"
    nml_doc_name += f"_{current_nA}"
    cell_doc_name += f"_{current_nA}"

    nml_doc_name = nml_doc_name.replace(".", "_").replace(" ", "_")
    cell_doc_name = cell_doc_name.replace(".", "_").replace(" ", "_") + ".cell.nml"

    nml_doc = component_factory(neuroml.NeuroMLDocument, id=nml_doc_name)
    cell_doc = read_neuroml2_file(f"{celldir}/{cellname}.cell.nml")

    # modify and store cell file
    cell = cell_doc.cells[0]  # type: neuroml.Cell
    create_modified_cell(cell, g_Ih_multiplier, g_Ca_LVAst_multiplier)
    # update channel file paths
    for inc in cell_doc.includes:
        inc.href = f"{get_relative_dir(celldir)}/{inc.href}"
    write_neuroml2_file(cell_doc, cell_doc_name)
    print("Written cell file to: " + cell_doc_name)

    # include cell file
    nml_doc.add("IncludeType", href=cell_doc_name)

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
    return nml_doc_name, cell_doc_name


def simulate_model(model_file_name: str, cellname: str, plot: bool = True,
                   skip_run=True, cwd=None):
    """Simulate the model, generating current plots if required.

    :param model_file_name: name of model file
    :type model_file_name: str
    :param cellname: name of cell
    :type cellname: str
    :param plot: toggle plotting
    :type plot: bool
    :returns: None
    """
    if cwd is not None:
        os.chdir(cwd)

    network_id = f"{cellname}_net"
    simulation_id = model_file_name.split(".")[0]
    simulation = LEMSSimulation(
        sim_id=simulation_id, duration=500, dt=0.01, simulation_seed=123
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


def runner(cellname, celldir, num_data_points, step_sim, if_curve, sim_current_na,
           ifcurve_custom_amps, scz=True, num_processes=None):
    """Common runner for experiment 1.

    This sets up and runs a simulation for each configuration, and then also
    runs the i-f curve generator.

    :param cellname: name of the cell
    :param celldir: name of the celldir
    :param num_data_points: number of ScZ data points to use
    :param step_sim: toggle step current simulation
    :param if_curve: toggle if curve simulation
    :param sim_current_na: current to provide for simulation
    :param ifcurve_custom_amps: list of currents for if curve generator
    :param scz: toggle if we're using ScZ data or just g values
    :returns: list of simulations
    """
    simlist = []
    celldir = get_abs_celldir(celldir)
    expdir = get_run_dir(cellname, "figure_01", scz)
    os.makedirs(expdir)
    os.chdir(expdir)

    if scz:
        data = [
            [1.0, 1.0],
            [data_means_cz[0], 1.0],
            [1.0, data_means_cz[1]],
            data_means_cz
        ]
        if len(scz_data) > 0:
            if num_data_points > len(scz_data):
                data.extend(scz_data)
            else:
                data.extend(random.choices(scz_data, k=num_data_points))
        print(f"Processing {len(data)} ScZ configurations")
    else:
        # for g
        # 2 is "normal"
        data = [0, 1.0, 2.0]
        print(f"Processing {len(data)} non-ScZ configurations")

    ctr = 0
    for d in data:
        if scz:
            mul_Ca_LVAst, mul_Ih = d
            model_file_name, cell_doc_name = create_model(
                cellname=cellname,
                celldir=celldir,
                current_nA=sim_current_na,
                g_Ih_multiplier=f"{mul_Ih}",
                g_Ca_LVAst_multiplier=f"{mul_Ca_LVAst}",
                cwd=f"{cellname}_{ctr}"
            )
        else:
            mul_Ih = d
            model_file_name, cell_doc_name = create_model(
                cellname=cellname,
                celldir=celldir,
                current_nA="0.5 nA",
                g_Ih_multiplier=f"{mul_Ih}",
                cwd=f"{cellname}_{ctr}"
            )
        simlist.append((model_file_name, cell_doc_name))
        ctr += 1

    # run simulations
    if step_sim:
        procs = []
        ctr = 0
        with Pool(processes=num_processes) as p:
            for model_file_name, cell_doc_name in simlist:
                proc = p.apply_async(
                    simulate_model,
                    args=(model_file_name, cellname),
                    kwd=dict(cwd=f"{model_file_name}_{ctr}")
                )
                procs.append(proc)
                ctr += 1
            for r in procs:
                r.wait()

    # do not use Pool here because we need to change the directory before we
    # can call the function
    if if_curve:
        procs = []
        ctr = 0
        proc_ctr = 0
        for model_file_name, cell_doc_name in simlist:
            os.chdir(f"{cellname}_{ctr}")
            proc = Process(
                target=generate_current_vs_frequency_curve,
                kwds=dict(nml2_file=cell_doc_name,
                          cell_id=cellname,
                          custom_amps_nA=ifcurve_custom_amps,
                          analysis_duration=2000,
                          analysis_delay=200,
                          temperature="34 degC",
                          simulator="jNeuroML_NEURON",
                          plot_if=True,
                          plot_iv=True,
                          pre_zero_pulse=300,
                          post_zero_pulse=300,
                          plot_voltage_traces=True,
                          save_iv_figure_to=f"{model_file_name}_iv.png",
                          save_iv_data_to=f"{model_file_name}_iv.dat",
                          save_if_figure_to=f"{model_file_name}_if.png",
                          save_if_data_to=f"{model_file_name}_if.dat",
                          save_voltage_traces_to=f"{model_file_name}_v.png",
                          show_plot_already=False,
                          )
            )
            proc.start()
            procs.append(proc)
            ctr += 1
            proc_ctr += 1
            os.chdir("..")
        # limit to num_processes, wait for these to finish
        if proc_ctr >= num_processes:
            for r in procs:
                r.join()
            proc_ctr = 0

    return simlist
