#!/usr/bin/env python3
"""
Post process and add biophysics to cells.

We make any updates to the morphology, and add biophysics.

File: postprocess_cell.py

Copyright 2022 Ankur Sinha
Author: Ankur Sinha <sanjay DOT ankur AT gmail DOT com>
"""


import numpy
import neuroml
from neuroml.loaders import read_neuroml2_file
from neuroml.writers import NeuroMLWriter
from pyneuroml.analysis import generate_current_vs_frequency_curve
from pyneuroml.pynml import write_neuroml2_file
from neuroml.neuro_lex_ids import neuro_lex_ids


def load_and_setup_cell(cellname: str):
    """Load a cell, and clean it to prepare it for further modifications.

    These operations are common for all cells.

    :param cellname: name of cell.
        the file containing the cell should then be <cell>.morph.cell.nml
    :returns: document with cell
    :rtype: neuroml.NeuroMLDocument

    """
    celldoc = read_neuroml2_file(f"{cellname}.morph.cell.nml")  # type: neuroml.NeuroMLDocument
    cell = celldoc.cells[0]  # type: neuroml.Cell
    celldoc.networks = []
    cell.id = cellname
    cell.notes = cell.notes.replace("L5PCtemplate_0_0", cellname)
    cell.notes += ". Reference: Hay E, Hill S, Schürmann F, Markram H, Segev I (2011) Models of Neocortical Layer 5b Pyramidal Cells Capturing a Wide Range of Dendritic and Perisomatic Active Properties. PLOS Computational Biology 7(7): e1002107. https://doi.org/10.1371/journal.pcbi.1002107"

    # create default groups if they don't exist
    [default_all_group, default_soma_group, default_dendrite_group, default_axon_group] = cell.setup_default_segment_groups(
        use_convention=True, default_groups=["all", "soma_group", "dendrite_group", "axon_group"]
    )

    # populate default groups
    for sg in cell.morphology.segment_groups:
        if "soma" in sg.id and sg.id != "soma_group":
            default_soma_group.add(neuroml.Include(segment_groups=sg.id))
        if "axon" in sg.id and sg.id != "axon_group":
            default_axon_group.add(neuroml.Include(segment_groups=sg.id))
        if "dend" in sg.id and sg.id != "dendrite_group":
            default_dendrite_group.add(neuroml.Include(segment_groups=sg.id))

    cell.optimise_segment_groups()

    return celldoc


def postprocess_L5PC():
    """Post process L5PC and add biophysics.

    """
    cellname = "L5PC"
    celldoc = load_and_setup_cell(cellname)
    cell = celldoc.cells[0]  # type: neuroml.Cell

    # apical dendrites are in groups called apic_
    # basal dendrites are in groups called dend_
    # populate the complete dendrite group, and new groups for all apical and
    # basal dendrites
    default_dendrite_group = cell.get_segment_group("dendrite_group")
    basal_group = cell.add_segment_group("basal_dendrite_group", neuro_lex_id=neuro_lex_ids["dend"], notes="Basal dendrites")
    apical_group = cell.add_segment_group("apical_dendrite_group", neuro_lex_id=neuro_lex_ids["dend"], notes="Apical dendrite_group")

    # remove extra groups
    segment_groups = [x for x in cell.morphology.segment_groups if "ModelViewParmSubset" not in x.id]
    cell.morphology.segment_groups = segment_groups

    for sg in cell.morphology.segment_groups:
        if "apic_" in sg.id:
            apical_group.add(neuroml.Include(segment_groups=sg.id))
        if "dend_" in sg.id:
            basal_group.add(neuroml.Include(segment_groups=sg.id))

    # optimise dendrite group
    default_dendrite_group.includes = []
    default_dendrite_group.includes.append(neuroml.Include(segment_groups=apical_group.id))
    default_dendrite_group.includes.append(neuroml.Include(segment_groups=basal_group.id))

    cell.optimise_segment_groups()
    cell.reorder_segment_groups()

    # include ca dynamics file
    celldoc.add("IncludeType", href="channels/CaDynamics_E2_NML2.nml", validate=False)
    celldoc.add("IncludeType", href="channels/CaDynamics_E2_NML2__decay460__gamma5_01Emin4.nml", validate=False)
    celldoc.add("IncludeType", href="channels/CaDynamics_E2_NML2__decay122__gamma5_09Emin4.nml", validate=False)

    # biophysics
    # all
    cell.set_resistivity("0.1 kohm_cm", group_id="all")
    cell.set_specific_capacitance("1 uF_per_cm2", group_id="all")
    cell.set_init_memb_potential("-80mV")

    # somatic
    soma_group = cell.get_segment_group("soma_group")
    sgid = soma_group.id
    print(f"Adding channels to {sgid}")
    # passive
    cell.add_channel_density(nml_cell_doc=celldoc,
                             cd_id="pas",
                             ion_channel="pas",
                             cond_density="0.0000338 S_per_cm2",
                             erev="-90 mV",
                             group_id=sgid,
                             ion="non_specific",
                             ion_chan_def_file="channels/pas.channel.nml")

    # K
    cell.add_channel_density(nml_cell_doc=celldoc,
                             cd_id="SK_E2_somatic",
                             ion_channel="SK_E2",
                             cond_density="0.441 S_per_cm2",
                             erev="-85 mV",
                             group_id=sgid,
                             ion="k",
                             ion_chan_def_file="channels/SK_E2.channel.nml")
    cell.add_channel_density(nml_cell_doc=celldoc,
                             cd_id="SKv3_1_somatic",
                             ion_channel="SKv3_1",
                             cond_density="0.693 S_per_cm2",
                             erev="-85 mV",
                             group_id=sgid,
                             ion="k",
                             ion_chan_def_file="channels/SKv3_1.channel.nml")
    cell.add_channel_density(nml_cell_doc=celldoc,
                             cd_id="K_Tst_somatic",
                             ion_channel="K_Tst",
                             cond_density="0.0812 S_per_cm2",
                             erev="-85 mV",
                             group_id=sgid,
                             ion="k",
                             ion_chan_def_file="channels/K_Tst.channel.nml")
    cell.add_channel_density(nml_cell_doc=celldoc,
                             cd_id="K_Pst_somatic",
                             ion_channel="K_Pst",
                             cond_density="0.00223 S_per_cm2",
                             erev="-85 mV",
                             group_id=sgid,
                             ion="k",
                             ion_chan_def_file="channels/K_Pst.channel.nml")
    cell.add_channel_density(nml_cell_doc=celldoc,
                             cd_id="Ih_somatic",
                             ion_channel="Ih",
                             cond_density="0.0002 S_per_cm2",
                             erev="-45 mV",
                             group_id=sgid,
                             ion="hcn",
                             ion_chan_def_file="channels/Ih.channel.nml")

    # Na
    cell.add_channel_density(nml_cell_doc=celldoc,
                             cd_id="NaTa_t_somatic",
                             ion_channel="NaTa_t",
                             cond_density="2.04 S_per_cm2",
                             erev="50 mV",
                             group_id=sgid,
                             ion="na",
                             ion_chan_def_file="channels/NaTa_t.channel.nml")
    cell.add_channel_density(nml_cell_doc=celldoc,
                             cd_id="Nap_Et2_somatic",
                             ion_channel="Nap_Et2",
                             cond_density="0.00172 S_per_cm2",
                             erev="50 mV",
                             group_id=sgid,
                             ion="na",
                             ion_chan_def_file="channels/Nap_Et2.channel.nml")
    # Ca
    # internal and external concentrations are set to defaults that NEURON
    # starts with
    cell.add_intracellular_property("Species", validate=False,
                                    id="ca",
                                    concentration_model="CaDynamics_E2_NML2__decay460__gamma5_01Emin4",
                                    ion="ca",
                                    initial_concentration="5.0E-11 mol_per_cm3",
                                    initial_ext_concentration="2.0E-6 mol_per_cm3",
                                    segment_groups=sgid)
    # https://www.neuron.yale.edu/neuron/static/new_doc/modelspec/programmatic/ions.html
    cell.add_channel_density_v(
        "ChannelDensityNernst",
        nml_cell_doc=celldoc,
        id="Ca_HVA_somatic",
        ion_channel="Ca_HVA",
        cond_density="0.000992 S_per_cm2",
        segment_groups=sgid,
        ion="ca",
        ion_chan_def_file="channels/Ca_HVA.channel.nml")
    cell.add_channel_density_v(
        "ChannelDensityNernst",
        nml_cell_doc=celldoc,
        id="Ca_LVAst_somatic",
        ion_channel="Ca_LVAst",
        cond_density="0.00343 S_per_cm2",
        segment_groups=sgid,
        ion="ca",
        ion_chan_def_file="channels/Ca_LVAst.channel.nml")

    # Apical
    sg = cell.get_segment_group("apical_dendrite_group")
    sgid = sg.id
    cell.set_specific_capacitance("2 uF_per_cm2",
                                  group_id=sgid)
    # passive
    cell.add_channel_density(nml_cell_doc=celldoc,
                             cd_id="pas",
                             ion_channel="pas",
                             cond_density="0.0000589 S_per_cm2",
                             erev="-90 mV",
                             group_id=sgid,
                             ion="non_specific",
                             ion_chan_def_file="channels/pas.channel.nml")

    # K
    cell.add_channel_density(nml_cell_doc=celldoc,
                             cd_id="SK_E2_apical",
                             ion_channel="SK_E2",
                             cond_density="0.0012 S_per_cm2",
                             erev="-85 mV",
                             group_id=sgid,
                             ion="k",
                             ion_chan_def_file="channels/SK_E2.channel.nml")
    cell.add_channel_density(nml_cell_doc=celldoc,
                             cd_id="SKv3_1_apical",
                             ion_channel="SKv3_1",
                             cond_density="0.000261 S_per_cm2",
                             erev="-85 mV",
                             group_id=sgid,
                             ion="k",
                             ion_chan_def_file="channels/SKv3_1.channel.nml")
    cell.add_channel_density(nml_cell_doc=celldoc,
                             cd_id="Im_apical",
                             ion_channel="Im",
                             cond_density="0.0000675 S_per_cm2",
                             erev="-85 mV",
                             group_id=sgid,
                             ion="k",
                             ion_chan_def_file="channels/Im.channel.nml")
    # Na
    cell.add_channel_density(nml_cell_doc=celldoc,
                             cd_id="NaTa_t_axonal",
                             ion_channel="NaTa_t",
                             cond_density="0.0213 S_per_cm2",
                             erev="50 mV",
                             group_id=sgid,
                             ion="na",
                             ion_chan_def_file="channels/NaTa_t.channel.nml")
    # Ca
    # internal and external concentrations are set to defaults that NEURON
    # starts with
    cell.add_intracellular_property("Species", validate=False,
                                    id="ca",
                                    concentration_model="CaDynamics_E2_NML2__decay122__gamma5_09Emin4",
                                    ion="ca",
                                    initial_concentration="5.0E-11 mol_per_cm3",
                                    initial_ext_concentration="2.0E-6 mol_per_cm3",
                                    segment_groups=sgid)
    # TODO: distribute CaHVA, CaLVA, Ih
    # Add parameter that we use to distribute Ih
    sg.add(
        "InhomogeneousParameter",
        id="PathLengthOverApicDends",
        variable="p",
        metric="Path Length from root",
        proximal=sg.component_factory(
            "ProximalDetails",
            translation_start="0")
    )
    # distribute Ih
    cdnonuniform_Ih = cell.add_channel_density_v(
        "ChannelDensityNonUniform",
        nml_cell_doc=celldoc,
        id="Ih_apical",
        ion_channel="Ih",
        ion="hcn",
        erev="-45 mV",
        validate=False
    )
    varparam_Ih = cdnonuniform_Ih.add(
        "VariableParameter",
        parameter="condDensity",
        segment_groups=sg.id,
        validate=False
    )
    # TODO: clarify unit conversions
    # 1300.5335 is the value of getLongestBranch("apic")
    # run test_L5PC_loader.hoc, and then run:
    # > access newcell.apic
    # > newcell.soma distance()
    # > newcell.getLongestBranch("apic")
    varparam_Ih.add(
        "InhomogeneousValue",
        inhomogeneous_parameters="PathLengthOverApicDends",
        value="1E9 * (0.2 * 1E-8) * ((2.087 * exp( 3.6161 * (p/1300.5335))) - 0.8696)"
    )

    cdnonuniform_ca_lva = cell.add_channel_density_v(
        "ChannelDensityNonUniformNernst",
        nml_cell_doc=celldoc,
        id="Ca_LVA_apic",
        ion_channel="Ca_LVAst",
        ion="ca",
        validate=False
    )
    varparam_ca_lva = cdnonuniform_ca_lva.add(
        "VariableParameter",
        parameter="condDensity",
        segment_groups=sg.id,
        validate=False
    )
    # TODO: clarify unit conversions
    # 1300.5335 is the value of getLongestBranch("apic")
    # run test_L5PC_loader.hoc, and then run:
    # > access newcell.apic
    # > newcell.soma distance()
    # > newcell.getLongestBranch("apic")
    # The conditional is implemented using a heaviside function:
    # https://github.com/NeuroML/org.neuroml.export/blob/master/src/main/java/org/neuroml/export/neuron/NRNUtils.java#L174
    # if both values in H are true, 0.09 is added to 0.01 = 1
    # otherwise, 0.01
    varparam_ca_lva.add(
        "InhomogeneousValue",
        inhomogeneous_parameters="PathLengthOverApicDends",
        value="1.0E9 * 1.0E-8 * 0.0187 * (0.01 + (0.09 * (H(p - 685) * H(885 - p))))"
    )

    cdnonuniform_ca_hva = cell.add_channel_density_v(
        "ChannelDensityNonUniformNernst",
        nml_cell_doc=celldoc,
        id="Ca_HVA_apic",
        ion_channel="Ca_HVA",
        ion="ca",
        validate=False
    )
    varparam_ca_hva = cdnonuniform_ca_hva.add(
        "VariableParameter",
        parameter="condDensity",
        segment_groups=sg.id,
        validate=False
    )
    # TODO: clarify unit conversions
    # The conditional is implemented using a heaviside function:
    # https://github.com/NeuroML/org.neuroml.export/blob/master/src/main/java/org/neuroml/export/neuron/NRNUtils.java#L174
    # if both values in H are true, 0.9 is added to 0.1 = 1
    # otherwise, 0.1
    varparam_ca_hva.add(
        "InhomogeneousValue",
        inhomogeneous_parameters="PathLengthOverApicDends",
        value="1.0E9 * 1.0E-8 * 0.000555 * (0.1 + (0.9 * (H(p - 685) * H(885 - p))))"
    )

    # basal
    cell.set_specific_capacitance("2 uF_per_cm2",
                                  group_id="basal_dendrite_group")
    # passive
    cell.add_channel_density(nml_cell_doc=celldoc,
                             cd_id="pas",
                             ion_channel="pas",
                             cond_density="0.0000467 S_per_cm2",
                             erev="-90 mV",
                             group_id=sgid,
                             ion="non_specific",
                             ion_chan_def_file="channels/pas.channel.nml")

    cell.add_channel_density(nml_cell_doc=celldoc,
                             cd_id="Ih_basal",
                             ion_channel="Ih",
                             cond_density="0.0002 S_per_cm2",
                             erev="-45 mV",
                             group_id=sgid,
                             ion="hcn",
                             ion_chan_def_file="channels/Ih.channel.nml")

    # axonal
    axon_group = cell.get_segment_group("axon_group")
    sgid = axon_group.id
    print(f"Adding channels to {sgid}")
    # passive
    cell.add_channel_density(nml_cell_doc=celldoc,
                             cd_id="pas",
                             ion_channel="pas",
                             cond_density="0.0000325 S_per_cm2",
                             erev="-90 mV",
                             group_id=sgid,
                             ion="non_specific",
                             ion_chan_def_file="channels/pas.channel.nml")

    # write to file
    write_neuroml2_file(celldoc, f"{cellname}.cell.nml")


def analyse_L5PC(hyperpolarising: bool = True, depolarising: bool = True):
    """Generate various curves for L5PC cells

    :returns: None

    """
    cellname = "L5PC"
    if hyperpolarising:
        # hyper-polarising inputs
        generate_current_vs_frequency_curve(
            nml2_file=f"{cellname}.cell.nml",
            cell_id=cellname,
            custom_amps_nA=list(numpy.arange(-0.05, -0.1, -0.01)),
            temperature="34 degC",
            pre_zero_pulse=200,
            post_zero_pulse=300,
            plot_voltage_traces=True,
            plot_iv=True,
            plot_if=False,
            simulator="jNeuroML_NEURON",
            analysis_delay=300.,
            analysis_duration=400.
        )

    if depolarising:
        # depolarising inputs
        generate_current_vs_frequency_curve(
            nml2_file=f"{cellname}.cell.nml",
            cell_id=cellname,
            plot_voltage_traces=True,
            spike_threshold_mV=-10.0,
            custom_amps_nA=list(numpy.arange(0, 0.3, 0.05)),
            temperature="34 degC",
            pre_zero_pulse=200,
            post_zero_pulse=300,
            plot_iv=True,
            simulator="jNeuroML_NEURON",
            analysis_delay=50.,
            analysis_duration=200.
        )


if __name__ == "__main__":
    # postprocess_L5PC()
    analyse_L5PC(False, True)
