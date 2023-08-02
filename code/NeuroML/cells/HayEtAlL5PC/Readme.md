This directory contains the L5PC cell from Hay et al. converted to NeuroML format.

- 00_swc_to_nml_morphology.py : loads in the cell to Neuron and exports the morphology to NeuroML in L5PC.morph.cell.nml.
- 01_post_process_nml.py: loads L5PC.morph.cell.nml and adds necessary segment groups; adds biophysics to create L5PC.cell.nml. Also includes some simple analysis code to generate if-curves for the converted cell.

Two additional scripts are included to compare the NEURON cell generated from the NeuroML conversion and the original NEURON cell.

- 02_get_info_neuron_from_nml_cell.py: after generating the NEURON cell from the NeuroML sources (by running the analysis functions in the post process script, for example), this is used to get information from NEURON about its representation of the cell.
- NEURON-tests/00_get_info_neuron_cell.py: same as above, but uses the test hoc to load the original NEURON cell.
