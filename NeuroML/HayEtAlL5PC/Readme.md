This directory contains the L5PC cell from Hay et al. converted to NeuroML format.

- cellmorph2nml.py: loads in the cell to Neuron and exports the morphology to NeuroML in L5PC.morph.cell.nml.
- postprocess_cell.py: loads L5PC.morph.cell.nml and adds necessary segment groups; adds biophysics to create L5PC.cell.nml. Also includes some simple analysis code to generate if-curves for the converted cell.
