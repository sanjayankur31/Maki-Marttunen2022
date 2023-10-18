# L5 ScZ modelling

Includes simulation code for:

Excitatory and inhibitory effects of HCN channel modulation on excitability of layer V pyramidal cells, PLOS Computational Biology, 2022. https://doi.org/10.1371/journal.pcbi.1010506

taken from the Model-DB submission:

http://modeldb.yale.edu/267293.

## Usage

This repository has two components:

- a public code component, hosted on GitHub: this includes all the source code
- a private data component, hosted in a private repository on G-node: this includes all the simulation generated data and analysis

You can clone the GitHub repository normally using Git.
The private repository on G-node is currently for collaborators only and will be made public on completion of the project.

## Using the private G-node repository

The G-node repository supports git-annex, which is used to store the simulation generated data.
The advantage of git-annex is that one can store large data files in the git-repository, but not necessarily keep a local copy of these on disk.
Required large files can be fetched on demand.
More information on git-annex can be found in its [documentation](https://git-annex.branchable.com/).

To use the G-node repository, clone it normally:

    git clone <url>

To also fetch the git-annex data, use:

    git annex sync --content

One can also just fetch individual files

    git annex get <path to file to fetch>


One can also use [datalad](https://handbook.datalad.org/en/latest/index.html) to work with this repository (which used git/git-annex under the hood).
