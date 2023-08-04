#!/usr/bin/env python3
"""
Common methods used by experiments

File: common/__init__.py

Copyright 2023 Ankur Sinha
Author: Ankur Sinha <sanjay DOT ankur AT gmail DOT com>
"""


import os
import datetime
import shutil


def get_abs_celldir(cellfolder_name):
    """Get the absolute cell dir

    We know where the cells are relative to this module, so we can calculate
    the absolute path

    :param cellfolder_name: name of the cell folder
    :type cellfolder_name: str
    :returns: absolute path to cell directory

    """
    celldir_rel = f"../../cells/{cellfolder_name}"
    celldir_abs = os.path.abspath(celldir_rel)
    print(f"Celldir is {celldir_abs}")
    return celldir_abs


def get_timestamp():
    """Get current time stamp"""
    return datetime.datetime.now().strftime("%Y%m%d%H%M%S")


def delete_neuron_special_dir():
    """Delete neuron special directory containing compiled mods
    """
    try:
        shutil.rmtree("x86_64")
    except FileNotFoundError:
        pass


def get_run_dir(cellname, experiment_name, root=None, timestamp=None):
    """Get a new directory to run the simulations:

    {root}/{experiment_name}/{timestamp}_{cellname}/

    :param cellname: name of cell
    :type cellname: str
    :param root: top level directory relative to which run dir is created
    :type root: str
    :param timestamp: timestamp string
    :type timestamp: str

    :returns: string of path
    """
    if timestamp is None:
        timestamp = get_timestamp()

    if root is None:
        root = "../../../../simdata/"

    root_abs = os.path.abspath(root)
    rootdir = f"{root_abs}/{experiment_name}/{timestamp}_{cellname}"
    print(f"Run dir is {rootdir}")
    return rootdir
