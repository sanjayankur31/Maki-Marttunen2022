#!/usr/bin/env python3
"""
Common methods used by experiments

File: common/__init__.py

Copyright 2023 Ankur Sinha
Author: Ankur Sinha <sanjay DOT ankur AT gmail DOT com>
"""


import datetime
import shutil


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
