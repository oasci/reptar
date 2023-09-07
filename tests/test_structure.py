# MIT License
#
# Copyright (c) 2022-2023, Alex M. Maldonado
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

"""Tests structure analysis"""

# pylint: skip-file

import sys
import os
import shutil
import pytest
import numpy as np
from reptar import File
from reptar.structure.bond import is_bonded_atomic_radii
from reptar.structure.topology import get_topology_change

sys.path.append("..")

# Ensures we execute from file directory (for relative paths).
os.chdir(os.path.dirname(os.path.realpath(__file__)))

# Source paths
DATA_DIR = "./data/"


def test_topology_change():
    rfile_source = File("./data/gfp-cro.zarr", mode="r")
    Z = rfile_source.get("atomic_numbers")
    R = rfile_source.get("geometry")

    is_bonded_kwargs = {"use_ray": False, "a0_scaling": 1.2}

    topo_changed = get_topology_change(
        Z,
        R[0],
        R[:2],
        is_bonded_atomic_radii,
        is_bonded_kwargs,
        use_ray=False,
    )
    assert np.allclose(np.array([False, False]), topo_changed)
