# MIT License
#
# Copyright (c) 2022-2023, OASCI
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

"""Tests writing ForceBalance files"""

# pylint: skip-file

import os
import sys
import pytest
import numpy as np
from reptar import File
from reptar.writers import write_qdata

sys.path.append("..")

# Ensures we execute from file directory (for relative paths).
os.chdir(os.path.dirname(os.path.realpath(__file__)))

# Source paths
DATA_DIR = "./data/"

# Writing paths
WRITING_DIR = "./tmp/writing/"
os.makedirs(WRITING_DIR, exist_ok=True)


def test_qdata_writer_1h2o():
    r"""Writing small schnetpack database"""

    exdir_path = os.path.join(DATA_DIR, "h2o-temelso.etal.exdir")
    write_path = os.path.join(WRITING_DIR, "qdata.txt")

    if os.path.exists(write_path):
        os.remove(write_path)

    rfile = File(exdir_path, mode="r")

    source_key = "6h2o/samples_1h2o/"

    Z = rfile.get(os.path.join(source_key, "atomic_numbers"))
    R = rfile.get(os.path.join(source_key, "geometry"))
    E = rfile.get(os.path.join(source_key, "energy_ele_mp2.def2tzvp_orca"))  # Eh
    G = rfile.get(os.path.join(source_key, "grads_mp2.def2tzvp_orca"))  # Eh/A
    F = np.negative(G[:])

    write_qdata(write_path, R, energy=E, forces=F)
