# MIT License
#
# Copyright (c) 2022, Alex M. Maldonado
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

"""Tests writing schnetpack databases"""

# Stuff for PyTest features like skip.
# pylint: disable=import-outside-toplevel, unused-import, duplicate-code

import os
import sys
import pytest
import numpy as np
from reptar import File
from reptar.writers import write_schnetpack_db

sys.path.append("..")
# pylint: disable-next=wrong-import-position
from .utils import _pdist

HARTREE_TO_EV = 27.21138602  # Psi4 v1.5

# Ensures we execute from file directory (for relative paths).
os.chdir(os.path.dirname(os.path.realpath(__file__)))

# Source paths
DATA_DIR = "./data/"

# Writing paths
WRITING_DIR = "./tmp/writing/"
os.makedirs(WRITING_DIR, exist_ok=True)


def test_schnetpack_db_writer_1h2o_120meoh_prod():
    r"""Writing small schnetpack database"""
    try:
        import schnetpack
    except ImportError:
        pytest.skip("schnetpack package not installed")

    exdir_path = os.path.join(DATA_DIR, "1h2o_120meoh_md.exdir")
    db_path = os.path.join(WRITING_DIR, "1h2o_120meoh_md_schnetpack.db")

    if os.path.exists(db_path):
        os.remove(db_path)

    i_test = 3

    rfile = File(exdir_path, mode="r")

    Z = rfile.get("eq_1/atomic_numbers")
    R = rfile.get("eq_1/geometry")[:10]
    E = rfile.get("eq_1/energy_pot")[:10]  # Hartree
    E *= HARTREE_TO_EV  # eV

    db = write_schnetpack_db(db_path, Z, R, energy=E, centering_function=None)
    row_atoms, row_props = db.get_properties(i_test)

    assert np.array_equal(Z, row_props["_atomic_numbers"])
    assert E[i_test] == row_props["energy"][0]
    assert np.array_equal(R[i_test], row_atoms.positions)  # Exact

    # Note that row_props['_positions'] is slightly different.
    # I think they are still centering the structure, so we check
    # pairwise distances.
    assert np.allclose(_pdist(R[i_test]), _pdist(row_props["_positions"]))
