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

"""Tests writing atomic simulation environment (ASE) databases"""

import pytest
import os
import numpy as np
from reptar import File
from reptar.writers import write_ase_db

import sys
sys.path.append("..")
from .paths import *

hartree2ev = 27.21138602  # Psi4 v1.5

# Ensures we execute from file directory (for relative paths).
os.chdir(os.path.dirname(os.path.realpath(__file__)))

# Source paths
xtb_dir = './tmp/xtb'

# Writing paths
writing_dir = './tmp/writing/'
os.makedirs(writing_dir, exist_ok=True)

@pytest.mark.dependency(
    depends=['tests/test_creator_xtb.py::test_1h2o_120meoh_md_exdir'],
    scope='session'
)
def test_ase_db_writer_1h2o_120meoh_prod():
    """Writing small ASE database"""
    try:
        import ase
    except (ModuleNotFoundError, ImportError):
        pytest.skip("ase package not installed")
    
    exdir_path = os.path.join(xtb_dir, '1h2o_120meoh_md.exdir')
    db_path = os.path.join(writing_dir, '1h2o_120meoh_md_ase.db')

    if os.path.exists(db_path):
        os.remove(db_path)

    i_test = 3

    rfile = File(exdir_path, mode='r')
    Z = rfile.get('prod_1/atomic_numbers')
    R = rfile.get('prod_1/geometry')[:10]
    E = rfile.get('prod_1/energy_pot')[:10]  # Hartree
    E *= hartree2ev  # eV

    db = write_ase_db(db_path, Z, R, energy=E)
    row = db.get(i_test+1)

    assert np.array_equal(Z, row['numbers'])
    assert np.array_equal(R[i_test], row['positions'])
    assert E[i_test] == row['energy']
    assert np.array_equal(np.array([False, False, False]), row['pbc'])
