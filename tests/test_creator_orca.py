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

"""Tests and example builders for ORCA calculations."""

import pytest
import os
import numpy as np
from reptar import creator
from reptar.utils import gen_entity_ids, gen_comp_ids

import sys
sys.path.append("..")
from .paths import *

# Ensures we execute from file directory (for relative paths).
os.chdir(os.path.dirname(os.path.realpath(__file__)))

# Writing paths
orca_dir = './tmp/orca'
os.makedirs(orca_dir, exist_ok=True)

@pytest.mark.order(0)
def test_6h2o_temelso_engrad_exdir():
    """
    """
    dir_path, out_path = get_6h2o_temelso_pr_engrad()
    exdir_path = os.path.join(orca_dir, '6h2o_temelso_pr_engrad.exdir')

    create = creator()
    create.load(exdir_path, mode='w', allow_remove=True)
    create.definitions(definitions=['qc', 'pes', 'molprop'])
    create.from_calc('/engrad', out_path=out_path)

    # Reading tests.
    create = creator()
    create.load(exdir_path, mode='r')
    assert create.rfile.get('/engrad/energy_scf') == -12419.360138637763
    assert np.array_equal(
        create.rfile.get('/engrad/dipole_moment'),
        np.array([2.5265211700823, -0.5830003327751, 1.2353903376292001])
    )

@pytest.mark.order(0)
def test_6h2o_temelso_engrad_json():
    """
    """
    dir_path, out_path = get_6h2o_temelso_pr_engrad()
    json_path = os.path.join(orca_dir, '6h2o_temelso_pr_engrad.json')

    # Writing tests.
    create = creator()
    create.load(json_path, mode='w')
    create.from_calc('/', out_path=out_path)
    create.rfile.save()

    # Reading tests.
    create = creator()
    create.load(json_path, mode='r')
    assert create.rfile.get('energy_scf') == -12419.360138637763
    assert np.array_equal(
        create.rfile.get('dipole_moment'),
        np.array([2.5265211700823, -0.5830003327751, 1.2353903376292001])
    )

@pytest.mark.order(0)
def test_6h2o_temelso_engrad_npz():
    """
    """
    dir_path, out_path = get_6h2o_temelso_pr_engrad()
    npz_path = os.path.join(orca_dir, '6h2o_temelso_pr_engrad.npz')

    create = creator()
    create.load(npz_path, mode='w')
    create.from_calc('/', out_path=out_path)
    create.rfile.save()

    # Reading tests.
    create = creator()
    create.load(npz_path, mode='r')
    assert create.rfile.get('energy_scf') == -12419.360138637763
    assert np.array_equal(
        create.rfile.get('dipole_moment'),
        np.array([2.5265211700823, -0.5830003327751, 1.2353903376292001])
    )
