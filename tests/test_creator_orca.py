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
from reptar import manager
from reptar.utils import get_md5, gen_entity_ids, gen_comp_ids

import sys
sys.path.append("..")
from .paths import *

# Ensures we execute from file directory (for relative paths).
os.chdir(os.path.dirname(os.path.realpath(__file__)))

@pytest.mark.order(0)
def test_6h2o_temelso_engrad_exdir():
    """
    """
    dir_path, out_path = get_6h2o_temelso_pr_engrad()
    exdir_path = dir_path + '.exdir'

    repman = manager()
    repman.load(exdir_path, mode='w', allow_remove=True)
    repman.add_definitions(definitions=['qc', 'pes', 'molprop'])
    repman.create_group('/engrad', out_path=out_path)

    # Reading tests.
    repman = manager()
    repman.load(exdir_path, mode='r')
    assert repman.data.get('/engrad/energy_scf') == -12419.360138637763
    assert repman.data.get('/engrad/md5') == "0345198c3a033c75c282744ab68969e0"
    assert np.array_equal(
        repman.data.get('/engrad/dipole_moment'),
        np.array([2.5265211700823, -0.5830003327751, 1.2353903376292001])
    )

@pytest.mark.order(0)
def test_6h2o_temelso_engrad_json():
    """
    """
    dir_path, out_path = get_6h2o_temelso_pr_engrad()
    json_path = dir_path + '.json'

    # Writing tests.
    repman = manager()
    repman.load(json_path, mode='w')
    repman.create_group('/', out_path=out_path)
    repman.data.save()

    # Reading tests.
    repman = manager()
    repman.load(json_path, mode='r')
    assert repman.data.get('energy_scf') == -12419.360138637763
    assert repman.data.get('md5') == "0470495abad0aa89590709c7137063b7"
    assert np.array_equal(
        repman.data.get('dipole_moment'),
        np.array([2.5265211700823, -0.5830003327751, 1.2353903376292001])
    )

@pytest.mark.order(0)
def test_6h2o_temelso_engrad_npz():
    """
    """
    dir_path, out_path = get_6h2o_temelso_pr_engrad()
    npz_path = dir_path + '.npz'

    repman = manager()
    repman.load(npz_path, mode='w')
    repman.create_group('/', out_path=out_path)
    repman.data.save()

    # Reading tests.
    repman = manager()
    repman.load(npz_path, mode='r')
    assert repman.data.get('energy_scf') == -12419.360138637763
    assert repman.data.get('md5') == "0470495abad0aa89590709c7137063b7"
    assert np.array_equal(
        repman.data.get('dipole_moment'),
        np.array([2.5265211700823, -0.5830003327751, 1.2353903376292001])
    )
