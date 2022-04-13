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

"""Tests sampling structures from exdir Groups"""

import pytest
import os
import shutil
import numpy as np
from reptar import manager, add_structures_to_group
import itertools

import sys
sys.path.append("..")
from .paths import *

# Ensures we execute from file directory (for relative paths).
os.chdir(os.path.dirname(os.path.realpath(__file__)))

def test_1h2o_120meoh_prod_sampling():
    """Sampling from xTB MD reptar file.
    """
    source_path = get_1h2o_120meoh_exdir_path()
    dest_path = './tmp/test_1h2o_120meoh_prod_sampling.exdir'

    source = manager()
    source.get_file(source_path, mode='r')
    grp_source = source.File['prod_1']

    dest = manager()
    dest.get_file(dest_path, mode='w', allow_remove=True)
    grp_dest = dest.create_group('wat.2met-pes')

    quantity = 100
    comp_labels = ('WAT', 'MET', 'MET')

    add_structures_to_group(
        grp_source, grp_dest, quantity, comp_labels,
        center_structures=False, copy_EF=False, write=True
    )
    assert np.array_equal(
        grp_dest['atomic_numbers'],
        np.array([8, 1, 1, 8, 1, 6, 1, 1, 1, 8, 1, 6, 1, 1, 1])
    )
    assert np.array_equal(
        grp_dest['entity_ids'],
        np.array([0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2])
    )
    assert grp_dest['geometry'].data.shape == (100, 15, 3)
    assert grp_dest['r_prov_specs'].data.shape == (100, 5)
    assert len(dict(grp_dest.attrs['r_prov_ids'])) == 1
    assert grp_dest.attrs['r_centered'] == False
    assert np.array_equal(
        grp_dest['comp_ids'].data, np.array(comp_labels)
    )

    add_structures_to_group(
        grp_source, grp_dest, quantity, comp_labels,
        center_structures=True, copy_EF=False, write=True
    )
    assert np.array_equal(
        grp_dest['atomic_numbers'],
        np.array([8, 1, 1, 8, 1, 6, 1, 1, 1, 8, 1, 6, 1, 1, 1])
    )
    assert np.array_equal(
        grp_dest['entity_ids'],
        np.array([0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2])
    )
    assert grp_dest['geometry'].data.shape == (200, 15, 3)
    assert grp_dest['r_prov_specs'].data.shape == (200, 5)
    assert len(dict(grp_dest.attrs['r_prov_ids'])) == 1
    assert grp_dest.attrs['r_centered'] == True
    assert np.array_equal(
        grp_dest['comp_ids'].data, np.array(comp_labels)
    )

def test_sampling_from_wat_2met_pes():
    """Sampling from a sampled exdir group.
    """
    src_path = './tmp/test_1h2o_120meoh_prod_sampling.exdir'

    source = manager()
    source.get_file(src_path, mode='a')
    grp_source = source.File['wat.2met-pes']
    grp_dest = source.create_group('wat.met-pes')

    quantity = 100
    comp_labels = ('WAT', 'MET')

    add_structures_to_group(
        grp_source, grp_dest, quantity, comp_labels,
        center_structures=True, copy_EF=False
    )

    assert np.array_equal(
        grp_dest['entity_ids'],
        np.array([0, 0, 0, 1, 1, 1, 1, 1, 1])
    )
    assert len(dict(grp_dest.attrs['r_prov_ids'])) == 1
    assert grp_dest.attrs['r_centered'] == True
    assert np.array_equal(
        grp_dest['comp_ids'].data, np.array(comp_labels)
    )
