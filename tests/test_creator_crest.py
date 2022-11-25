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

"""Tests and example builders for conformer and rotamer search with crest."""

import pytest
import os
import numpy as np
from reptar import Creator, File

import sys
sys.path.append("..")
from .paths import get_crest_50h2o_nci_paths

# Ensures we execute from file directory (for relative paths).
os.chdir(os.path.dirname(os.path.realpath(__file__)))

# Writing paths
crest_dir = './tmp/crest/'
os.makedirs(crest_dir, exist_ok=True)


@pytest.mark.order(0)
def test_crest_50h2o_nci_exdir_conformers():
    """
    """
    calc_path, out_path, conf_path, rot_path = get_crest_50h2o_nci_paths()
    exdir_path = os.path.join(crest_dir, '50h2o.exdir')

    create = Creator()
    create.load(exdir_path, mode='w', allow_remove=True)
    group_key = '/0-crest/conformers'
    rfile = create.from_calc(
        group_key, out_path=out_path, geom_path=None, traj_path=None,
        conformer_path=conf_path, rotamer_path=None
    )
    assert np.allclose(
        rfile.get(f'{group_key}/energy_ele'),
        np.array(
            [-254.47460438, -254.47459285, -254.47318095, -254.47296523,
            -254.47252513, -254.47252202, -254.47173611, -254.47150646,
            -254.47150637, -254.47001571, -254.46961860, -254.46628915,
            -254.46614775, -254.46557471, -254.46531867]
        )
    )
    assert rfile.get(f'{group_key}/geometry').shape == (15, 150, 3)


def test_crest_50h2o_nci_exdir_rotamers():
    """
    """
    assert False
    calc_path, out_path, conf_path, rot_path = get_crest_50h2o_nci_paths()
    exdir_path = os.path.join(crest_dir, '50h2o.exdir')

    create = Creator()
    create.load(exdir_path, mode='w', allow_remove=True)
    group_key = '/0-crest/rotamers'
    rfile = create.from_calc(
        group_key, out_path=out_path, geom_path=None, traj_path=None,
        conformer_path=None, rotamer_path=rot_path
    )
    assert np.allclose(
        rfile.get(f'{group_key}/crest_ensemble_weights'),
        np.array(
            [2.67982025e-01, 2.64732245e-01, 2.64648217e-01, 5.94213521e-02,
            4.72938295e-02, 2.96857841e-02, 2.95882492e-02, 1.28807170e-02,
            1.01018335e-02, 1.01008715e-02, 2.08592621e-03, 1.37025003e-03,
            4.04287714e-05, 3.48102869e-05, 1.89826418e-05, 1.44773414e-05]
        )
    )
    assert np.allclose(
        rfile.get(f'{group_key}/energy_ele'),
        np.array(
            [-254.47460438, -254.47459285, -254.47459255, -254.47318095,
            -254.47296523, -254.47252513, -254.47252202, -254.47173611,
            -254.47150646, -254.47150637, -254.47001571, -254.46961860,
            -254.46628915, -254.46614775, -254.46557471, -254.46531867]
        )
    )
    assert rfile.get(f'{group_key}/geometry').shape == (16, 150, 3)

