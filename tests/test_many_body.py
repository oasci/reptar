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

"""Tests many-body utilities"""

import pytest
import os
import shutil
import numpy as np
from reptar import File
from reptar.mb import mb_contributions

import sys
sys.path.append("..")
from .paths import *

# Ensures we execute from file directory (for relative paths).
os.chdir(os.path.dirname(os.path.realpath(__file__)))

def test_2h2o_mb_from_data():
    """Computing 2-body energy and gradient contributions.
    """
    exdir_path = get_140h2o_samples_path()

    exdir_file = File(exdir_path, mode='r')

    # Computing reptar values
    E = np.array(exdir_file.get('2h2o/energy_ele_rimp2.ccpvtz'))[:10]  # Eh
    G = np.array(exdir_file.get('2h2o/grads_rimp2.ccpvtz'))[:10]  # Eh/A
    entity_ids = np.array(exdir_file.get('2h2o/entity_ids'))
    r_prov_ids = exdir_file.get('2h2o/r_prov_ids')
    r_prov_specs = np.array(exdir_file.get('2h2o/r_prov_specs'))[:10]

    E_lower = np.array(exdir_file.get('1h2o/energy_ele_rimp2.ccpvtz')) # Eh
    G_lower = np.array(exdir_file.get('1h2o/grads_rimp2.ccpvtz'))  # Eh/A
    entity_ids_lower = np.array(exdir_file.get('1h2o/entity_ids'))
    r_prov_ids_lower = exdir_file.get('1h2o/r_prov_ids')
    r_prov_specs_lower = np.array(exdir_file.get('1h2o/r_prov_specs'))

    E_mb, G_mb = mb_contributions(
        E, G, entity_ids, r_prov_ids, r_prov_specs, E_lower, G_lower,
        entity_ids_lower, r_prov_ids_lower, r_prov_specs_lower,
        n_workers=1
    )

    assert E_mb[5] == -0.0025468531256365168
    assert G_mb[5][1][2] == -0.02043695304313046
