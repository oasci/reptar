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

"""Tests and example builders for molecular dynamics with xtb."""

import pytest
import os
import numpy as np
from reptar import manager
from reptar.utils import get_md5, gen_entity_ids, gen_comp_ids

import sys
sys.path.append("..")
from .paths import get_1h2o_120meoh_eq_paths, get_1h2o_120meoh_prod_paths

# Ensures we execute from file directory (for relative paths).
os.chdir(os.path.dirname(os.path.realpath(__file__)))

@pytest.mark.order(0)
def test_1h2o_120meoh_md():
    """
    """
    dir_path, _, out_path_eq, geom_path_eq, traj_path_eq = get_1h2o_120meoh_eq_paths()
    _, _, out_path_prod, geom_path_prod, traj_path_prod = get_1h2o_120meoh_prod_paths()
    exdir_path = dir_path

    num_waters = 1
    atoms_per_water = 3
    num_methanols = 120
    atoms_per_methanol = 6

    # For water molecules.
    entity_ids = gen_entity_ids(atoms_per_water, num_waters)
    comp_ids = gen_comp_ids('WAT', num_waters, entity_ids)
    # For methanol molecules.
    entity_ids = gen_entity_ids(
        atoms_per_methanol, num_methanols, starting_idx=np.max(entity_ids)+1,
        add_to=entity_ids
    )
    comp_ids = gen_comp_ids(
        'MET', num_methanols, entity_ids, add_to=comp_ids
    )

    repman = manager()
    repman.get_file(exdir_path, mode='w', allow_remove=True)
    repman.add_definitions(definitions=['md', 'xtb', 'pes'])
    repman.create_group(
        'eq_1', parent=repman.File, out_path=out_path_eq,
        geom_path=geom_path_eq, traj_path=traj_path_eq
    )
    repman.add_ids(repman.File['eq_1'], entity_ids, comp_ids)
    repman.create_group(
        'prod_1', parent=repman.File, out_path=out_path_prod,
        geom_path=geom_path_prod, traj_path=traj_path_prod
    )
    repman.add_ids(repman.File['prod_1'], entity_ids, comp_ids)
    repman.close()
    
    repman = manager()
    repman.get_file(
        exdir_path, mode='r'
    )

    assert repman.File['prod_1']['geometry'].shape == (1001, 723, 3)
    assert repman.File['prod_1']['geometry'][0][3][1] == -2.64576119977354
    assert repman.File['prod_1']['geometry'][-1][-1][0] == -1.62420092727186
    assert repman.File['prod_1']['energy_pot'].shape == (1001,)
    assert repman.File['prod_1']['energy_pot'][0] == -991.881189902216
    assert repman.File['prod_1']['energy_pot'][-1] == -991.818996146108


