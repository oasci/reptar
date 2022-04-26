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
from reptar import manager, data
from reptar.utils import get_md5, gen_entity_ids, gen_comp_ids

import sys
sys.path.append("..")
from .paths import get_1h2o_120meoh_eq_paths, get_1h2o_120meoh_prod_paths

# Ensures we execute from file directory (for relative paths).
os.chdir(os.path.dirname(os.path.realpath(__file__)))

@pytest.mark.order(0)
def test_1h2o_120meoh_md_exdir():
    """
    """
    dir_path, _, out_path_eq, geom_path_eq, traj_path_eq = get_1h2o_120meoh_eq_paths()
    _, _, out_path_prod, geom_path_prod, traj_path_prod = get_1h2o_120meoh_prod_paths()
    exdir_path = dir_path + '.exdir'

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
    repman.load(exdir_path, mode='w', allow_remove=True)
    repman.add_definitions(definitions=['md', 'xtb', 'pes'])
    repman.create_group(
        '/eq_1', out_path=out_path_eq,
        geom_path=geom_path_eq, traj_path=traj_path_eq
    )
    repman.add_ids('/eq_1', entity_ids, comp_ids)
    repman.create_group(
        'prod_1', out_path=out_path_prod,
        geom_path=geom_path_prod, traj_path=traj_path_prod
    )
    repman.add_ids('prod_1', entity_ids, comp_ids)
    
    repman = manager()
    repman.load(
        exdir_path, mode='r'
    )

    assert repman.data.get('prod_1/geometry').shape == (1001, 723, 3)
    assert repman.data.get('prod_1/geometry')[0][3][1] == -2.64576119977354
    assert repman.data.get('prod_1/geometry')[-1][-1][0] == -1.62420092727186
    assert repman.data.get('prod_1/energy_pot').shape == (1001,)
    assert repman.data.get('prod_1/energy_pot')[0] == -991.881189902216
    assert repman.data.get('prod_1/energy_pot')[-1] == -991.818996146108
    assert repman.data.get('prod_1/wall_potential')[0]['sphere_radius'] == 12.500003

@pytest.mark.order(0)
def test_1h2o_120meoh_md_json():
    """
    """
    dir_path, _, out_path_eq, geom_path_eq, traj_path_eq = get_1h2o_120meoh_eq_paths()
    _, _, out_path_prod, geom_path_prod, traj_path_prod = get_1h2o_120meoh_prod_paths()
    json_path = dir_path + '.json'

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
    repman.load(json_path, mode='w')
    repman.create_group(
        '/eq_1', out_path=out_path_eq,
        geom_path=geom_path_eq, traj_path=traj_path_eq
    )
    repman.add_ids('/eq_1', entity_ids, comp_ids)
    repman.create_group(
        'prod_1', out_path=out_path_prod,
        geom_path=geom_path_prod, traj_path=traj_path_prod
    )
    repman.add_ids('prod_1', entity_ids, comp_ids)

    assert repman.data.get('prod_1/geometry').shape == (1001, 723, 3)
    assert repman.data.get('prod_1/geometry')[0][3][1] == -2.64576119977354
    assert repman.data.get('prod_1/geometry')[-1][-1][0] == -1.62420092727186
    assert repman.data.get('prod_1/energy_pot').shape == (1001,)
    assert repman.data.get('prod_1/energy_pot')[0] == -991.881189902216
    assert repman.data.get('prod_1/energy_pot')[-1] == -991.818996146108
    assert repman.data.get('prod_1/wall_potential')[0]['sphere_radius'] == 12.500003

    repman.data.save(json_prettify=True)

    repman = manager()
    repman.load(json_path, mode='r')

    assert repman.data.get('prod_1/geometry').shape == (1001, 723, 3)
    assert repman.data.get('prod_1/geometry')[0][3][1] == -2.64576119977354
    assert repman.data.get('prod_1/geometry')[-1][-1][0] == -1.62420092727186
    assert repman.data.get('prod_1/energy_pot').shape == (1001,)
    assert repman.data.get('prod_1/energy_pot')[0] == -991.881189902216
    assert repman.data.get('prod_1/energy_pot')[-1] == -991.818996146108
    assert repman.data.get('prod_1/wall_potential')[0]['sphere_radius'] == 12.500003

@pytest.mark.order(0)
def test_1h2o_120meoh_md_prod_exdir_to_npz():
    """
    """
    dir_path, _, out_path_eq, geom_path_eq, traj_path_eq = get_1h2o_120meoh_eq_paths()
    exdir_path = dir_path + '.exdir'
    npz_path = dir_path + '-prod.npz'

    repman_exdir = manager()
    repman_exdir.load(exdir_path, mode='r')
    prod_dict = repman_exdir.data.as_dict('prod_1')
    
    data_npz = data(npz_path, mode='w', from_dict=prod_dict)
    data_npz.save()

    assert data_npz.get('geometry').shape == (1001, 723, 3)
    assert data_npz.get('geometry')[0][3][1] == -2.64576119977354
    assert data_npz.get('geometry')[-1][-1][0] == -1.62420092727186
    assert data_npz.get('energy_pot').shape == (1001,)
    assert data_npz.get('energy_pot')[0] == -991.881189902216
    assert data_npz.get('energy_pot')[-1] == -991.818996146108
    assert data_npz.get('wall_potential')[0]['sphere_radius'] == 12.500003
