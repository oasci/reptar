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

"""Tests and example builders for molecular dynamics with ASE."""

# pylint: skip-file

import os
import sys

import numpy as np
import pytest

from reptar import Creator, File
from reptar.utils import gen_comp_ids, gen_entity_ids

sys.path.append("..")
from .paths import get_1h2o_57h2o_pm_periodic_paths

# Ensures we execute from file directory (for relative paths).
os.chdir(os.path.dirname(os.path.realpath(__file__)))

# Writing paths
ASE_DIR = "./tmp/ase/"
os.makedirs(ASE_DIR, exist_ok=True)


@pytest.mark.order(0)
def test_1h2o_57h2o_ase_exdir():
    _, traj_path = get_1h2o_57h2o_pm_periodic_paths()
    exdir_path = os.path.join(ASE_DIR, "1h2o_57h2o_ase.exdir")

    num_waters = 58
    atoms_per_water = 3

    # For water molecules.
    entity_ids = gen_entity_ids(atoms_per_water, num_waters)
    comp_ids = gen_comp_ids("h2o", num_waters)

    group_key = "nvt_run_0"
    create = Creator()
    create.load(exdir_path, mode="w")
    create.from_calc(group_key, traj_path=traj_path)
    create.ids(group_key, entity_ids, comp_ids)

    create = Creator()
    create.load(exdir_path, mode="r")

    assert create.rfile.get(f"{group_key}/geometry").shape == (201, 174, 3)
    assert create.rfile.get(f"{group_key}/geometry")[0][3][1] == pytest.approx(
        1.5628686905646898, abs=1e-10
    )
    assert create.rfile.get(f"{group_key}/geometry")[-1][-1][0] == pytest.approx(
        6.727016221639617, abs=1e-10
    )
    assert create.rfile.get(f"{group_key}/energy_pot").shape == (201,)
    assert create.rfile.get(f"{group_key}/energy_pot")[0] == pytest.approx(
        -4427.15748237676, abs=1e-10
    )
    assert create.rfile.get(f"{group_key}/energy_pot")[-1] == pytest.approx(
        -4427.009820153819, abs=1e-10
    )
    assert create.rfile.get(f"{group_key}/temp")[24] == pytest.approx(
        289.7385847153916, abs=1e-10
    )
    assert create.rfile.get(f"{group_key}/velcs")[24][12][1] == pytest.approx(
        -0.003541032926883468, abs=1e-10
    )


def test_1h2o_57h2o_ase_zarr():
    _, traj_path = get_1h2o_57h2o_pm_periodic_paths()
    zarr_path = os.path.join(ASE_DIR, "1h2o_57h2o_ase.zarr")

    num_waters = 58
    atoms_per_water = 3

    # For water molecules.
    entity_ids = gen_entity_ids(atoms_per_water, num_waters)
    comp_ids = gen_comp_ids("h2o", num_waters)

    group_key = "nvt_run_0"
    create = Creator()
    create.load(zarr_path, mode="w")
    create.from_calc(group_key, traj_path=traj_path)
    create.ids(group_key, entity_ids, comp_ids)

    create = Creator()
    create.load(zarr_path, mode="r")

    assert create.rfile.get(f"{group_key}/geometry").shape == (201, 174, 3)
    assert create.rfile.get(f"{group_key}/geometry")[0][3][1] == pytest.approx(
        1.5628686905646898, abs=1e-10
    )
    assert create.rfile.get(f"{group_key}/geometry")[-1][-1][0] == pytest.approx(
        6.727016221639617, abs=1e-10
    )
    assert create.rfile.get(f"{group_key}/energy_pot").shape == (201,)
    assert create.rfile.get(f"{group_key}/energy_pot")[0] == pytest.approx(
        -4427.15748237676, abs=1e-10
    )
    assert create.rfile.get(f"{group_key}/energy_pot")[-1] == pytest.approx(
        -4427.009820153819, abs=1e-10
    )
    assert create.rfile.get(f"{group_key}/temp")[24] == pytest.approx(
        289.7385847153916, abs=1e-10
    )
    assert create.rfile.get(f"{group_key}/velcs")[24][12][1] == pytest.approx(
        -0.003541032926883468, abs=1e-10
    )
