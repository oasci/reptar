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

"""Tests and example builders for ORCA calculations."""

# pylint: skip-file

import os
import sys

import numpy as np
import pytest

from reptar import Creator
from reptar.utils import gen_comp_ids, gen_entity_ids

sys.path.append("..")
from .paths import get_6h2o_temelso_pr_engrad

# Ensures we execute from file directory (for relative paths).
os.chdir(os.path.dirname(os.path.realpath(__file__)))

# Writing paths
ORCA_DIR = "./tmp/orca"
os.makedirs(ORCA_DIR, exist_ok=True)


@pytest.mark.order(0)
def test_6h2o_temelso_engrad_exdir():
    _, out_path = get_6h2o_temelso_pr_engrad()
    exdir_path = os.path.join(ORCA_DIR, "6h2o_temelso_pr_engrad.exdir")

    create = Creator()
    create.load(exdir_path, mode="w")
    create.definitions(definitions=["qc", "pes", "molprop"])
    create.from_calc("/engrad", out_path=out_path)

    # Reading tests.
    create = Creator()
    create.load(exdir_path, mode="r")
    assert create.rfile.get("/engrad/energy_scf") == -456.40308701
    assert create.rfile.get("/engrad/energy_ele") == -457.961331933491
    assert np.array_equal(
        create.rfile.get("/engrad/dipole_moment"),
        np.array([2.5265211700823, -0.5830003327751, 1.2353903376292001]),
    )


def test_6h2o_temelso_engrad_zarr():
    _, out_path = get_6h2o_temelso_pr_engrad()
    zarr_path = os.path.join(ORCA_DIR, "6h2o_temelso_pr_engrad.zarr")

    create = Creator()
    create.load(zarr_path, mode="w")
    create.from_calc("/engrad", out_path=out_path)

    # Reading tests.
    create = Creator()
    create.load(zarr_path, mode="r")
    assert create.rfile.get("/engrad/energy_scf") == -456.40308701
    assert create.rfile.get("/engrad/energy_ele") == -457.961331933491
    assert np.array_equal(
        create.rfile.get("/engrad/dipole_moment"),
        np.array([2.5265211700823, -0.5830003327751, 1.2353903376292001]),
    )


def test_6h2o_temelso_engrad_json():
    _, out_path = get_6h2o_temelso_pr_engrad()
    json_path = os.path.join(ORCA_DIR, "6h2o_temelso_pr_engrad.json")

    # Writing tests.
    create = Creator()
    create.load(json_path, mode="w")
    create.from_calc("/", out_path=out_path)
    create.rfile.save()

    # Reading tests.
    create = Creator()
    create.load(json_path, mode="r")
    assert create.rfile.get("energy_scf") == -456.40308701
    assert create.rfile.get("energy_ele") == -457.961331933491
    assert np.array_equal(
        create.rfile.get("dipole_moment"),
        np.array([2.5265211700823, -0.5830003327751, 1.2353903376292001]),
    )


def test_6h2o_temelso_engrad_npz():
    _, out_path = get_6h2o_temelso_pr_engrad()
    npz_path = os.path.join(ORCA_DIR, "6h2o_temelso_pr_engrad.npz")

    create = Creator()
    create.load(npz_path, mode="w")
    create.from_calc("/", out_path=out_path)
    create.rfile.save()

    # Reading tests.
    create = Creator()
    create.load(npz_path, mode="r")
    assert create.rfile.get("energy_scf") == -456.40308701
    assert create.rfile.get("energy_ele") == -457.961331933491
    assert np.array_equal(
        create.rfile.get("dipole_moment"),
        np.array([2.5265211700823, -0.5830003327751, 1.2353903376292001]),
    )
