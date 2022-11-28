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
from reptar import Creator, Sampler
from reptar.sampling import r_from_entities
import itertools

import sys

sys.path.append("..")
from .paths import *

# Ensures we execute from file directory (for relative paths).
os.chdir(os.path.dirname(os.path.realpath(__file__)))

# Source paths
data_dir = "./data/"

# Writing paths
sampling_dir = "./tmp/sampling/"
os.makedirs(sampling_dir, exist_ok=True)


def test_1h2o_120meoh_prod_sampler():
    """Sampling from xTB MD reptar file."""
    source_path = os.path.join(data_dir, "1h2o_120meoh_md.exdir")
    dest_path = f"{sampling_dir}/test_1h2o_120meoh_prod_sampling.exdir"

    source = Creator()
    source.load(source_path, mode="r")
    source_key = "/eq_1"

    dest = Creator()
    dest.load(dest_path, mode="w", allow_remove=True)
    dest_key = "/wat.2met-pes"
    dest.rfile.create_group(dest_key)

    quantity = 20
    comp_labels = ("WAT", "MET", "MET")

    sampler = Sampler(
        source.rfile,
        source_key,
        dest.rfile,
        dest_key,
        criteria=None,
        center_structures=False,
        E_key=None,
        G_key=None,
        dry_run=False,
        all_init_size=50000,
        use_ray=False,
        n_workers=2,
        ray_address="auto",
    )
    sampler.sample(comp_labels, quantity, R_source_idxs=None, specific_entities=None)

    assert np.array_equal(
        dest.rfile.get(f"{dest_key}/atomic_numbers"),
        np.array([8, 1, 1, 8, 1, 6, 1, 1, 1, 8, 1, 6, 1, 1, 1]),
    )
    assert np.array_equal(
        dest.rfile.get(f"{dest_key}/entity_ids"),
        np.array([0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2]),
    )
    assert dest.rfile.get(f"{dest_key}/geometry").shape == (quantity, 15, 3)
    assert dest.rfile.get(f"{dest_key}/r_prov_specs").shape == (quantity, 5)
    assert len(dest.rfile.get(f"{dest_key}/r_prov_ids")) == 1
    assert dest.rfile.get(f"{dest_key}/r_centered") == False
    assert np.array_equal(dest.rfile.get(f"{dest_key}/comp_ids"), np.array(comp_labels))

    # Checks that the correct geometries are sampled and in the correct order.
    R_source = source.rfile.get(f"{source_key}/geometry")
    entity_ids_source = source.rfile.get(f"{source_key}/entity_ids")

    R_sampled = dest.rfile.get(f"{dest_key}/geometry")
    r_prov_specs = dest.rfile.get(f"{dest_key}/r_prov_specs")

    for i in range(len(r_prov_specs)):
        r_sampled = R_sampled[i]
        r_prov_spec = r_prov_specs[i]
        source_idx = r_prov_spec[1]

        idx_atom_check = 0
        for entity_id in r_prov_spec[2:]:
            entity_mask = entity_ids_source == entity_id
            r_source_frag = R_source[source_idx][entity_mask]
            assert r_sampled[idx_atom_check][0] == r_source_frag[0][0]
            idx_atom_check += int(r_source_frag.shape[0])

    # Test additional sampling
    # Now we center structures to test.
    sampler.center_structures = True
    sampler.sample(comp_labels, quantity, R_source_idxs=None, specific_entities=None)
    assert np.array_equal(
        dest.rfile.get(f"{dest_key}/atomic_numbers"),
        np.array([8, 1, 1, 8, 1, 6, 1, 1, 1, 8, 1, 6, 1, 1, 1]),
    )
    assert np.array_equal(
        dest.rfile.get(f"{dest_key}/entity_ids"),
        np.array([0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2]),
    )
    assert dest.rfile.get(f"{dest_key}/geometry").shape == (quantity * 2, 15, 3)
    assert dest.rfile.get(f"{dest_key}/r_prov_specs").shape == (quantity * 2, 5)
    assert len(dest.rfile.get(f"{dest_key}/r_prov_ids")) == 1
    assert dest.rfile.get(f"{dest_key}/r_centered") == True
    assert np.array_equal(dest.rfile.get(f"{dest_key}/comp_ids"), np.array(comp_labels))


def test_ray_1h2o_120meoh_prod_sampler():
    """Sampling from xTB MD reptar file with ray."""
    try:
        import ray
    except (ModuleNotFoundError, ImportError):
        pytest.skip("ray package not installed")

    source_path = os.path.join(data_dir, "1h2o_120meoh_md.exdir")
    dest_path = f"{sampling_dir}/test_1h2o_120meoh_prod_sampling_ray.exdir"

    source = Creator()
    source.load(source_path, mode="r")
    source_key = "/eq_1"

    dest = Creator()
    dest.load(dest_path, mode="w", allow_remove=True)
    dest_key = "/wat.2met-pes"
    dest.rfile.create_group(dest_key)

    quantity = 20
    comp_labels = ("WAT", "MET", "MET")

    sampler = Sampler(
        source.rfile,
        source_key,
        dest.rfile,
        dest_key,
        criteria=None,
        center_structures=False,
        E_key=None,
        G_key=None,
        dry_run=False,
        all_init_size=50000,
        use_ray=True,
        n_workers=2,
        ray_address="auto",
    )
    sampler.sample(comp_labels, quantity, R_source_idxs=None, specific_entities=None)

    assert np.array_equal(
        dest.rfile.get(f"{dest_key}/atomic_numbers"),
        np.array([8, 1, 1, 8, 1, 6, 1, 1, 1, 8, 1, 6, 1, 1, 1]),
    )
    assert np.array_equal(
        dest.rfile.get(f"{dest_key}/entity_ids"),
        np.array([0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2]),
    )
    assert dest.rfile.get(f"{dest_key}/geometry").shape == (quantity, 15, 3)
    assert dest.rfile.get(f"{dest_key}/r_prov_specs").shape == (quantity, 5)
    assert len(dest.rfile.get(f"{dest_key}/r_prov_ids")) == 1
    assert dest.rfile.get(f"{dest_key}/r_centered") == False
    assert np.array_equal(dest.rfile.get(f"{dest_key}/comp_ids"), np.array(comp_labels))

    # Checks that the correct geometries are sampled and in the correct order.
    R_source = source.rfile.get(f"{source_key}/geometry")
    entity_ids_source = source.rfile.get(f"{source_key}/entity_ids")

    R_sampled = dest.rfile.get(f"{dest_key}/geometry")
    r_prov_specs = dest.rfile.get(f"{dest_key}/r_prov_specs")

    for i in range(len(r_prov_specs)):
        r_sampled = R_sampled[i]
        r_prov_spec = r_prov_specs[i]
        source_idx = r_prov_spec[1]

        idx_atom_check = 0
        for entity_id in r_prov_spec[2:]:
            entity_mask = entity_ids_source == entity_id
            r_source_frag = R_source[source_idx][entity_mask]
            assert r_sampled[idx_atom_check][0] == r_source_frag[0][0]
            idx_atom_check += int(r_source_frag.shape[0])

    # Test additional sampling
    # Now we center structures to test.
    sampler.center_structures = True
    sampler.sample(comp_labels, quantity, R_source_idxs=None, specific_entities=None)
    assert np.array_equal(
        dest.rfile.get(f"{dest_key}/atomic_numbers"),
        np.array([8, 1, 1, 8, 1, 6, 1, 1, 1, 8, 1, 6, 1, 1, 1]),
    )
    assert np.array_equal(
        dest.rfile.get(f"{dest_key}/entity_ids"),
        np.array([0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2]),
    )
    assert dest.rfile.get(f"{dest_key}/geometry").shape == (quantity * 2, 15, 3)
    assert dest.rfile.get(f"{dest_key}/r_prov_specs").shape == (quantity * 2, 5)
    assert len(dest.rfile.get(f"{dest_key}/r_prov_ids")) == 1
    assert dest.rfile.get(f"{dest_key}/r_centered") == True
    assert np.array_equal(dest.rfile.get(f"{dest_key}/comp_ids"), np.array(comp_labels))


def test_sampling_from_wat_2met_pes():
    """Sampling from a sampled exdir group."""
    src_path = f"{sampling_dir}/test_1h2o_120meoh_prod_sampling.exdir"

    source = Creator()
    source.load(src_path, mode="a")
    source_key = "/wat.2met-pes"
    dest_key = "/wat.met-pes"
    source.rfile.create_group(dest_key)

    quantity = "all"
    comp_labels = ("WAT", "MET")

    sampler = Sampler(
        source.rfile,
        source_key,
        source.rfile,
        dest_key,
        criteria=None,
        center_structures=False,
        E_key=None,
        G_key=None,
        dry_run=False,
        all_init_size=50000,
        use_ray=False,
        n_workers=2,
        ray_address="auto",
    )
    sampler.sample(comp_labels, quantity, R_source_idxs=None, specific_entities=None)

    entity_ids_source = source.rfile.get(f"{source_key}/entity_ids")
    entity_ids_dest = source.rfile.get(f"{dest_key}/entity_ids")
    assert np.array_equal(entity_ids_dest, np.array([0, 0, 0, 1, 1, 1, 1, 1, 1]))
    assert len(source.rfile.get(f"{dest_key}/r_prov_ids")) == 1
    assert source.rfile.get(f"{dest_key}/r_centered") == False
    assert np.array_equal(
        source.rfile.get(f"{dest_key}/comp_ids"), np.array(comp_labels)
    )

    # We know that the first structure in destination is the first possible
    # dimer from the source. This is only because quantity is 'all'.
    R_source = source.rfile.get(f"{source_key}/geometry")
    R_dest = source.rfile.get(f"{dest_key}/geometry")
    R_ref = r_from_entities(R_source[0], entity_ids_source, (0, 1))
    assert np.allclose(R_ref, R_dest[0])

    assert np.array_equal(
        source.rfile.get(f"{source_key}/r_prov_specs")[0][:4],
        source.rfile.get(f"{dest_key}/r_prov_specs")[0],
    )

    # We also know the last structure is the last possible dimer.
    assert np.array_equal(
        source.rfile.get(f"{source_key}/r_prov_specs")[-1][np.array([0, 1, 2, 4])],
        source.rfile.get(f"{dest_key}/r_prov_specs")[-1],
    )


def test_ray_sampling_from_wat_2met_pes():
    """Sampling from a sampled exdir group with ray."""
    try:
        import ray
    except (ModuleNotFoundError, ImportError):
        pytest.skip("ray package not installed")

    src_path = f"{sampling_dir}/test_1h2o_120meoh_prod_sampling_ray.exdir"

    source = Creator()
    source.load(src_path, mode="a")
    source_key = "/wat.2met-pes"
    dest_key = "/wat.met-pes"
    source.rfile.create_group(dest_key)

    quantity = "all"
    comp_labels = ("WAT", "MET")

    sampler = Sampler(
        source.rfile,
        source_key,
        source.rfile,
        dest_key,
        criteria=None,
        center_structures=False,
        E_key=None,
        G_key=None,
        dry_run=False,
        all_init_size=50000,
        use_ray=True,
        n_workers=2,
        ray_address="auto",
    )
    sampler.sample(comp_labels, quantity, R_source_idxs=None, specific_entities=None)

    entity_ids_source = source.rfile.get(f"{source_key}/entity_ids")
    entity_ids_dest = source.rfile.get(f"{dest_key}/entity_ids")
    assert np.array_equal(entity_ids_dest, np.array([0, 0, 0, 1, 1, 1, 1, 1, 1]))
    assert len(source.rfile.get(f"{dest_key}/r_prov_ids")) == 1
    assert source.rfile.get(f"{dest_key}/r_centered") == False
    assert np.array_equal(
        source.rfile.get(f"{dest_key}/comp_ids"), np.array(comp_labels)
    )

    # We know that the first structure in destination is the first possible
    # dimer from the source. This is only because quantity is 'all'.
    R_source = source.rfile.get(f"{source_key}/geometry")
    R_dest = source.rfile.get(f"{dest_key}/geometry")
    R_ref = r_from_entities(R_source[0], entity_ids_source, (0, 1))
    assert np.allclose(R_ref, R_dest[0])

    assert np.array_equal(
        source.rfile.get(f"{source_key}/r_prov_specs")[0][:4],
        source.rfile.get(f"{dest_key}/r_prov_specs")[0],
    )

    # We also know the last structure is the last possible dimer.
    assert np.array_equal(
        source.rfile.get(f"{source_key}/r_prov_specs")[-1][np.array([0, 1, 2, 4])],
        source.rfile.get(f"{dest_key}/r_prov_specs")[-1],
    )
