# MIT License
#
# Copyright (c) 2022-2023, Alex M. Maldonado
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

"""Tests for Psi4 calculations"""

# pylint: skip-file

import sys
import shutil
import os
import pytest
import numpy as np
from reptar import File, Saver
from reptar.calculators.drivers import DriverEnergy, DriverEnGrad
from reptar.calculators.psi4_workers import psi4_energy, psi4_engrad


sys.path.append("..")
from .paths import get_140h2o_samples_path

# Ensures we execute from file directory (for relative paths).
os.chdir(os.path.dirname(os.path.realpath(__file__)))

# Writing paths
WRITING_DIR = "./tmp/calculators/"
os.makedirs(WRITING_DIR, exist_ok=True)


def test_calculator_psi4_1h2o_engrad():
    try:
        import psi4
    except ImportError:
        pytest.skip("psi4 package not installed")

    exdir_path_source = get_140h2o_samples_path()
    rfile_source = File(exdir_path_source, mode="r")

    exdir_path_dest = os.path.join(WRITING_DIR, "1h2o-psi4.exdir")

    if os.path.exists(exdir_path_dest):
        shutil.rmtree(exdir_path_dest)
    rfile = File(exdir_path_dest, mode="w")

    # Copy over a few structures for calculations.
    group_key = "1h2o"
    start_slice = None
    end_slice = 5

    R_ref = np.array(
        [
            [
                [-0.060401859445, -0.008161713319, 0.038695740786],
                [0.081656670824, 0.256790001513, -0.919880479138],
                [0.877197077102, -0.127226284815, 0.305602103897],
            ],
            [
                [-0.029215828214, -0.054854435349, -0.009617450002],
                [0.564775863493, 0.357610588759, -0.602273276481],
                [-0.100987055674, 0.513180168492, 0.754946191924],
            ],
            [
                [-0.028530376949, -0.014602259575, -0.051906786482],
                [0.797151206373, 0.499067262665, -0.076776383384],
                [-0.344243644933, -0.267262622045, 0.900774472492],
            ],
            [
                [-0.026402272283, -0.048446860866, 0.021616226734],
                [0.919026404127, -0.099176175566, -0.267588806956],
                [-0.499901596297, 0.868249421724, -0.075559569775],
            ],
            [
                [0.018370719066, -0.050323835823, 0.047003163594],
                [-0.619142950603, 0.648350631239, 0.10363099464],
                [0.327515623389, 0.150518791537, -0.849786162559],
            ],
        ]
    )
    Z = rfile_source.get(f"{group_key}/atomic_numbers")
    R = rfile_source.get(f"{group_key}/geometry")[start_slice:end_slice]
    assert np.allclose(R, R_ref)

    rfile.put(f"{group_key}/atomic_numbers", Z)
    rfile.put(f"{group_key}/geometry", R)

    # Setup energy and gradient arrays
    method_label = "df.mp2.def2tzvp"

    E_key = f"{group_key}/energy_ele_{method_label}"
    E = np.empty(R.shape[0], dtype=np.float64)
    E[:] = np.nan
    rfile.put(E_key, E)

    G_key = f"{group_key}/grads_{method_label}"
    G = np.empty(R.shape, dtype=np.float64)
    G[:] = np.nan
    rfile.put(G_key, G)

    driver_kwargs = {
        "use_ray": False,
        "n_workers": 1,
        "n_cpus_per_worker": 1,
        "chunk_size": 1,
        "start_slice": None,
        "end_slice": None,
    }

    worker_kwargs = {
        "charge": 0,
        "mult": 1,
        "method": "mp2",
        "threads": 1,
        "options": {
            "basis": "def2-tzvppd",
            "df_basis_scf": "def2-universal-jkfit",
            "df_basis_mp2": "def2-tzvppd-ri",
            "reference": "rhf",
            "e_convergence": 10,
            "d_convergence": 10,
            "scf_type": "df",
            "mp2_type": "df",
            "qc_module": "dfmp2",
            "print": 2,
        },
    }

    saver = Saver(exdir_path_dest, (E_key, G_key))

    driver = DriverEnGrad(psi4_engrad, worker_kwargs, **driver_kwargs)

    saver.save(E, G)

    E_ref = np.array(
        [
            -76.343645506584,
            -76.345480677691,
            -76.338126861147,
            -76.335569277436,
            -76.344566833326,
        ]
    )
    G_ref = np.array(
        [
            [
                [-0.056453332701, -0.016034305011, 0.065841463133],
                [0.020135069212, 0.017627192921, -0.065522078855],
                [0.036318263489, -0.00159288791, -0.000319384278],
            ],
            [
                [0.035116948314, 0.0326620761, -0.025755143029],
                [-0.031198527186, -0.026192952531, 0.026049191432],
                [-0.003918421127, -0.006469123568, -0.000294048403],
            ],
            [
                [0.035245100851, 0.026091664752, -0.070900241198],
                [0.025238572836, 0.017529978829, -0.031422241428],
                [-0.060483673688, -0.043621643581, 0.102322482627],
            ],
            [
                [0.026718369939, -0.059122634819, 0.007491615688],
                [0.060197098403, -0.038574411766, -0.008816867362],
                [-0.086915468342, 0.097697046585, 0.001325251674],
            ],
            [
                [-0.016459922509, -0.002869394785, 0.035143984283],
                [0.027302614579, -0.016913072538, -0.023381485441],
                [-0.01084269207, 0.019782467323, -0.011762498842],
            ],
        ]
    )
    E, G = driver.run(Z, R, E, G, saver=saver)

    assert np.allclose(E, E_ref)
    assert np.allclose(G, G_ref)


def test_ray_calculator_psi4_1h2o_engrad():
    try:
        import psi4
    except ImportError:
        pytest.skip("psi4 package not installed")
    try:
        import ray
    except ImportError:
        pytest.skip("ray package not installed")

    exdir_path_source = get_140h2o_samples_path()
    rfile_source = File(exdir_path_source, mode="r")

    exdir_path_dest = os.path.join(WRITING_DIR, "1h2o-psi4.exdir")

    if os.path.exists(exdir_path_dest):
        shutil.rmtree(exdir_path_dest)
    rfile = File(exdir_path_dest, mode="w")

    # Copy over a few structures for calculations.
    group_key = "1h2o"
    start_slice = None
    end_slice = 5

    R_ref = np.array(
        [
            [
                [-0.060401859445, -0.008161713319, 0.038695740786],
                [0.081656670824, 0.256790001513, -0.919880479138],
                [0.877197077102, -0.127226284815, 0.305602103897],
            ],
            [
                [-0.029215828214, -0.054854435349, -0.009617450002],
                [0.564775863493, 0.357610588759, -0.602273276481],
                [-0.100987055674, 0.513180168492, 0.754946191924],
            ],
            [
                [-0.028530376949, -0.014602259575, -0.051906786482],
                [0.797151206373, 0.499067262665, -0.076776383384],
                [-0.344243644933, -0.267262622045, 0.900774472492],
            ],
            [
                [-0.026402272283, -0.048446860866, 0.021616226734],
                [0.919026404127, -0.099176175566, -0.267588806956],
                [-0.499901596297, 0.868249421724, -0.075559569775],
            ],
            [
                [0.018370719066, -0.050323835823, 0.047003163594],
                [-0.619142950603, 0.648350631239, 0.10363099464],
                [0.327515623389, 0.150518791537, -0.849786162559],
            ],
        ]
    )
    Z = rfile_source.get(f"{group_key}/atomic_numbers")
    R = rfile_source.get(f"{group_key}/geometry")[start_slice:end_slice]
    assert np.allclose(R, R_ref)

    rfile.put(f"{group_key}/atomic_numbers", Z)
    rfile.put(f"{group_key}/geometry", R)

    # Setup energy and gradient arrays
    method_label = "df.mp2.def2tzvp"

    E_key = f"{group_key}/energy_ele_{method_label}"
    E = np.empty(R.shape[0], dtype=np.float64)
    E[:] = np.nan
    rfile.put(E_key, E)

    G_key = f"{group_key}/grads_{method_label}"
    G = np.empty(R.shape, dtype=np.float64)
    G[:] = np.nan
    rfile.put(G_key, G)

    driver_kwargs = {
        "use_ray": True,
        "n_workers": 2,
        "n_cpus_per_worker": 1,
        "chunk_size": 2,
        "start_slice": None,
        "end_slice": None,
    }

    worker_kwargs = {
        "charge": 0,
        "mult": 1,
        "method": "mp2",
        "threads": 1,
        "options": {
            "basis": "def2-tzvppd",
            "df_basis_scf": "def2-universal-jkfit",
            "df_basis_mp2": "def2-tzvppd-ri",
            "reference": "rhf",
            "e_convergence": 10,
            "d_convergence": 10,
            "scf_type": "df",
            "mp2_type": "df",
            "qc_module": "dfmp2",
            "print": 2,
        },
    }

    saver = Saver(exdir_path_dest, (E_key, G_key))

    driver = DriverEnGrad(psi4_engrad, worker_kwargs, **driver_kwargs)

    saver.save(E, G)

    E_ref = np.array(
        [
            -76.343645506584,
            -76.345480677691,
            -76.338126861147,
            -76.335569277436,
            -76.344566833326,
        ]
    )
    G_ref = np.array(
        [
            [
                [-0.056453332701, -0.016034305011, 0.065841463133],
                [0.020135069212, 0.017627192921, -0.065522078855],
                [0.036318263489, -0.00159288791, -0.000319384278],
            ],
            [
                [0.035116948314, 0.0326620761, -0.025755143029],
                [-0.031198527186, -0.026192952531, 0.026049191432],
                [-0.003918421127, -0.006469123568, -0.000294048403],
            ],
            [
                [0.035245100851, 0.026091664752, -0.070900241198],
                [0.025238572836, 0.017529978829, -0.031422241428],
                [-0.060483673688, -0.043621643581, 0.102322482627],
            ],
            [
                [0.026718369939, -0.059122634819, 0.007491615688],
                [0.060197098403, -0.038574411766, -0.008816867362],
                [-0.086915468342, 0.097697046585, 0.001325251674],
            ],
            [
                [-0.016459922509, -0.002869394785, 0.035143984283],
                [0.027302614579, -0.016913072538, -0.023381485441],
                [-0.01084269207, 0.019782467323, -0.011762498842],
            ],
        ]
    )
    E, G = driver.run(Z, R, E, G, saver=saver)

    assert np.allclose(E, E_ref)
    assert np.allclose(G, G_ref)


def test_calculator_psi4_1h2o_energy():
    try:
        import psi4
    except ImportError:
        pytest.skip("psi4 package not installed")

    exdir_path_source = get_140h2o_samples_path()
    rfile_source = File(exdir_path_source, mode="r")

    exdir_path_dest = os.path.join(WRITING_DIR, "1h2o-psi4.exdir")

    if os.path.exists(exdir_path_dest):
        shutil.rmtree(exdir_path_dest)
    rfile = File(exdir_path_dest, mode="w")

    # Copy over a few structures for calculations.
    group_key = "1h2o"
    start_slice = None
    end_slice = 5

    R_ref = np.array(
        [
            [
                [-0.060401859445, -0.008161713319, 0.038695740786],
                [0.081656670824, 0.256790001513, -0.919880479138],
                [0.877197077102, -0.127226284815, 0.305602103897],
            ],
            [
                [-0.029215828214, -0.054854435349, -0.009617450002],
                [0.564775863493, 0.357610588759, -0.602273276481],
                [-0.100987055674, 0.513180168492, 0.754946191924],
            ],
            [
                [-0.028530376949, -0.014602259575, -0.051906786482],
                [0.797151206373, 0.499067262665, -0.076776383384],
                [-0.344243644933, -0.267262622045, 0.900774472492],
            ],
            [
                [-0.026402272283, -0.048446860866, 0.021616226734],
                [0.919026404127, -0.099176175566, -0.267588806956],
                [-0.499901596297, 0.868249421724, -0.075559569775],
            ],
            [
                [0.018370719066, -0.050323835823, 0.047003163594],
                [-0.619142950603, 0.648350631239, 0.10363099464],
                [0.327515623389, 0.150518791537, -0.849786162559],
            ],
        ]
    )
    Z = rfile_source.get(f"{group_key}/atomic_numbers")
    R = rfile_source.get(f"{group_key}/geometry")[start_slice:end_slice]
    assert np.allclose(R, R_ref)

    rfile.put(f"{group_key}/atomic_numbers", Z)
    rfile.put(f"{group_key}/geometry", R)

    # Setup energy and gradient arrays
    method_label = "df.mp2.def2tzvp"

    E_key = f"{group_key}/energy_ele_{method_label}"
    E = np.empty(R.shape[0], dtype=np.float64)
    E[:] = np.nan
    rfile.put(E_key, E)

    driver_kwargs = {
        "use_ray": False,
        "n_workers": 1,
        "n_cpus_per_worker": 1,
        "chunk_size": 1,
        "start_slice": None,
        "end_slice": None,
    }

    worker_kwargs = {
        "charge": 0,
        "mult": 1,
        "method": "mp2",
        "threads": 1,
        "options": {
            "basis": "def2-tzvppd",
            "df_basis_scf": "def2-universal-jkfit",
            "df_basis_mp2": "def2-tzvppd-ri",
            "reference": "rhf",
            "e_convergence": 10,
            "d_convergence": 10,
            "scf_type": "df",
            "mp2_type": "df",
            "qc_module": "dfmp2",
            "print": 2,
        },
    }

    saver = Saver(exdir_path_dest, (E_key,))

    driver = DriverEnergy(psi4_energy, worker_kwargs, **driver_kwargs)

    saver.save(E)

    E_ref = np.array(
        [
            -76.343645506584,
            -76.345480677691,
            -76.338126861147,
            -76.335569277436,
            -76.344566833326,
        ]
    )

    E = driver.run(Z, R, E, saver=saver)

    assert np.allclose(E, E_ref)


def test_ray_calculator_psi4_1h2o_energy():
    try:
        import psi4
    except ImportError:
        pytest.skip("psi4 package not installed")
    try:
        import ray
    except ImportError:
        pytest.skip("ray package not installed")

    exdir_path_source = get_140h2o_samples_path()
    rfile_source = File(exdir_path_source, mode="r")

    exdir_path_dest = os.path.join(WRITING_DIR, "1h2o-psi4.exdir")

    if os.path.exists(exdir_path_dest):
        shutil.rmtree(exdir_path_dest)
    rfile = File(exdir_path_dest, mode="w")

    # Copy over a few structures for calculations.
    group_key = "1h2o"
    start_slice = None
    end_slice = 5

    R_ref = np.array(
        [
            [
                [-0.060401859445, -0.008161713319, 0.038695740786],
                [0.081656670824, 0.256790001513, -0.919880479138],
                [0.877197077102, -0.127226284815, 0.305602103897],
            ],
            [
                [-0.029215828214, -0.054854435349, -0.009617450002],
                [0.564775863493, 0.357610588759, -0.602273276481],
                [-0.100987055674, 0.513180168492, 0.754946191924],
            ],
            [
                [-0.028530376949, -0.014602259575, -0.051906786482],
                [0.797151206373, 0.499067262665, -0.076776383384],
                [-0.344243644933, -0.267262622045, 0.900774472492],
            ],
            [
                [-0.026402272283, -0.048446860866, 0.021616226734],
                [0.919026404127, -0.099176175566, -0.267588806956],
                [-0.499901596297, 0.868249421724, -0.075559569775],
            ],
            [
                [0.018370719066, -0.050323835823, 0.047003163594],
                [-0.619142950603, 0.648350631239, 0.10363099464],
                [0.327515623389, 0.150518791537, -0.849786162559],
            ],
        ]
    )
    Z = rfile_source.get(f"{group_key}/atomic_numbers")
    R = rfile_source.get(f"{group_key}/geometry")[start_slice:end_slice]
    assert np.allclose(R, R_ref)

    rfile.put(f"{group_key}/atomic_numbers", Z)
    rfile.put(f"{group_key}/geometry", R)

    # Setup energy and gradient arrays
    method_label = "df.mp2.def2tzvp"

    E_key = f"{group_key}/energy_ele_{method_label}"
    E = np.empty(R.shape[0], dtype=np.float64)
    E[:] = np.nan
    rfile.put(E_key, E)

    driver_kwargs = {
        "use_ray": True,
        "n_workers": 2,
        "n_cpus_per_worker": 1,
        "chunk_size": 2,
        "start_slice": None,
        "end_slice": None,
    }

    worker_kwargs = {
        "charge": 0,
        "mult": 1,
        "method": "mp2",
        "threads": 1,
        "options": {
            "basis": "def2-tzvppd",
            "df_basis_scf": "def2-universal-jkfit",
            "df_basis_mp2": "def2-tzvppd-ri",
            "reference": "rhf",
            "e_convergence": 10,
            "d_convergence": 10,
            "scf_type": "df",
            "mp2_type": "df",
            "qc_module": "dfmp2",
            "print": 2,
        },
    }

    saver = Saver(exdir_path_dest, (E_key,))

    driver = DriverEnergy(psi4_energy, worker_kwargs, **driver_kwargs)

    saver.save(E)

    E_ref = np.array(
        [
            -76.343645506584,
            -76.345480677691,
            -76.338126861147,
            -76.335569277436,
            -76.344566833326,
        ]
    )

    E = driver.run(Z, R, E, saver=saver)

    assert np.allclose(E, E_ref)
