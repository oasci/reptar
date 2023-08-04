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

"""Tests for reptar-calc script"""

# pylint: skip-file

import subprocess
import sys
import shutil
import os
import pytest
import numpy as np
from reptar import File
import qcelemental as qcel


sys.path.append("..")
from .paths import get_140h2o_samples_path

# Ensures we execute from file directory (for relative paths).
os.chdir(os.path.dirname(os.path.realpath(__file__)))

# Writing paths
WRITING_DIR = "./tmp/calculators/"
os.makedirs(WRITING_DIR, exist_ok=True)


def test_script_psi4_engrad():
    exdir_path_source = get_140h2o_samples_path()
    rfile_source = File(exdir_path_source, mode="r")

    path_dest = os.path.join(WRITING_DIR, "1h2o-psi4-script.zarr")

    if os.path.exists(path_dest):
        shutil.rmtree(path_dest)
    rfile = File(path_dest, mode="w")

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

    E_ref = np.array(
        [
            -76.324720546651,
            -76.326430526027,
            -76.319203316807,
            -76.316664051316,
            -76.325563963740,
        ]
    )
    G_ref = np.array(
        [
            [
                [-0.055221236465, -0.015923938795, 0.065250246288],
                [0.019954697536, 0.017399784508, -0.064689768889],
                [0.035266538929, -0.001475845713, -0.000560477399],
            ],
            [
                [0.036079133182, 0.034255674828, -0.025677243030],
                [-0.032235782259, -0.026943315044, 0.027050358924],
                [-0.003843350923, -0.007312359783, -0.001373115895],
            ],
            [
                [0.036081127644, 0.026572729767, -0.070270813882],
                [0.024223984734, 0.016900786434, -0.031425197594],
                [-0.060305112378, -0.043473516201, 0.101696011476],
            ],
            [
                [0.027419791458, -0.058540364025, 0.007108724011],
                [0.059224212124, -0.038558554903, -0.008509393861],
                [-0.086644003582, 0.097098918928, 0.001400669850],
            ],
            [
                [-0.017082675549, -0.001472896446, 0.034049099197],
                [0.028273021961, -0.018011563022, -0.023411323608],
                [-0.011190346412, 0.019484459468, -0.010637775589],
            ],
        ]
    )

    command = [
        "reptar-calc",
        "./data/script-configs/config-psi4-engrads.yaml",
        "--log_level",
        "debug",
    ]
    subprocess.run(command, check=True, shell=True)

    assert np.allclose(rfile.get(f"{group_key}/energy_ele_mp2.def2tzvp"), E_ref)
    assert np.allclose(rfile.get(f"{group_key}/grads_mp2.def2tzvp"), G_ref)


def test_script_xtb_opt():
    if shutil.which("xtb") is None:
        pytest.skip("xtb package not installed")

    exdir_path_source = get_140h2o_samples_path()
    rfile_source = File(exdir_path_source, mode="r")

    path_dest = os.path.join(WRITING_DIR, "1h2o-xtb-script.zarr")

    if os.path.exists(path_dest):
        shutil.rmtree(path_dest)
    rfile = File(path_dest, mode="w")

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
    rfile.put(f"{group_key}/geometry_init", R)

    command = [
        "reptar-calc",
        "./data/script-configs/config-xtb-opt.yaml",
        "--log_level",
        "debug",
    ]
    subprocess.run(command, check=True, shell=True)

    R_opt_ref = np.array(
        [
            [
                [-2.17748066e-02, -6.00011964e-04, 5.68107361e-03],
                [5.08637699e-02, 2.56960978e-01, -9.15445558e-01],
                [8.69362925e-01, -1.34958963e-01, 3.34181850e-01],
            ],
            [
                [-3.45410880e-02, -5.82001057e-02, -3.90920794e-03],
                [5.74136373e-01, 3.60327538e-01, -6.15855369e-01],
                [-1.05022305e-01, 5.13808890e-01, 7.62820042e-01],
            ],
            [
                [-5.14412257e-02, -2.87345922e-02, -5.32433476e-02],
                [7.64062111e-01, 4.75688796e-01, -2.89445751e-02],
                [-2.88243701e-01, -2.29751823e-01, 8.54279225e-01],
            ],
            [
                [-5.22537078e-02, -6.08576268e-02, 3.32707800e-02],
                [8.59051639e-01, -3.74425799e-02, -2.65131775e-01],
                [-4.14075395e-01, 8.18926592e-01, -8.96711551e-02],
            ],
            [
                [1.14998840e-02, -2.09476581e-02, 1.24209049e-02],
                [-6.53408474e-01, 6.55319370e-01, 1.55943370e-01],
                [3.68651982e-01, 1.14173875e-01, -8.67516279e-01],
            ],
        ]
    )
    E_ref = np.array([-5.07054445, -5.07054444, -5.07054445, -5.07054445, -5.07054445])
    conv_opt_ref = np.array([True, True, True, True, True])

    assert np.allclose(rfile.get(f"{group_key}/geometry"), R_opt_ref)
    assert np.allclose(rfile.get(f"{group_key}/energy_ele_gfn2"), E_ref)
    assert np.allclose(rfile.get(f"{group_key}/conv_opt"), conv_opt_ref)
