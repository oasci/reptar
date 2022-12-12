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

"""Tests for xTB calculations"""

# Stuff for PyTest features like skip.
# pylint: disable=import-outside-toplevel, unused-import, duplicate-code

import sys
import shutil
import os
import pytest
import numpy as np
from reptar import File, Saver
from reptar.calculators.drivers import DriverEnergy, DriverEnGrad
from reptar.calculators.xtb_workers import xtb_engrad


sys.path.append("..")
# pylint: disable-next=wrong-import-position
from .paths import get_140h2o_samples_path

# Ensures we execute from file directory (for relative paths).
os.chdir(os.path.dirname(os.path.realpath(__file__)))

# Writing paths
WRITING_DIR = "./tmp/calculators/"
os.makedirs(WRITING_DIR, exist_ok=True)


def test_calculator_xtb_1h2o_engrad():
    try:
        import xtb
    except ImportError:
        pytest.skip("xtb package not installed")

    exdir_path_source = get_140h2o_samples_path()
    rfile_source = File(exdir_path_source, mode="r")

    exdir_path_dest = os.path.join(WRITING_DIR, "1h2o-xtb.exdir")

    if os.path.exists(exdir_path_dest):
        shutil.rmtree(exdir_path_dest)
    rfile = File(exdir_path_dest, mode="w")

    # Copy over a few structures for calculations.
    group_key = "1h2o"
    start_slice = None
    end_slice = 5

    # pylint: disable-next=invalid-name
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
    method_label = "gfn2xtb"

    E_key = f"{group_key}/energy_ele_{method_label}"  # pylint: disable=invalid-name
    E = np.empty(R.shape[0], dtype=np.float64)
    E[:] = np.nan
    rfile.put(E_key, E)

    G_key = f"{group_key}/grads_{method_label}"  # pylint: disable=invalid-name
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
        "calc_acc": 0.001,
        "max_iterations": 300,
        "params": None,
    }

    saver = Saver(exdir_path_dest, (E_key, G_key))

    driver = DriverEnGrad(xtb_engrad, worker_kwargs, **driver_kwargs)

    saver.save(E, G)

    # pylint: disable-next=invalid-name
    E_ref = np.array(
        [
            -3.656060950664,
            -3.324654135118,
            -3.732034784342,
            -3.796356937974,
            -3.409646976024,
        ]
    )
    # pylint: disable-next=invalid-name
    G_ref = np.array(
        [
            [
                [6.672041620663, 0.729849613603, -3.668229063218],
                [-0.636198759233, -1.551648795581, 5.581724697274],
                [-6.035842861431, 0.821799181978, -1.913495634056],
            ],
            [
                [4.252787805953, 7.416829680566, 0.762905027301],
                [-4.892364467061, -3.28486681912, 5.007384223041],
                [0.639576661108, -4.131962861446, -5.770289250341],
            ],
            [
                [3.950069236123, 2.18864836491, 4.387309005897],
                [-5.592459441233, -3.485759723786, 0.279153980821],
                [1.64239020511, 1.297111358875, -4.666462986718],
            ],
            [
                [3.469120306905, 4.047174550132, -2.210702433228],
                [-5.865794145448, 0.407529699691, 1.769141641574],
                [2.396673838544, -4.454704249823, 0.441560791654],
            ],
            [
                [-2.720951568168, 6.590718234792, -5.571701096746],
                [5.026376064981, -5.372855856569, -0.665128770194],
                [-2.305424496813, -1.217862378223, 6.23682986694],
            ],
        ]
    )
    E, G = driver.run(Z, R, E, G, saver=saver)

    assert np.allclose(E, E_ref)
    assert np.allclose(G, G_ref)
