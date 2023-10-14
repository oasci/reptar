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

"""Tests for reptar-dihedral-scan script"""

# pylint: skip-file

import os
import shutil
import subprocess
import sys

import numpy as np
import pytest
import qcelemental as qcel

from reptar import File
from reptar.structure.scan import create_grid

sys.path.append("..")

# Ensures we execute from file directory (for relative paths).
os.chdir(os.path.dirname(os.path.realpath(__file__)))

# Writing paths
WRITING_DIR = "./tmp/calculators/"
os.makedirs(WRITING_DIR, exist_ok=True)


def test_script_gfp_cro_geom_scan_gen():
    R = np.array(
        [
            [-2.59295103e00, -1.94121745e00, 8.36369574e-03],
            [-1.48845914e00, -1.85186140e00, 9.44147787e-01],
            [-4.89849235e-01, -3.03550749e00, 7.67343803e-01],
            [7.46438026e-04, -3.15170571e00, -6.75484274e-01],
            [5.65624393e-01, -2.91874754e00, 1.67914251e00],
            [-7.99963347e-01, -5.27998697e-01, 8.25626215e-01],
            [4.61118601e-01, -3.74641488e-01, 5.85088394e-01],
            [-1.45353553e00, 6.55528233e-01, 1.08958984e00],
            [-5.46945494e-01, 1.68168514e00, 8.85358976e-01],
            [-8.20335842e-01, 2.86764396e00, 9.41999289e-01],
            [7.14110524e-01, 9.87997019e-01, 5.87110785e-01],
            [-2.82237937e00, 8.75897051e-01, 1.47334554e00],
            [-3.73690252e00, 1.13094286e00, 2.64604161e-01],
            [-4.70211380e00, 4.27802776e-01, 3.04580150e-02],
            [1.88302142e00, 1.63102318e00, 3.51698517e-01],
            [3.16387715e00, 1.06859449e00, 3.67808251e-02],
            [3.35438537e00, -3.04475191e-01, -1.63061881e-01],
            [4.27570003e00, 1.91613875e00, -8.24700605e-02],
            [4.60116934e00, -8.06814681e-01, -4.60383379e-01],
            [5.52290574e00, 1.41741871e00, -3.76143425e-01],
            [5.69826434e00, 4.63073120e-02, -5.66636900e-01],
            [6.94110919e00, -4.07370060e-01, -8.53874581e-01],
            [-4.69652245e00, -1.18215632e00, -1.60329717e00],
            [-3.31104389e00, -2.61741929e00, 2.10994563e-01],
            [6.92882903e00, -1.36673353e00, -9.61497180e-01],
            [-1.89807405e00, -1.91622964e00, 1.96361587e00],
            [-1.01338207e00, -3.95868097e00, 1.04939132e00],
            [7.86479111e-01, -3.90169629e00, -7.21939197e-01],
            [-8.09062001e-01, -3.45285552e00, -1.33418887e00],
            [3.92511013e-01, -2.19718286e00, -1.01470154e00],
            [1.04314316e00, -2.10460179e00, 1.44894793e00],
            [-3.20170334e00, 1.56363030e-02, 2.02767098e00],
            [-2.83977198e00, 1.76413456e00, 2.11506881e00],
            [-3.65786270e00, -1.24034015e00, -3.03643767e00],
            [1.83048128e00, 2.71126707e00, 4.16792970e-01],
            [2.50430342e00, -9.61497792e-01, -8.99433241e-02],
            [4.14813080e00, 2.97855469e00, 6.39151161e-02],
            [4.73498415e00, -1.86978528e00, -6.13996584e-01],
            [6.38089941e00, 2.06452960e00, -4.65641353e-01],
            [-2.62900722e00, -1.27800676e00, -1.17256863e00],
            [-3.80914149e00, -1.63938902e00, -2.03795016e00],
            [-1.79552975e00, -4.48740078e-01, -1.47581210e00],
            [-3.94750376e00, -2.71639583e00, -2.08590155e00],
            [-3.37267388e00, 2.19911754e00, -4.49874047e-01],
            [-2.54399485e00, 2.71442767e00, -1.76711878e-01],
            [-4.01374312e00, 2.56461462e00, -1.68795141e00],
            [-4.26264018e00, 3.62595486e00, -1.68154464e00],
            [-3.35792282e00, 2.35188886e00, -2.53470205e00],
            [-4.92439255e00, 1.97526101e00, -1.78180512e00],
        ]
    )
    R = np.tile(R, (8, 1, 1))
    Z = np.array(
        [
            7,
            6,
            6,
            6,
            8,
            6,
            7,
            7,
            6,
            8,
            6,
            6,
            6,
            8,
            6,
            6,
            6,
            6,
            6,
            6,
            6,
            8,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            1,
            6,
            6,
            8,
            1,
            7,
            1,
            6,
            1,
            1,
            1,
        ],
        dtype=np.uint8,
    )

    scan_info = {
        "angle": [
            {
                "atoms": [30, 4, 2],
                "fragment": [30],
                "range": {"start": 100, "stop": 200, "step": 100},
            }
        ],
        "dihedral": [
            {
                "atoms": [39, 0, 1, 5],
                "fragment": [39, 41, 40, 22, 33, 42],
                "range": {"start": -180, "stop": 0, "step": 180},
            }
        ],
    }

    R_rotated = create_grid(Z, R, scan_info, use_ray=True, n_workers=2)
    assert R_rotated.shape == (32, 49, 3)
