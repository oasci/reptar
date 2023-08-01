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

"""Tests for cube utilities"""

# pylint: skip-file

import sys
import os
import pytest
import numpy as np
import qcelemental as qcel
from reptar import File
from reptar.calculators.cube import initialize_grid_arrays



sys.path.append("..")
from .paths import get_140h2o_samples_path

# Ensures we execute from file directory (for relative paths).
os.chdir(os.path.dirname(os.path.realpath(__file__)))

# Writing paths
WRITING_DIR = "./tmp/calculators/"
os.makedirs(WRITING_DIR, exist_ok=True)


def test_ray_cube_max_points():
    exdir_path_source = get_140h2o_samples_path()
    rfile_source = File(exdir_path_source, mode="r")

    # Copy over a few structures for calculations.
    group_key = "1h2o"
    start_slice = None
    end_slice = 3

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
        ]
    )
    R = rfile_source.get(f"{group_key}/geometry")[start_slice:end_slice]
    assert np.allclose(R, R_ref)

    spacing = np.array([0.2, 0.2, 0.2], dtype=np.float64)  # Bohr
    overage = np.array([4.0, 4.0, 4.0], dtype=np.float64)  # Bohr
    grid_R, grid_V = initialize_grid_arrays(
        R / qcel.constants.bohr2angstroms, overage, spacing
    )
    assert grid_R.shape == (3, 129948, 3)
    assert grid_V.shape == (3, 129948)
