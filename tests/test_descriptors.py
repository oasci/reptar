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

"""Tests computing structural descriptors"""

import pytest
import os
import numpy as np
from reptar import File
from reptar import descriptors

import sys

sys.path.append("..")
from .paths import *

# Ensures we execute from file directory (for relative paths).
os.chdir(os.path.dirname(os.path.realpath(__file__)))


def test_max_atom_pair_dist_3h2o():
    """Computing 2-body energy and gradient contributions."""
    exdir_path = get_140h2o_samples_path()
    exdir_file = File(exdir_path, mode="r")

    group_key = "3h2o"

    i_test = 482
    Z = exdir_file.get(f"{group_key}/atomic_numbers")
    R = exdir_file.get(f"{group_key}/geometry")

    criteria = descriptors.Criteria(descriptors.max_atom_pair_dist, {}, 9.0)
    accept, v_desc = criteria.accept(Z, R)

    assert accept[i_test] == True
    assert v_desc[i_test] == 8.96801858001878


def test_com_distance_sum_3h2o():
    """Computing 2-body energy and gradient contributions."""
    exdir_path = get_140h2o_samples_path()
    exdir_file = File(exdir_path, mode="r")

    group_key = "3h2o"

    i_test = 482
    Z = exdir_file.get(f"{group_key}/atomic_numbers")
    R = exdir_file.get(f"{group_key}/geometry")
    entity_ids = exdir_file.get(f"{group_key}/entity_ids")

    criteria = descriptors.Criteria(
        descriptors.com_distance_sum, {"entity_ids": entity_ids}, 9.0
    )
    accept, v_desc = criteria.accept(Z, R)

    assert accept[i_test] == False
    assert v_desc[i_test] == 9.020819657050316
