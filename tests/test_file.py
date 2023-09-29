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

"""Tests for File"""

# pylint: skip-file

import sys
import shutil
import os
import pytest
import numpy as np
from reptar import File

sys.path.append("..")

# Ensures we execute from file directory (for relative paths).
os.chdir(os.path.dirname(os.path.realpath(__file__)))


def nesting_dicts():
    nested_dict = {
        "keep": [281, "a"],
        "start": {
            "continue": {
                "keep": [9192],
                "remove": 1928,
                "remove_too": np.array([8, 1, 1], dtype=np.uint8),
            },
            "keep": 2891,
        },
    }
    rfile_dict_true = {
        "keep": np.array([281, "a"], dtype="<U21"),
        "start": {
            "continue": {
                "keep": 9192,
                "remove": 1928,
                "remove_too": np.array([8, 1, 1], dtype=np.uint8),
            },
            "keep": 2891,
        },
    }
    rfile_dict_true_deleted = {
        "keep": np.array([281, "a"], dtype="<U21"),
        "start": {
            "continue": {
                "keep": 9192,
            },
            "keep": 2891,
        },
    }
    return nested_dict, rfile_dict_true, rfile_dict_true_deleted


def test_file_del_nested_json():
    r"""Delete a key nested inside a json."""
    nested_dict, rfile_dict_true, rfile_dict_true_deleted = nesting_dicts()
    rfile = File("./tmp/test.json", mode="w", from_dict=nested_dict)
    np.testing.assert_equal(rfile.File_, rfile_dict_true)

    rfile.remove("/start/continue/remove/")
    rfile.remove("/start/continue/remove_too")
    np.testing.assert_equal(rfile.File_, rfile_dict_true_deleted)


def test_file_del_nested_exdir():
    r"""Delete a key nested inside an exdir file."""
    nested_dict, _, _ = nesting_dicts()

    # Remove exdir file if it already exists
    test_exdir_path = "./tmp/test.exdir"
    if os.path.exists(test_exdir_path):
        shutil.rmtree(test_exdir_path)
    rfile = File(test_exdir_path, mode="w", from_dict=nested_dict)
    assert rfile.get("start/continue/remove", missing_is_none=True) == 1928

    rfile.remove("/start/continue/remove_too/")
    rfile.remove("/start/continue/remove")
    assert rfile.get("start/continue/remove", missing_is_none=True) is None
    assert rfile.get("start/continue/remove_too", missing_is_none=True) is None


def test_file_del_nested_zarr():
    r"""Delete a key nested inside an zarr file."""
    nested_dict, _, _ = nesting_dicts()

    # Remove exdir file if it already exists
    test_zarr_path = "./tmp/test.zarr"
    if os.path.exists(test_zarr_path):
        shutil.rmtree(test_zarr_path)
    rfile = File(test_zarr_path, mode="w", from_dict=nested_dict)
    assert rfile.get("start/continue/remove", missing_is_none=True) == 1928

    rfile.remove("/start/continue/remove_too/")
    rfile.remove("/start/continue/remove")
    assert rfile.get("start/continue/remove", missing_is_none=True) is None
    assert rfile.get("start/continue/remove_too", missing_is_none=True) is None
