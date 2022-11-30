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

"""Tests writing PDB files from group"""

# Stuff for PyTest features like skip.
# pylint: disable=import-outside-toplevel, unused-import, duplicate-code

import sys
import os
import pytest
import numpy as np
from reptar import File
from reptar.writers import write_pdb

sys.path.append("..")

# Ensures we execute from file directory (for relative paths).
os.chdir(os.path.dirname(os.path.realpath(__file__)))

# Source paths
DATA_DIR = "./data/"

# Writing paths
WRITING_DIR = "./tmp/writing/"
os.makedirs(WRITING_DIR, exist_ok=True)


def test_pdb_writer_1h2o_120meoh_prod():
    r"""Writing short PDB file from exdir file"""
    exdir_path = os.path.join(DATA_DIR, "1h2o_120meoh_md.exdir")
    pdb_path = os.path.join(WRITING_DIR, "1h2o_120meoh_md_eq_1.pdb")

    rfile = File(exdir_path, mode="r")

    Z = rfile.get("eq_1/atomic_numbers")
    R = rfile.get("eq_1/geometry")[:5]
    entity_ids = rfile.get("eq_1/entity_ids")
    comp_ids = rfile.get("eq_1/comp_ids")
    write_pdb(pdb_path, Z, R, entity_ids, comp_ids)

    # TODO: Write tests
