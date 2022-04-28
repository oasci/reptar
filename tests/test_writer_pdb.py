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
import numpy as np
from reptar import creator
from reptar.writers import pdbWriter

import sys
sys.path.append("..")
from .paths import *

# Ensures we execute from file directory (for relative paths).
os.chdir(os.path.dirname(os.path.realpath(__file__)))


def test_pdb_writer_1h2o_120meoh_prod():
    """Writing short PDB file from exdir file.
    """
    exdir_path = '../examples/xtb/v6.4.1/md/1h2o.120meoh.pm.gfn2.exdir'

    create = creator()
    create.load(exdir_path, mode='r')

    writer = pdbWriter()
    writer.write(
        create.data, '/prod_1', file_name='prod_1',
        save_dir='./tmp', R_limits=(None, 5)
    )

    # TODO: Write any tests?

    os.remove('./tmp/prod_1.pdb')


