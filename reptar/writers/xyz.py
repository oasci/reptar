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

import numpy as np
from .writing_utils import string_xyz_arrays
    
def write_xyz(
    xyz_path, Z, R, comments=None, data_precision=10
):
    """Write standard XYZ file.
    
    Parameters
    ----------
    xyz_path : :obj:`str`
        Path to XYZ file to write.
    Z : :obj:`numpy.ndarray`, ndim: ``1``
        Atomic numbers of all atoms in the system.
    R : :obj:`numpy.ndarray`, ndim: ``3``
        Cartesian coordinates of all structures in the same order as ``Z``.
    comments : :obj:`list`, default: ``None``
        Comment lines for each XYZ structure.
    data_precision : :obj:`int`, default: ``10``
        Number of decimal points for printing array data.

    """
    n_atoms = len(Z)
    with open(xyz_path, 'w') as f:
        for i in range(len(R)):
            f.write(f'{n_atoms}\n')
            if comments is not None:
                comment = comments[i]
                if comment[-2:] != '\n':
                    comment += '\n'
            else:
                comment = '\n'
            f.write(comment)
            f.write(string_xyz_arrays(Z, R[i], precision=data_precision))
