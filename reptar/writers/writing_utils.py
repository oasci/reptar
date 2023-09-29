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

import numpy as np
from qcelemental import periodictable as ptable


def string_xyz_arrays(Z, R, *args, precision=10):
    r"""Create string of array data in XYZ format for a single structure.

    Parameters
    ----------
    Z : :obj:`numpy.ndarray`, ndim=1
        Atomic numbers of all atoms in the system.
    R : :obj:`numpy.ndarray`, ndim=2
        Cartesian coordinates of all atoms in the same order as ``Z``.
    args
        Other :obj:`numpy.ndarray` (ndim >= 1) to add where it's assumed the
        zero axis is with respect to ``R``. For example, if we have atomic
        forces the array shape would be ``(n, 3)`` where ``n`` is the number of
        atoms in the structure.
    precision : :obj:`int`, default: ``10``
        Number of decimal points for printing array data.

    Returns
    -------
    :obj:`str`
        The XYZ string for a single structure. This does not include the number
        of atoms
    """
    struct_string = ""
    for i, z in enumerate(Z):
        atom_string = ptable.to_symbol(z)
        for arr in (R, *args):
            if arr is not None:
                atom_string += "    "
                atom_string += np.array2string(
                    arr[i],
                    suppress_small=True,
                    separator="    ",
                    formatter={"float_kind": lambda x: f"%.{precision}f" % x},
                )[1:-1]
        atom_string = atom_string.replace(" -", "-")
        atom_string += "\n"
        struct_string += atom_string
    return struct_string
