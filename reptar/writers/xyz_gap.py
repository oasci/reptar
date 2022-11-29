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


def write_xyz_gap(
    xyz_path, lattice, Z, R, E, F=None, lattice_precision=3, data_precision=10
):
    """Write an extended XYZ file in the GAP format.

    Parameters
    ----------
    xyz_path : :obj:`str`
        Path to the XYZ file to write.
    lattice : :obj:`numpy.ndarray`, ndim: ``2``
        The three cartesian lattice vectors describing the periodic cell (in
        Angstroms).

        All structures need a this lattice even structures are not supposed to
        be periodic. Just use a lattice vector that is larger than twice the
        cutoff of any potential you plan to create or use.
    Z : :obj:`numpy.ndarray`, ndim: ``1``
        Atomic numbers of all atoms in the system.
    R : :obj:`numpy.ndarray`, ndim: ``3``
        Cartesian coordinates of all structures in the same order as ``Z``.
    E : :obj:`numpy.ndarray`, ndim: ``1``
        Energy of all structures in ``R`` in units of eV.
    F : :obj:`numpy.ndarray`, ndim: ``3``, default: ``None``
        Atomic forces of all structures in ``R`` in units of eV/A.
    lattice_precision : :obj:`int`, default: ``3``
        Number of decimal points to print for lattice dimensions.
    data_precision : :obj:`int`, default: ``10``
        Number of decimal points for printing array data.

    """
    lat_str = np.array2string(
        lattice.flatten(),
        formatter={"float_kind": lambda x: f"%.{lattice_precision}f" % x},
    )
    # TODO: is pbc always T T T?
    lat_str = 'pbc="T T T" Lattice="' + lat_str[1:-1] + '"'

    prop_line = "Properties=species:S:1:pos:R:3"
    if F is not None:
        prop_line += ":forces:R:3"

    e_formatter = lambda x: f"%.{data_precision}f" % x

    n_atoms = len(Z)
    F_arr = None
    with open(xyz_path, "w", encoding="utf-8") as f:
        for i in range(len(R)):
            f.write(f"{n_atoms}\n")
            f.write(
                " ".join([f"energy={e_formatter(E[i])}", prop_line, lat_str]) + "\n"
            )
            if F is not None:
                F_arr = F[i]
            f.write(string_xyz_arrays(Z, R[i], F_arr, precision=data_precision))
