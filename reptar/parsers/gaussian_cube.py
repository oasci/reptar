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

from ..calculators.cube import get_grid_points


def _get_cube_line(line):
    return (float(i) for i in line.strip().split())


def _get_cube_coords(f_obj):
    r"""Reads header information to generate Cartesian coordinates of cube data.

    Parameters
    ----------
    f_obj
        File object

    Returns
    -------
    :obj:`int`
        Total number of atoms.
    :obj:`numpy.ndarray`
        Cartesian coordinates of points in space.
    """
    n_atoms, origin_x, origin_y, origin_z = _get_cube_line(f_obj.readline())
    n_x, spacing_x, _, _ = _get_cube_line(f_obj.readline())
    n_y, _, spacing_y, _ = _get_cube_line(f_obj.readline())
    n_z, _, _, spacing_z = _get_cube_line(f_obj.readline())

    origin = np.array([origin_x, origin_y, origin_z], dtype=np.float64)
    n_points = np.array([int(n_x), int(n_y), int(n_z)], dtype=np.uint32)
    spacing = np.array([spacing_x, spacing_y, spacing_z], dtype=np.float64)

    coords = get_grid_points(origin, n_points, spacing)
    return int(n_atoms), coords


def parse_cube(file_path):
    r"""Parse Gaussian cube file into NumPy arrays.

    Parameters
    ----------
    file_path : :obj:`str`
        Path to cube file.

    Returns
    -------
    :obj:`numpy.ndarray`
        Cartesian coordinates of points where a property is probed.
    :obj:`numpy.ndarray`
        Property values at the same Cartesian coordinates.
    """
    with open(file_path, encoding="utf-8") as f_cube:
        f_cube.readline()
        f_cube.readline()
        n_atoms, cube_R = _get_cube_coords(f_cube)  # pylint: disable=invalid-name
        # Skip over atom positions.
        for _ in range(n_atoms):
            f_cube.readline()
        # pylint: disable-next=invalid-name
        cube_V = np.zeros((cube_R.shape[0]), dtype=np.float64)
        idx = 0
        for line in f_cube:
            for val in line.strip().split():
                cube_V[idx] = float(val)
                idx += 1
    return cube_R, cube_V
