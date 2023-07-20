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

import numpy as np


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
    n_x = int(n_x)
    n_y = int(n_y)
    n_z = int(n_z)
    n_points = n_x * n_y * n_z
    coords = np.empty((n_points, 3), dtype=np.float64)
    x_coords = np.linspace(
        origin_x, origin_x + (n_x * spacing_x), num=n_x, endpoint=False
    )
    y_coords = np.linspace(
        origin_y, origin_y + (n_y * spacing_y), num=n_y, endpoint=False
    )
    z_coords = np.linspace(
        origin_z, origin_z + (n_z * spacing_z), num=n_z, endpoint=False
    )
    x_coords, y_coords, z_coords = np.meshgrid(
        x_coords, y_coords, z_coords, indexing="ij"
    )
    coords[:, 0] = x_coords.flatten()
    coords[:, 1] = y_coords.flatten()
    coords[:, 2] = z_coords.flatten()
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
    with open(file_path, "r", encoding="utf-8") as f_cube:
        f_cube.readline()
        f_cube.readline()
        n_atoms, cube_R = _get_cube_coords(f_cube)  # pylint: disable=invalid-name
        # Skip over atom positions.
        for _ in range(n_atoms):
            f_cube.readline()
        cube_V = np.zeros((cube_R.shape[0]), dtype=np.float64)
        idx = 0
        for line in f_cube:
            for val in line.strip().split():
                cube_V[idx] = float(val)
                idx += 1
    return cube_R, cube_V
