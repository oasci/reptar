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


def _get_cube_coords(f_obj):
    r"""Reads header information to generate Cartesian coordinates of cube data.

    Parameters
    ----------
    f_obj
        File object

    Returns
    -------
    :obj:`numpy.ndarray`
        Cartesian coordinates of points in space.
    """
    _, origin_x, origin_y, origin_z = f_obj.readline().strip().split()
    n_x, spacing_x, _, _ = f_obj.readline().strip().split()
    n_y, _, spacing_y, _ = f_obj.readline().strip().split()
    n_z, _, _, spacing_z = f_obj.readline().strip().split()
    coords = np.empty((n_x, n_y, n_z), dtype=np.float64)
    coords[:, 0] = np.linspace(origin_x, spacing_x, num=n_x)
    coords[:, 1] = np.linspace(origin_y, spacing_y, num=n_y)
    coords[:, 2] = np.linspace(origin_z, spacing_z, num=n_z)
    return coords


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
        cube_R = _get_cube_coords(f_cube)
        n_x, n_y, n_z = cube_R.shape
        data = np.zeros((n_x * n_y * n_z))
        idx = 0
        for line in f_cube:
            for val in line.strip().split():
                data[idx] = float(val)
                idx += 1
    data = np.reshape(data, cube_R.shape)
    return cube_R, data
