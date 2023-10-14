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

from __future__ import annotations

import numpy as np


def get_R_span(R: np.ndarray) -> np.ndarray:  # pylint: disable=invalid-name
    """Span of ``x``, `y``, and ``z`` coordinates of all structures.

    Parameters
    ----------
    R
        Cartesian coordinates of all structures to consider.

    Returns
    -------
    :obj:`numpy.ndarray`
        Span of structure coordinates.
    """
    return np.subtract(np.max(R, axis=-2), np.min(R, axis=-2))


# pylint: disable=invalid-name
def get_n_grid_points(
    spacing: np.ndarray, R_span: np.ndarray, overage: np.ndarray
) -> np.ndarray:
    """Compute the number of cubic grid points in the ``x``, ``y``, ``z`` directions.

    Parameters
    ----------
    spacing
        Spacing between grid points for all ``x``, ``y``, and ``z`` directions.
    R_span
        The span of ``x``, ``y``, and ``z`` coordinates for one or many structures.
    overage
        Additional space in the minimum or maximum ``x``, ``y``, and ``z`` directions.
        Must be in the same units as ``R_span``.

    Returns
    -------
    :obj:`numpy.ndarray`
        Total number of points in the cubic grid.

    Notes
    -----
    This function is based off of `Psi4's grid generation code <https://github.com/
    psi4/psi4/blob/2cd33eda01b7018a23739d00c1cdd51ca87faa64/psi4/src/psi4/libcubeprop
    /csg.cc#L94-L223>`__.
    """
    if not isinstance(overage, np.ndarray):
        overage = np.array(overage)
    if not isinstance(spacing, np.ndarray):
        spacing = np.array(spacing)

    cube_span = R_span + 2.0 * overage
    n_points = np.divide(cube_span, spacing).astype(np.uint32)
    # Adjustment on the next line is from here
    # https://github.com/psi4/psi4/blob/2cd33eda01b7018a23739d00c1cdd51ca87faa64/psi4/src/psi4/libcubeprop/csg.cc#L133
    n_points += np.array(spacing * n_points < cube_span)
    # Adjustment on the next line is from here
    # https://github.com/psi4/psi4/blob/2cd33eda01b7018a23739d00c1cdd51ca87faa64/psi4/src/psi4/libcubeprop/csg.cc#L208-L210
    n_points += 1
    return n_points


# pylint: disable-next=invalid-name
def get_total_grid_points(
    spacing: np.ndarray, R_span: np.ndarray, overage: np.ndarray
) -> np.ndarray:
    """Compute the total number of cubic grid points.

    Parameters
    ----------
    spacing
        Spacing between grid points for all ``x``, ``y``, and ``z`` directions.
    R_span
        The span of ``x``, ``y``, and ``z`` coordinates for one or many structures.
    overage
        Additional space in the minimum or maximum ``x``, ``y``, and ``z`` directions.
        Must be in the same units as ``R_span``.

    Returns
    -------
    :obj:`numpy.ndarray`
        Total number of points in the cubic grid.
    """
    return np.prod(get_n_grid_points(spacing, R_span, overage), axis=-1)


# pylint: disable-next=invalid-name
def get_grid_origin(
    R_min: np.ndarray,
    R_max: np.ndarray,
    spacing: np.ndarray | float,
    n_points: np.ndarray | float,
) -> np.ndarray:
    """Determine origin of a cubic grid.

    Parameters
    ----------
    R_min : :obj:`numpy.ndarray`
        Minimum atomic coordinate in ``x``, ``y``, and ``z`` axis for each structure.
    R_max : :obj:`numpy.ndarray`
        Maximum atomic coordinate in ``x``, ``y``, and ``z`` axis for each structure.
    spacing : :obj:`float`
        Space between grid points for all ``x``, ``y``, and ``z`` directions.
    """
    return R_min - (spacing * (n_points - 1) - (R_max - R_min)) / 2.0


def get_grid_points(
    origin: np.ndarray, n_points: np.ndarray, spacing: np.ndarray | float
) -> np.ndarray:
    """Generate grid coordinates from coordinate vectors.

    Parameters
    ----------
    origin
        Origin of the cubic grid (with the minimum coordinates).
    n_points
        Number of points in the ``x``, ``y``, and ``z`` direction.
    spacing
        Space between grid points for all ``x``, ``y``, and ``z`` directions.

    Returns
    -------
    :obj:`numpy.ndarray`
        Cartesian coordinates of all cubic grid points.
    """
    if not isinstance(origin, np.ndarray):
        origin = np.array(origin)
    if not isinstance(n_points, np.ndarray):
        n_points = np.array(n_points)
    if not isinstance(spacing, np.ndarray):
        spacing = np.array(spacing)

    assert origin.ndim == 1
    assert n_points.ndim == 1
    assert spacing.ndim == 1

    total_n_points = np.prod(n_points)

    point_coords = np.empty((total_n_points, 3), dtype=np.float64)
    x_coords = np.linspace(
        origin[0],
        origin[0] + (n_points[0] * spacing[0]),
        num=n_points[0],
        endpoint=False,
    )
    y_coords = np.linspace(
        origin[1],
        origin[1] + (n_points[1] * spacing[1]),
        num=n_points[1],
        endpoint=False,
    )
    z_coords = np.linspace(
        origin[2],
        origin[2] + (n_points[2] * spacing[2]),
        num=n_points[2],
        endpoint=False,
    )
    x_coords, y_coords, z_coords = np.meshgrid(
        x_coords, y_coords, z_coords, indexing="ij"
    )
    point_coords[:, 0] = x_coords.flatten()
    point_coords[:, 1] = y_coords.flatten()
    point_coords[:, 2] = z_coords.flatten()
    return point_coords


# pylint: disable-next=invalid-name
def get_max_grid_points(
    R_span: np.ndarray, overage: np.ndarray | float, spacing: np.ndarray | float
) -> np.ndarray:
    """Determine cubic grid spacing for a set of structures that will keep the number
    of points the same.

    Parameters
    ----------
    R_span
        Span of structure coordinates.
    overage
        Additional space in the minimum and maximum ``x``, ``y``, and ``z`` directions.
        Must be in the same units as ``R``.
    spacing
        The space between grid points in ``x``, ``y``, and ``z`` directions. If only
        :obj:`float` is provided, we broadcast it in all dimensions.
    """
    total_points = get_total_grid_points(spacing, R_span, overage)
    return np.max(total_points)


def initialize_grid_arrays(R, overage=None, spacing=None, max_points=None):
    """Create arrays for grid points and values based on the number of points of the
    largest grid.

    You must specify either ``overage`` and ``spacing`` or ``max_points``.

    Parameters
    ----------
    R : :obj:`numpy.ndarray`
        Cartesian coordinates of all structures to consider.
    overage : :obj:`numpy.ndarray`, ndim: ``1``
        Additional space in the minimum and maximum ``x``, ``y``, and ``z`` directions.
        Must be in the same units as ``R``.
    spacing : :obj:`float` or :obj:`numpy.ndarray`
        The space between grid points in ``x``, ``y``, and ``z`` directions. If only
        :obj:`float` is provided, we broadcast it in all dimensions. Must be in the
        same units as ``R``.
    max_points : :obj:`int`, default: ``None``
        Manually specify the number of maximum points to use.

    Returns
    -------
    :obj:`numpy.ndarray`
        Cartesian coordinates of points where a property is probed.
    :obj:`numpy.ndarray`
        Property values at the same Cartesian coordinates.
    """
    # pylint: disable=invalid-name
    if R.ndim == 2:
        R = R[None, ...]
    if max_points is None:
        R_span = get_R_span(R)
        max_points = get_max_grid_points(R_span, overage, spacing)
    cube_R = np.full((R.shape[0], max_points, 3), np.nan, dtype=np.float64)
    cube_V = np.full((R.shape[0], max_points), np.nan, dtype=np.float64)
    return cube_R, cube_V
