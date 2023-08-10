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

"""Interact with bonding of structures."""

from __future__ import annotations
from collections.abc import Iterable
from mendeleev import element
from scipy.spatial.distance import cdist
import numpy as np
from ray.util.multiprocessing import Pool
from ..logger import ReptarLogger

log = ReptarLogger(__name__)


def get_atomic_radii(Z: Iterable[int]) -> np.ndarray:
    """Obtain atomic radius using
    `mendeleev <https://mendeleev.readthedocs.io/en/stable/index.html>`__.

    Parameters
    ----------
    Z
        Atomic number

    Returns
    -------

        Atomic radii in Angstroms.
    """
    return np.array([element(int(z)).atomic_radius * 0.01 for z in Z], dtype=np.float64)


def is_bonded_atomic_radii(
    Z: np.ndarray, R: np.ndarray, a0_scaling: float = 1.2, use_ray: bool = False
) -> np.ndarray:
    r"""Two atoms are defined as bonded if the distance between atoms, :math:`d`, is
    less than :math:`c (r_0 + r_1)` where :math:`r_0` and :math:`r_1` are the atomic
    radii of the two atoms and :math:`c` is the chosen constant ``a0_scaling``.

    Parameters
    ----------
    Z
        Atomic numbers.
    R
        Cartesian coordinates.
    a0_scaling
        Scaling factor for atomic radii.
    use_ray
        Parallelize computations. Usually ``True`` is only faster when you have
        10,000+ structures.

    Returns
    -------

        Adjacency matrices of many structures with ``True`` indicated that the atoms
        are bonded.
    """
    atomic_radii = np.tile(get_atomic_radii(Z), (Z.shape[0], 1))
    sum_pair_atomic_radii = atomic_radii + atomic_radii.T

    def _cdist_self(a):
        if np.isnan(a).any():
            n_atoms = a.shape[0]
            return np.full((n_atoms, n_atoms), np.nan, dtype=np.float64)
        return cdist(a, a, metric="euclidean")

    if use_ray:
        with Pool() as p:
            pairwise_dist = np.array(p.map(_cdist_self, R), dtype=np.float64)
    else:
        pairwise_dist = np.empty((R.shape[0], R.shape[1], R.shape[1]))
        for i, r in enumerate(R):
            pairwise_dist[i] = _cdist_self(r)

    return pairwise_dist < a0_scaling * sum_pair_atomic_radii


def identify_bonds(
    Z: np.ndarray, R: np.ndarray, f_is_bonded: callable, f_kwargs: dict
) -> np.ndarray:
    r"""Identify covalent bonds of a structure from atomic coordinates.

    Parameters
    ----------
    Z
        Atomic numbers.
    R
        Cartesian coordinates where ``ndim`` is ``3``.
    f_is_bonded
        A function that computes an adjacency matrix for multiple structures. See
        :func:`~reptar.structure.bond.is_bonded_atomic_radii` for an example.
    f_kwargs
        Keyword arguments for ``f_is_bonded``.

    Returns
    -------

        Array of atom indices (column ``0`` and ``1``) and bond lengths (column ``2``).
    """
    if R.ndim == 2:
        R = R[None, ...]
    is_bonded = f_is_bonded(Z=Z, R=R, **f_kwargs)
    if is_bonded.shape[0] == 1:
        is_bonded = is_bonded[0]
    return is_bonded
