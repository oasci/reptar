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

"""Interact with dihedral angles of structures."""

from __future__ import annotations
from collections.abc import Iterable
from ase import Atoms
import numpy as np
from ..logger import ReptarLogger

log = ReptarLogger(__name__)


def get_angle_mask(n_atoms: int, fragment_indices: Iterable[int]):
    """Specify boolean mask for which atoms to rotate.

    Parameters
    ----------
    n_atoms
        Total number of atoms in the structure.
    fragment_indices
        Atom indices of one fragment involved in rotating the dihedral.

    Returns
    -------
    :obj:`numpy.ndarray`
        Mask of which atoms to rotate for the dihedral. ``0`` for cap and ``1`` for
        all other atoms. We keep the cap stationary so that the likely larger section
        of the residue is rotated to minimize chance of clashes.
    """
    mask = np.ones(n_atoms, dtype=np.uint64)
    mask[np.array(fragment_indices)] = 0
    return mask


def set_angles(
    Z: np.ndarray,
    R: np.ndarray,
    angle_types: Iterable[str],
    atoms: Iterable[Iterable[int]],
    angles: Iterable[Iterable[float]],
    masks: Iterable[Iterable[int]] = None,
) -> np.ndarray:
    """Set any number of dihedral angles of a single structure.

    Parameters
    ----------
    Z
        Atomic numbers that are shared with every structure in ``R``.
    R
        Cartesian coordinates of one or more structures.
    atoms
        Atomic indices of the four atoms that make up the dihedral angle to change.
    angles
        Desired angle (degrees) to set the dihedral to.
    masks
        Describes the two subgroups to move based on the angle. This reduces the chances
        of a clash or greatly distorting the geometry.

    Returns
    -------

        Cartesian coordinates with new dihedral angles.

    Notes
    -----
    For more information, see
    https://wiki.fysik.dtu.dk/ase/ase/atoms.html#ase.Atoms.set_dihedral
    """
    if R.ndim != 2:
        raise ValueError(f"R must be two dimensional, but contains {R.ndim}")
    ase_atoms = Atoms(Z, positions=R)
    for angle_type, atom_indices, angle, mask in zip(angle_types, atoms, angles, masks):
        angle_type = angle_type.strip().lower()
        if angle_type == "dihedral":
            ase_atoms.set_dihedral(*atom_indices, angle, mask=mask)
        elif angle_type == "angle":
            ase_atoms.set_angle(*atom_indices, angle, mask=mask)
        else:
            raise ValueError(
                f"angle_type must be 'dihedral' or 'angle' but is {angle_type}"
            )
    return ase_atoms.get_positions()


def sample_angles(
    Z: np.ndarray,
    R: np.ndarray,
    angle_types: Iterable[str],
    atoms: Iterable[Iterable[int]],
    angles: Iterable[Iterable[float]],
    masks: Iterable[Iterable[int]] = None,
) -> np.ndarray:
    r"""Generate structures with different sets of dihedral angels. No optimizations or
    energy calculations are performed.

    Parameters
    ----------
    Z
        Atomic numbers.
    R
        Cartesian coordinates of structures. If multiple are given, then the angles
        are applied to each structure.
    angle_types
        Type of angle. Can be either ``angle`` or ``dihedral``.
    atoms
        Atom indices in continuous order along the dihedral. Each row defines the four
        atoms that make up that dihedral.
    angles
        Angles to set for all dihedrals. Columns correspond to a single dihedral where
        each row is a unique set of angles.
    masks
        Describes the two subgroups to move based on the angle. This reduces the chances
        of a clash or greatly distorting the geometry.

    Returns
    -------

        Coordinates of each dihedral set.
    """
    if R.ndim == 2:
        R = R[None, ...]

    n_angle_sets = len(angles)
    R_rotated = np.tile(R, (n_angle_sets, 1, 1))  # pylint: disable=invalid-name
    i_angle = 0
    for i, r in enumerate(R_rotated):
        # Reset angle counter if we are setting angles to multiple structures.
        if i_angle == n_angle_sets:
            i_angle = 0
        R_rotated[i] = set_angles(Z, r, angle_types, atoms, angles[i_angle], masks)
        i_angle += 1
    return R_rotated
