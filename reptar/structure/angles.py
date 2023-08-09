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
import itertools
from ase import Atoms
import numpy as np
import ray
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
        Desired angle (degrees) to set.
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


class SetAngles:
    r"""Set angles of structures"""

    def __init__(
        self,
        Z,
        angle_types: Iterable[str],
        atoms: Iterable[Iterable[int]],
        masks: Iterable[Iterable[int]],
    ) -> None:
        r"""Set any number of dihedral angles of a single structure.

        Parameters
        ----------
        Z
            Atomic numbers that are shared with every structure in ``R``.
        angle_types
            Type of angle. Can be either ``angle`` or ``dihedral``.
        atoms
            Atomic indices of the four atoms that make up the dihedral angle to change.
        masks
            Describes the two subgroups to move based on the angle. This reduces the
            chances of a clash or greatly distorting the geometry.

        Returns
        -------

            Cartesian coordinates with new dihedral angles.

        Notes
        -----
        For more information, see
        https://wiki.fysik.dtu.dk/ase/ase/atoms.html#ase.Atoms.set_dihedral
        """
        self.Z = Z
        self.angle_types = angle_types
        self.atoms = atoms
        self.masks = masks

    def __call__(self, batch: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
        r"""Rotate structures on a batch basis.

        Parameters
        ----------
        batch
            Contains ``{"data": R, "data_1": angle_values}``.

        Returns
        -------

            Dictionary containing rotated structures under ``"data"`` key.
        """
        R, angle_values = batch["data"], batch["data_1"]
        if R.ndim == 2:
            R = R[None, ...]
        gen_angles = itertools.cycle(angle_values)
        for i, r in enumerate(R):
            R[i] = set_angles(
                self.Z,
                r,
                self.angle_types,
                self.atoms,
                next(gen_angles),
                self.masks,
            )
        return {"data": R}


def sample_angles(
    Z: np.ndarray,
    R: np.ndarray,
    angle_types: Iterable[str],
    atoms: Iterable[Iterable[int]],
    angles: np.ndarray,
    masks: Iterable[Iterable[int]],
    use_ray: bool = False,
    n_workers: int = 4,
    batch_size: int = 100,
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
    use_ray
        Parallelize geometry modifications using ray.

    Returns
    -------

        Coordinates of each dihedral set.
    """
    # pylint: disable=invalid-name,no-member
    if R.ndim == 2:
        R = R[None, ...]

    n_angle_sets = angles.shape[0]
    R_rotated = np.tile(R, (n_angle_sets, 1, 1))  # pylint: disable=invalid-name

    if not use_ray:
        sa = SetAngles(Z, angle_types, atoms, masks)
        R_rotated = sa({"data": R_rotated, "data_1": angles})["data"]
    else:
        ds = ray.data.from_numpy(np.array_split(R_rotated, n_workers))
        ds = ds.zip(
            ray.data.from_numpy(np.array_split(angles, n_workers))
        ).materialize()
        ds = ds.map_batches(
            SetAngles,
            batch_size=batch_size,
            batch_format="numpy",
            # pylint: disable-next=unexpected-keyword-arg
            compute=ray.data.ActorPoolStrategy(size=n_workers),
            zero_copy_batch=True,
            fn_constructor_args=(Z, angle_types, atoms, masks),
            num_cpus=1,
        )
        results = ds.materialize()
        R_rotated = ray.get(results.to_numpy_refs()[0])["data"]

    return R_rotated
