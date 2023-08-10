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
from typing import Any
from collections.abc import Callable
import numpy as np
import ray
from .bond import identify_bonds
from ..logger import ReptarLogger

log = ReptarLogger(__name__)


class TopologyCheck:
    r"""Check changes in topology with respect to some reference."""

    def __init__(
        self,
        Z: np.ndarray,
        R_ref: np.ndarray,  # pylint: disable=invalid-name
        f_is_bonded: Callable,
        f_kwargs: dict[str, Any] | None = None,
    ) -> None:
        r"""Set any number of dihedral angles of a single structure.

        Parameters
        ----------
        Z
            Atomic numbers that are shared with every structure in ``R``.
        R_ref
            Reference Cartesian coordinates to compare bonding to.
        f_is_bonded
            Function that, with ``Z``, ``R``, and ``f_kwargs``, will compute a
            :obj:`bool` adjacency matrix for multiple structures.
        f_kwargs
            Keyword arguments for ``f_is_bonded``.

        Example
        -------
        >>> topo_change = TopologyCheck(Z, R_ref, is_bonded_atomic_radii)
        >>> did_topo_change = topo_change({"data": R})["data"]
        array([False, False, True, False])
        """
        self.Z = Z
        self.f_is_bonded = f_is_bonded
        if f_kwargs is None:
            f_kwargs = {}
        self.f_kwargs = f_kwargs
        self.bonding_ref = self._get_bonding(R_ref)

    def _get_bonding(self, R):
        r"""Wrapper around :func:`~reptar.structure.bond.identify_bonds`."""
        return identify_bonds(self.Z, R, self.f_is_bonded, self.f_kwargs)

    def _compare_bonding(self, bonding: np.ndarray) -> np.ndarray:
        r"""Compare bonding to reference.

        Parameters
        ----------
        bonding
            Adjacency matrices for all structures to compare.

        Returns
        -------

            If bonding has changed with respect to the reference. ``True`` means that
            covalent bonding has changed.
        """
        return ~np.all(bonding == self.bonding_ref, axis=(2, 1))

    def __call__(self, batch: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
        r"""Check if topology has changed.

        Parameters
        ----------
        batch
            Contains ``{"data": R}``.

        Returns
        -------

            Dictionary containing ``data`` key.
        """
        R = batch["data"]
        if R.ndim == 2:
            R = R[None, ...]
        bonding = self._get_bonding(R)
        return {"data": self._compare_bonding(bonding)}


def get_topology_change(
    Z: np.ndarray,
    R_ref: np.ndarray,  # pylint: disable=invalid-name
    R: np.ndarray,
    f_is_bonded: Callable,
    f_kwargs: dict[str, Any],
    use_ray: bool = False,
    n_workers: int = 2,
) -> np.ndarray:
    r"""Checks is the topology of a system changes based on changes in bonding.

    Parameters
    ----------
    Z
        Atomic numbers
    R_ref
        Reference structure to compare bonding to.
    R
        Structures to check if the topology has changed.
    f_is_bonded
        A function that computes an adjacency matrix for multiple structures. See
        :func:`~reptar.structure.bond.is_bonded_atomic_radii` for an example.
    f_kwargs
        Keyword arguments for ``f_is_bonded``.
    use_ray
        Parallelize computations using ray.
    n_workers
        Number of parallel workers.

    Returns
    -------

        If bonding has changed with respect to the reference. ``True`` means that
        covalent bonding has changed.
    """
    # pylint: disable=invalid-name, no-member
    if R.ndim == 2:
        R = R[None, ...]

    if not use_ray:
        topo_check = TopologyCheck(Z, R_ref, f_is_bonded, f_kwargs)
        did_topo_change = topo_check({"data": R})["data"]
    else:
        ds = ray.data.from_numpy(np.array_split(R, n_workers))

        # Ensure that the results are kept in the same order as R
        context = ray.data.DataContext.get_current()
        context.execution_options = ray.data.ExecutionOptions(preserve_order=True)

        results = ds.map_batches(
            TopologyCheck,
            batch_size=None,
            batch_format="numpy",
            # pylint: disable-next=unexpected-keyword-arg
            compute=ray.data.ActorPoolStrategy(size=n_workers),
            zero_copy_batch=True,
            fn_constructor_args=(Z, R_ref, f_is_bonded, f_kwargs),
            num_cpus=1,
        )
        results = results.materialize()
        block_refs = results.to_numpy_refs()
        did_topo_change = np.hstack(
            [ray.get(block_ref)["data"] for block_ref in block_refs]
        )

    return did_topo_change
