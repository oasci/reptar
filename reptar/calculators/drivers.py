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
from ..utils import chunk_iterable
from ..logger import ReptarLogger

log = ReptarLogger(__name__)

try:
    import ray
except ImportError:
    pass


class DriverEnergy:
    r"""Supervisor of energy+gradient workers.

    Creates and manages ray tasks using specified worker.
    """

    def __init__(
        self,
        worker,
        worker_kwargs,
        use_ray=False,
        n_workers=1,
        n_cpus_per_worker=1,
        chunk_size=50,
        start_slice=None,
        end_slice=None,
        ray_address="auto",
    ):
        r"""
        Parameters
        ----------
        worker : ``callable``
            The desired worker function to compute energy and gradients. It should
            be the same as any previous calculations. The ``ray.remote`` decorator
            will be applied automatically.
        worker_kwargs : :obj:`tuple`, ndim: ``1``
            The other keyword arguments for the worker function after
            ``R_idxs``, ``Z``, and ``R``.
        use_ray : :obj:`bool`, default: ``False``
            Use ray to parallelize calculations. If ``False``, calculations are
            done serially. ``False`` can be useful when running locally or only
            a few calculations are needed. ``True`` is useful for tons of
            calculations.
        n_workers : :obj:`int`, default: ``1``
            Number of parallel workers to use if ``use_ray`` is ``True``.
        n_cpus_per_worker : :obj:`int`, default: ``1``
            Number of CPU cores to provide each worker.
        chunk_size : :obj:`int`, default: ``50``
            Number of calculations per task to do. This should be enough to make
            the ray task overhead significantly less than calculations.
        start_slice : :obj:`int`, default: ``None``
            Trims ``R`` to start at this index.
        end_slice : :obj:`int`, default: ``None``
            Trims ``R`` to end at this index.
        ray_address : :obj:`str`, default: ``'auto'``
            Ray cluster address to connect to.
        """
        self.worker = worker
        self.worker_kwargs = worker_kwargs
        self.start_slice = start_slice
        self.end_slice = end_slice

        self.n_workers = n_workers
        self.n_cpus_per_worker = n_cpus_per_worker
        self.chunk_size = chunk_size

        self.use_ray = use_ray
        if use_ray:
            if not ray.is_initialized():
                ray.init(address=ray_address)

    @staticmethod
    def _idx_todo(E, start_slice=None, end_slice=None):
        r"""Indices of missing energies (calculations to do).

        Determines this by finding all ``NaN`` elements.

        Returns
        -------
        :obj:`numpy.ndarray`
            Indices for ``R`` that are missing energie values.
        """
        # pylint: disable-next=invalid-name
        E = E[start_slice:end_slice]
        idx_todo = np.argwhere(np.isnan(E))[:, 0]
        return idx_todo

    def run(self, Z, R, E, saver=None):
        r"""Run the calculations.

        Parameters
        ----------
        Z : :obj:`numpy.ndarray`, ndim: ``1``
            Atomic numbers of the atoms with respect to ``R``. Should have shape
            ``(i,)`` where ``i`` is the number of atoms in the system.
        R : :obj:`numpy.ndarray`, ndim: ``3``
            Cartesian coordinates of all structures in group. Should have shape
            ``(j, i, 3)`` where ``j`` is the number of structures. Units are in
            Angstroms.
        E : :obj:`numpy.ndarray`, ndim: ``1``
            Total electronic energies of all structures in ``R``. Energies that need
            to be calculated should have ``NaN`` in the element corresponding
            to the same index in ``R``. Should have shape ``(j,)``. Units are in
            Hartrees.
        saver : :obj:`reptar.Saver`, optional
            Save data after every worker finishes. Passes ``E`` into
            :meth:`~reptar.Saver.save`.

        Returns
        -------
        :obj:`numpy.ndarray`
            The total electronic energy array, ``E``, after all computations.
        """
        # Check arrays (obsessively)
        assert R.ndim == 3
        assert Z.shape[0] == R.shape[1]
        assert R.shape[0] == E.shape[0]
        assert R.shape[2] == 3

        worker = self.worker
        idxs_todo = self._idx_todo(E, self.start_slice, self.end_slice)
        chunker = chunk_iterable(idxs_todo, self.chunk_size)

        if not self.use_ray:
            for idx in idxs_todo:
                # pylint: disable-next=invalid-name
                _, E_done = worker([idx], Z, R, **self.worker_kwargs)
                print(E_done)
                E[idx] = E_done

                if saver is not None:
                    saver.save(E)
        else:
            worker = ray.remote(worker)
            Z = ray.put(Z)
            R = ray.put(R)

            # Initialize ray workers
            workers = []

            def add_worker(workers, chunker):
                try:
                    chunk = list(next(chunker))
                    workers.append(
                        worker.options(num_cpus=self.n_cpus_per_worker).remote(
                            chunk, Z, R, **self.worker_kwargs
                        )
                    )
                except StopIteration:
                    pass

            for _ in range(self.n_workers):
                add_worker(workers, chunker)

            # Start calculations
            while len(workers) != 0:
                done_id, workers = ray.wait(workers)

                idx_done, E_done = ray.get(done_id)[0]  # pylint: disable=invalid-name
                E[idx_done] = E_done

                if saver is not None:
                    saver.save(E)

                add_worker(workers, chunker)

            # Cleanup ray object store
            del Z
            del R

        return E


class DriverEnGrad:
    r"""Supervisor of energy+gradient workers.

    Creates and manages ray tasks using specified worker.
    """

    def __init__(
        self,
        worker,
        worker_kwargs,
        use_ray=False,
        n_workers=1,
        n_cpus_per_worker=1,
        chunk_size=50,
        start_slice=None,
        end_slice=None,
        ray_address="auto",
    ):
        r"""
        Parameters
        ----------
        worker : ``callable``
            The desired worker function to compute energy and gradients. It should
            be the same as any previous calculations. The ``ray.remote`` decorator
            will be applied automatically.
        worker_kwargs : :obj:`tuple`, ndim: ``1``
            The other keyword arguments for the worker function after
            ``R_idxs``, ``Z``, and ``R``.
        use_ray : :obj:`bool`, default: ``False``
            Use ray to parallelize calculations. If ``False``, calculations are
            done serially. ``False`` can be useful when running locally or only
            a few calculations are needed. ``True`` is useful for tons of
            calculations.
        n_workers : :obj:`int`, default: ``1``
            Number of parallel workers to use if ``use_ray`` is ``True``.
        n_cpus_per_worker : :obj:`int`, default: ``1``
            Number of CPU cores to provide each worker.
        chunk_size : :obj:`int`, default: ``50``
            Number of calculations per task to do. This should be enough to make
            the ray task overhead significantly less than calculations.
        start_slice : :obj:`int`, default: ``None``
            Trims ``R`` to start at this index.
        end_slice : :obj:`int`, default: ``None``
            Trims ``R`` to end at this index.
        ray_address : :obj:`str`, default: ``'auto'``
            Ray cluster address to connect to.
        """
        # Storing arrays and other information
        self.worker = worker
        self.worker_kwargs = worker_kwargs
        self.start_slice = start_slice
        self.end_slice = end_slice

        self.n_workers = n_workers
        self.n_cpus_per_worker = n_cpus_per_worker
        self.chunk_size = chunk_size

        self.use_ray = use_ray
        if use_ray:
            if not ray.is_initialized():
                ray.init(address=ray_address)

    @staticmethod
    def _idx_todo(E, start_slice=None, end_slice=None):
        r"""Indices of missing energies (calculations to do).

        Determines this by finding all ``NaN`` elements.

        Returns
        -------
        :obj:`numpy.ndarray`
            Indices for ``R`` that are missing energie values.
        """
        # pylint: disable-next=invalid-name
        E = E[start_slice:end_slice]
        idx_todo = np.argwhere(np.isnan(E))[:, 0]
        return idx_todo

    def run(self, Z, R, E, G, saver=None):
        r"""Run the calculations.

        Parameters
        ----------
        Z : :obj:`numpy.ndarray`, ndim: ``1``
            Atomic numbers of the atoms with respect to ``R``. Should have shape
            ``(i,)`` where ``i`` is the number of atoms in the system.
        R : :obj:`numpy.ndarray`, ndim: ``3``
            Cartesian coordinates of all structures in group. Should have shape
            ``(j, i, 3)`` where ``j`` is the number of structures. Units are in
            Angstroms.
        E : :obj:`numpy.ndarray`, ndim: ``1``
            Total electronic energies of all structures in ``R``. Energies that need
            to be calculated should have ``NaN`` in the element corresponding
            to the same index in ``R``. Should have shape ``(j,)``. Units are in
            Hartrees.
        G : :obj:`numpy.ndarray`, ndim: ``3``
            Atomic gradients of all structures in ``R``. Gradients that need to be
            calculated should have ``NaN`` in the corresponding elements. Should
            have the same shape as ``R``. Units are in Hartrees/Angstrom.
        saver : :obj:`reptar.Saver`, optional
            Save data after every worker finishes. Passes ``E, G`` into
            :meth:`~reptar.Saver.save`.

        Returns
        -------
        :obj:`numpy.ndarray`
            The total electronic energy array, ``E``, after all computations.
        :obj:`numpy.ndarray`
            The atomic gradients array, ``E``, after all computations.
        """
        # Check arrays (obsessively)
        assert R.ndim == 3
        assert Z.shape[0] == R.shape[1]
        assert Z.shape[0] == G.shape[1]
        assert R.shape[0] == E.shape[0]
        assert R.shape[0] == G.shape[0]
        assert R.shape[2] == 3
        assert G.shape[2] == 3

        worker = self.worker
        idxs_todo = self._idx_todo(E, self.start_slice, self.end_slice)
        chunker = chunk_iterable(idxs_todo, self.chunk_size)

        if not self.use_ray:
            for idx in idxs_todo:
                # pylint: disable-next=invalid-name
                _, E_done, G_done = worker([idx], Z, R, **self.worker_kwargs)
                E[idx] = E_done
                G[idx] = G_done

                if saver is not None:
                    saver.save((E, G))
        else:
            worker = ray.remote(worker)
            Z = ray.put(Z)
            R = ray.put(R)

            # Initialize ray workers
            workers = []

            def add_worker(workers, chunker):
                try:
                    chunk = list(next(chunker))
                    workers.append(
                        worker.options(num_cpus=self.n_cpus_per_worker).remote(
                            chunk, Z, R, **self.worker_kwargs
                        )
                    )
                except StopIteration:
                    pass

            for _ in range(self.n_workers):
                add_worker(workers, chunker)

            # Start calculations
            while len(workers) != 0:
                done_id, workers = ray.wait(workers)

                # pylint: disable-next=invalid-name
                idx_done, E_done, G_done = ray.get(done_id)[0]
                E[idx_done] = E_done
                G[idx_done] = G_done

                if saver is not None:
                    saver.save(E, G)

                add_worker(workers, chunker)

            # Cleanup ray object store
            del Z
            del R

        return E, G


class DriverOpt:
    r"""Supervisor of optimization workers.

    Creates and manages ray tasks using specified worker.
    """

    def __init__(
        self,
        worker,
        worker_kwargs,
        use_ray=False,
        n_workers=1,
        n_cpus_per_worker=1,
        chunk_size=1,
        start_slice=None,
        end_slice=None,
        ray_address="auto",
    ):
        r"""
        Parameters
        ----------
        worker : ``callable``
            The desired worker function to compute energy and gradients. It should
            be the same as any previous calculations. The ``ray.remote`` decorator
            will be applied automatically.
        worker_kwargs : :obj:`tuple`, ndim: ``1``
            The other keyword arguments for the worker function after
            ``R_idxs``, ``Z``, and ``R``.
        use_ray : :obj:`bool`, default: ``False``
            Use ray to parallelize calculations. If ``False``, calculations are
            done serially. ``False`` can be useful when running locally or only
            a few calculations are needed. ``True`` is useful for tons of
            calculations.
        n_workers : :obj:`int`, default: ``1``
            Number of parallel workers to use if ``use_ray`` is ``True``.
        n_cpus_per_worker : :obj:`int`, default: ``1``
            Number of CPU cores to provide each worker.
        chunk_size : :obj:`int`, default: ``1``
            Number of calculations per task to do. This should be enough to make
            the ray task overhead significantly less than calculations.
        start_slice : :obj:`int`, default: ``None``
            Trims ``R`` to start at this index.
        end_slice : :obj:`int`, default: ``None``
            Trims ``R`` to end at this index.
        ray_address : :obj:`str`, default: ``'auto'``
            Ray cluster address to connect to. Passes ``opt_conv, R_opt, E, G`` into
            :meth:`~reptar.Saver.save`.
        """
        # Storing arrays and other information
        self.worker = worker
        self.worker_kwargs = worker_kwargs
        self.start_slice = start_slice
        self.end_slice = end_slice

        self.n_workers = n_workers
        self.n_cpus_per_worker = n_cpus_per_worker
        self.chunk_size = chunk_size

        self.use_ray = use_ray
        if use_ray:
            if not ray.is_initialized():
                ray.init(address=ray_address)

    @staticmethod
    def _idx_todo(opt_conv, start_slice=None, end_slice=None):
        r"""Indices of NaN geometries (calculations to do).

        Returns
        -------
        :obj:`numpy.ndarray`
            Indices for ``R`` that are missing energie values.
        """
        return np.where(~opt_conv[start_slice:end_slice])[0]

    # pylint: disable-next=invalid-name
    def run(self, Z, R, R_opt, E, G, saver=None):
        r"""Run the calculations.

        Parameters
        ----------
        Z : :obj:`numpy.ndarray`, ndim: ``1``
            Atomic numbers of the atoms with respect to ``R``. Should have shape
            ``(i,)`` where ``i`` is the number of atoms in the system.
        R : :obj:`numpy.ndarray`, ndim: ``3``
            Cartesian coordinates of all structures in group. Should have shape
            ``(j, i, 3)`` where ``j`` is the number of structures. Units are in
            Angstroms.
        R_opt : :obj:`numpy.ndarray`, ndim: ``3``
            Optimized Cartesian coordinates of all structures in group. This
            is used to identifying which optimizations need to be done.
        E : :obj:`numpy.ndarray`, ndim: ``1``
            Total electronic energies of all structures in ``R``. Energies that need
            to be calculated should have ``NaN`` in the element corresponding
            to the same index in ``R``. Should have shape ``(j,)``. Units are in
            Hartrees.
        G : :obj:`numpy.ndarray`, ndim: ``3``
            Atomic gradients of all structures in ``R``. Gradients that need to be
            calculated should have ``NaN`` in the corresponding elements. Should
            have the same shape as ``R``. Units are in Hartrees/Angstrom.
        saver : :obj:`reptar.Saver`, optional
            Save data after every worker finishes.

        Returns
        -------
        :obj:`numpy.ndarray`
            The total electronic energy array, ``E``, after all computations.
        :obj:`numpy.ndarray`
            The atomic gradients array, ``E``, after all computations.
        """
        # Check arrays (obsessively)
        assert R.ndim == 3
        assert Z.shape[0] == R.shape[1]
        assert R.shape[2] == 3

        opt_conv = ~np.isnan(R_opt[:, 0, 0])

        worker = self.worker
        idxs_todo = self._idx_todo(opt_conv, self.start_slice, self.end_slice)
        chunker = chunk_iterable(idxs_todo, self.chunk_size)

        if not self.use_ray:
            for idx in idxs_todo:
                # pylint: disable-next=invalid-name
                _, opt_conv_done, R_opt_done, E_done, G_done = worker(
                    [idx], Z, R, **self.worker_kwargs
                )
                opt_conv[idx] = opt_conv_done
                R_opt[idx] = R_opt_done[0]
                E[idx] = E_done[0]
                G[idx] = G_done[0]

                if saver is not None:
                    saver.save(opt_conv, R_opt, E, G)
        else:
            worker = ray.remote(worker)
            Z = ray.put(Z)
            R = ray.put(R)

            # Initialize ray workers
            workers = []

            def add_worker(workers, chunker):
                try:
                    chunk = list(next(chunker))
                    workers.append(
                        worker.options(num_cpus=self.n_cpus_per_worker).remote(
                            chunk, Z, R, **self.worker_kwargs
                        )
                    )
                except StopIteration:
                    pass

            for _ in range(self.n_workers):
                add_worker(workers, chunker)

            # Start calculations
            while len(workers) != 0:
                done_id, workers = ray.wait(workers)

                # pylint: disable-next=invalid-name
                idxs_done, opt_conv_done, R_opt_done, E_done, G_done = ray.get(done_id)[
                    0
                ]
                opt_conv[idxs_done] = opt_conv_done
                R_opt[idxs_done] = R_opt_done
                E[idxs_done] = E_done
                G[idxs_done] = G_done

                if saver is not None:
                    saver.save(opt_conv, R_opt, E, G)

                add_worker(workers, chunker)

        return opt_conv, R_opt, E, G
