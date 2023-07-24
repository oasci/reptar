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

from abc import ABC, abstractmethod
import numpy as np
from ..utils import chunk_iterable
from ..logger import ReptarLogger

log = ReptarLogger(__name__)

try:
    import ray
except ImportError:
    pass


def add_worker(
    workers_list, worker, chunker, worker_args, worker_kwargs, n_cpus_per_worker=1
):
    """Add a ray worker to a running list.

    Parameters
    ----------
    workers_list : :obj:`list`
        List of all running ray workers.
    worker : ``callable``
        Ray worker.
    chunker : ``iterable``
        Provides a chunk of indices assigned to this worker.
    worker_args : :obj:`list`
        Positional arguments for worker with the exception of ``idxs`` (provided by
        ``chunker``).
    worker_kwargs : :obj:`dict`
        Keyword arguments for worker.
    n_cpus_per_worker : :obj:`int`, default: ``1``
        Number of CPU cores to assign the worker.

    Returns
    -------
    :obj:`list`
        Updated ``workers_list``.
    """
    try:
        chunk = list(next(chunker))
        workers_list.append(
            worker.options(num_cpus=n_cpus_per_worker).remote(
                chunk, *worker_args, **worker_kwargs
            )
        )
    except StopIteration:
        pass
    return workers_list


class Driver(ABC):
    r"""Template driver class"""

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
    @abstractmethod
    def idx_todo(args, start_slice=None, end_slice=None):
        r"""Determine which calculations need to be done.

        Parameters
        ----------
        args : :obj:`tuple`
            ``run()`` arguments determined by the specific driver.
        start_slice : :obj:`int`, default: ``None``
            Slice arrays in ``args`` starting at this index.
        end_slice : :obj:`int`, default: ``None``
            Slice arrays` in ``args`` stopping at this index.

        Returns
        -------
        :obj:`numpy.ndarray`
            Indices of source data where calculations need to be ran.
        """

    @abstractmethod
    def check_run_args(self, args):
        r"""Perform any checks on ``run()`` arguments."""

    @abstractmethod
    def prep_worker_args(self, args):
        r"""Prepare worker arguments."""

    @abstractmethod
    def setup_results(self, args):
        r"""Setup attribute containing the results of ``run()``."""

    @abstractmethod
    def process_worker_returns(self, returns):
        r"""Process worker returns and store in results."""

    @abstractmethod
    def prepare_run_returns(self):
        r"""Prepare and process final ``run()`` returns"""

    def _run(self, *args, saver=None):
        self.check_run_args(args)
        self.saver = saver

        worker = self.worker
        idxs_todo = self.idx_todo(args, self.start_slice, self.end_slice)
        chunker = chunk_iterable(idxs_todo, self.chunk_size)
        worker_args = self.prep_worker_args(args)
        self.setup_results(args)

        if not self.use_ray:
            for idx in idxs_todo:
                # pylint: disable-next=invalid-name
                worker_returns = worker([idx], *worker_args, **self.worker_kwargs)
                self.process_worker_returns(worker_returns)
        else:
            worker = ray.remote(worker)

            # Initialize ray workers
            workers = []

            for _ in range(self.n_workers):
                workers = add_worker(
                    workers,
                    worker,
                    chunker,
                    worker_args,
                    self.worker_kwargs,
                    self.n_cpus_per_worker,
                )

            # Start calculations
            while len(workers) != 0:
                done_id, workers = ray.wait(workers)

                worker_returns = ray.get(done_id)[0]  # pylint: disable=invalid-name
                self.process_worker_returns(worker_returns)

                workers = add_worker(
                    workers,
                    worker,
                    chunker,
                    worker_args,
                    self.worker_kwargs,
                    self.n_cpus_per_worker,
                )

        return self.prepare_run_returns()


class DriverEnergy(Driver):
    r"""Supervisor of energy workers."""

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
        Driver.__init__(
            self,
            worker,
            worker_kwargs,
            use_ray,
            n_workers,
            n_cpus_per_worker,
            chunk_size,
            start_slice,
            end_slice,
            ray_address,
        )

    @staticmethod
    def idx_todo(args, start_slice=None, end_slice=None):
        # args = (Z, R, E)
        idx_todo = np.argwhere(np.isnan(args[2][start_slice:end_slice]))[:, 0]
        return idx_todo

    def check_run_args(self, args):
        # args = (Z, R, E)
        assert args[1].ndim == 3
        assert args[0].shape[0] == args[1].shape[1]
        assert args[1].shape[0] == args[2].shape[0]
        assert args[1].shape[2] == 3

    def prep_worker_args(self, args):
        # args = (Z, R, E)
        if self.use_ray:
            return ray.put(args[0]), ray.put(args[1])
        return args[0], args[1]

    def setup_results(self, args):
        # args = (Z, R, E)
        self.results = {"E": args[2]}

    def process_worker_returns(self, returns):
        # returns = (idx_done, E_done)
        self.results["E"][returns[0]] = returns[1]

        if self.saver is not None:
            self.saver.save(self.results["E"])

    def prepare_run_returns(self):
        return self.results["E"]

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
        return Driver._run(self, Z, R, E, saver=saver)


class DriverEnGrad(Driver):
    r"""Supervisor of energy+gradient workers."""

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
        Driver.__init__(
            self,
            worker,
            worker_kwargs,
            use_ray,
            n_workers,
            n_cpus_per_worker,
            chunk_size,
            start_slice,
            end_slice,
            ray_address,
        )

    @staticmethod
    def idx_todo(args, start_slice=None, end_slice=None):
        # args = (Z, R, E, G)
        idx_todo = np.argwhere(np.isnan(args[2][start_slice:end_slice]))[:, 0]
        return idx_todo

    def check_run_args(self, args):
        # args = (Z, R, E, G)
        assert args[1].ndim == 3
        assert args[0].shape[0] == args[1].shape[1]
        assert args[0].shape[0] == args[3].shape[1]
        assert args[1].shape[0] == args[3].shape[0]
        assert args[1].shape[0] == args[3].shape[0]
        assert args[1].shape[2] == 3
        assert args[3].shape[2] == 3

    def prep_worker_args(self, args):
        # args = (Z, R, E, G)

        if self.use_ray:
            return ray.put(args[0]), ray.put(args[1])
        return args[0], args[1]

    def setup_results(self, args):
        # args = (Z, R, E, G)
        self.results = {"E": args[2], "G": args[3]}

    def process_worker_returns(self, returns):
        # returns = (idx_done, E_done, G_done)
        self.results["E"][returns[0]] = returns[1]
        self.results["G"][returns[0]] = returns[2]

        if self.saver is not None:
            self.saver.save(self.results["E"], self.results["G"])

    def prepare_run_returns(self):
        return self.results["E"], self.results["G"]

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
            The atomic gradients array, ``G``, after all computations.
        """
        return Driver._run(self, Z, R, E, G, saver=saver)


class DriverOpt(Driver):
    r"""Supervisor of optimization workers."""

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
        Driver.__init__(
            self,
            worker,
            worker_kwargs,
            use_ray,
            n_workers,
            n_cpus_per_worker,
            chunk_size,
            start_slice,
            end_slice,
            ray_address,
        )

    @staticmethod
    def idx_todo(args, start_slice=None, end_slice=None):
        # args = (Z, R, conv_opt, R_opt, E_opt)
        return np.where(~args[2][start_slice:end_slice])[0]

    def check_run_args(self, args):
        # args = (Z, R, conv_opt, R_opt, E_opt)
        assert args[1].ndim == 3
        assert args[0].shape[0] == args[1].shape[1]
        assert args[1].shape[2] == 3

    def prep_worker_args(self, args):
        # args = (Z, R, conv_opt, R_opt, E_opt)
        if self.use_ray:
            return ray.put(args[0]), ray.put(args[1])
        return args[0], args[1]

    def setup_results(self, args):
        # args = (Z, R, conv_opt, R_opt, E_opt)
        self.results = {"conv_opt": args[2], "R_opt": args[3], "E_opt": args[4]}

    def process_worker_returns(self, returns):
        # returns = (idx_done, conv_opt_done, R_opt_done, E_opt_one)
        self.results["conv_opt"][returns[0]] = returns[1]
        self.results["R_opt"][returns[0]] = returns[2]
        self.results["E_opt"][returns[0]] = returns[3]

        if self.saver is not None:
            self.saver.save(
                self.results["conv_opt"], self.results["R_opt"], self.results["E_opt"]
            )

    def prepare_run_returns(self):
        return self.results["conv_opt"], self.results["R_opt"], self.results["E_opt"]

    # pylint: disable-next=invalid-name
    def run(self, Z, R, conv_opt, R_opt, E_opt, saver=None):
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
        conv_opt : :obj:`numpy.ndarray`, ndim: ``1``
            Boolean flags for if the geometry optimization was successful.
        R_opt : :obj:`numpy.ndarray`, ndim: ``3``
            Optimized Cartesian coordinates of all structures in group. This
            is used to identifying which optimizations need to be done.
        E_opt : :obj:`numpy.ndarray`, ndim: ``1``
            Total electronic energies of all structures in ``R``. Energies that need
            to be calculated should have ``NaN`` in the element corresponding
            to the same index in ``R``. Should have shape ``(j,)``. Units are in
            Hartrees.
        saver : :obj:`reptar.Saver`, optional
            Save data after every worker finishes.

        Returns
        -------
        :obj:`numpy.ndarray`
            If the geometry optimization converged, ``conv_opt``.
        :obj:`numpy.ndarray`
            The last step of geometry optimization coordinates, ``R_opt``.
        :obj:`numpy.ndarray`
            The last electronic energy of the geometry optimization, ``E_opt``.
        """
        return Driver._run(self, Z, R, conv_opt, R_opt, E_opt, saver=saver)


class DriverCube(Driver):
    r"""Supervisor of cube workers."""

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
        Driver.__init__(
            self,
            worker,
            worker_kwargs,
            use_ray,
            n_workers,
            n_cpus_per_worker,
            chunk_size,
            start_slice,
            end_slice,
            ray_address,
        )

    @staticmethod
    def idx_todo(args, start_slice=None, end_slice=None):
        # args = (Z, R, cube_R, cube_V)
        # Checks the first element per structure in cube_V for NaN
        return np.argwhere(np.isnan(args[3][start_slice:end_slice][:, 0]))[:, 0]

    def check_run_args(self, args):
        # args = (Z, R, cube_R, cube_V)
        assert args[1].ndim == 3
        assert args[0].shape[0] == args[1].shape[1]
        assert args[1].shape[2] == 3

    def prep_worker_args(self, args):
        # args = (Z, R, cube_R, cube_V)
        if self.use_ray:
            return ray.put(args[0]), ray.put(args[1]), args[3].shape[-1]
        return args[0], args[1], args[3].shape[-1]

    def setup_results(self, args):
        # args = (Z, R, cube_R, cube_V)
        self.results = {"cube_R": args[2], "cube_V": args[3]}

    def process_worker_returns(self, returns):
        # returns = (idx_done, cube_R, cube_V)
        self.results["cube_R"][returns[0]] = returns[1]
        self.results["cube_V"][returns[0]] = returns[2]

        if self.saver is not None:
            self.saver.save(self.results["cube_R"], self.results["cube_V"])

    def prepare_run_returns(self):
        return self.results["cube_R"], self.results["cube_V"]

    # pylint: disable-next=invalid-name
    def run(self, Z, R, cube_R, cube_V, saver=None):
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
        cube_R : :obj:`numpy.ndarray`, ndim: ``3``
            Cartesian coordinates of points where a property is probed.
        cube_V : :obj:`numpy.ndarray`, ndim: ``2``
            Property values at the same Cartesian coordinates as ``cube_R``.
        saver : :obj:`reptar.Saver`, optional
            Save data after every worker finishes.

        Returns
        -------
        :obj:`numpy.ndarray`
            Cartesian coordinates of points where a property is probed.
        :obj:`numpy.ndarray`
            Property values at the same Cartesian coordinates.
        """
        return Driver._run(self, Z, R, cube_R, cube_V, saver=saver)
