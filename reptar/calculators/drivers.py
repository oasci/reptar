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

from __future__ import annotations
from collections.abc import Callable, Iterator
from . import Data
from ..utils import chunk_iterable, common_elements
from ..logger import ReptarLogger

log = ReptarLogger(__name__)

try:
    import ray
except ImportError:
    pass


def add_worker(
    workers_list: list[Callable],
    worker: Callable[["..."], Data],
    chunker: Iterator[list[int]],
    worker_args: tuple,
    worker_kwargs: dict,
    n_cpus_per_worker: int = 1,
):
    """Add a ray worker to a running list.

    Parameters
    ----------
    workers_list
        List of all running ray workers.
    worker
        Ray worker.
    chunker
        Provides a chunk of indices assigned to this worker.
    worker_args
        Positional arguments for worker with the exception of ``idxs`` (provided by
        ``chunker``).
    worker_kwargs
        Keyword arguments for worker.
    n_cpus_per_worker
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


class Driver:
    r"""Template driver class"""

    def __init__(
        self,
        worker: Callable[["..."], Data],
        worker_kwargs: dict,
        use_ray: bool = False,
        n_workers: int = 1,
        n_cpus_per_worker: int = 1,
        chunk_size: int = 50,
        start_slice: int = None,
        end_slice: int = None,
        ray_address: str = "auto",
    ) -> None:
        r"""
        Parameters
        ----------
        worker
            The desired worker function to compute energy and gradients. It should
            be the same as any previous calculations. The ``ray.remote`` decorator
            will be applied automatically.
        worker_kwargs
            The other keyword arguments for the worker function after
            ``R_idxs``, ``Z``, and ``R``.
        use_ray
            Use ray to parallelize calculations. If ``False``, calculations are
            done serially. ``False`` can be useful when running locally or only
            a few calculations are needed. ``True`` is useful for tons of
            calculations.
        n_workers
            Number of parallel workers to use if ``use_ray`` is ``True``.
        n_cpus_per_worker
            Number of CPU cores to provide each worker.
        chunk_size
            Number of calculations per task to do. This should be enough to make
            the ray task overhead significantly less than calculations.
        start_slice
            Trims ``R`` to start at this index.
        end_slice
            Trims ``R`` to end at this index.
        ray_address
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

    def run(self, data: Data, tasks: Iterator[str]) -> Data:
        r"""

        Parameters
        ----------
        data
            Parent data from current reptar file with initialized arrays required
            for ``worker`` arguments.
        tasks
            Types of computations to perform with the provided worker. For example,

            - ``E`` for electronic energy,
            - ``G`` for atomic gradients,
            - ``opt`` for geometry optimizations,
            - ``cube`` for properties on a cubic grid.

            .. note::

                The above labels are used to keep track of calculation completeness.
                :obj:`numpy.NaN` is used to flag incomplete/unfinished calculations.
                Other data need to be initialized for the calculation. For example,
                ``R_opt`` is needed for geometry optimizations and ``cube_R`` for
                cube property.

        Notes
        -----
        If ``data`` is starting calculations from scratch, it needs to have initialized
        arrays with :obj:`numpy.NaN` in the correct shape. Each ``job_type`` requires
        different data listed below. ``Z`` and ``R`` can never contain :obj:`numpy.NaN`.

        - ``E``: ``Z``, ``R``, ``E``.
        - ``G``: ``Z``, ``R``, ``E``, ``G``.
        - ``opt``: ``Z``, ``R``, ``E``, ``conv_opt``, ``R_opt``.
        - ``cube``: ``Z``, ``R``, ``cube_R``, ``cube_V``.
        """
        worker = self.worker
        data_todo = data.idxs_todo()

        # We make a list of indices that are missing all requested calculations.
        idxs_todo = data_todo[tasks[0]]
        if len(tasks) > 1:
            for i in range(1, len(tasks)):
                idxs_todo = common_elements(idxs_todo, data_todo[i])
        chunker = chunk_iterable(idxs_todo, self.chunk_size)

        if not self.use_ray:
            for idx in idxs_todo:
                data_worker = worker([idx], tasks, data, **self.worker_kwargs)
                data.update(data_worker)
                data.save()
        else:
            worker = ray.remote(worker)
            worker_args = (ray.put(tasks), ray.put(data))

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

                data_worker = ray.get(done_id)[0]
                data.update(data_worker)
                data.save()

                workers = add_worker(
                    workers,
                    worker,
                    chunker,
                    worker_args,
                    self.worker_kwargs,
                    self.n_cpus_per_worker,
                )

        return data
