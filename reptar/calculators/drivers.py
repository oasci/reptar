# MIT License
# 
# Copyright (c) 2022, Alex M. Maldonado
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

from ..utils import chunk_iterable
import math
import numpy as np
try:
    import ray
except ImportError:
    pass

class driverENGRAD:
    """Supervisor of energy+gradient workers.

    Creates and manages ray tasks using specified worker.
    """
    def __init__(
        self, Z, R, E, G, worker, worker_kwargs, n_cpus, n_cpus_worker=1,
        chunk_size=50, start_slice=None, end_slice=None, ray_address='auto'
    ):
        """
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
        worker : ``callable``
            The desired worker function to compute energy and gradients. It should
            be the same as any previous calculations. The ``ray.remote`` decorator
            will be applied automatically.
        worker_kwargs : :obj:`tuple`, ndim: ``1``
            The other keyword arguments for the worker function after ``Z``,
            ``R``, and ``R_idxs``.
        n_cpus : :obj:`int`
            Total number of CPUs we can use for ray tasks.
        n_cpus_worker : :obj:`int`, default: ``1``
            Number of CPUs to dedicate to each task.
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
        # Check arrays (obsessively)
        assert R.ndim == 3
        assert Z.shape[0] == R.shape[1]
        assert Z.shape[0] == G.shape[1]
        assert R.shape[0] == E.shape[0]
        assert R.shape[0] == G.shape[0]
        assert R.shape[2] == 3
        assert G.shape[2] == 3

        # Storing arrays and other information
        self.Z = ray.put(Z)
        self.R = ray.put(R)
        self.E = E
        self.G = G
        self.worker = worker
        self.worker_kwargs = worker_kwargs
        self.start_slice = start_slice
        self.end_slice = end_slice

        self.n_cpus = n_cpus
        self.n_cpus_worker = n_cpus_worker
        self.chunk_size = chunk_size
        self.n_workers = math.floor(n_cpus/n_cpus_worker)

        if not ray.is_initialized():
            ray.init(address=ray_address)
    
    def _idx_todo(self):
        """Indices of missing energies (calculations to do).

        Determines this by finding all ``NaN`` elements.
        
        Returns
        -------
        :obj:`numpy.ndarray`
            Indices for ``R`` that are missing energie values.
        """
        E = self.E[self.start_slice:self.end_slice]
        idx_todo = np.argwhere(np.isnan(E))[:,0]
        return idx_todo
    
    def run(self, saver=None):
        """Run the calculations.

        Parameters
        ----------
        saver : :obj:`reptar.calculators.save.Saver`, optional
            Save data after every worker finishes.

        Returns
        -------
        :obj:`numpy.ndarray`
            The total electronic energy array, ``E``, after all computations.
        :obj:`numpy.ndarray`
            The atomic gradients array, ``E``, after all computations.
        """
        worker = ray.remote(self.worker)
        idxs_todo = self._idx_todo()
        chunker = chunk_iterable(idxs_todo, self.chunk_size)

        # Initialize ray workers
        workers = []
        def add_worker(workers, chunker):
            try:
                chunk = list(next(chunker))
                workers.append(
                    worker.options(num_cpus=self.n_cpus_worker).remote(
                        chunk, self.Z, self.R, **self.worker_kwargs
                    )
                )
            except StopIteration:
                pass
        for _ in range(self.n_workers):
            add_worker(workers, chunker)
        
        # Start calculations
        while len(workers) != 0:
            done_id, workers = ray.wait(workers)
            
            idx_done, E_done, G_done = ray.get(done_id)[0]
            self.E[idx_done] = E_done
            self.G[idx_done] = G_done

            if saver is not None:
                saver.save((self.E, self.G))
            
            add_worker(workers, chunker)
        
        return self.E, self.G
