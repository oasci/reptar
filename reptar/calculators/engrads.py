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

import math
import numpy as np
import ray

class driverENGRAD:
    """Supervisor of energy+gradient workers.

    Creates and manages ray tasks using specified worker.

    Parameters
    ----------
    Z : :obj:`numpy.ndarray`
        Atomic numbers of the atoms with repsect to ``R``. Should have shape
        ``(i,)`` where ``i`` is the number of atoms in the system.
    R : :obj:`numpy.ndarray`
        Cartesian coordinates of all structures in group. Should have shape
        ``(j, i, 3)`` where ``j`` is the number of structures. Units are in
        Angstroms.
    E : :obj:`numpy.ndarray`
        Total electronic energies of all structures in ``R``. Energies that need
        to be calculated should have ``NaN`` in the element corresponding
        to the same index in ``R``. Should have shape ``(j,)``. Units are in
        Hartrees.
    G : :obj:`numpy.ndarray`
        Atomic gradients of all structures in ``R``. Gradients that need to be
        calculated should have ``NaN`` in the corresponding elements. Should
        have the same shape as ``R``. Units are in Hartrees/Angstroms
    worker : ``function``
        The desired worker function to compute energy and gradients. It should
        be the same as any previous calculations. The ``ray.remote`` decorator
        will be applied automatically.
    worker_args : :obj:`tuple`
        The other required arguments for the worker function after ``Z``,
        ``R``, and ``R_idxs``.
    n_cpus : :obj:`int`
        Total number of CPUs to use for all workers (same as ``num_cpus``
        when initalizing ray).
    n_cpus_worker : :obj:`int`, optional
        Number of CPUs to dedicate to each worker. Defaults to ``4``.
    start_slice : :obj:`int`, optional
        Trims ``R`` to start at this index.
    end_slice : :obj:`int`, optional
        Trims ``R`` to end at this index.
    """
    def __init__(
        self, Z, R, E, G, worker, worker_args, n_cpus, n_cpus_worker=4,
        start_slice=None, end_slice=None
    ):
        assert ray.is_initialized() == True
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
        self.worker_args = worker_args
        self.start_slice = start_slice
        self.end_slice = end_slice
        self.n_parallel_workers = math.floor(n_cpus/n_cpus_worker)
        
        # Initializes all ray tasks.
        worker = ray.remote(worker)
        self.worker_list = [
            worker.options(num_cpus=n_cpus_worker).remote(
                self.Z, self.R, R_idxs, *worker_args
            ) for R_idxs in self._task_R_idxs()
        ]
    
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
    
    def _task_R_idxs(self):
        """Generates the assigned structures for each worker.

        Yields
        ------
        :obj:`numpy.ndarray`
            ``R`` indices for the worker.
        """
        if not hasattr(self, 'task_size'):
            self.task_size = self._task_size()
            if self.task_size == 0:
                e = 'No calculations to run (task_size is zero)'
                raise ValueError(e)
        task_size = self.task_size
        idx_todo = self._idx_todo()
        for i in range(0, len(idx_todo), task_size):
            r_idxs = idx_todo[i:i + task_size]
            yield r_idxs
    
    def _task_size(self):
        """Determines the number of calculations per worker.
        
        Returns
        -------
        :obj:`int`
            Desired number of calculations per work.
        """
        n_todo = len(self._idx_todo())
        task_size = math.ceil(n_todo/self.n_parallel_workers)
        return task_size
    
    def run(self):
        """Run the workers.

        Returns
        -------
        :obj:`numpy.ndarray`
            The total electronic energy array, ``E``, after all computations.
        :obj:`numpy.ndarray`
            The atomic gradients array, ``E``, after all computations.
        """
        worker_list = self.worker_list
        while len(worker_list) != 0:
            done_id, worker_list = ray.wait(worker_list)
            
            R_idx_done, E_done, G_done = ray.get(done_id)[0]
            self.E[R_idx_done] = E_done
            self.G[R_idx_done] = G_done
        
        return self.E, self.G
