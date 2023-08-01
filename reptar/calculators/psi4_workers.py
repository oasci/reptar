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
import os
from collections.abc import Iterator
from tempfile import TemporaryDirectory
import numpy as np
from . import Data
from .utils import initialize_worker_data
from ..logger import ReptarLogger
from ..parsers.gaussian_cube import parse_cube

log = ReptarLogger(__name__)

try:
    import psi4
    from optking.exceptions import AlgError
except ImportError:
    pass


def _do_psi4_energy(
    method: str, mol: psi4.core.Molecule
) -> tuple[float, psi4.core.Wavefunction]:
    r"""Perform energy calculation in Psi4."""
    e, wfn = psi4.energy(method, molecule=mol, return_wfn=True)
    return e, wfn


def _do_psi4_grad(
    method: str, mol: psi4.core.Molecule
) -> tuple[float, np.ndarray, psi4.core.Wavefunction]:
    r"""Perform gradient calculation in Psi4."""
    g, wfn = psi4.gradient(method, molecule=mol, return_wfn=True)
    g = g.to_array()
    g /= psi4.constants.bohr2angstroms
    e = wfn.energy()
    return e, g, wfn


def _do_psi4_cube(
    method: str,
    mol: psi4.core.Molecule,
    cube_path: str,
    wfn: psi4.core.Wavefunction | None = None,
) -> tuple[np.ndarray, np.ndarray, psi4.core.Wavefunction]:
    r"""Perform cube property calculation in Psi4."""
    if wfn is None:
        _, wfn = psi4.energy(method, molecule=mol, return_wfn=True)
    psi4.cubeprop(wfn)
    cube_r, cube_v = parse_cube(cube_path)

    return cube_r, cube_v, wfn


def _do_psi4_opt(
    method: str, mol: psi4.core.Molecule
) -> tuple[psi4.core.Wavefunction, bool, np.ndarray, float]:
    r"""Perform geometry optimization in Psi4."""
    try:
        e, wfn = psi4.opt(method, molecule=mol, return_wfn=True)
        r_opt = np.asarray(wfn.molecule().geometry())
        r_conv_opt = True
    except (psi4.OptimizationConvergenceError, AlgError) as ex:
        r_conv_opt = False
        r_opt = np.asarray(ex.wfn.molecule().geometry())
        e = ex.wfn.energy()
    finally:
        r_opt *= psi4.constants.bohr2angstroms  # pylint: disable=undefined-variable
    return r_conv_opt, r_opt, e, wfn


def _do_psi4_task(
    task: str,
    data_worker: Data,
    idx: int,
    method: str,
    mol: psi4.core.Molecule,
    options: dict[str, "Any"],
    cube_file_name: str | None,
    wfn: psi4.core.Wavefunction | None,
) -> tuple[Data, psi4.core.Wavefunction]:
    worker_idx = np.array([idx], dtype=np.uint64)
    if task == "opt":
        r_conv_opt, r_opt, e, wfn = _do_psi4_opt(method, mol)
        data_worker.add_subset("R_opt", worker_idx, r_opt)
        data_worker.add_subset("conv_opt", worker_idx, r_conv_opt)
        data_worker.add_subset("E", worker_idx, e)
    elif task == "G":
        e, g, wfn = _do_psi4_grad(method, mol)
        data_worker.add_subset("E", worker_idx, e)
        data_worker.add_subset("G", worker_idx, g)
    elif task == "E":
        e, wfn = _do_psi4_energy(method, mol)
        data_worker.add_subset("E", worker_idx, e)
    elif task == "cube":
        with TemporaryDirectory() as temp_dir:
            options["cubeprop_filepath"] = temp_dir
            psi4.set_options(options)
            cube_path = os.path.join(temp_dir, cube_file_name)
            cube_r, cube_v, wfn = _do_psi4_cube(method, mol, cube_path, wfn)

            cube_r_slice = (
                slice(None, None, None),
                slice(None, cube_r.shape[0], None),
                slice(None, None, None),
            )
            cube_v_slice = (slice(None, None, None), slice(None, cube_v.shape[0], None))
            data_worker.add_subset("cube_R", cube_r_slice, cube_r)
            data_worker.add_subset("cube_V", cube_v_slice, cube_v)
    return data_worker, wfn


def psi4_worker(
    idxs: Iterator[int],
    tasks: Iterator[str],
    data: Data,
    charge: int = 0,
    mult: int = 1,
    method: str = "mp2",
    options: dict[str, "Any"] = None,
    threads: int = 1,
    mem: str | float | int = "1 GB",
    total_grid_points: int | None = None,
) -> Data:
    r"""Worker function for optimizations using Psi4.

    Parameters
    ----------
    idxs
        Indices of the structures from ``R`` to compute energies and gradients
        for.
    tasks
        Calculations this worker needs to run in the order specified here. In general,
        we recommend using this ordering: ``opt``, ``E``, ``G``, ``cube``.
    data
        All required data required for computations such as ``Z``, ``R``, ``conv_opt``,
        etc.
    charge
        Total system charge.
    mult
        Total system multiplicity.
    method
        Specifies the Psi4 method used for all calculations. For more information,
        please see `the Psi4 documentation <https://psicode.org/psi4manual/master
        /methods.html>`__ for your specific version.
    options
        `Psi4 control keywords <https://psicode.org/psi4manual/master/
        psithoninput.html#job-control-keywords>`__
        using the PsiAPI format. Some common ones are ``basis``,
        ``e_convergence``, ``d_convergence``, and ``reference``.
    threads
        Number of threads for Psi4. This is almost always the number of cores
        being used for the worker. For more information, see the
        `documentation <https://psicode.org/psi4manual/master/api/
        psi4.core.set_num_threads.html#psi4.core.set_num_threads>`__.
    mem
        The amount of memory available. For more information, see the
        `documentation <https://psicode.org/psi4manual/master/api/
        psi4.driver.set_memory.html#psi4.driver.set_memory>`__.
    total_grid_points
        Number of grid points to initialize cube arrays with.

    Returns
    -------
    :obj:`reptar.calculators.Data`
        Data object from this worker.

    Notes
    -----
    Psi4 uses QCElemental to build molecules from arrays.
    There is some postprocessing in the `from_arrays routine
    <http://docs.qcarchive.molssi.org/projects/QCElemental/en/latest/api/
    qcelemental.molparse.from_arrays.html#from-arrays>`__.
    QCElemental will translate and rotate molecules which can affect computed
    properties. Here, we set ``fix_com`` and ``fix_orientation`` to ``True`` to avoid
    this.
    """
    implemented_tasks = ["E", "G", "opt", "cube"]
    for task in tasks:
        if task not in implemented_tasks:
            raise ValueError(f"Task ({task}) is not implemented in this worker")

    psi4.core.set_num_threads(threads)
    psi4.set_memory(mem)

    # Obtain all required input data
    Z = data.Z
    R = data.R[idxs]

    # Initialize data object with resulting arrays
    data_worker = initialize_worker_data(
        tasks, Data(), R, total_grid_points=total_grid_points
    )
    data_worker.idxs_source = idxs

    cube_file_name = None
    if "cube" in tasks:
        assert "cubeprop_tasks" in options
        if options["cubeprop_tasks"][0].lower() == "esp":
            cube_file_name = "ESP.cube"
        elif options["cubeprop_tasks"][0].lower() == "density":
            cube_file_name = "Dt.cube"
        with TemporaryDirectory() as temp_dir:
            options["cubeprop_filepath"] = temp_dir

    if options is not None:
        psi4.set_options(options)

    for i in range(R.shape[0]):
        mol = psi4.core.Molecule.from_arrays(
            geom=R[i],
            elez=Z,
            molecular_charge=charge,
            molecular_multiplicity=mult,
            fix_com=True,
            fix_orientation=True,
        )
        wfn = None
        for task in tasks:
            data_worker, wfn = _do_psi4_task(
                task, data_worker, i, method, mol, options, cube_file_name, wfn
            )
    return data_worker
