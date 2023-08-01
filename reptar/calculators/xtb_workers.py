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
from collections.abc import Iterator
import os
import subprocess
from tempfile import TemporaryDirectory
import numpy as np
from . import Data
from .utils import initialize_worker_data
from ..writers.xyz import write_xyz
from ..utils import parse_xyz
from ..logger import ReptarLogger

log = ReptarLogger(__name__)

try:
    from xtb.interface import Calculator, Param

    _HAS_XTB = True
except ImportError:
    _HAS_XTB = False


def xtb_python_worker(
    idxs: Iterator[int],
    tasks: Iterator[str],
    data: Data,
    charge: int = 0,
    mult: int = 1,
    acc: float = 0.1,
    max_iterations: int = 300,
    params: None = None,
) -> Data:
    r"""Ray remote function for computing total electronic energy and atomic
    gradients using xtb.

    Parameters
    ----------
    idxs
        Indices of the structures from ``R`` to compute energies and gradients
        for.
    tasks
        Calculations this worker needs to run in the order specified here. In general,
        we recommend using this ordering: ``E`` or ``G``.
    Z
        Atomic numbers of the atoms with respect to ``R``.
    R
        Cartesian coordinates of all structures in group. This includes
        unused structures.
    charge
        Total molecular charge.
    mult
        Total multiplicity.
    acc
        Numerical accuracy for calculation. For more information, visit the
        `documentation <https://xtb-python.readthedocs.io/en/latest/
        general-api.html#xtb.interface.Calculator.set_accuracy>`__.
    max_iterations
        Maximum number of iterations for self-consistent charge methods. If the
        calculations fails to converge in a given number of cycles, no error
        is necessarily shown.
    params
        xTB parameters. Defaults to ``xtb.interface.Param.GFN2xTB`` if ``None``.

    Returns
    -------
    :obj:`reptar.calculators.Data`
        Calculation results from this worker
    """
    if not _HAS_XTB:
        raise EnvironmentError("xtb-python not installed for this worker")

    implemented_tasks = ["E", "G"]
    for task in tasks:
        if task not in implemented_tasks:
            raise ValueError(f"Task ({task}) is not implemented in this worker")

    n_upair_ele = int(mult - 1)
    Z = data.Z
    R = data.R[idxs]

    data_worker = initialize_worker_data(tasks, Data(), R, None)
    data_worker.idxs_source = idxs

    if params is None:
        params = Param.GFN2xTB
    for i in range(R.shape[0]):
        calc = Calculator(params, Z, R[i], charge, n_upair_ele)
        calc.set_accuracy(acc)
        calc.set_max_iterations(max_iterations)
        res = calc.singlepoint()
        g = res.get_gradient()
        g /= 0.52917721067  # psi4.constants.bohr2angstroms
        data_worker.add_subset("E", np.array([i]), res.get_energy())
        data_worker.add_subset("G", np.array([i]), g)
    return data_worker


def _do_xtb_task(
    task: str, data_worker: Data, idx: int, xtb_command: list[str], output_path: str
) -> Data:
    worker_idx = np.array([idx], dtype=np.uint64)
    if task == "opt":
        xtb_command.append("--opt")
        with open(output_path, "w", encoding="utf-8") as f_out:
            subprocess.run(xtb_command, check=False, shell=False, stdout=f_out)

        _, comments, r_opt = parse_xyz("xtbopt.xyz")
        e = float(comments[0].split()[1])

        r_conv_opt = False
        with open(output_path, "r", encoding="utf-8") as f_out:
            for line in reversed(list(f_out)):
                if "*** GEOMETRY OPTIMIZATION CONVERGED AFTER" in line:
                    r_conv_opt = True
                    break
        r_opt = np.array(r_opt)[0]

        data_worker.add_subset("R_opt", worker_idx, r_opt)
        data_worker.add_subset("conv_opt", worker_idx, r_conv_opt)
        data_worker.add_subset("E", worker_idx, e)
    return data_worker


def xtb_worker(
    idxs: Iterator[int],
    tasks: Iterator[str],
    data: Data,
    input_lines: Iterator[str],
    acc: float = 0.1,
    n_cores: int = 1,
    xtb_path: str = "xtb",
    log_dir: str | None = None,
):
    r"""Ray remote function for calculations with xTB.

    Parameters
    ----------
    idxs
        Indices of the structures from ``R`` to compute energies and gradients
        for.
    tasks
        Calculations this worker needs to run in the order specified here. In general,
        we recommend using this ordering: ``opt``.
    data
        All required data required for computations such as ``Z``, ``R``, ``conv_opt``,
        etc.
    input_lines
        Lines for xTB input file. Must at least include ``$chrg`` and ``$spin`` blocks.
    acc
        Numerical accuracy for calculation. For more information, visit the
        `documentation <https://xtb-docs.readthedocs.io/en/latest/sp.html#accuracy>`__.
    n_cores
        Number of cores to use for xTB calculation.
    xtb_path
        Path to xtb executable to use. Defaults to assuming ``xtb`` is in your path.
    log_dir
        Work directory for the xtb calculations. If nothing is specified, no logs are
        stored.

    Returns
    -------
    :obj:`reptar.calculators.Data`
        Calculation results from this worker
    """
    # pylint: disable=consider-using-with
    implemented_tasks = ["opt"]
    for task in tasks:
        if task not in implemented_tasks:
            raise ValueError(f"Task ({task}) is not implemented in this worker")

    Z = data.Z
    R = data.R[idxs]

    # Initialize data object with resulting arrays
    data_worker = initialize_worker_data(tasks, Data(), R)
    data_worker.idxs_source = idxs

    log.debug("Setting up work and log directories")
    work_dir = TemporaryDirectory()
    if log_dir is not None:
        log_dir = os.path.abspath(log_dir)
        os.makedirs(log_dir, exist_ok=True)
    else:
        log_dir = work_dir.name
    cwd_path = os.getcwd()
    os.chdir(work_dir.name)

    xtb_input_path = "xtb-opt.in"
    with open(xtb_input_path, "w", encoding="utf-8") as f:
        f.writelines(input_lines)

    xtb_command = [
        xtb_path,
        "--input",
        xtb_input_path,
        "",  # xyz_path
        "--acc",
        str(acc),
        "--parallel",
        str(n_cores),
    ]

    log.debug("Starting xTB computations")
    for i in range(R.shape[0]):
        xyz_path = f"R_{i}.xyz"
        xtb_command[3] = xyz_path
        # Write temporary input file for coordinates
        write_xyz(xyz_path, Z, R[i])

        output_path = os.path.join(log_dir, f"{idxs[i]}.out")

        for task in tasks:
            data_worker = _do_xtb_task(task, data_worker, i, xtb_command, output_path)
    os.chdir(cwd_path)
    return data_worker
