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

import os
import subprocess
from tempfile import NamedTemporaryFile, TemporaryDirectory
import numpy as np
from ..writers.xyz import write_xyz
from ..utils import parse_xyz
from ..logger import ReptarLogger

log = ReptarLogger(__name__)

try:
    from xtb.interface import Calculator, Param

    _HAS_XTB = True
except ImportError:
    _HAS_XTB = False


def xtb_python_engrad(
    idxs, Z, R, charge=0, mult=1, acc=0.1, max_iterations=300, params=None
):
    r"""Ray remote function for computing total electronic energy and atomic
    gradients using xtb.

    Parameters
    ----------
    idxs : :obj:`numpy.ndarray`, ndim: ``1``
        Indices of the structures from ``R`` to compute energies and gradients
        for.
    Z : :obj:`numpy.ndarray`, ndim: ``1``
        Atomic numbers of the atoms with respect to ``R``.
    R : :obj:`numpy.ndarray`, ndim: ``3``
        Cartesian coordinates of all structures in group. This includes
        unused structures.
    charge : :obj:`int`, default: ``0``
        Total molecular charge.
    mult : :obj:`int`, default: ``1``
        Total multiplicity.
    acc : :obj:`int`, default: ``0.1``
        Numerical accuracy for calculation. For more information, visit the
        `documentation <https://xtb-python.readthedocs.io/en/latest/
        general-api.html#xtb.interface.Calculator.set_accuracy>`__.
    max_iterations : :obj:`int`, default: ``300``
        Maximum number of iterations for self-consistent charge methods. If the
        calculations fails to converge in a given number of cycles, no error
        is necessarily shown.
    params : default: ``None``
        xTB parameters. Defaults to ``xtb.interface.Param.GFN2xTB`` if ``None``.

    Returns
    -------
    :obj:`numpy.ndarray`
        ``idxs``
    :obj:`numpy.ndarray`
        Total electronic energy of computed structures in the same order as
        ``idxs``. Units of Hartree.
    :obj:`numpy.ndarray`
        Atomic gradients of computed structures in the same order as ``idxs``.
        Units of Hartree/Angstrom.
    """
    assert _HAS_XTB

    n_upair_ele = int(mult - 1)
    R = R[idxs]
    G = np.zeros(R.shape)
    E = np.zeros(R.shape[0])
    if params is None:
        params = Param.GFN2xTB
    for i, r in enumerate(R):
        calc = Calculator(params, Z, r, charge, n_upair_ele)
        calc.set_accuracy(acc)
        calc.set_max_iterations(max_iterations)
        res = calc.singlepoint()
        g = res.get_gradient()
        g /= 0.52917721067  # psi4.constants.bohr2angstroms
        G[i] = g
        E[i] = res.get_energy()
    return idxs, E, G


def xtb_opt(
    idxs,
    Z,
    R,
    input_lines,
    acc=0.1,
    n_cores=1,
    xtb_path="xtb",
    work_dir=None,
):
    r"""Ray remote function for computing total electronic energy and atomic
    gradients using xtb.

    Parameters
    ----------
    idxs : :obj:`numpy.ndarray`, ndim: ``1``
        Indices of the structures from ``R`` to compute energies and gradients
        for.
    Z : :obj:`numpy.ndarray`, ndim: ``1``
        Atomic numbers of the atoms with respect to ``R``.
    R : :obj:`numpy.ndarray`, ndim: ``3``
        Cartesian coordinates of all structures in group. This includes
        unused structures.
    input_lines : :obj:`list`
        Lines for xTB input file. Must at least include ``$chrg`` and ``$spin`` blocks.
    acc : :obj:`float`, default: ``0.1``
        Numerical accuracy for calculation. For more information, visit the
        `documentation <https://xtb-docs.readthedocs.io/en/latest/sp.html#accuracy>`__.
    n_cores : :obj:`int`, default ``1``
        Number of cores to use for xTB calculation.
    xtb_path : :obj:`str`, default: ``"xtb"``
        Path to xtb executable to use. Defaults to assuming ``xtb`` is in your path.
    work_dir : :obj:`str`, default: ``None``
        Work directory for the xtb calculations. If nothing is specified, the
        temporary directory is used.

    Returns
    -------
    :obj:`numpy.ndarray`
        ``idxs``
    :obj:`numpy.ndarray`
        If the optimizations converged or not.
    :obj:`numpy.ndarray`
        Optimized geometries.
    :obj:`numpy.ndarray`
        Total electronic energies of optimized structures. Units of Hartree.
    """
    # pylint: disable=consider-using-with
    log.debug("Initializing optimization arrays")
    log.debug("R array indices to do:")
    log.log_array(idxs, level=10)
    R = R[idxs]
    opt_conv = np.full(R.shape[0], False, dtype=np.bool8)
    R_opt = np.empty(R.shape, dtype=np.float64)  # pylint: disable=invalid-name
    E_opt = np.zeros(R.shape[0])  # pylint: disable=invalid-name

    log.debug("Setting up work directory and input files")
    if work_dir is None:
        temp_dir = TemporaryDirectory()
        work_dir = temp_dir.name
    cwd_path = os.getcwd()
    os.makedirs(work_dir, exist_ok=True)
    os.chdir(work_dir)

    xtb_input = NamedTemporaryFile(suffix=".in")
    with open(xtb_input.name, "w", encoding="utf-8") as f:
        f.writelines(input_lines)

    log.debug("Starting xTB computations")
    for i, r in enumerate(R):
        # Write temporary input file for coordinates
        xyz_input = NamedTemporaryFile(suffix=".xyz")
        write_xyz(xyz_input.name, Z, r)

        output_path = f"{idxs[i]}.out"

        xtb_command = [
            xtb_path,
            xyz_input.name,
            "--acc",
            str(acc),
            "--input",
            xtb_input.name,
            "--opt",
            "--parallel",
            str(n_cores),
        ]

        with open(output_path, "w", encoding="utf-8") as f_out:
            subprocess.run(xtb_command, check=False, shell=False, stdout=f_out)

        _, comments, r_opt = parse_xyz("xtbopt.xyz")
        e = float(comments[0].split()[1])

        r_opt_conv = False
        with open(output_path, "r", encoding="utf-8") as f_out:
            for line in reversed(list(f_out)):
                if "*** GEOMETRY OPTIMIZATION CONVERGED AFTER" in line:
                    r_opt_conv = True
                    break
        r_opt = np.array(r_opt)[0]

        opt_conv[i] = r_opt_conv  # pylint: disable=used-before-assignment
        R_opt[i] = r_opt
        E_opt[i] = e  # pylint: disable=used-before-assignment
    for tmp_file in [
        "charges",
        "wbo",
        "xtbopt.log",
        "xtbopt.xyz",
        "xtbrestart",
        "xtbtopo.mol",
        ".xtboptok",
    ]:
        os.remove(tmp_file)
    os.chdir(cwd_path)
    return idxs, opt_conv, R_opt, E_opt
