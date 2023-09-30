# MIT License
#
# Copyright (c) 2022-2023, OASCI
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

from typing import Any

import os
from collections.abc import Iterable

import numpy as np

from ..logger import ReptarLogger
from . import Data
from .cube import initialize_grid_arrays

log = ReptarLogger(__name__)


def _prep_xtb_opt_block(block: dict[str, Any]) -> list[str]:
    """Prepare xTB constrain block lines for input file.

    Parameters
    ----------
    block
        Keywords for the optimization block in xTB.

    Returns
    -------
    :obj:`list`
        Lines for the input file. Does not contain newline characters.
    """
    lines = []
    lines.append("$opt")
    for key, value in block.items():
        lines.append(f"    {key}={value}")
    lines.append("$end")
    return lines


def _prep_xtb_constrain_block(constraints):
    """Prepare xTB constrain block lines for input file.

    Parameters
    ----------
    constraints : :obj:`tuple`
        Internal constraints in a nested tuple to add to the xtb input file. The first
        element is the label (e.g., ``"distance"`` or ``"dihedral"``) and the second
        element is a list containing the values after  ``:`` used to control xtb.
        Atom indices should start from ``0``. For more information, see
        https://xtb-docs.readthedocs.io/en/latest/xcontrol.html#constraining-potentials

    Returns
    -------
    :obj:`list`
        Lines for the input file. Does not contain newline characters.
    """
    lines = []
    lines.append("$constrain")
    for constraint in constraints:
        constraint_name = constraint[0]
        if constraint_name == "force constant":
            constraint_line = f"    force constant={constraint[1]}"
        else:
            constraint_line = f"    {constraint[0]}: "
            # xtb indices start at 1
            constraint_line += ", ".join([str(i + 1) for i in constraint[1][:-1]])
            constraint_line += f", {constraint[1][-1]}"
        lines.append(constraint_line)
    lines.append("$end")
    return lines


def prep_xtb_input_lines(
    charge: int,
    multiplicity: int,
    opt_block: dict | None = None,
    constraints: tuple[str] | None = None,
    save_traj: bool = False,
) -> list[str]:
    """Prepare lines for xtb input file.

    Parameters
    ----------
    charge : :obj:`int`
        Total charge of the system.
    multiplicity : :obj:`int`
        Spin state multiplicity of the system.
    opt_block : :obj:`dict`
        Optimization block keywords and values.
    constraints : :obj:`tuple`
        Internal constraints in a nested tuple to add to the xtb input file. The first
        element is the label (e.g., ``"distance"`` or ``"dihedral"``) and the second
        element is a list containing the values after  ``:`` used to control xtb.
        Atom indices should start from ``0``. For more information, see
        https://xtb-docs.readthedocs.io/en/latest/xcontrol.html#constraining-potentials

    Returns
    -------
    :obj:`list`
        Lines of xtb input file.
    """
    spin = int(multiplicity - 1)
    xtb_input_lines = [f"$chrg {charge}", f"$spin {spin}"]
    if opt_block is not None:
        opt_lines = _prep_xtb_opt_block(opt_block)
        xtb_input_lines.extend(opt_lines)
    if constraints is not None:
        constraint_lines = _prep_xtb_constrain_block(constraints)
        xtb_input_lines.extend(constraint_lines)
    if save_traj:
        xtb_input_lines.extend(["$opt", "    logfile=xtbopt.trj", "$end"])
    xtb_input_lines = [i + "\n" for i in xtb_input_lines]
    return xtb_input_lines


def cleanup_xtb_calc(work_dir: str = "./") -> None:
    """Remove xTB files that are not commonly needed.

    Parameters
    ----------
    work_dir : :obj:`str`, default: ``./``
        Work directory to look for files to remove.
    """
    for tmp_file in [
        "charges",
        "wbo",
        "xtbopt.log",
        "xtbopt.xyz",
        "xtbrestart",
        "xtbtopo.mol",
        ".xtboptok",
    ]:
        try:
            os.remove(os.path.join(work_dir, tmp_file))
        except FileNotFoundError:
            pass


def initialize_worker_data(
    tasks: Iterable[str],
    data: Data,
    R: np.ndarray,
    total_grid_points: int | None = None,
) -> Data:
    if ("E" in tasks) or ("G" in tasks) or ("opt" in tasks):
        data.E = np.full(R.shape[0], np.nan, dtype=np.float64)
    if "G" in tasks:
        data.G = np.full(R.shape, np.nan, dtype=np.float64)
    if "opt" in tasks:
        data.conv_opt = np.full(R.shape[0], False, dtype=np.bool8)
        data.R_opt = np.full(R.shape, np.nan, dtype=np.float64)
    if "cube" in tasks:
        data.cube_R, data.cube_V = initialize_grid_arrays(
            R, max_points=total_grid_points
        )
    return data
