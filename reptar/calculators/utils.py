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
from typing import Any
from collections.abc import Iterable
import numpy as np
from . import Data
from .cube import initialize_grid_arrays
from ..logger import ReptarLogger

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


# pylint: disable=invalid-name
def prep_calc_data(
    tasks: Iterable[str],
    rfile: "File",
    source_key: str,
    source_labels: dict[str, str],
    dest_key: str,
    dest_labels: dict[str, str],
) -> Data:
    """Prepare group and data for reptar calculations

    Parameters
    ----------
    tasks
        Reptar calculations that will be ran with :class:`~reptar.calculators.Driver`.
    rfile
        File to prepare a group for optimization-like data.
    source_key
        Key to group containing data sources.
    source_labels
        Labels of data from group ``source_key`` to populate the data object. No
        ``_key`` properties in :class:`~reptar.calculators.Data` will be added from
        this ``dict``.
    dest_key
        Key to store calculation results.
    dest_labels
        Labels of data in group ``dest_key`` to populate the data object. Any data
        that is not provided will be initialized here based on ``tasks``.
        ``_key`` properties in :class:`~reptar.calculators.Data` will be added from
        this ``dict``.

    Returns
    -------
    :obj:`reptar.calculators.Data`
        Data for calculations.
    """
    data = Data()
    data.rfile = rfile

    log.debug("Retrieving source data")
    for data_attr, label in source_labels.items():
        data_key = os.path.join(source_key, label)
        value = rfile.get(data_key)
        setattr(data, data_attr, value)
    data.validate(None)  # Checks for Z and R

    log.debug("Checking calculation data")
    try:
        rfile.get(dest_key)
    except RuntimeError as e:
        if " does not exist" in str(e):
            rfile.create_group(dest_key)
        else:
            raise RuntimeError from e

    for data_attr, label in dest_labels.items():
        data_key = os.path.join(dest_key, label)
        setattr(data, data_attr + "_key", data_key)

        # If included from source, then we may have previous data already.
        # Check before initializing.
        if getattr(data, data_attr) is None:
            data.initialize_array(data_attr)
    data.validate(tasks)
    data.save()
    return data


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
