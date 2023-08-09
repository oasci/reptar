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

"""Interact with dihedral angles of structures."""

from __future__ import annotations
import argparse
import itertools
import os
import numpy as np
from .angles import get_angle_mask, sample_angles
from .. import File
from ..writers import write_xyz
from ..calculators import Data
from ..calculators.run import run_calcs
from ..utils import _load_config
from ..logger import ReptarLogger, set_log_level

log = ReptarLogger(__name__)


def gen_scan_range(range_info: dict[str, float]):
    range_scan = range(
        range_info["start"],
        range_info["stop"] + range_info["step"],
        range_info["step"],
    )
    return range_scan


def get_dihedral_info(
    info, n_atoms, angle_type, angle_types, angle_atoms, angle_ranges, angle_masks
):
    for i in info:
        angle_types.append(angle_type)
        angle_ranges.append(gen_scan_range(i["range"]))
        angle_atoms.append(i["atoms"])
        angle_masks.append(get_angle_mask(n_atoms, i["fragment"]))
    return angle_types, angle_atoms, angle_ranges, angle_masks


def _prep_angle_constraints(angle_type, config, constrain=None):
    if constrain is None:
        constrain = []
    angle_type = angle_type.strip().lower()
    if config["worker"]["path"] == "reptar.calculators.xtb_workers.xtb_worker":
        constrain.extend(
            [[angle_type, [*v["atoms"], "auto"]] for v in config["scan"][angle_type]]
        )
    return constrain


# pylint: disable=too-many-statements
def geometry_scan(config: dict, ray_address: str = "") -> None:

    log.info("Opening file")
    rfile_path = os.path.abspath(config["rfile"]["path"])
    rfile = File(rfile_path, mode="a")

    log.info("Getting source data")
    data = Data(rfile)
    source_info = config["rfile"]["source"]
    data.Z = rfile.get(os.path.join(source_info["key"], source_info["labels"]["Z"]))
    data.R = rfile.get(os.path.join(source_info["key"], source_info["labels"]["R"]))
    if data.R.ndim == 2:
        R = R[None, ...]
    data.R = data.R[source_info["R_slice"]][None, ...]  # Slice a single structure
    if data.R.shape[0] != 1:
        raise ValueError(
            f"R can only contain one structure, but contains {data.R.shape[0]}"
        )

    log.info("Defining grid")
    scan_info = config["scan"]

    angle_types = []
    angle_atoms = []
    angle_ranges = []
    angle_masks = []
    n_atoms = data.R.shape[1]
    constraints = config["worker"]["constrain"]
    for angle_type in ["dihedral", "angle"]:
        if angle_type in scan_info.keys():
            angle_types, angle_atoms, angle_ranges, angle_masks = get_dihedral_info(
                scan_info[angle_type],
                n_atoms,
                angle_type,
                angle_types,
                angle_atoms,
                angle_ranges,
                angle_masks,
            )
            constraints = _prep_angle_constraints(angle_type, config, constraints)
    angle_values = np.array(tuple(itertools.product(*angle_ranges)), dtype=np.float64)
    if len(angle_types) != 0:
        log.info("Generating %i geometries", int(len(angle_values)))
        t_start = log.t_start()
        data.R = sample_angles(
            data.Z,
            data.R,
            angle_types,
            angle_atoms,
            angle_values,
            angle_masks,
            use_ray=config["driver"]["kwargs"]["use_ray"],
            n_workers=int(
                config["driver"]["kwargs"]["n_workers"]
                * config["driver"]["kwargs"]["n_cpus_per_worker"]
            ),
        )
        log.t_stop(t_start, precision=2)

    log.info("Saving structures to file")
    dest_info = config["rfile"]["destination"]
    data.Z_key = os.path.join(dest_info["key"], dest_info["labels"]["Z"])
    data.R_key = os.path.join(dest_info["key"], dest_info["labels"]["R_opt"])
    data.save()
    del data

    # Need to update values
    if config["do_optimizations"]:
        log.info("Setting up geometry optimizations configuration")

        config_calc = config.copy()
        config_calc["tasks"] = ["opt"]
        config_calc["rfile"]["source"]["key"] = dest_info["key"]
        del config_calc["rfile"]["source"]["R_slice"]
        del config_calc["scan"]

        config_calc["worker"]["constrain"] = constraints

        data = run_calcs(config_calc, ray_address)
    else:
        log.info("No geometry optimizations will be performed")

    if isinstance(config["write_xyz"], str):
        if config["do_optimizations"]:
            R = data.R_opt
            comments = [str(e) for e in data.E]
        else:
            R = data.R
            comments = None
        write_xyz(
            os.path.abspath(config["write_xyz"]),
            data.Z,
            R,
            comments=comments,
        )


# pylint: disable=duplicate-code
def main():
    parser = argparse.ArgumentParser(
        description="Scan dihedral angles with an optional geometry optimization"
    )
    parser.add_argument(
        "config_path",
        type=str,
        nargs="?",
        help="Path to YAML configuration file",
    )
    parser.add_argument(
        "--ray_address",
        type=str,
        nargs="?",
        default="",
        help="Desired ray address (will override config file)",
    )
    parser.add_argument(
        "--log_level",
        type=str,
        nargs="?",
        default="info",
        help="Desired logging level",
    )

    args = parser.parse_args()
    set_log_level(args.log_level.upper())
    config = _load_config(args.config_path)
    geometry_scan(config, ray_address=args.ray_address)


if __name__ == "__main__":
    main()
