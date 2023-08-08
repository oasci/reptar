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

import argparse
import os
import numpy as np
from .writing_utils import string_xyz_arrays
from .. import File
from ..utils import parse_xyz, atoms_by_number
from ..logger import ReptarLogger

log = ReptarLogger(__name__)


def write_xyz(xyz_path, Z, R, comments=None, data_precision=10):
    r"""Write standard XYZ file.

    Parameters
    ----------
    xyz_path : :obj:`str`
        Path to XYZ file to write.
    Z : :obj:`numpy.ndarray`, ndim: ``1``
        Atomic numbers of all atoms in the system.
    R : :obj:`numpy.ndarray`, ndim: ``3``
        Cartesian coordinates of all structures in the same order as ``Z``.
    comments : :obj:`list`, default: ``None``
        Comment lines for each XYZ structure.
    data_precision : :obj:`int`, default: ``10``
        Number of decimal points for printing array data.

    """
    if R.ndim == 2:
        R = R[None, ...]

    n_atoms = len(Z)
    with open(xyz_path, "w", encoding="utf-8") as f:
        for i, r in enumerate(R):
            f.write(f"{n_atoms}\n")
            if comments is not None:
                comment = comments[i]
                if comment[-2:] != "\n":
                    comment += "\n"
            else:
                comment = "\n"
            f.write(comment)
            f.write(string_xyz_arrays(Z, r, precision=data_precision))


def run_xyz_to_file(xyz_path, save_path, group_key, overwrite=False):
    log.info("Opening data file")
    if os.path.exists(save_path):
        if overwrite:
            rfile = File(save_path, mode="a")
        else:
            raise RuntimeError(f"{save_path} exists and overwrite is False")
    else:
        rfile = File(save_path, mode="w")

    log.info("Parsing XYZ file")
    Z, _, R = parse_xyz(os.path.abspath(xyz_path))
    Z = atoms_by_number(Z[0])
    Z = np.array(Z, dtype="uint8")
    R = np.array(R, dtype="float64")

    log.info("Storing coordinates")
    rfile.put(os.path.join(group_key, "atomic_numbers"), Z)
    rfile.put(os.path.join(group_key, "geometry"), R)


# pylint: disable=duplicate-code
def main_xyz_to_file():
    parser = argparse.ArgumentParser(
        description="Make a reptar file from XYZ",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "xyz_path",
        type=str,
        nargs="?",
        help="Path to XYZ file",
    )
    parser.add_argument(
        "--save_path",
        type=str,
        nargs="?",
        default="data.zarr",
        help="File name",
    )
    parser.add_argument(
        "--group_key",
        type=str,
        nargs="?",
        default="./",
        help="Group key to store in file",
    )
    parser.add_argument(
        "--overwrite",
        default=False,
        action="store_true",
        help="If file exists, overwrite all data",
    )
    parser.add_argument(
        "--log_level",
        type=str,
        nargs="?",
        default="info",
        help="Desired logging level",
    )
    args = parser.parse_args()
    run_xyz_to_file(args.xyz_path, args.save_path, args.group_key, args.overwrite)


def run_write_xyz(group_path, group_key, comment_key, save_dir):
    group_path = os.path.abspath(group_path)

    if ".exdir" in group_path:
        idx_split = group_path.rfind(".exdir") + 6
    elif ".zarr" in group_path:
        idx_split = group_path.rfind(".zarr") + 5
    elif ".json" in group_path:
        idx_split = group_path.rfind(".json") + 5
    elif ".npz" in group_path:
        idx_split = group_path.rfind(".npz") + 4
    else:
        raise ValueError("File type is not recognized.")

    if group_key == "":
        # Need to find parse group_key
        file_path, group_key = group_path[:idx_split], group_path[idx_split:]
    else:
        file_path = group_path

    print("Collecting data")
    rfile = File(file_path, mode="r")
    Z = rfile.get(os.path.join(group_key, "atomic_numbers"))
    R = rfile.get(os.path.join(group_key, "geometry"))
    if comment_key != "":
        comments = rfile.get(os.path.join(group_key, comment_key))
    else:
        comments = None

    if isinstance(comments, np.ndarray):
        print("Handling comment data")
        if comments.ndim > 1:
            raise ValueError(
                f"Comments has {comments.ndim} dimensions, but it must be 1."
            )

        comments = [str(i) for i in comments]
    print("Writing XYZ file")
    write_xyz(os.path.join(save_dir, "data.xyz"), Z, R, comments=comments)


def main_write_xyz():

    parser = argparse.ArgumentParser(
        description="Write XYZ file from a reptar supported file type."
    )
    parser.add_argument(
        "group_path",
        type=str,
        nargs="?",
        help="Path to group path. This can be to a file or exdir/zarr group.",
    )
    parser.add_argument(
        "--group_key",
        type=str,
        nargs="?",
        default="",
        help="Manually specify group key for file.",
    )
    parser.add_argument(
        "--comment_key",
        type=str,
        nargs="?",
        default="",
        help="Label of data to include as xyz comment.",
    )
    parser.add_argument(
        "--save_dir",
        type=str,
        nargs="?",
        default=".",
        help="Directory to save XYZ file.",
    )

    args = parser.parse_args()
    run_write_xyz(args.group_path, args.group_key, args.comment_key, args.save_dir)
